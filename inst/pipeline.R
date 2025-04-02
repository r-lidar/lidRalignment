library(lidR)
library(lidRalignment)

# ==== USER CONFIGURATION ====

fref = "als_file.las"
fmov = "mls_file.las"

radius = 20
overlap = NA

ref_is_tls = FALSE
mov_is_tls = TRUE

cc = find_cloudcompare()

# ==== AUTOMATIC CONFIGURATION ====

set_lidr_threads(0)

reader_ref = readALS
filter_ref = ""
csf_ref = csf(rigidness = 2)

reader_mov = readTLS
filter_mov = "-keep_random_fraction 0.2"
csf_mov = csf(rigidness = 3, class_threshold = 0.05, cloth_resolution = 0.25)

strategy = "chm-dtm"
res = 2

if (!ref_is_tls & is.na(overlap)) overlap = 90
if (ref_is_tls & is.na(overlap)) overlap = 30
if (ref_is_tls) { reader_ref = readTLS ; csf_ref = csf_mov ; res = 0.5}
if (ref_is_tls & mov_is_tls) { filter_ref = filter_mov = filter = "-keep_random_fraction 0.1" ; strategy = "mls-tls" }

# ==== LOAD THE DATA ====

las_ref = reader_ref(fref, select = "", filter = filter_ref)
las_mov = reader_mov(fmov, select = "", filter = filter_mov)

# ==== DATA PREPARATION ====

global_shift_x = mean(las_ref$X)
global_shift_y = mean(las_ref$Y)

center_x = mean(las_mov$X)
center_y = mean(las_mov$Y)

# Clip a 20 m radius. This is sufficient and ensures both datasets have the same size.
las_ref = clip_circle(las_ref, global_shift_x, global_shift_y, radius)
las_mov = clip_circle(las_mov, center_x, center_y, radius)

# Remove outlier noise to ensure that the CHM and DTM are meaningful.
las_ref = classify_noise(las_ref, ivf(1))
las_mov = classify_noise(las_mov, ivf(1))
las_ref = remove_noise(las_ref)
las_mov = remove_noise(las_mov)

# Classify ground points to compute a DTM.
las_ref = classify_ground(las_ref, csf_ref, last_returns = F)
las_mov = classify_ground(las_mov, csf_mov, last_returns = F)

# Now that noise has been removed and the ground is classified, compute the Z offset to align the point clouds along the Z-axis.
# We use the minimum Z value of the ground points (more robust than the absolute minimum Z).
global_shift_z = min(filter_ground(las_ref)$Z)
center_z = min(filter_ground(las_mov)$Z)

# We cannot align using all points; we need to extract alignable features.
ref = extract_features(las_ref, strategy = strategy)
mov = extract_features(las_mov, strategy = strategy)

# Translate the point clouds to center them at (0,0,0) regardless of the original coordinate system.
ref = translate_las(ref, -global_shift_x, -global_shift_y, -global_shift_z)
mov = translate_las(mov, -center_x, -center_y, -center_z)

# Display the alignment.
show_alignment(ref, mov, size = 3)

# Notify when this part of the script is complete.
if (interactive()) beepr::beep(1)

# ==== ALIGNMENT PROCESS ====

# The first step is a coarse brute-force alignment.
# ICP cannot align the point clouds if they are too misaligned.
# The brute-force alignment rotates along Z and translates along XY.
# This is sufficient to achieve a rough initial registration.
M0 = brute_force_registration(ref, mov, res = res)
if (interactive()) show_alignment(ref, mov, M0, size = 3)

# Apply the initial transformation to the moving point cloud.
mov2 = transform_las(mov, M0)

# Perform a finer alignment using ICP with CloudCompare.
M1 = icp(ref, mov2, overlap = overlap, cc = cc)
if (interactive()) show_alignment(ref, mov2, M1, size = 3)

# Perform a final fine Z registration on ground points
Mz = diag(4)
if (!ref_is_tls)
{
  M = combine_transformations(M0, M1)
  ref_gnd = filter_ground(ref)
  mov_gnd = filter_ground(mov)
  mov_gnd = transform_las(mov_gnd, M)
  if (interactive()) show_alignment(ref_gnd, mov_gnd, size = 3)

  Mz = icp(ref_gnd, mov_gnd, overlap = overlap, skip_txy = TRUE, rot = "NONE", cc = cc)
  if (interactive()) show_alignment(ref_gnd, mov_gnd, Mz, size = 3)
}

# ==== FINAL REGISTRATION ====

# The moving point cloud was transformed as follows:
#  1. Initially centered at (0,0,0) in its own coordinate system.
#  2. Then aligned with the reference, which itself was centered at (0,0,0) in the global coordinate system.
# We need to combine these transformations and apply a final one to register the point cloud
# in the global coordinate system.
Mlocal  = translation_matrix(-center_x, -center_y, -center_z)
Mglobal = translation_matrix(global_shift_x, global_shift_y, global_shift_z)
Mregistration = combine_transformations(Mlocal, M, Mz, Mglobal)

# Apply the final transformation to register the full point cloud.
ofile = transform_las(fmov, Mregistration, sf::st_crs(ref))
