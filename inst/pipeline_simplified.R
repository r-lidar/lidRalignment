library(lidR)
library(lidRalignment)
options(lidR.progress = FALSE)

fref = "als_file.las"
fmov = "mls_file.las"

# Setup the point clouds to align
alignment = AlignmentScene$new(fref, fmov, verbose = TRUE)
alignment$set_ref_is_ground_based(TRUE)
alignment$set_mov_is_ground_based(TRUE)
alignment$set_radius(20)


alignment$prepare()
#alignment$plot("raw")

alignment$coarse_align(debug = T)
#alignment$plot("coarse")

alignment$fine_align()
#alignment$plot("fine")

alignment$extra_fine_align()
alignment$plot("extra", compare_to = "fine")

M = alignment$get_registration_matrix()
M

writeLAS(alignment$trunks_mov, "~/Téléchargements/trunk_mov.las")
writeLAS(alignment$trunks_ref, "~/Téléchargements/trunk_ref.las")

mov = readLAS("~/Téléchargements/trunk_mov.las")
ref = readLAS("~/Téléchargements/trunk_ref.las")

lidRalignment:::show_alignment(ref, mov, size = 3)

M = lidRalignment:::cc_icp(ref, mov, overlap = 30)
M = lidRalignment:::icp(ref, mov, overlap = 30)

lidRalignment:::show_alignment(ref, mov, M, size = 3)

# Apply the final transformation to register the entire point cloud.
crs = sf::st_crs(readLASheader(fref))
ofile = transform_las(fmov, M, crs)


# Setup the pipeline. It is important to tell the object
# what we are aligning in order to perform or not
# the last extra fine alignment
alignment = AlignmentScene$new(fref, fmov)
alignment$set_ref_is_ground_based(TRUE)
alignment$set_mov_is_ground_based(TRUE)

# Run the alignment pipeline
alignment$align()

# Visualize the different level of alignment
alignment$plot("raw")
alignment$plot("coarse")
alignment$plot("fine")
alignment$plot("extra", compare_to = "fine")
