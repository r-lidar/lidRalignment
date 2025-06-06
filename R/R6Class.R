#' AlignmentScene R6 Class for Point Cloud Alignment
#'
#' @description
#' The `AlignmentScene` R6 class provides methods for aligning two 3D point clouds plots, typically
#' from terrestrial (MLS/TLS) or aerial (ALS/UAV) LiDAR sources. It implements a multi-step alignment
#' process with coarse, fine, and extra-fine alignment stages.
#'
#' @export
#' @include cloudcompare.R
AlignmentScene <- R6::R6Class("AlignmentScene",
  public = list(

    #cc = find_cloudcompare(),

    fref = NULL,
    fmov = NULL,

    full_ref = NULL,
    full_mov = NULL,
    chmdtm_ref = NULL,
    chmdtm_mov = NULL,
    trunks_ref = NULL,
    trunks_mov = NULL,

    M0 = diag(4),
    M1 = diag(4),
    Mz = diag(4),
    Mex = diag(4),
    Mlocal = diag(4),
    Mglobal = diag(4),

    radius = 20,
    ref_is_ground_based = TRUE,
    mov_is_ground_based = TRUE,

    prepare_done = FALSE,
    coarse_done = FALSE,
    fine_done = FALSE,
    extra_done = FALSE,

    verbose = TRUE,

    #' @description
    #' Create a new AlignmentScene object.
    #' @param fref,fmov character, path to the las laz file with the point cloud to align. fref is
    #' the reference point cloud, fmov is the point cloud to be aligned
    #' @param radius clip radius. The routine does not align the full point cloud but clip before
    #' to compute the registration matrix
    #' @param verbose logical
    #' @return A new `AlignmentScene` object.
    initialize = function(fref, fmov, radius = 20, verbose = TRUE)
    {
      self$verbose = TRUE

      href = lidR::readLASheader(fref)
      hmov = lidR::readLASheader(fmov)

      dref = lidR::density(href)
      dmov = lidR::density(hmov)

      if (dref < 200)
        self$set_ref_is_ground_based(FALSE)
      else if (dref > 4000)
        self$set_ref_is_ground_based(TRUE)
      else
      {
        self$set_ref_is_ground_based(TRUE)
        warning("The system was not able to guess if the reference file was ground-based (MLS/TLS) or air-based (ALS/UAV). It assumed ground-based. You can use set_ref_is_ground_based() to assign manually")
      }

      if (dmov < 200)
        self$set_mov_is_ground_based(FALSE)
      else if (dmov > 4000)
        self$set_mov_is_ground_based(TRUE)
      else
      {
        self$set_mov_is_ground_based(TRUE)
        warning("The system was not able to guess if the moving file was ground-based (MLS/TLS) or air-based (ALS/UAV). It assumed ground-based. You can use set_mov_is_ground_based() to assign manually")
      }

      if (self$ref_is_ground_based)
        message("Reference point cloud detected as ground-based (MLS/TLS) based on point density. If not correct, use set_ref_is_ground_based().")
      else
        message("Reference point cloud detected as air-based (ALS/ULS) based on point density. If not correct, use set_ref_is_ground_based().")

      if (self$mov_is_ground_based)
        message("Moving point cloud detected as ground-based (MLS/TLS) based on point density. If not correct, use set_mov_is_ground_based().")
      else
        message("Moving point cloud detected as air-based (ALS/ULS) based on point density. If not correct, use set_mov_is_ground_based().")

      self$fref = fref
      self$fmov = fmov
      self$radius = radius
    },

    #' @description
    #' print
    print = function() {
      cat("Files:\n")
      cat(" - Reference: ", self$fref, "\n")
      cat(" - Moving: ", self$fmov, "\n")

      cat("Lidar:\n")
      if (self$ref_is_ground_based)
        cat(" - Reference: TLS/MLS\n")
      else
        cat(" - Reference: ALS/ULS\n")

      if (self$mov_is_ground_based)
        cat(" - Moving: TLS/MLS\n")
      else
        cat(" - Moving: ALS/ULS\n")

      cat("Process:\n")
      if (self$prepare_done)
        cat(" 1. prepare: done\n")
      else
        cat(" 1. prepare: waiting\n")

      if (self$coarse_done)
        cat(" 2. coarse alignment: done\n")
      else
        cat(" 2. coarse alignment: waiting\n")

      if (self$fine_done)
        cat(" 3. fine alignment: done\n")
      else
        cat(" 3. fine alignment: waiting\n")

      if (self$extra_done)
        cat(" 4. extra fine alignment: done\n")
      else
        cat(" 4. extra fine alignment: waiting\n")
    },

    #' @description
    #' Set clip radius
    #' @param radius numeric.
    set_radius = function(radius)
    {
      if (self$prepare_done)
        stop("radius cannot be assigned if the function prepare() as already been called")

      self$radius = radius
    },

    #' @description
    #' Tell the pipeline if the reference point cloud is ground-based (TLS, MLS)
    #' @param val logical
    set_ref_is_ground_based = function(val = TRUE)
    {
      self$ref_is_ground_based = val
    },

    #' @description
    #' Tell the pipeline if the moving point cloud is ground-based (TLS, MLS)
    #' @param val logical
    set_mov_is_ground_based = function(val = TRUE)
    {
      self$mov_is_ground_based = val
    },

    #' @description
    #' First function to run. It read the point cloud, extract features, and perform raw alignment
    prepare = function() {

      filter_ref = ""
      select = "c"
      csf_ref = lidR::csf(rigidness = 2)

      filter_mov = "-keep_random_fraction 0.2"
      csf_mov = lidR::csf(rigidness = 3, class_threshold = 0.05, cloth_resolution = 0.25)

      if (self$ref_is_ground_based) { csf_ref = csf_mov }
      if (self$ref_is_ground_based & self$mov_is_ground_based) { filter_ref = filter_mov = filter = "-keep_random_fraction 0.1"}

      # ==== LOAD THE DATA ====

      if(self$verbose) cat("Reading reference point cloud...\n")
      full_ref = lidR::readTLS(self$fref, select = select, filter = filter_ref)

      if(self$verbose) cat("Reading point cloud to align...\n")
      full_mov = lidR::readTLS(self$fmov, select = select, filter = filter_mov)

      # ==== DATA PREPARATION ====

      if(self$verbose) cat("Checking for pre-classified ground points...\n")

      ref_has_gnd = sum(full_ref$Classification == 2L) > 0L
      mov_has_gnd = sum(full_mov$Classification == 2L) > 0L
      dtm_ref = NULL
      dtm_mov = NULL

      if (ref_has_gnd)
      {
        if(self$verbose) cat("  Reference point cloud has ground points...\n")
        dtm_ref = lidR::rasterize_terrain(full_ref, 0.5, lidR::tin())
        dtm_ref = lidR:::raster_as_las(dtm_ref)
        dtm_ref@data$Classification = 2L
        full_ref = suppressWarnings(rbind(full_ref, dtm_ref))
      }

      if (mov_has_gnd)
      {
        if(self$verbose) cat("  Moving point cloud has ground points...\n")
        dtm_mov = lidR::rasterize_terrain(full_mov, 1, lidR::tin())
        dtm_mov = lidR:::raster_as_las(dtm_mov)
        dtm_mov@data$Classification = 2L
        full_mov = suppressWarnings(rbind(full_mov, dtm_mov))
      }

      global_shift_x = mean(full_ref$X)
      global_shift_y = mean(full_ref$Y)

      center_x = mean(full_mov$X)
      center_y = mean(full_mov$Y)

      # Clip a 20 m radius. This is sufficient and ensures both datasets have the same size.
      if(self$verbose) cat("Clipping point clouds with a buffer...\n")
      full_ref = lidR::clip_circle(full_ref, global_shift_x, global_shift_y, self$radius+2) # +5m buffer for ground classification
      full_mov = lidR::clip_circle(full_mov, center_x, center_y, self$radius+2)

      # Remove outlier noise to ensure that the CHM and DTM are meaningful.
      if(self$verbose) cat("Cleaning noise points...\n")
      full_ref = lidR::classify_noise(full_ref, lidR::ivf(1))
      full_mov = lidR::classify_noise(full_mov, lidR::ivf(1))
      full_ref = lidR::remove_noise(full_ref)
      full_mov = lidR::remove_noise(full_mov)

      # Classify ground points to compute a DTM.
      if(self$verbose) cat("Classifying ground points...\n")

      if (!ref_has_gnd)
        full_ref = lidR::classify_ground(full_ref, csf_ref, last_returns = FALSE)
      else
        if(self$verbose) cat(" Skip reference point cloud.\n")

      if (!mov_has_gnd)
        full_mov = lidR::classify_ground(full_mov, csf_mov, last_returns = FALSE)
      else
        if(self$verbose) cat(" Skip moving point cloud.\n")

      if(self$verbose) cat("Clipping point clouds...\n")
      full_ref = lidR::clip_circle(full_ref, global_shift_x, global_shift_y, self$radius)
      full_mov = lidR::clip_circle(full_mov, center_x, center_y, self$radius)

      # Now that noise has been removed and the ground is classified, compute the Z offset to align the point clouds along the Z-axis.
      # We use the minimum Z value of the ground points (more robust than the absolute minimum Z).
      global_shift_z = min(lidR::filter_ground(full_ref)$Z)
      center_z = min(lidR::filter_ground(full_mov)$Z)

      # We cannot align using all points; we need to extract alignable features.
      cat("Extracting features to align from reference point cloud...\n")
      chmdtm_ref = extract_features(full_ref, strategy = "chm-dtm", verbose = self$verbose)
      cat("Extracting features to align from moving point cloud...\n")
      chmdtm_mov = extract_features(full_mov, strategy = "chm-dtm", verbose = self$verbose)

      # Translate the point clouds to center them at (0,0,0) regardless of the original coordinate system.
      chmdtm_ref = translate_las(chmdtm_ref, -global_shift_x, -global_shift_y, -global_shift_z)
      chmdtm_mov = translate_las(chmdtm_mov, -center_x, -center_y, -center_z)

      self$full_ref = full_ref
      self$full_mov = full_mov

      self$chmdtm_ref = chmdtm_ref
      self$chmdtm_mov = chmdtm_mov

      self$Mlocal  = translation_matrix(-center_x, -center_y, -center_z)
      self$Mglobal = translation_matrix(-global_shift_x, -global_shift_y, -global_shift_z)

      self$prepare_done = TRUE
    },

    #' @description
    #' Second function to run. It perform a brute force alignment
    #' @param res numeric. Subsampling resolution. Keep it as is.
    #' @param max_offset numeric. Maximum translation possible. Increase if the point cloud are
    #' strongly missaligned on XY.
    #' @param debug logical.
    coarse_align = function(res = 2, max_offset = 10, debug = FALSE)
    {
      if (!self$prepare_done)
        stop("coarse_align() must be called after prepare()")

      self$M0 = brute_force_registration(self$chmdtm_ref, self$chmdtm_mov, res = res, max_offset = max_offset, debug = debug, verbose = FALSE)
      self$coarse_done = TRUE
    },

    #' @description
    #' Third function to run. It perform an ICP alignment
    #' @param overlap numeric trimmed ICP overlap. Can be 10, 20, 30, 40 to 100.
    #' @param use_cc bool. Use CloudCompare ICP instead of native ICP.
    fine_align = function(overlap = "auto", use_cc = FALSE)
    {
      if (!self$coarse_done)
        stop("fine_align() must be called after coarse_align()")

      if(self$verbose) cat("  Iterative closest point fine alignment\n")

      ref = self$chmdtm_ref
      mov = self$chmdtm_mov

      if (overlap == "auto")
        overlap = adjust_overlap(90, self$chmdtm_ref, self$chmdtm_mov, self$M0)
      else
        stopifnot(overlap %in% c(1:10*10))

      if(self$verbose) cat("   overlap =", overlap, "\n")

      # Aligning only the CHM
      ref_chm = ref[ref$Classification != 2L]
      mov_chm = mov[mov$Classification != 2L]
      mov_chm = transform_las(mov_chm, self$M0)

      if (!use_cc)
        self$M1 = icp(ref_chm, mov_chm, rz_only = FALSE, overlap = overlap)
      else
        self$M1 = cc_icp(ref_chm, mov_chm, rot = "XYZ", overlap = overlap)

      if (!use_cc)
      {
        rmsi = attr(self$M1, "RMSi")
        rmsf = attr(self$M1, "RMSf")
        cat("    RMS initial", round(rmsi, 6), "\n")
        cat("    RMS final", round(rmsf, 6), "\n")
      }

      M = combine_transformations(self$M0, self$M1)

      if(self$verbose) cat("  Iterative closest point fine Z alignment\n")

      ref_gnd = lidR::filter_ground(ref)
      mov_gnd = lidR::filter_ground(mov)
      mov_gnd = transform_las(mov_gnd, M)

      if (!use_cc)
        self$Mz = icp(ref_gnd, mov_gnd, overlap = overlap, tz_only = TRUE)
      else
        self$Mz = cc_icp(ref_gnd, mov_gnd, overlap = overlap, skip_txy = TRUE, rot = "NONE", verbose = FALSE)

      if (!use_cc)
      {
        rmsi = attr(self$Mz, "RMSi")
        rmsf = attr(self$Mz, "RMSf")
        cat("    RMS initial", round(rmsi, 6), "\n")
        cat("    RMS final", round(rmsf, 6), "\n")
      }

      self$fine_done = TRUE
    },

    #' @description
    #' Fourth function to run. It perform an ICP alignment on the trunks
    extra_fine_align = function()
    {
      if (!self$fine_done)
        stop("extra_fine_align() must be called after fine_align()")

      if (!self$ref_is_ground_based | !self$mov_is_ground_based)
      {
        message("extra_fine_align() applies only on two ground-based point clouds")
        self$extra_done = FALSE
        return(invisible())
      }

      M = combine_transformations(self$M0, self$M1, self$Mz)

      cat("Extracting features to align from reference point cloud...\n")
      trunks_ref = extract_features(self$full_ref, strategy = "trunks")
      cat("Extracting features to align from moving point cloud...\n")
      trunks_mov = extract_features(self$full_mov, strategy = "trunks")

      trunks_ref = transform_las(trunks_ref, self$Mglobal)
      trunks_mov = transform_las(trunks_mov, self$Mlocal)

      mov2 = transform_las(trunks_mov, M)
      overlap = adjust_overlap(30, self$chmdtm_ref, self$chmdtm_mov, M)

      cat("  Iterative closest point extra fine alignment...\n")

      #self$Mex = cc_icp(trunks_ref, mov2, overlap = overlap, cc = self$cc, verbose = FALSE)
      self$Mex = icp(trunks_ref, mov2, overlap = overlap)
      rmsi = attr(self$Mex, "RMSi")
      rmsf = attr(self$Mex, "RMSf")
      cat("    RMS initial", round(rmsi, 6), "\n")
      cat("    RMS final", round(rmsf, 6), "\n")

      self$trunks_ref = trunks_ref
      self$trunks_mov = trunks_mov
      self$extra_done = TRUE
    },

    #' @description
    #' Full pipeline
    #' @param res numeric. Subsampling resolution. Keep it as is.
    #' @param max_offset numeric. Maximum translation possible. Increase if the point cloud are
    #' strongly missaligned on XY.
    align = function(res = 2, max_offset = 8)
    {
      self$prepare()
      cat("Coarse alignment...\n")
      self$coarse_align(res, max_offset)
      cat("Fine alignment...\n")
      self$fine_align()
      if (self$ref_is_ground_based & self$mov_is_ground_based)
        self$extra_fine_align()

      cat("Alignment completed!")
    },

    #' @description
    #' Plot alignment
    #' @param which string. The alignment level to plot
    #' @param compare_to string. The alignment level to plot and to compare
    plot = function(which = c("raw", "coarse", "fine", "extra"), compare_to = which)
    {
      which = match.arg(which, c("raw", "coarse", "fine", "extra"))
      compare_to = match.arg(compare_to, c("raw", "coarse", "fine", "extra"))

      if (which == "coarse" & !self$coarse_done) stop("Corse alignment not performed yet")
      if (which == "fine" & !self$fine_done) stop("Fine alignment not performed yet")
      if (which == "extra" & !self$extra_done) stop("Extra fine alignment not performed yet")

      if (compare_to != which)
      {
        M = diag(4)

        if (compare_to == "coarse") M = self$M0
        if (compare_to == "fine")   M = combine_transformations(self$M0, self$M1, self$Mz)
        if (compare_to == "extra")  M = combine_transformations(self$M0, self$M1, self$Mz, self$Mex)

        if (which == "extra")
          show_alignment(self$trunks_ref, self$trunks_mov, M, size = 3)
        else
          show_alignment(self$chmdtm_ref, self$chmdtm_mov, M, size = 3)
      }

      M = diag(4)

      if (which == "coarse") M = self$M0
      if (which == "fine")   M = combine_transformations(self$M0, self$M1, self$Mz)
      if (which == "extra")  M = combine_transformations(self$M0, self$M1, self$Mz, self$Mex)

      if (which == "extra")
        show_alignment(self$trunks_ref, self$trunks_mov, M, size = 3)
      else
        show_alignment(self$chmdtm_ref, self$chmdtm_mov, M, size = 3)
    },

    #' @description
    #' Get the registration matrix. This matrix is more or less accurate as a function of the level
    #' of registration performed
    get_registration_matrix = function()
    {
      M = combine_transformations(self$Mlocal, self$M0, self$M1, self$Mz, self$Mex, solve(self$Mglobal))
      return(M)
    }
  ))
