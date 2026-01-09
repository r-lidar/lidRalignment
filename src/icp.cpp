#include <RcppEigen.h>
#include <vector>
#include <limits>
#include <iostream>
#include <numeric>
#include "nanoflann/nanoflann.h"


struct PointCloudAdaptor
{
  const std::vector<Eigen::Vector3d>& pts;
  PointCloudAdaptor(const std::vector<Eigen::Vector3d>& points) : pts(points) {}
  inline size_t kdtree_get_point_count() const { return pts.size(); }
  inline double kdtree_get_pt(const size_t idx, int dim) const { return pts[idx](dim); }
  template<class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudAdaptor>, PointCloudAdaptor, 3>;

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_icp(
    Eigen::MatrixXd& source_mat,
    Eigen::MatrixXd& target_mat,
    bool tz_only = false,
    bool rz_only = true,
   int max_iterations = 100,
   int overlap = 100,
   double tolerance = 1e-5)
{
  // 1. Convert Input Matrices to std::vector<Eigen::Vector3d>
  std::vector<Eigen::Vector3d> source, target;
  source.reserve(source_mat.rows());
  target.reserve(target_mat.rows());

  for (int i = 0; i < source_mat.rows(); ++i)
    source.emplace_back(source_mat(i, 0), source_mat(i, 1), source_mat(i, 2));
  for (int i = 0; i < target_mat.rows(); ++i)
    target.emplace_back(target_mat(i, 0), target_mat(i, 1), target_mat(i, 2));

  // 2. Initialize KDTree for Target
  PointCloudAdaptor adaptor(target);
  KDTree index(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  // 3. Initialize Transformation
  Eigen::Matrix4d T_global = Eigen::Matrix4d::Identity();
  std::vector<Eigen::Vector3d> aligned = source;

  double previous_rms = std::numeric_limits<double>::max();
  double rmsi = 0.0;
  double rms_final = 0.0;

  struct Match { size_t src_idx; size_t tgt_idx; double dist_sq; };

  // 4. ICP Loop
  for (int iter = 0; iter < max_iterations; ++iter)
  {
    // --- A. Find Nearest Neighbors ---
    std::vector<Match> matches;
    matches.reserve(aligned.size());

    for (size_t i = 0; i < aligned.size(); ++i)
    {
      size_t idx;
      double dist_sq;
      nanoflann::KNNResultSet<double> result(1);
      result.init(&idx, &dist_sq);
      index.findNeighbors(result, aligned[i].data(), nanoflann::SearchParameters(1));
      matches.push_back({i, idx, dist_sq});
    }

    // --- B. Filter by Overlap ---
    // Sort by distance to find the best % of matches
    std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) {
      return a.dist_sq < b.dist_sq;
    });

    size_t n = matches.size() * overlap / 100.0;
    if (n < 3) n = 3; // Minimum points for SVD

    // --- C. Compute Centroids ---
    Eigen::Vector3d c_src = Eigen::Vector3d::Zero();
    Eigen::Vector3d c_dst = Eigen::Vector3d::Zero();

    for (size_t i = 0; i < n; ++i)
    {
      c_src += aligned[matches[i].src_idx];
      c_dst += target[matches[i].tgt_idx];
    }
    c_src /= static_cast<double>(n);
    c_dst /= static_cast<double>(n);

    // --- D. Compute Cross-Covariance ---
    Eigen::MatrixXd src_centered(3, n);
    Eigen::MatrixXd dst_centered(3, n);

    for (size_t i = 0; i < n; ++i)
    {
      src_centered.col(i) = aligned[matches[i].src_idx] - c_src;
      dst_centered.col(i) = target[matches[i].tgt_idx] - c_dst;
    }

    Eigen::Matrix3d H = src_centered * dst_centered.transpose();

    // --- E. SVD for Rotation ---
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R_optimal = svd.matrixV() * svd.matrixU().transpose();

    // Fix reflection case (determinant must be 1, not -1)
    if (R_optimal.determinant() < 0)
    {
      Eigen::Matrix3d V = svd.matrixV();
      V.col(2) *= -1;
      R_optimal = V * svd.matrixU().transpose();
    }

    // --- F. Apply Logic Constraints ---
    Eigen::Matrix3d R_step = R_optimal;

    // Constraint 1: Rotation Constraints
    if (rz_only)
    {
      // Project full 3D rotation to 2D (Yaw) only
      double theta = atan2(R_optimal(1, 0), R_optimal(0, 0));
      R_step = Eigen::Matrix3d::Identity();
      R_step(0, 0) =  cos(theta);
      R_step(0, 1) = -sin(theta);
      R_step(1, 0) =  sin(theta);
      R_step(1, 1) =  cos(theta);
    }
    else if (tz_only)
    {
      // If Translation-Z only is requested AND Rotation-Z is NOT requested,
      // we assume strict vertical lift (no rotation).
      R_step = Eigen::Matrix3d::Identity();
    }

    // Constraint 2: Calculate Translation based on the DECIDED Rotation
    // This must happen AFTER R_step is finalized.
    Eigen::Vector3d t_step = c_dst - R_step * c_src;

    // Constraint 3: Translation Constraints
    if (tz_only)
    {
      // Force X and Y movement to zero
      t_step.x() = 0.0;
      t_step.y() = 0.0;
    }

    // --- G. Update Transformations ---

    // Update Global Matrix
    Eigen::Matrix4d T_step_mat = Eigen::Matrix4d::Identity();
    T_step_mat.block<3, 3>(0, 0) = R_step;
    T_step_mat.block<3, 1>(0, 3) = t_step;

    T_global = T_step_mat * T_global;

    // Update Points
    for (auto& pt : aligned)
    {
      pt = R_step * pt + t_step;
    }

    // --- H. Check Convergence (RMS) ---
    double current_rms = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
      current_rms += std::sqrt(matches[i].dist_sq); // matches store squared dist
    }
    current_rms /= static_cast<double>(n);

    if (iter == 0) rmsi = current_rms;
    rms_final = current_rms;

    if (std::abs(previous_rms - current_rms) < tolerance)
    {
      break;
    }
    previous_rms = current_rms;
  }

  // 5. Output Results
  Rcpp::NumericMatrix M = Rcpp::wrap(T_global);
  M.attr("RMSi") = rmsi;
  M.attr("RMSf") = rms_final;

  return M;
}