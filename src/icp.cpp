// File: icp_overlap_cloudcompare.cpp

#include <RcppEigen.h>
#include <vector>
#include <limits>
#include <iostream>
#include <numeric>
#include "nanoflann.hpp"

struct PointCloudAdaptor
{
  const std::vector<Eigen::Vector3d>& pts;
  PointCloudAdaptor(const std::vector<Eigen::Vector3d>& points) : pts(points) {}

  inline size_t kdtree_get_point_count() const { return pts.size(); }
  inline double kdtree_get_pt(const size_t idx, int dim) const { return pts[idx](dim); }

  template<class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudAdaptor>, PointCloudAdaptor, 3>;

// [[Rcpp::export]]
Eigen::MatrixXd cpp_icp(Eigen::MatrixXd& source_mat, Eigen::MatrixXd& target_mat,
                        bool tz_only = false, bool rz_only = true,
                        int max_iterations = 100, int overlap = 100, double tolerance = 1e-5)
{
  std::vector<Eigen::Vector3d> source, target;
  for (int i = 0; i < source_mat.rows(); ++i) source.emplace_back(source_mat(i, 0), source_mat(i, 1), source_mat(i, 2));
  for (int i = 0; i < target_mat.rows(); ++i) target.emplace_back(target_mat(i, 0), target_mat(i, 1), target_mat(i, 2));

  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  std::vector<Eigen::Vector3d> aligned = source;

  PointCloudAdaptor adaptor(target);
  KDTree index(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  Eigen::Matrix4d previous_T = T;
  double previous_error = std::numeric_limits<double>::max();

  for (int iter = 0; iter < max_iterations; ++iter)
  {
    struct Match { size_t src_idx; size_t tgt_idx; double dist; };
    std::vector<Match> matches;

    for (size_t i = 0; i < aligned.size(); ++i)
    {
      size_t idx;
      double dist;
      nanoflann::KNNResultSet<double> result(1);
      result.init(&idx, &dist);
      index.findNeighbors(result, aligned[i].data(), nanoflann::SearchParameters(1));
      matches.push_back({i, idx, dist});
    }

    std::sort(matches.begin(), matches.end(), [](const Match& a, const Match& b) { return a.dist < b.dist; });
    size_t n = matches.size() * overlap / 100.f;
    if (n == 0) n = 1;

    Eigen::Vector3d c_src = Eigen::Vector3d::Zero();
    Eigen::Vector3d c_dst = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < n; ++i)
    {
      c_src += aligned[matches[i].src_idx];
      c_dst += target[matches[i].tgt_idx];
    }
    c_src /= n;
    c_dst /= n;

    Eigen::MatrixXd src_mat(3, n);
    Eigen::MatrixXd dst_mat(3, n);
    for (size_t i = 0; i < n; ++i)
    {
      src_mat.col(i) = aligned[matches[i].src_idx] - c_src;
      dst_mat.col(i) = target[matches[i].tgt_idx] - c_dst;
    }

    Eigen::Matrix3d H = src_mat * dst_mat.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixV() * svd.matrixU().transpose();

    if (R.determinant() < 0)
    {
      Eigen::Matrix3d V = svd.matrixV();
      V.col(2) *= -1;
      R = V * svd.matrixU().transpose();
    }

    Eigen::Vector3d t = c_dst - R * c_src;

    Eigen::Vector3d t_locked = t;
    Eigen::Matrix3d R_locked = R;

    if (tz_only)
    {
      t_locked.x() = 0.0;
      t_locked.y() = 0.0;
      R_locked = Eigen::Matrix3d::Identity();
    }

    if (rz_only)
    {
      double theta = atan2(R(1, 0), R(0, 0));
      R_locked = Eigen::Matrix3d::Identity();
      R_locked(0, 0) =  cos(theta);
      R_locked(0, 1) = -sin(theta);
      R_locked(1, 0) =  sin(theta);
      R_locked(1, 1) =  cos(theta);
      R_locked(2, 2) =  1.0;
    }

    Eigen::Matrix4d T_step = Eigen::Matrix4d::Identity();
    T_step.block<3, 3>(0, 0) = R_locked;
    T_step.block<3, 1>(0, 3) = t_locked;
    T = T_step * T;

    for (size_t i = 0; i < aligned.size(); ++i)
    {
      aligned[i] = (T_step.block<3,3>(0,0) * aligned[i] + T_step.block<3,1>(0,3));
    }

    double total_error = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
      total_error += matches[i].dist;
    }
    total_error /= n;

    //Rprintf("Iteration %d : n = %lu, error %lf\n", iter, n, total_error);

    if (std::abs(previous_error - total_error) < tolerance)
    {
      break;
    }
    previous_error = total_error;
  }

  return T;
}
