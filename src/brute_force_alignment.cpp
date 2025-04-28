#include <Rcpp.h>
#include <omp.h>
#include <vector>
#include <cmath>
#include "nanoflann.hpp"

// Simple KD-tree wrapper for nanoflann
struct PointCloud
{
  const float* data;
  size_t n_pts;

  inline size_t kdtree_get_point_count() const { return n_pts; }
  inline float kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    return data[idx * 3 + dim];
  }
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3> KDtree;

float pseudo_rms(const KDtree& tree, const float* mov_data, size_t n_pts)
{
  std::vector<float> Dsq;

  for (size_t i = 0; i < n_pts; ++i)
  {
    float query_pt[3] = {mov_data[i*3+0], mov_data[i*3+1], mov_data[i*3+2]};
    size_t ret_index;
    float out_dist_sqr;
    nanoflann::KNNResultSet<float> resultSet(1);
    resultSet.init(&ret_index, &out_dist_sqr);
    tree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(1));
    Dsq.push_back(out_dist_sqr);
  }

  // Sort distances to get the best 50% (smallest distances)
  std::sort(Dsq.begin(), Dsq.end());

  // Calculate RMS based on the best 50%
  size_t half_n_pts = n_pts / 2;
  float rms = 0.0f;
  for (size_t i = 0; i < half_n_pts; ++i)
  {
    rms += Dsq[i]; // Use the smallest half distances
  }

  return std::sqrt(rms / half_n_pts);
}

void transform(float* points, size_t n_pts, float angle, float dx, float dy, float dz)
{
  float cos_a = std::cos(angle);
  float sin_a = std::sin(angle);

  for (size_t i = 0; i < n_pts; ++i)
  {
    float x = points[i*3+0];
    float y = points[i*3+1];

    points[i*3+0] = cos_a * x - sin_a * y + dx;
    points[i*3+1] = sin_a * x + cos_a * y + dy;
    points[i*3+2] += dz;
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame rms_scan_grid(Rcpp::NumericMatrix ref, Rcpp::NumericMatrix mov, Rcpp::NumericMatrix param_grid)
{
  const size_t nref = ref.nrow();
  const size_t nmov = mov.nrow();
  const size_t np = param_grid.nrow();

  const size_t ref_ncol = ref.ncol();
  const size_t mov_ncol = mov.ncol();

  std::vector<float> ref_pts;
  ref_pts.reserve(nref * ref_ncol);

  for (size_t i = 0; i < nref; ++i) {
    for (size_t j = 0; j < ref_ncol; ++j) {
      ref_pts.push_back(static_cast<float>(ref(i, j)));
    }
  }

  std::vector<float> mov_pts;
  mov_pts.reserve(nmov * mov_ncol);

  for (size_t i = 0; i < nmov; ++i) {
    for (size_t j = 0; j < mov_ncol; ++j) {
      mov_pts.push_back(static_cast<float>(mov(i, j)));
    }
  }

  // Build the KDTree
  PointCloud ref_cloud{ref_pts.data(), nref};
  KDtree tree(3, ref_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree.buildIndex();

  std::vector<float> angles(np), dxs(np), dys(np), dzs(np), rmss(np);

  #pragma omp parallel for
  for (int i = 0; i < static_cast<int>(np); ++i)
  {
    float angle = param_grid(i, 0);
    float dx    = param_grid(i, 1);
    float dy    = param_grid(i, 2);
    float dz    = param_grid(i, 3);

    // Copy and transform
    std::vector<float> temp_pts(mov_pts);

    // Transform points once
    transform(temp_pts.data(), nmov, angle, dx, dy, dz);

    // Compute the RMS based on nearest neighbors
    float rms = pseudo_rms(tree, temp_pts.data(), nmov);

    // Store results
    angles[i] = angle;
    dxs[i] = dx;
    dys[i] = dy;
    dzs[i] = dz;
    rmss[i] = rms;  // Since we only computed one RMS value per grid

    // Optional: Debugging output every 1000 steps
    // if (i % 1000 == 0) printf("%d/%d\n", i, (int)np);
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("angle") = angles,
    Rcpp::Named("dx") = dxs,
    Rcpp::Named("dy") = dys,
    Rcpp::Named("dz") = dzs,
    Rcpp::Named("rms") = rmss
  );
}
