#include <Rcpp.h>
#include <omp.h>
#include <vector>
#include <cmath>
#include "nanoflann/nanoflann.h"

#include <vector>
#include <limits>

constexpr float INF = std::numeric_limits<float>::max();

struct Vec3
{
  float x, y, z;
};

struct BoundingBox
{
  Vec3 low;
  Vec3 high;
};

struct PointCloud
{
  std::vector<Vec3> points;
  BoundingBox bbox;

  PointCloud(const std::vector<Vec3>& pts) : points(pts)
  {
    compute_bbox();
  }

  void compute_bbox()
  {
    if (points.empty())
    {
      bbox.low  = {0, 0, 0};
      bbox.high = {0, 0, 0};
      return;
    }

    bbox.low  = {INF, INF, INF};
    bbox.high = { -INF, INF, -INF};

    for (const auto& p : points)
    {
      if (p.x < bbox.low.x)  bbox.low.x = p.x;
      if (p.x > bbox.high.x) bbox.high.x = p.x;

      if (p.y < bbox.low.y)  bbox.low.y = p.y;
      if (p.y > bbox.high.y) bbox.high.y = p.y;

      if (p.z < bbox.low.z)  bbox.low.z = p.z;
      if (p.z > bbox.high.z) bbox.high.z = p.z;
    }
  }

  inline size_t kdtree_get_point_count() const { return points.size(); }

  inline float kdtree_get_pt(const size_t idx, const size_t dim) const
  {
    const float* ptr = &points[idx].x;
    return ptr[dim];
  }

  template <class BBOX> bool kdtree_get_bbox(BBOX& bb) const
  {
    bb[0].low  = bbox.low.x;
    bb[0].high = bbox.high.x;
    bb[1].low  = bbox.low.y;
    bb[1].high = bbox.high.y;
    bb[2].low  = bbox.low.z;
    bb[2].high = bbox.high.z;
    return true; // bbox is provided
  }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3> KDtree;

float pseudo_rms(const KDtree& tree, const std::vector<Vec3>& mov_pts)
{
  std::vector<float> Dsq;
  Dsq.reserve(mov_pts.size());

  for (const auto& pt : mov_pts)
  {
    float query_pt[3] = {pt.x, pt.y, pt.z};
    size_t ret_index;
    float out_dist_sqr;
    nanoflann::KNNResultSet<float> resultSet(1);
    resultSet.init(&ret_index, &out_dist_sqr);
    tree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters(1));
    Dsq.push_back(out_dist_sqr);
  }

  std::sort(Dsq.begin(), Dsq.end());

  size_t half_n_pts = mov_pts.size() / 2;
  float rms = 0.0f;
  for (size_t i = 0; i < half_n_pts; ++i)
  {
    rms += Dsq[i];
  }

  return std::sqrt(rms / half_n_pts);
}

void transform(std::vector<Vec3>& points, float angle, float dx, float dy, float dz)
{
  float cos_a = std::cos(angle);
  float sin_a = std::sin(angle);

  for (auto& pt : points)
  {
    float x_new = cos_a * pt.x - sin_a * pt.y + dx;
    float y_new = sin_a * pt.x + cos_a * pt.y + dy;
    float z_new = pt.z + dz;
    pt.x = x_new;
    pt.y = y_new;
    pt.z = z_new;
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame rms_scan_grid(Rcpp::NumericMatrix ref, Rcpp::NumericMatrix mov, Rcpp::NumericMatrix param_grid)
{
  const size_t nref = ref.nrow();
  const size_t nmov = mov.nrow();
  const size_t np = param_grid.nrow();

  std::vector<Vec3> ref_pts(nref);
  for (size_t i = 0; i < nref; ++i)
  {
    ref_pts[i] = Vec3{
      static_cast<float>(ref(i, 0)),
      static_cast<float>(ref(i, 1)),
      static_cast<float>(ref(i, 2))
    };
  }

  std::vector<Vec3> mov_pts(nmov);
  for (size_t i = 0; i < nmov; ++i)
  {
    mov_pts[i] = Vec3{
      static_cast<float>(mov(i, 0)),
      static_cast<float>(mov(i, 1)),
      static_cast<float>(mov(i, 2))
    };
  }

  PointCloud ref_cloud{ref_pts};
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
    std::vector<Vec3> temp_pts = mov_pts;
    transform(temp_pts, angle, dx, dy, dz);

    // Compute RMS
    float rms = pseudo_rms(tree, temp_pts);

    // Store results
    angles[i] = angle;
    dxs[i] = dx;
    dys[i] = dy;
    dzs[i] = dz;
    rmss[i] = rms;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("angle") = angles,
    Rcpp::Named("dx") = dxs,
    Rcpp::Named("dy") = dys,
    Rcpp::Named("dz") = dzs,
    Rcpp::Named("rms") = rmss
  );
}
