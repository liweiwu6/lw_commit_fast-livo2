/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#ifndef LIVO_POINT_H_
#define LIVO_POINT_H_

#include <boost/noncopyable.hpp>
#include "common_lib.h"
#include "frame.h"

class Feature;

/// A visual map point on the surface of the scene.
class VisualPoint : boost::noncopyable
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Vector3d pos_;                //!< 3d pos of the point in the world coordinate frame. 世界坐标系中点的三维位置
  Vector3d normal_;             //!< Surface normal at point. 点处的表面法线
  Matrix3d normal_information_; //!< Inverse covariance matrix of normal estimation. 正态估计的逆协方差矩阵
  Vector3d previous_normal_;    //!< Last updated normal vector. 上次更新的法向量
  list<Feature *> obs_;         //!< Reference patches which observe the point. 观察点的参考补丁
  Eigen::Matrix3d covariance_;  //!< Covariance of the point. 点的协方差
  bool is_converged_;           //!< True if the point is converged. 如果点收敛，则为true
  bool is_normal_initialized_;  //!< True if the normal is initialized. 法线初始化成功则为true
  bool has_ref_patch_;          //!< True if the point has a reference patch. 点有参考补丁则为true
  Feature *ref_patch;           //!< Reference patch of the point. 点的参考补丁

  VisualPoint(const Vector3d &pos);
  ~VisualPoint();
  void findMinScoreFeature(const Vector3d &framepos, Feature *&ftr) const;
  void deleteNonRefPatchFeatures();
  void deleteFeatureRef(Feature *ftr);
  void addFrameRef(Feature *ftr);
  bool getCloseViewObs(const Vector3d &pos, Feature *&obs, const Vector2d &cur_px) const;
};

#endif // LIVO_POINT_H_
