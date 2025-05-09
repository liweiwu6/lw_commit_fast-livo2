/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "voxel_map.h"

void calcBodyCov(Eigen::Vector3d &pb, const float range_inc, const float degree_inc, Eigen::Matrix3d &cov)
{
  if (pb[2] == 0) pb[2] = 0.0001;
  float range = sqrt(pb[0] * pb[0] + pb[1] * pb[1] + pb[2] * pb[2]);//计算到原点的距离
  float range_var = range_inc * range_inc;
  Eigen::Matrix2d direction_var;
  direction_var << pow(sin(DEG2RAD(degree_inc)), 2), 0, 0, pow(sin(DEG2RAD(degree_inc)), 2);//二维方向方差矩阵
  Eigen::Vector3d direction(pb);
  direction.normalize();//归一化
  Eigen::Matrix3d direction_hat;//反对称矩阵 用于表示方向向量的叉乘操作
  direction_hat << 0, -direction(2), direction(1), direction(2), 0, -direction(0), -direction(1), direction(0), 0;
  Eigen::Vector3d base_vector1(1, 1, -(direction(0) + direction(1)) / direction(2));//基向量1
  base_vector1.normalize();
  Eigen::Vector3d base_vector2 = base_vector1.cross(direction);//基向量2
  base_vector2.normalize();
  Eigen::Matrix<double, 3, 2> N;
  N << base_vector1(0), base_vector2(0), base_vector1(1), base_vector2(1), base_vector1(2), base_vector2(2);
  Eigen::Matrix<double, 3, 2> A = range * direction_hat * N;
  cov = direction * range_var * direction.transpose() + A * direction_var * A.transpose();
}

void loadVoxelConfig(ros::NodeHandle &nh, VoxelMapConfig &voxel_config)//读取LIO参数
{
  nh.param<bool>("publish/pub_plane_en", voxel_config.is_pub_plane_map_, false);
  
  nh.param<int>("lio/max_layer", voxel_config.max_layer_, 1);
  nh.param<double>("lio/voxel_size", voxel_config.max_voxel_size_, 0.5);
  nh.param<double>("lio/min_eigen_value", voxel_config.planner_threshold_, 0.01);
  nh.param<double>("lio/sigma_num", voxel_config.sigma_num_, 3);
  nh.param<double>("lio/beam_err", voxel_config.beam_err_, 0.02);
  nh.param<double>("lio/dept_err", voxel_config.dept_err_, 0.05);
  nh.param<vector<int>>("lio/layer_init_num", voxel_config.layer_init_num_, vector<int>{5,5,5,5,5});
  nh.param<int>("lio/max_points_num", voxel_config.max_points_num_, 50);
  nh.param<int>("lio/max_iterations", voxel_config.max_iterations_, 5);

  nh.param<bool>("local_map/map_sliding_en", voxel_config.map_sliding_en, false);
  nh.param<int>("local_map/half_map_size", voxel_config.half_map_size, 100);
  nh.param<double>("local_map/sliding_thresh", voxel_config.sliding_thresh, 8);
}

void VoxelOctoTree::init_plane(const std::vector<pointWithVar> &points, VoxelPlane *plane)//todo 平面拟合
{
  plane->plane_var_ = Eigen::Matrix<double, 6, 6>::Zero();//平面协方差矩阵
  plane->covariance_ = Eigen::Matrix3d::Zero();//协方差矩阵
  plane->center_ = Eigen::Vector3d::Zero();//平面中心
  plane->normal_ = Eigen::Vector3d::Zero();//法向量
  plane->points_size_ = points.size();//点云数量
  plane->radius_ = 0;//半径
  for (auto pv : points)
  {
    plane->covariance_ += pv.point_w * pv.point_w.transpose();//累加外积，用于后续协方差计算
    plane->center_ += pv.point_w;//累加坐标，用于计算中心
  }
  plane->center_ = plane->center_ / plane->points_size_;//点云中心
  plane->covariance_ = plane->covariance_ / plane->points_size_ - plane->center_ * plane->center_.transpose();//协方差计算
  Eigen::EigenSolver<Eigen::Matrix3d> es(plane->covariance_);
  Eigen::Matrix3cd evecs = es.eigenvectors();//特征向量
  Eigen::Vector3cd evals = es.eigenvalues();//特征值
  Eigen::Vector3d evalsReal;
  evalsReal = evals.real();//取出特征值实部
  Eigen::Matrix3f::Index evalsMin, evalsMax;
  evalsReal.rowwise().sum().minCoeff(&evalsMin);//最小特征值索引
  evalsReal.rowwise().sum().maxCoeff(&evalsMax);//最大特征值索引
  int evalsMid = 3 - evalsMin - evalsMax;//中间特征值索引
  Eigen::Vector3d evecMin = evecs.real().col(evalsMin);//最小特征值对应的特征向量
  Eigen::Vector3d evecMid = evecs.real().col(evalsMid);//中间特征值对应的特征向量
  Eigen::Vector3d evecMax = evecs.real().col(evalsMax);//最大特征值对应的特征向量
  Eigen::Matrix3d J_Q;
  J_Q << 1.0 / plane->points_size_, 0, 0, 0, 1.0 / plane->points_size_, 0, 0, 0, 1.0 / plane->points_size_;//构造一个与点数量相关的对角矩阵，用于后续雅可比计算
  // && evalsReal(evalsMid) > 0.05
  //&& evalsReal(evalsMid) > 0.01
  if (evalsReal(evalsMin) < planer_threshold_)//判断最小特征值是否小于planer_threshold_，认为点云近似共面
  {
    for (int i = 0; i < points.size(); i++)//* 遍历每一个临时点
    {
      Eigen::Matrix<double, 6, 3> J;//雅可比矩阵
      Eigen::Matrix3d F;//中间矩阵
      for (int m = 0; m < 3; m++)
      {
        if (m != (int)evalsMin)//不是最小特征值
        {//* 结合了点到中心的距离、点数量、特征值差以及特征向量的组合，反映了该点在当前特征方向上的变化对整体平面参数的影响
          Eigen::Matrix<double, 1, 3> F_m =
              (points[i].point_w - plane->center_).transpose() / ((plane->points_size_) * (evalsReal[evalsMin] - evalsReal[m])) *
              (evecs.real().col(m) * evecs.real().col(evalsMin).transpose() + evecs.real().col(evalsMin) * evecs.real().col(m).transpose());
          F.row(m) = F_m;
        }
        else
        {
          Eigen::Matrix<double, 1, 3> F_m;
          F_m << 0, 0, 0;// 最小特征值对应行置零
          F.row(m) = F_m;
        }
      }
      J.block<3, 3>(0, 0) = evecs.real() * F;//雅可比矩阵前三行
      J.block<3, 3>(3, 0) = J_Q;//后三行
      plane->plane_var_ += J * points[i].var * J.transpose();//利用雅可比矩阵J、当前点的协方差points[i].var，通过三者的乘积累加到平面参数的协方差
    }

    plane->normal_ << evecs.real()(0, evalsMin), evecs.real()(1, evalsMin), evecs.real()(2, evalsMin);//协方差矩阵最小特征值对应的特征向量(通常作为平面的法向量)
    plane->y_normal_ << evecs.real()(0, evalsMid), evecs.real()(1, evalsMid), evecs.real()(2, evalsMid);//协方差矩阵中间特征值对应的特征向量
    plane->x_normal_ << evecs.real()(0, evalsMax), evecs.real()(1, evalsMax), evecs.real()(2, evalsMax);//协方差矩阵最大特征值对应的特征向量
    plane->min_eigen_value_ = evalsReal(evalsMin);//保存特征值，用于衡量点云在各主方向上的离散程度
    plane->mid_eigen_value_ = evalsReal(evalsMid);//保存特征值，用于衡量点云在各主方向上的离散程度
    plane->max_eigen_value_ = evalsReal(evalsMax);//保存特征值，用于衡量点云在各主方向上的离散程度
    plane->radius_ = sqrt(evalsReal(evalsMax));//通过对最大特征值开方得到，反映了点云在最大主方向上的空间尺度
    // 根据平面法向量和平面中心点计算，表示平面方程中的常数项，使得平面方程为 normal_.dot(x) + d_ = 0
    plane->d_ = -(plane->normal_(0) * plane->center_(0) + plane->normal_(1) * plane->center_(1) + plane->normal_(2) * plane->center_(2));
    plane->is_plane_ = true;
    plane->is_update_ = true;
    if (!plane->is_init_)//初始化就是给平面分配一个ID
    {
      plane->id_ = voxel_plane_id;//在头文件中定义，默认值为0
      voxel_plane_id++;
      plane->is_init_ = true;
    }
  }
  else//没找到平面
  {
    plane->is_update_ = true;
    plane->is_plane_ = false;
  }
}

void VoxelOctoTree::init_octo_tree()//初始化八叉树
{
  if (temp_points_.size() > points_size_threshold_)//如果点云数量大于阈值
  {
    init_plane(temp_points_, plane_ptr_);//todo 检查是否可以拟合成一个平面
    if (plane_ptr_->is_plane_ == true)
    {
      octo_state_ = 0;//无需细分
      // new added
      if (temp_points_.size() > max_points_num_)//如果点云数量大于最大值50
      {
        update_enable_ = false;//不再更新
        std::vector<pointWithVar>().swap(temp_points_);//!清空临时点云集合，体素只保留平面参数，不保留所有点，节省内存，提高效率
        new_points_ = 0;
      }
    }
    else
    {
      octo_state_ = 1;//需要细分
      cut_octo_tree();//todo进一步细分八叉树（不断细分，直到找到平面或者达到最大层数）
    }
    init_octo_ = true;//已经初始化
    new_points_ = 0;//* 初始化完成会把新点设置为0
  }
}

void VoxelOctoTree::cut_octo_tree()
{
  if (layer_ >= max_layer_)//超过最大层数，停止
  {
    octo_state_ = 0;//无需细分
    return;
  }
  for (size_t i = 0; i < temp_points_.size(); i++)//* 将临时点分配到八叉树的子节点
  {
    int xyz[3] = {0, 0, 0};
    if (temp_points_[i].point_w[0] > voxel_center_[0]) { xyz[0] = 1; }
    if (temp_points_[i].point_w[1] > voxel_center_[1]) { xyz[1] = 1; }
    if (temp_points_[i].point_w[2] > voxel_center_[2]) { xyz[2] = 1; }
    int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];//确定八叉树的子节点
        // 8个外框示意图
    // clang-format off
    // 第一层：左上1 右上2 左下3 右下4
    // 第二层：左上5 右上6 左下7 右下8
    //     ---> x    /-------/-------/|
    //    /|        /-------/-------/||
    //   / |       /-------/-------/ ||
    //  y  |z      |       |       | /|
    //             |_______|_______|/|/
    //             |       |       | /
    //             |_______|_______|/
    if (leaves_[leafnum] == nullptr)//初始化子结点
    {
      leaves_[leafnum] = new VoxelOctoTree(max_layer_, layer_ + 1, layer_init_num_[layer_ + 1], max_points_num_, planer_threshold_);
      leaves_[leafnum]->layer_init_num_ = layer_init_num_;
      leaves_[leafnum]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
      leaves_[leafnum]->quater_length_ = quater_length_ / 2;
    }
    leaves_[leafnum]->temp_points_.push_back(temp_points_[i]);
    leaves_[leafnum]->new_points_++;
  }
  for (uint i = 0; i < 8; i++)//遍历8个子节点
  {
    if (leaves_[i] != nullptr)
    {
      if (leaves_[i]->temp_points_.size() > leaves_[i]->points_size_threshold_)
      {
        init_plane(leaves_[i]->temp_points_, leaves_[i]->plane_ptr_);//拟合平面
        if (leaves_[i]->plane_ptr_->is_plane_)
        {
          leaves_[i]->octo_state_ = 0;//找到平面，无需细分
          // new added
          if (leaves_[i]->temp_points_.size() > leaves_[i]->max_points_num_)
          {
            leaves_[i]->update_enable_ = false;//有足够点云，停止更新
            std::vector<pointWithVar>().swap(leaves_[i]->temp_points_);//清空临时点云集合
            new_points_ = 0;
          }
        }
        else
        {
          leaves_[i]->octo_state_ = 1;
          leaves_[i]->cut_octo_tree();
        }
        leaves_[i]->init_octo_ = true;
        leaves_[i]->new_points_ = 0;
      }
    }
  }
}

void VoxelOctoTree::UpdateOctoTree(const pointWithVar &pv)//todo 更新八叉树
{ 
  if (!init_octo_)//八叉树未初始化
  {
    new_points_++;
    temp_points_.push_back(pv);
    if (temp_points_.size() > points_size_threshold_) { init_octo_tree(); }
  }
  else//八叉树已经初始化
  {
    if (plane_ptr_->is_plane_)//存在平面
    {
      if (update_enable_)//可以更新
      {
        new_points_++;
        temp_points_.push_back(pv);
        if (new_points_ > update_size_threshold_)
        {
          init_plane(temp_points_, plane_ptr_);
          new_points_ = 0;
        }
        if (temp_points_.size() >= max_points_num_)
        {
          update_enable_ = false;
          std::vector<pointWithVar>().swap(temp_points_);
          new_points_ = 0;
        }
      }
    }
    else//不存在平面
    {
      if (layer_ < max_layer_)//小于最大层数
      {
        int xyz[3] = {0, 0, 0};
        if (pv.point_w[0] > voxel_center_[0]) { xyz[0] = 1; }
        if (pv.point_w[1] > voxel_center_[1]) { xyz[1] = 1; }
        if (pv.point_w[2] > voxel_center_[2]) { xyz[2] = 1; }
        int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];
        if (leaves_[leafnum] != nullptr) { leaves_[leafnum]->UpdateOctoTree(pv); }
        else
        {
          leaves_[leafnum] = new VoxelOctoTree(max_layer_, layer_ + 1, layer_init_num_[layer_ + 1], max_points_num_, planer_threshold_);
          leaves_[leafnum]->layer_init_num_ = layer_init_num_;
          leaves_[leafnum]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_;
          leaves_[leafnum]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
          leaves_[leafnum]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
          leaves_[leafnum]->quater_length_ = quater_length_ / 2;
          leaves_[leafnum]->UpdateOctoTree(pv);
        }
      }
      else//大于最大层数
      {
        if (update_enable_)//可以更新
        {
          new_points_++;
          temp_points_.push_back(pv);
          if (new_points_ > update_size_threshold_)
          {
            init_plane(temp_points_, plane_ptr_);
            new_points_ = 0;
          }
          if (temp_points_.size() > max_points_num_)
          {
            update_enable_ = false;
            std::vector<pointWithVar>().swap(temp_points_);
            new_points_ = 0;
          }
        }
      }
    }
  }
}

VoxelOctoTree *VoxelOctoTree::find_correspond(Eigen::Vector3d pw)
{
  if (!init_octo_ || plane_ptr_->is_plane_ || (layer_ >= max_layer_)) return this;

  int xyz[3] = {0, 0, 0};
  xyz[0] = pw[0] > voxel_center_[0] ? 1 : 0;
  xyz[1] = pw[1] > voxel_center_[1] ? 1 : 0;
  xyz[2] = pw[2] > voxel_center_[2] ? 1 : 0;
  int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];

  // printf("leafnum: %d. \n", leafnum);

  return (leaves_[leafnum] != nullptr) ? leaves_[leafnum]->find_correspond(pw) : this;
}

VoxelOctoTree *VoxelOctoTree::Insert(const pointWithVar &pv)
{
  if ((!init_octo_) || (init_octo_ && plane_ptr_->is_plane_) || (init_octo_ && (!plane_ptr_->is_plane_) && (layer_ >= max_layer_)))
  {
    new_points_++;
    temp_points_.push_back(pv);
    return this;
  }

  if (init_octo_ && (!plane_ptr_->is_plane_) && (layer_ < max_layer_))
  {
    int xyz[3] = {0, 0, 0};
    xyz[0] = pv.point_w[0] > voxel_center_[0] ? 1 : 0;
    xyz[1] = pv.point_w[1] > voxel_center_[1] ? 1 : 0;
    xyz[2] = pv.point_w[2] > voxel_center_[2] ? 1 : 0;
    int leafnum = 4 * xyz[0] + 2 * xyz[1] + xyz[2];
    if (leaves_[leafnum] != nullptr) { return leaves_[leafnum]->Insert(pv); }
    else
    {
      leaves_[leafnum] = new VoxelOctoTree(max_layer_, layer_ + 1, layer_init_num_[layer_ + 1], max_points_num_, planer_threshold_);
      leaves_[leafnum]->layer_init_num_ = layer_init_num_;
      leaves_[leafnum]->voxel_center_[0] = voxel_center_[0] + (2 * xyz[0] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[1] = voxel_center_[1] + (2 * xyz[1] - 1) * quater_length_;
      leaves_[leafnum]->voxel_center_[2] = voxel_center_[2] + (2 * xyz[2] - 1) * quater_length_;
      leaves_[leafnum]->quater_length_ = quater_length_ / 2;
      return leaves_[leafnum]->Insert(pv);
    }
  }
  return nullptr;
}

void VoxelMapManager::StateEstimation(StatesGroup &state_propagat)
{
  cross_mat_list_.clear();
  cross_mat_list_.reserve(feats_down_size_);//特征点的叉乘矩阵
  body_cov_list_.clear();
  body_cov_list_.reserve(feats_down_size_);//特征点的协方差矩阵

  // build_residual_time = 0.0;
  // ekf_time = 0.0;
  // double t0 = omp_get_wtime();

  for (size_t i = 0; i < feats_down_body_->size(); i++)//* 计算所有特征点的协方差矩阵和叉乘矩阵
  {
    V3D point_this(feats_down_body_->points[i].x, feats_down_body_->points[i].y, feats_down_body_->points[i].z);//* lidar坐标系 
    if (point_this[2] == 0) { point_this[2] = 0.001; }//防止除零
    M3D var;
    calcBodyCov(point_this, config_setting_.dept_err_, config_setting_.beam_err_, var);//计算特征点的协方差矩阵 这里是lidar坐标系
    body_cov_list_.push_back(var);
    point_this = extR_ * point_this + extT_;// * lidar->imu extR_为lidar与imu外参
    M3D point_crossmat;
    point_crossmat << SKEW_SYM_MATRX(point_this);
    cross_mat_list_.push_back(point_crossmat); //叉乘矩阵
  }

  vector<pointWithVar>().swap(pv_list_);//清空点云集合，释放内存
  pv_list_.resize(feats_down_size_);//重新分配内存

  int rematch_num = 0;
  MD(DIM_STATE, DIM_STATE) G, H_T_H, I_STATE;//eigen矩阵
  G.setZero();//清零
  H_T_H.setZero();//清零
  I_STATE.setIdentity();//初始化为单位矩阵

  bool flg_EKF_inited, flg_EKF_converged, EKF_stop_flg = 0;
  for (int iterCount = 0; iterCount < config_setting_.max_iterations_; iterCount++)//* 迭代求解 max_iterations_最大迭代次数为5
  {
    double total_residual = 0.0;
    pcl::PointCloud<pcl::PointXYZI>::Ptr world_lidar(new pcl::PointCloud<pcl::PointXYZI>);
    TransformLidar(state_.rot_end, state_.pos_end, feats_down_body_, world_lidar);//转换雷达坐标系 lidar->world 
    M3D rot_var = state_.cov.block<3, 3>(0, 0);
    M3D t_var = state_.cov.block<3, 3>(3, 3);
    for (size_t i = 0; i < feats_down_body_->size(); i++)//* 计算每个点的协方差矩阵，写入到pv_list_中
    {
      pointWithVar &pv = pv_list_[i];//这里应该都是空的 进行绑定 填充pv_list_
      pv.point_b << feats_down_body_->points[i].x, feats_down_body_->points[i].y, feats_down_body_->points[i].z;//lidar坐标系下的点
      pv.point_w << world_lidar->points[i].x, world_lidar->points[i].y, world_lidar->points[i].z;//世界坐标系下的点

      M3D cov = body_cov_list_[i];
      M3D point_crossmat = cross_mat_list_[i];
      // * //////////从激光雷达坐标系转换到世界坐标系///////////////////传播误差//////////////////////////////////////////////////平移不确定性
      cov = state_.rot_end * cov * state_.rot_end.transpose() + (-point_crossmat) * rot_var * (-point_crossmat.transpose()) + t_var;
      pv.var = cov;//世界坐标系下的点的协方差矩阵
      pv.body_var = body_cov_list_[i];//激光雷达坐标系下的点的协方差矩阵
    }
    ptpl_list_.clear();//清空残差列表

    // double t1 = omp_get_wtime();

    BuildResidualListOMP(pv_list_, ptpl_list_);//todo 构建残差列表 这里对pv_list_

    // build_residual_time += omp_get_wtime() - t1;

    for (int i = 0; i < ptpl_list_.size(); i++)//计算残差总和
    {
      total_residual += fabs(ptpl_list_[i].dis_to_plane_);
    }
    effct_feat_num_ = ptpl_list_.size();//有效特征点数量
    cout << "[ LIO ] Raw feature num: " << feats_undistort_->size() << ", downsampled feature num:" << feats_down_size_ 
         << " effective feature num: " << effct_feat_num_ << " average residual: " << total_residual / effct_feat_num_ << endl;

    /*** Computation of Measuremnt Jacobian matrix H and measurents covarience  计算测量雅可比矩阵 H 和测量协方差
     * ***/ //todo
    MatrixXd Hsub(effct_feat_num_, 6);//测量雅可比矩阵的一部分
    MatrixXd Hsub_T_R_inv(6, effct_feat_num_);//雅可比矩阵的转置乘以测量协方差的逆
    VectorXd R_inv(effct_feat_num_);//测量协方差的逆对角元素
    VectorXd meas_vec(effct_feat_num_);//测量残差向量
    meas_vec.setZero();//初始化为0
    for (int i = 0; i < effct_feat_num_; i++)//*计算每个点的测量协方差和雅可比矩阵H
    {
      auto &ptpl = ptpl_list_[i];
      V3D point_this(ptpl.point_b_);//lidar坐标系
      point_this = extR_ * point_this + extT_;//IMU坐标系
      V3D point_body(ptpl.point_b_);//lidar坐标系
      M3D point_crossmat;//叉乘矩阵
      point_crossmat << SKEW_SYM_MATRX(point_this);//调用宏，生成叉乘矩阵（反对称矩阵）

      /***  get the normal vector of closest surface/corner   获取最近表面/角点的法向量       ***/

      V3D point_world = state_propagat.rot_end * point_this + state_propagat.pos_end;//世界坐标系
      Eigen::Matrix<double, 1, 6> J_nq;//雅可比矩阵
      J_nq.block<1, 3>(0, 0) = point_world - ptpl_list_[i].center_;//点到平面中心的向量
      J_nq.block<1, 3>(0, 3) = -ptpl_list_[i].normal_;//平面法向量
      M3D var;//协方差矩阵
      // V3D normal_b = state_.rot_end.inverse() * ptpl_list_[i].normal_;
      // V3D point_b = ptpl_list_[i].point_b_;
      // double cos_theta = fabs(normal_b.dot(point_b) / point_b.norm());
      // ptpl_list_[i].body_cov_ = ptpl_list_[i].body_cov_ * (1.0 / cos_theta) * (1.0 / cos_theta);

      // point_w cov
      // var = state_propagat.rot_end * extR_ * ptpl_list_[i].body_cov_ * (state_propagat.rot_end * extR_).transpose() +
      //       state_propagat.cov.block<3, 3>(3, 3) + (-point_crossmat) * state_propagat.cov.block<3, 3>(0, 0) * (-point_crossmat).transpose();

      // point_w cov (another_version)
      // var = state_propagat.rot_end * extR_ * ptpl_list_[i].body_cov_ * (state_propagat.rot_end * extR_).transpose() +
      //       state_propagat.cov.block<3, 3>(3, 3) - point_crossmat * state_propagat.cov.block<3, 3>(0, 0) * point_crossmat;

      // point_body cov 
      var = state_propagat.rot_end * extR_ * ptpl_list_[i].body_cov_ * (state_propagat.rot_end * extR_).transpose();//将激光雷达坐标系下的协方差矩阵转换到世界坐标系下

      double sigma_l = J_nq * ptpl_list_[i].plane_var_ * J_nq.transpose();//测量噪声

      R_inv(i) = 1.0 / (0.001 + sigma_l + ptpl_list_[i].normal_.transpose() * var * ptpl_list_[i].normal_);//协方差的逆
      // R_inv(i) = 1.0 / (sigma_l + ptpl_list_[i].normal_.transpose() * var * ptpl_list_[i].normal_);

      /*** calculate the Measuremnt Jacobian matrix H  计算测量雅可比矩阵 H   ***/
      V3D A(point_crossmat * state_.rot_end.transpose() * ptpl_list_[i].normal_);
      Hsub.row(i) << VEC_FROM_ARRAY(A), ptpl_list_[i].normal_[0], ptpl_list_[i].normal_[1], ptpl_list_[i].normal_[2];
      Hsub_T_R_inv.col(i) << A[0] * R_inv(i), A[1] * R_inv(i), A[2] * R_inv(i), ptpl_list_[i].normal_[0] * R_inv(i),
          ptpl_list_[i].normal_[1] * R_inv(i), ptpl_list_[i].normal_[2] * R_inv(i);
      meas_vec(i) = -ptpl_list_[i].dis_to_plane_;
    }
    EKF_stop_flg = false;//停止标志
    flg_EKF_converged = false;//收敛标志
    /*** Iterative Kalman Filter Update    迭代卡尔曼滤波更新   ***/  // todo
    MatrixXd K(DIM_STATE, effct_feat_num_);// ? 好像没有用到这个变量
    // auto &&Hsub_T = Hsub.transpose();
    auto &&HTz = Hsub_T_R_inv * meas_vec;
    // fout_dbg<<"HTz: "<<HTz<<endl;
    H_T_H.block<6, 6>(0, 0) = Hsub_T_R_inv * Hsub;
    // EigenSolver<Matrix<double, 6, 6>> es(H_T_H.block<6,6>(0,0));
    MD(DIM_STATE, DIM_STATE) &&K_1 = (H_T_H.block<DIM_STATE, DIM_STATE>(0, 0) + state_.cov.block<DIM_STATE, DIM_STATE>(0, 0).inverse()).inverse();
    G.block<DIM_STATE, 6>(0, 0) = K_1.block<DIM_STATE, 6>(0, 0) * H_T_H.block<6, 6>(0, 0);
    auto vec = state_propagat - state_;
    VD(DIM_STATE)
    solution = K_1.block<DIM_STATE, 6>(0, 0) * HTz + vec.block<DIM_STATE, 1>(0, 0) - G.block<DIM_STATE, 6>(0, 0) * vec.block<6, 1>(0, 0);
    int minRow, minCol;
    state_ += solution;
    auto rot_add = solution.block<3, 1>(0, 0);
    auto t_add = solution.block<3, 1>(3, 0);
    if ((rot_add.norm() * 57.3 < 0.01) && (t_add.norm() * 100 < 0.015)) { flg_EKF_converged = true; }//收敛判断
    V3D euler_cur = state_.rot_end.eulerAngles(2, 1, 0);

    /*** Rematch Judgement  重新匹配判断 ***/

    if (flg_EKF_converged || ((rematch_num == 0) && (iterCount == (config_setting_.max_iterations_ - 2)))) { rematch_num++; }

    /*** Convergence Judgements and Covariance Update   收敛判断和协方差更新  ***/
    if (!EKF_stop_flg && (rematch_num >= 2 || (iterCount == config_setting_.max_iterations_ - 1)))
    {
      /*** Covariance Update 协方差更新   ***/
      // _state.cov = (I_STATE - G) * _state.cov;
      state_.cov.block<DIM_STATE, DIM_STATE>(0, 0) =
          (I_STATE.block<DIM_STATE, DIM_STATE>(0, 0) - G.block<DIM_STATE, DIM_STATE>(0, 0)) * state_.cov.block<DIM_STATE, DIM_STATE>(0, 0);
      // total_distance += (_state.pos_end - position_last).norm();
      position_last_ = state_.pos_end;//用于判断地图是否需要滑动调整
      geoQuat_ = tf::createQuaternionMsgFromRollPitchYaw(euler_cur(0), euler_cur(1), euler_cur(2));

      // VD(DIM_STATE) K_sum  = K.rowwise().sum();
      // VD(DIM_STATE) P_diag = _state.cov.diagonal();
      EKF_stop_flg = true;
    }
    if (EKF_stop_flg) break;
  }

  // double t2 = omp_get_wtime();
  // scan_count++;
  // ekf_time = t2 - t0 - build_residual_time;

  // ave_build_residual_time = ave_build_residual_time * (scan_count - 1) / scan_count + build_residual_time / scan_count;
  // ave_ekf_time = ave_ekf_time * (scan_count - 1) / scan_count + ekf_time / scan_count;

  // cout << "[ Mapping ] ekf_time: " << ekf_time << "s, build_residual_time: " << build_residual_time << "s" << endl;
  // cout << "[ Mapping ] ave_ekf_time: " << ave_ekf_time << "s, ave_build_residual_time: " << ave_build_residual_time << "s" << endl;
}

void VoxelMapManager::TransformLidar(const Eigen::Matrix3d rot, const Eigen::Vector3d t, const PointCloudXYZI::Ptr &input_cloud,
                                     pcl::PointCloud<pcl::PointXYZI>::Ptr &trans_cloud)
{
  pcl::PointCloud<pcl::PointXYZI>().swap(*trans_cloud);//清空
  trans_cloud->reserve(input_cloud->size());//重新分配内存
  for (size_t i = 0; i < input_cloud->size(); i++)
  {
    pcl::PointXYZINormal p_c = input_cloud->points[i];
    Eigen::Vector3d p(p_c.x, p_c.y, p_c.z);
    p = (rot * (extR_ * p + extT_) + t);//转换到世界坐标系
    pcl::PointXYZI pi;
    pi.x = p(0);
    pi.y = p(1);
    pi.z = p(2);
    pi.intensity = p_c.intensity;
    trans_cloud->points.push_back(pi);
  }
}

void VoxelMapManager::BuildVoxelMap()
{
  //读取地图参数
  float voxel_size = config_setting_.max_voxel_size_;//体素大小
  float planer_threshold = config_setting_.planner_threshold_;//!平面阈值,没找到怎么赋值的
  int max_layer = config_setting_.max_layer_;//最大层数
  int max_points_num = config_setting_.max_points_num_;//最大点数
  std::vector<int> layer_init_num = config_setting_.layer_init_num_;//每层初始化点数

  std::vector<pointWithVar> input_points;//存储点云数据

  for (size_t i = 0; i < feats_down_world_->size(); i++)//计算每个点的协方差
  {
    pointWithVar pv;
    pv.point_w << feats_down_world_->points[i].x, feats_down_world_->points[i].y, feats_down_world_->points[i].z;//世界坐标系点
    V3D point_this(feats_down_body_->points[i].x, feats_down_body_->points[i].y, feats_down_body_->points[i].z);//雷达坐标系点
    M3D var;//协方差矩阵
    calcBodyCov(point_this, config_setting_.dept_err_, config_setting_.beam_err_, var);//计算this的协方差矩阵
    M3D point_crossmat;//叉乘矩阵
    point_crossmat << SKEW_SYM_MATRX(point_this);//反对称矩阵
    var = (state_.rot_end * extR_) * var * (state_.rot_end * extR_).transpose() +
          (-point_crossmat) * state_.cov.block<3, 3>(0, 0) * (-point_crossmat).transpose() + state_.cov.block<3, 3>(3, 3);//计算点的协方差
    pv.var = var;//更新协方差
    input_points.push_back(pv);
  }

  uint plsize = input_points.size();
  for (uint i = 0; i < plsize; i++)//将点加入到地图中
  {
    const pointWithVar p_v = input_points[i];
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)//计算点在体素地图中的位置
    {
      loc_xyz[j] = p_v.point_w[j] / voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
    auto iter = voxel_map_.find(position);
    if (iter != voxel_map_.end())//已经存在对应的体素，将点加入体素中
    {
      voxel_map_[position]->temp_points_.push_back(p_v);//临时集合
      voxel_map_[position]->new_points_++;//更新新点计数
    }
    else
    {
      VoxelOctoTree *octo_tree = new VoxelOctoTree(max_layer, 0, layer_init_num[0], max_points_num, planer_threshold);
      voxel_map_[position] = octo_tree;//八叉树
      voxel_map_[position]->quater_length_ = voxel_size / 4;//四分之一体素长度（计算子体素中心点在世界系的坐标，所以需要子体素的半边长，即1/4父体素边长）
      voxel_map_[position]->voxel_center_[0] = (0.5 + position.x) * voxel_size;//体素中心
      voxel_map_[position]->voxel_center_[1] = (0.5 + position.y) * voxel_size;//体素中心
      voxel_map_[position]->voxel_center_[2] = (0.5 + position.z) * voxel_size;//体素中心
      voxel_map_[position]->temp_points_.push_back(p_v);//临时集合
      voxel_map_[position]->new_points_++;//更新新点计数
      voxel_map_[position]->layer_init_num_ = layer_init_num;//每层初始化点数
    }
  }
  for (auto iter = voxel_map_.begin(); iter != voxel_map_.end(); ++iter)
  {
    iter->second->init_octo_tree();//todo初始化八叉树
  }
}

V3F VoxelMapManager::RGBFromVoxel(const V3D &input_point)
{
  int64_t loc_xyz[3];
  for (int j = 0; j < 3; j++)
  {
    loc_xyz[j] = floor(input_point[j] / config_setting_.max_voxel_size_);
  }

  VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
  int64_t ind = loc_xyz[0] + loc_xyz[1] + loc_xyz[2];
  uint k((ind + 100000) % 3);
  V3F RGB((k == 0) * 255.0, (k == 1) * 255.0, (k == 2) * 255.0);
  // cout<<"RGB: "<<RGB.transpose()<<endl;
  return RGB;
}

void VoxelMapManager::UpdateVoxelMap(const std::vector<pointWithVar> &input_points)//更新地图
{
  float voxel_size = config_setting_.max_voxel_size_;//体素大小
  float planer_threshold = config_setting_.planner_threshold_;//平面阈值
  int max_layer = config_setting_.max_layer_;//最大层数
  int max_points_num = config_setting_.max_points_num_;//最大点数
  std::vector<int> layer_init_num = config_setting_.layer_init_num_;//每层初始化点数
  uint plsize = input_points.size();
  for (uint i = 0; i < plsize; i++)//* 遍历所有点
  {
    const pointWithVar p_v = input_points[i];
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)//计算点在体素地图中的位置
    {
      loc_xyz[j] = p_v.point_w[j] / voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);
    auto iter = voxel_map_.find(position);
    if (iter != voxel_map_.end()) { voxel_map_[position]->UpdateOctoTree(p_v); }//体素存在，更新
    else//新建体素
    {
      VoxelOctoTree *octo_tree = new VoxelOctoTree(max_layer, 0, layer_init_num[0], max_points_num, planer_threshold);
      voxel_map_[position] = octo_tree;
      voxel_map_[position]->layer_init_num_ = layer_init_num;
      voxel_map_[position]->quater_length_ = voxel_size / 4;
      voxel_map_[position]->voxel_center_[0] = (0.5 + position.x) * voxel_size;
      voxel_map_[position]->voxel_center_[1] = (0.5 + position.y) * voxel_size;
      voxel_map_[position]->voxel_center_[2] = (0.5 + position.z) * voxel_size;
      voxel_map_[position]->UpdateOctoTree(p_v);
    }
  }
}

void VoxelMapManager::BuildResidualListOMP(std::vector<pointWithVar> &pv_list, std::vector<PointToPlane> &ptpl_list)//输入，输出
{
  int max_layer = config_setting_.max_layer_;//最大层数
  double voxel_size = config_setting_.max_voxel_size_;//体素大小
  double sigma_num = config_setting_.sigma_num_;//用于残差判断的标准差倍数
  std::mutex mylock;
  ptpl_list.clear();//清空
  std::vector<PointToPlane> all_ptpl_list(pv_list.size());
  std::vector<bool> useful_ptpl(pv_list.size());
  std::vector<size_t> index(pv_list.size());
  for (size_t i = 0; i < index.size(); ++i)//初始化一些参数
  {
    index[i] = i;
    useful_ptpl[i] = false;
  }
  #ifdef MP_EN
    omp_set_num_threads(MP_PROC_NUM);//多线程计算
    #pragma omp parallel for
  #endif
  for (int i = 0; i < index.size(); i++)//迭代计算扫描中的所有LiDAR点与地图中各平面之间的残差
  {
    pointWithVar &pv = pv_list[i];//引用
    float loc_xyz[3];
    for (int j = 0; j < 3; j++)//计算点在体素地图中的位置
    {
      loc_xyz[j] = pv.point_w[j] / voxel_size;
      if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
    }
    VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);//体素位置
    auto iter = voxel_map_.find(position);
    if (iter != voxel_map_.end())//找到对应的体素
    {
      VoxelOctoTree *current_octo = iter->second;
      PointToPlane single_ptpl;
      bool is_sucess = false;
      double prob = 0;
      build_single_residual(pv, current_octo, 0, is_sucess, prob, single_ptpl);//todo 递归函数，用于计算单个点与体素平面的残差 同时把平面法向量赋值给点
      if (!is_sucess)//如果没有找到对应的体素,查找附近的体素
      {
        VOXEL_LOCATION near_position = position;
        if (loc_xyz[0] > (current_octo->voxel_center_[0] + current_octo->quater_length_)) { near_position.x = near_position.x + 1; }
        else if (loc_xyz[0] < (current_octo->voxel_center_[0] - current_octo->quater_length_)) { near_position.x = near_position.x - 1; }
        if (loc_xyz[1] > (current_octo->voxel_center_[1] + current_octo->quater_length_)) { near_position.y = near_position.y + 1; }
        else if (loc_xyz[1] < (current_octo->voxel_center_[1] - current_octo->quater_length_)) { near_position.y = near_position.y - 1; }
        if (loc_xyz[2] > (current_octo->voxel_center_[2] + current_octo->quater_length_)) { near_position.z = near_position.z + 1; }
        else if (loc_xyz[2] < (current_octo->voxel_center_[2] - current_octo->quater_length_)) { near_position.z = near_position.z - 1; }
        auto iter_near = voxel_map_.find(near_position);
        if (iter_near != voxel_map_.end()) { build_single_residual(pv, iter_near->second, 0, is_sucess, prob, single_ptpl); }
      }
      if (is_sucess)
      {
        mylock.lock();
        useful_ptpl[i] = true;
        all_ptpl_list[i] = single_ptpl;
        mylock.unlock();
      }
      else
      {
        mylock.lock();
        useful_ptpl[i] = false;
        mylock.unlock();
      }
    }
  }
  for (size_t i = 0; i < useful_ptpl.size(); i++)//取出有效的残差
  {
    if (useful_ptpl[i]) { ptpl_list.push_back(all_ptpl_list[i]); }
  }
}

void VoxelMapManager::build_single_residual(pointWithVar &pv, const VoxelOctoTree *current_octo, const int current_layer, bool &is_sucess,
                                            double &prob, PointToPlane &single_ptpl) 
{
  int max_layer = config_setting_.max_layer_;
  double sigma_num = config_setting_.sigma_num_;

  double radius_k = 3;
  Eigen::Vector3d p_w = pv.point_w;//点的世界坐标
  if (current_octo->plane_ptr_->is_plane_)//检查是否含有平面信息
  {
    VoxelPlane &plane = *current_octo->plane_ptr_;//取出平面信息
    Eigen::Vector3d p_world_to_center = p_w - plane.center_;//计算点到平面中心的向量
    float dis_to_plane = fabs(plane.normal_(0) * p_w(0) + plane.normal_(1) * p_w(1) + plane.normal_(2) * p_w(2) + plane.d_);//点到平面的距离
    float dis_to_center = (plane.center_(0) - p_w(0)) * (plane.center_(0) - p_w(0)) + (plane.center_(1) - p_w(1)) * (plane.center_(1) - p_w(1)) +
                          (plane.center_(2) - p_w(2)) * (plane.center_(2) - p_w(2));//点到平面中心的距离
    float range_dis = sqrt(dis_to_center - dis_to_plane * dis_to_plane);//计算点到平面中心的水平距离（投影距离）

    if (range_dis <= radius_k * plane.radius_)//判断点是否在平面的有效范围内
    {
      Eigen::Matrix<double, 1, 6> J_nq;//雅可比矩阵
      J_nq.block<1, 3>(0, 0) = p_w - plane.center_;//点到平面中心的向量
      J_nq.block<1, 3>(0, 3) = -plane.normal_;//平面法向量的负值
      double sigma_l = J_nq * plane.plane_var_ * J_nq.transpose();//平面协方差  
      sigma_l += plane.normal_.transpose() * pv.var * plane.normal_;//点的协方差
      if (dis_to_plane < sigma_num * sqrt(sigma_l))//匹配阈值
      {
        is_sucess = true;
        double this_prob = 1.0 / (sqrt(sigma_l)) * exp(-0.5 * dis_to_plane * dis_to_plane / sigma_l);//更新匹配概率相关数据
        if (this_prob > prob)//更新最高概率
        {
          prob = this_prob;
          pv.normal = plane.normal_;//!将平面法向量赋值给点
          single_ptpl.body_cov_ = pv.body_var;
          single_ptpl.point_b_ = pv.point_b; 
          single_ptpl.point_w_ = pv.point_w;
          single_ptpl.plane_var_ = plane.plane_var_;
          single_ptpl.normal_ = plane.normal_;
          single_ptpl.center_ = plane.center_;
          single_ptpl.d_ = plane.d_;
          single_ptpl.layer_ = current_layer;
          single_ptpl.dis_to_plane_ = plane.normal_(0) * p_w(0) + plane.normal_(1) * p_w(1) + plane.normal_(2) * p_w(2) + plane.d_;//*计算点 p_w 到平面的有符号距离
        }
        return;
      }
      else
      {
        // is_sucess = false;
        return;
      }
    }
    else
    {
      // is_sucess = false;
      return;
    }
  }
  else//无平面信息，递归查找子结点
  {
    if (current_layer < max_layer)
    {
      for (size_t leafnum = 0; leafnum < 8; leafnum++)
      {
        if (current_octo->leaves_[leafnum] != nullptr)
        {

          VoxelOctoTree *leaf_octo = current_octo->leaves_[leafnum];
          build_single_residual(pv, leaf_octo, current_layer + 1, is_sucess, prob, single_ptpl);
        }
      }
      return;
    }
    else { return; }//达到最大层级且仍未找到匹配，则返回
  }
}

void VoxelMapManager::pubVoxelMap()//发布体素地图
{
  double max_trace = 0.25;// 最大协方差迹
  double pow_num = 0.2;// 颜色映射参数
  ros::Rate loop(500);//控制发布频率
  float use_alpha = 0.8;//透明度
  visualization_msgs::MarkerArray voxel_plane;//存放要发布的平面
  voxel_plane.markers.reserve(1000000);//预留空间
  std::vector<VoxelPlane> pub_plane_list;//存放要发布的平面列表
  for (auto iter = voxel_map_.begin(); iter != voxel_map_.end(); iter++)//遍历地图
  {
    GetUpdatePlane(iter->second, config_setting_.max_layer_, pub_plane_list);//获取更新的平面列表
  }
  for (size_t i = 0; i < pub_plane_list.size(); i++)//遍历平面列表
  {
    V3D plane_cov = pub_plane_list[i].plane_var_.block<3, 3>(0, 0).diagonal();//获取平面的协方差矩阵
    double trace = plane_cov.sum();//计算协方差矩阵的迹
    if (trace >= max_trace) { trace = max_trace; }
    trace = trace * (1.0 / max_trace);
    trace = pow(trace, pow_num);
    uint8_t r, g, b;
    mapJet(trace, 0, 1, r, g, b);//映射为RGB颜色。trace越大，颜色越偏向红色，trace越小，颜色越偏向蓝色，这样可以直观反映平面拟合的不确定性。（迹越小，代表平面越可靠）
    Eigen::Vector3d plane_rgb(r / 256.0, g / 256.0, b / 256.0);
    double alpha;
    if (pub_plane_list[i].is_plane_) { alpha = use_alpha; }//透明度 alpha 则根据平面是否有效（is_plane_）决定 这里只显示平面的体素
    else { alpha = 0; }//0为透明
    pubSinglePlane(voxel_plane, "plane", pub_plane_list[i], alpha, plane_rgb);//封装为一个ROS Marker
  }
  voxel_map_pub_.publish(voxel_plane);
  loop.sleep();
}

void VoxelMapManager::GetUpdatePlane(const VoxelOctoTree *current_octo, const int pub_max_voxel_layer, std::vector<VoxelPlane> &plane_list)// 获取更新的平面列表
{
  if (current_octo->layer_ > pub_max_voxel_layer) { return; }//检查当前八叉树节点的层级是否超过了最大发布层级
  if (current_octo->plane_ptr_->is_update_) { plane_list.push_back(*current_octo->plane_ptr_); }//检查是否以更新，如果更新了就加入到发布列表
  if (current_octo->layer_ < current_octo->max_layer_)
  {
    if (!current_octo->plane_ptr_->is_plane_)
    {
      for (size_t i = 0; i < 8; i++)
      {
        if (current_octo->leaves_[i] != nullptr) { GetUpdatePlane(current_octo->leaves_[i], pub_max_voxel_layer, plane_list); }//递归查找，获取更新的平面
      }
    }
  }
  return;
}

void VoxelMapManager::pubSinglePlane(visualization_msgs::MarkerArray &plane_pub, const std::string plane_ns, const VoxelPlane &single_plane,
                                     const float alpha, const Eigen::Vector3d rgb)//封装为一个ROS Marker
{
  visualization_msgs::Marker plane;
  plane.header.frame_id = "camera_init";
  plane.header.stamp = ros::Time();
  plane.ns = plane_ns;
  plane.id = single_plane.id_;
  plane.type = visualization_msgs::Marker::CYLINDER;
  plane.action = visualization_msgs::Marker::ADD;
  plane.pose.position.x = single_plane.center_[0];//平面中心坐标
  plane.pose.position.y = single_plane.center_[1];//平面中心坐标
  plane.pose.position.z = single_plane.center_[2];//平面中心坐标
  geometry_msgs::Quaternion q;
  CalcVectQuation(single_plane.x_normal_, single_plane.y_normal_, single_plane.normal_, q);//计算得到四元数
  plane.pose.orientation = q;
  plane.scale.x = 3 * sqrt(single_plane.max_eigen_value_);//特征值 表示X，Y方向的拓展
  plane.scale.y = 3 * sqrt(single_plane.mid_eigen_value_);
  plane.scale.z = 2 * sqrt(single_plane.min_eigen_value_);//法线方向（厚度）上的弥散
  plane.color.a = alpha;//透明度
  plane.color.r = rgb(0);
  plane.color.g = rgb(1);
  plane.color.b = rgb(2);
  plane.lifetime = ros::Duration();
  plane_pub.markers.push_back(plane);
}

void VoxelMapManager::CalcVectQuation(const Eigen::Vector3d &x_vec, const Eigen::Vector3d &y_vec, const Eigen::Vector3d &z_vec,
                                      geometry_msgs::Quaternion &q)//计算得到四元数 (将一组三维正交向量（描述空间方向）转换为四元数)
{
  Eigen::Matrix3d rot;
  rot << x_vec(0), x_vec(1), x_vec(2), y_vec(0), y_vec(1), y_vec(2), z_vec(0), z_vec(1), z_vec(2);
  Eigen::Matrix3d rotation = rot.transpose();
  Eigen::Quaterniond eq(rotation);
  q.w = eq.w();
  q.x = eq.x();
  q.y = eq.y();
  q.z = eq.z();
}

void VoxelMapManager::mapJet(double v, double vmin, double vmax, uint8_t &r, uint8_t &g, uint8_t &b)//根据trace计算颜色
{
  r = 255;
  g = 255;
  b = 255;

  if (v < vmin) { v = vmin; }

  if (v > vmax) { v = vmax; }

  double dr, dg, db;

  if (v < 0.1242)
  {
    db = 0.504 + ((1. - 0.504) / 0.1242) * v;
    dg = dr = 0.;
  }
  else if (v < 0.3747)
  {
    db = 1.;
    dr = 0.;
    dg = (v - 0.1242) * (1. / (0.3747 - 0.1242));
  }
  else if (v < 0.6253)
  {
    db = (0.6253 - v) * (1. / (0.6253 - 0.3747));
    dg = 1.;
    dr = (v - 0.3747) * (1. / (0.6253 - 0.3747));
  }
  else if (v < 0.8758)
  {
    db = 0.;
    dr = 1.;
    dg = (0.8758 - v) * (1. / (0.8758 - 0.6253));
  }
  else
  {
    db = 0.;
    dg = 0.;
    dr = 1. - (v - 0.8758) * ((1. - 0.504) / (1. - 0.8758));
  }

  r = (uint8_t)(255 * dr);
  g = (uint8_t)(255 * dg);
  b = (uint8_t)(255 * db);
}

void VoxelMapManager::mapSliding()//地图滑动调整
{
  if((position_last_ - last_slide_position).norm() < config_setting_.sliding_thresh)
  {
    std::cout<<RED<<"[DEBUG]: Last sliding length "<<(position_last_ - last_slide_position).norm()<<RESET<<"\n";
    return;
  }

  //get global id now
  last_slide_position = position_last_;
  double t_sliding_start = omp_get_wtime();
  float loc_xyz[3];
  for (int j = 0; j < 3; j++)
  {
    loc_xyz[j] = position_last_[j] / config_setting_.max_voxel_size_;
    if (loc_xyz[j] < 0) { loc_xyz[j] -= 1.0; }
  }
  // VOXEL_LOCATION position((int64_t)loc_xyz[0], (int64_t)loc_xyz[1], (int64_t)loc_xyz[2]);//discrete global
  clearMemOutOfMap((int64_t)loc_xyz[0] + config_setting_.half_map_size, (int64_t)loc_xyz[0] - config_setting_.half_map_size,
                    (int64_t)loc_xyz[1] + config_setting_.half_map_size, (int64_t)loc_xyz[1] - config_setting_.half_map_size,
                    (int64_t)loc_xyz[2] + config_setting_.half_map_size, (int64_t)loc_xyz[2] - config_setting_.half_map_size);
  double t_sliding_end = omp_get_wtime();
  std::cout<<RED<<"[DEBUG]: Map sliding using "<<t_sliding_end - t_sliding_start<<" secs"<<RESET<<"\n";
  return;
}

void VoxelMapManager::clearMemOutOfMap(const int& x_max,const int& x_min,const int& y_max,const int& y_min,const int& z_max,const int& z_min )//清除超出地图范围的体素
{
  int delete_voxel_cout = 0;
  // double delete_time = 0;
  // double last_delete_time = 0;
  for (auto it = voxel_map_.begin(); it != voxel_map_.end(); )
  {
    const VOXEL_LOCATION& loc = it->first;
    bool should_remove = loc.x > x_max || loc.x < x_min || loc.y > y_max || loc.y < y_min || loc.z > z_max || loc.z < z_min;
    if (should_remove){
      // last_delete_time = omp_get_wtime();
      delete it->second;
      it = voxel_map_.erase(it);
      // delete_time += omp_get_wtime() - last_delete_time;
      delete_voxel_cout++;
    } else {
      ++it;
    }
  }
  std::cout<<RED<<"[DEBUG]: Delete "<<delete_voxel_cout<<" root voxels"<<RESET<<"\n";
  // std::cout<<RED<<"[DEBUG]: Delete "<<delete_voxel_cout<<" voxels using "<<delete_time<<" s"<<RESET<<"\n";
}