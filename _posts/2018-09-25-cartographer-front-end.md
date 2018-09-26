---
layout:     post   				    # 使用的布局（不需要改）
title:      Cartographer 的前端算法思路				# 标题 
subtitle:   #副标题
date:       2018-09-25 				# 时间
author:     Edward Liu 						# 作者
header-img: #这篇文章标题背景图片
catalog: true 						# 是否归档
tags:								#标签
    - SLAM
    - Cartographer
---

## Cartographer 的前端算法思路

前一篇博客里面提到的是 Cartographer 前端实现中非常小的一个部分的算法思路，参照了《Real time correlative scan matching》里的算法实现了一部分实时scan match 的功能，不过这并不是Cartographer中前端的全部，甚至是可以通过参数disable的一部分功能。   
在 Cartographer 对应的论文《Real-Time Loop Closure in 2D LIDAR SLAM》中提到的前端算法中只有Ceres scan matching，其实就是基于ceres solver实现的非线性优化模型，今天我们看一下具体的算法模型。同样，在读懂算法和代码之前需要一些基础知识：

> [Ceres Solver tutorial](http://ceres-solver.org/index.html) & 最小二乘求解非线性优化问题  
> [Cartographer Git](https://github.com/googlecartographer/cartographer)  
> [双三次插值](https://www.wikiwand.com/zh-hans/%E5%8F%8C%E4%B8%89%E6%AC%A1%E6%8F%92%E5%80%BC)  

### 论文内容简单翻译
这里我们主要看论文的"IV. Local 2D SLAM" —— "C. ceres scan matching"部分：
> Prior to inserting a scan into a submap, the scan pose ξ is optimized relative to the current local submap using a Ceres-based [14] scan matcher. The scan matcher is responsible for finding a scan pose that maximizes the probabilities at the scan points in the submap. We cast this as a nonlinear least squares problem  
> **我们用一个基于 ceres solver 的 scan matcher 优化获得当前的 scan pose $\xi$，这个 scan matcher 负责找到一个位姿使得 scan 中的所有点在当前 local map 中的概率和最大，于是我们定义一下最小二乘问题：**  
> $$ \underset{\xi}{\arg\min}\sum_{k=1}^K(1-M_{smooth}(T_{\xi}h_k))^2 \tag{CS} $$  
> where $T_{\xi}$ transforms $h_k$ from the scan frame to the submap frame according to the scan pose. The function $ M_{smooth} : \mathbb{R}^2 → \mathbb{R} $ is a smooth version of the probability values in the local submap. We use **bicubic interpolation**. As a result, values outside the interval [0, 1] can occur but are considered harmless.  
> ***式中 $T_{\xi}$ 将 $h_k$ 中的 scan 点全部转化到 local map 坐标系下，$ M_{smooth} : \mathbb{R}^2 → \mathbb{R} $  则是一个将 local map 中的各点概率进行一个平滑处理的函数，这里我们用双三次插值。这样可能会出现概率小于0或者大于1的情况，不过这种情况并不会产生错误。***  
> Mathematical optimization of this smooth function usually gives better precision than the resolution of the grid. Since this is a local optimization, good initial estimates are required. An IMU capable of measuring angular velocities can be used to estimate the rotational component $\theta$ of the pose between scan matches. A higher frequency of scan matches or a pixel-accurate scan matching approach, although more computationally intensive, can be used in the absence of an IMU.  
> ***数学优化问题通常会提供一个比网格地图的分辨率精度更高的优化结果。由于这是一个实时的局部优化，需要一个好的初始位姿估计。我们可以用 IMU 来估计 scan match 中的旋转角度 $\theta$，当然如果 scan matching 的频率很高，是可以不使用IMU的。***

### 算法与代码分析
Cartographer 的 ceres scan matcher 将上面的最小二成问题分解成了三个 Cost Function：
- Translation delta cost function
- Rotational drlta cost function
- ★★★ Occupied space cost function ★★★

#### Ceres Scan Matcher
由于代码量很小而且不难懂，我们直接从代码来看算法可能比较直观（用注释来简单解释算法结构，下面算法的Cost functions的解释也都放在**代码注释**里）：
```cpp
/*
 input:
 1.上一个 scan 的位姿 previous_pose  
 2.当前的 scan 的位姿的初始估计 initial_pose_estimate  
 3.当前 scan 点云（2D）point_cloud  
 4.local map 概率分布栅格图 probability_grid  
 output  
 1. 计算得到的位姿估计 pose_estimate  
 2. ceres solver 计算的总结 summary  
*/  
void CeresScanMatcher::Match(const transform::Rigid2d& previous_pose,
                             const transform::Rigid2d& initial_pose_estimate,
                             const sensor::PointCloud& point_cloud,
                             const ProbabilityGrid& probability_grid,
                             transform::Rigid2d* const pose_estimate,
                             ceres::Solver::Summary* const summary) const 
{
    // 优化的变量（3个）
    double ceres_pose_estimate[3] = {initial_pose_estimate.translation().x(),
                                     initial_pose_estimate.translation().y(),
                                     initial_pose_estimate.rotation().angle()};
    ceres::Problem problem;
    CHECK_GT(options_.occupied_space_weight(), 0.);
	
    // 下面分别加入了三个 Cost Function  
    // 这里的 ceres 相关的只是需要读者自行阅读 ceres solver的教程，教程写的很详细也很好理解  
    problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<OccupiedSpaceCostFunctor, ceres::DYNAMIC, 3>(
        new OccupiedSpaceCostFunctor(
          options_.occupied_space_weight() / std::sqrt(static_cast<double>(point_cloud.size())),
            point_cloud, probability_grid),
          point_cloud.size()),
        nullptr, 
        ceres_pose_estimate);
    CHECK_GT(options_.translation_weight(), 0.);
	
    problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<TranslationDeltaCostFunctor, 2, 3>(
        new TranslationDeltaCostFunctor(options_.translation_weight(),previous_pose)),
        nullptr,
        ceres_pose_estimate);
	CHECK_GT(options_.rotation_weight(), 0.);
	
    problem.AddResidualBlock(
      new ceres::AutoDiffCostFunction<RotationDeltaCostFunctor, 1, 3>(
        new RotationDeltaCostFunctor(options_.rotation_weight(),ceres_pose_estimate[2])),
      nullptr,
      ceres_pose_estimate);

    ceres::Solve(ceres_solver_options_, &problem, summary);
    *pose_estimate = transform::Rigid2d({ceres_pose_estimate[0], ceres_pose_estimate[1]}, ceres_pose_estimate[2]);
}
```
#### Translation & Rotational Cost Function
下面是两个Cost Function 的实现，就是保证 translation 和 rotation 最小的 Cost Function，很简单易懂不多做解释：
```cpp
class TranslationDeltaCostFunctor {
 public:
  // Constructs a new TranslationDeltaCostFunctor from the given
  // 'initial_pose_estimate' (x, y, theta).
  explicit TranslationDeltaCostFunctor(
      const double scaling_factor,
      const transform::Rigid2d& initial_pose_estimate)
      : scaling_factor_(scaling_factor),
        x_(initial_pose_estimate.translation().x()),
        y_(initial_pose_estimate.translation().y()) {}

  TranslationDeltaCostFunctor(const TranslationDeltaCostFunctor&) = delete;
  TranslationDeltaCostFunctor& operator=(const TranslationDeltaCostFunctor&) =
      delete;

  template <typename T>
  bool operator()(const T* const pose, T* residual) const {
	// 获得一个 2 维的残差，即 x,y 方向上的位移
    residual[0] = scaling_factor_ * (pose[0] - x_);
    residual[1] = scaling_factor_ * (pose[1] - y_);
    return true;
  }

 private:
  const double scaling_factor_;
  const double x_;
  const double y_;
};

class RotationDeltaCostFunctor {
 public:
  // Constructs a new RotationDeltaCostFunctor for the given 'angle'.
  explicit RotationDeltaCostFunctor(const double scaling_factor,
                                    const double angle)
      : scaling_factor_(scaling_factor), angle_(angle) {}

  RotationDeltaCostFunctor(const RotationDeltaCostFunctor&) = delete;
  RotationDeltaCostFunctor& operator=(const RotationDeltaCostFunctor&) = delete;

  template <typename T>
  bool operator()(const T* const pose, T* residual) const {
    //　获得一个 1 维的残差，即旋转角度
    residual[0] = scaling_factor_ * (pose[2] - angle_);
    return true;
  }

 private:
  const double scaling_factor_;
  const double angle_;
};
```
#### Occupied Space Cost Funtion
```cpp
class OccupiedSpaceCostFunctor {
 public:
  // Creates an OccupiedSpaceCostFunctor using the specified map, resolution
  // level, and point cloud.
  OccupiedSpaceCostFunctor(const double scaling_factor,
                           const sensor::PointCloud& point_cloud,
                           const ProbabilityGrid& probability_grid)
      : scaling_factor_(scaling_factor),
        point_cloud_(point_cloud),
        probability_grid_(probability_grid) {}

  OccupiedSpaceCostFunctor(const OccupiedSpaceCostFunctor&) = delete;
  OccupiedSpaceCostFunctor& operator=(const OccupiedSpaceCostFunctor&) = delete;

  template <typename T>
  bool operator()(const T* const pose, T* residual) const {
    // 第一步，将 pose 转换成一个 3 × 3 的变换矩阵
    Eigen::Matrix<T, 2, 1> translation(pose[0], pose[1]);
    Eigen::Rotation2D<T> rotation(pose[2]);
    Eigen::Matrix<T, 2, 2> rotation_matrix = rotation.toRotationMatrix();
    Eigen::Matrix<T, 3, 3> transform;
    transform << rotation_matrix, translation, T(0.), T(0.), T(1.);

    // 这里将构造时传入的概率栅格图（local submap）加载到一个双三次插值器中
    // GridArrayAdapter 的实现这里省略了，想了解具体实现的可以在 Cartographer 的代码里找到
    // 功能主要是在概率栅格图对应的 index 中取出相应的概率值
    const GridArrayAdapter adapter(probability_grid_);
    ceres::BiCubicInterpolator<GridArrayAdapter> interpolator(adapter);
    const MapLimits& limits = probability_grid_.limits();
    
    // 遍历 point_cloud_（当前 scan）中的所有点，用变换矩阵将其变换到 local map 的坐标系中
    // 取出每个点对应的栅格地图概率（双三次插值之后的）p
    // 我们要求的是这个 p 最大的情况，也就是 1-p 最小的情况，所以最后残差是 factor*(1-p)
    // 这个残差的纬度就是 point_cloud_ 中的点数
    for (size_t i = 0; i < point_cloud_.size(); ++i) {
      // Note that this is a 2D point. The third component is a scaling factor.
      const Eigen::Matrix<T, 3, 1> point((T(point_cloud_[i].x())),
                                         (T(point_cloud_[i].y())), T(1.));
      const Eigen::Matrix<T, 3, 1> world = transform * point;
      interpolator.Evaluate(
          (limits.max().x() - world[0]) / limits.resolution() - 0.5 +
              T(kPadding),
          (limits.max().y() - world[1]) / limits.resolution() - 0.5 +
              T(kPadding),
          &residual[i]);
      residual[i] = scaling_factor_ * (1. - residual[i]);
    }
    return true;
  }
  ...
  };
```
其实不难发现，这个Occupied Space Cost Function 的模型和上一篇博客中的 Real time correlative scan matching 的思路基本上一致，只是求解方法变成了最小二乘问题的求解，这样求解过程我们无需关心，只需要关心建模本身。

### 参考资料
[1] [Ceres Solver tutorial](http://ceres-solver.org/index.html)
[2] [Real-Time Loop Closure in 2D LIDAR SLAM](https://static.googleusercontent.com/media/research.google.com/zh-CN//pubs/archive/45466.pdf)