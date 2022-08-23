#pragma once 
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <deque>
#include <map>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>
using namespace Eigen;
using namespace std;



struct SFMFeature
{
    bool state;
    int id;  // 特征点的ID
    vector<pair<int,Vector2d>> observation;  // 帧号以及去畸变的相机系归一化坐标的x和y
    double position[3];
    double depth;
};

struct ReprojectionError3D
{
	ReprojectionError3D(double observed_u, double observed_v)
		:observed_u(observed_u), observed_v(observed_v)
		{}

    /**
     * notes:
     *      1. 最后一维必须是残差结果
     *      2. 参数与AddResidualBlock函数中kernel function之后的参数对应
     */
	template <typename T>
	bool operator()(const T* const camera_R, const T* const camera_T, const T* point, T* residuals) const
	{
		T p[3];
		ceres::QuaternionRotatePoint(camera_R, point, p);	// 旋转这个点
		p[0] += camera_T[0]; p[1] += camera_T[1]; p[2] += camera_T[2];	// 这其实就是Rcw * pw + tcw
		// 得到该相机坐标系下的3d坐标
		T xp = p[0] / p[2];
    	T yp = p[1] / p[2];	// 归一化处理
    	// 跟现有观测形成残差：注意这里使用的是归一化相机系下面的残差
    	/**
    	 * u = fx * x + cx, v = fy * y + cy
    	 *
    	 * (u - u0)^2 + (v - v0)^2 = fx^2 * (x - x0)^2 + fy^2 * (y - y0)^2 ！= (x - x0)^2 + (y - y0)^2
    	 * 如果fx和fy并不接近，那么优化归一化相机系的残差相当于认为给重投影误差加了权重，这样做并不符合重投影误差的定义，其与重投影误差并不等效
    	 */
    	residuals[0] = xp - T(observed_u);
    	residuals[1] = yp - T(observed_v);
    	return true;
	}

	static ceres::CostFunction* Create(const double observed_x,
	                                   const double observed_y) 
	{
	    /**
	     * notes:
	     *      2-residuals dimension
	     *      4-rotation dimension
	     *      3-translation dimension
	     *      3-point dimension
	     */
	  return (new ceres::AutoDiffCostFunction<  // 自动求导需要重载一下括号运算符
	          ReprojectionError3D, 2, 4, 3, 3>(
	          	new ReprojectionError3D(observed_x,observed_y)));
	}

	double observed_u;
	double observed_v;
};

class GlobalSFM
{
public:
	GlobalSFM();
	bool construct(int frame_num, Quaterniond* q, Vector3d* T, int l,
			  const Matrix3d relative_R, const Vector3d relative_T,
			  vector<SFMFeature> &sfm_f, map<int, Vector3d> &sfm_tracked_points);

private:
	bool solveFrameByPnP(Matrix3d &R_initial, Vector3d &P_initial, int i, vector<SFMFeature> &sfm_f);

	void triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
							Vector2d &point0, Vector2d &point1, Vector3d &point_3d);
	void triangulateTwoFrames(int frame0, Eigen::Matrix<double, 3, 4> &Pose0, 
							  int frame1, Eigen::Matrix<double, 3, 4> &Pose1,
							  vector<SFMFeature> &sfm_f);

	int feature_num;
};