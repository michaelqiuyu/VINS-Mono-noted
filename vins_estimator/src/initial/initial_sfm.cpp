#include "initial_sfm.h"

GlobalSFM::GlobalSFM(){}

/**
 * @brief 对特征点三角化
 * 
 * @param[in] Pose0 两帧位姿
 * @param[in] Pose1 
 * @param[in] point0 特征点在两帧下的观测
 * @param[in] point1 
 * @param[out] point_3d 三角化结果
 */

void GlobalSFM::triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
						Vector2d &point0, Vector2d &point1, Vector3d &point_3d)
{
	// 通过奇异值分解求解一个Ax = 0得到
	Matrix4d design_matrix = Matrix4d::Zero();
	design_matrix.row(0) = point0[0] * Pose0.row(2) - Pose0.row(0);
	design_matrix.row(1) = point0[1] * Pose0.row(2) - Pose0.row(1);
	design_matrix.row(2) = point1[0] * Pose1.row(2) - Pose1.row(0);
	design_matrix.row(3) = point1[1] * Pose1.row(2) - Pose1.row(1);
	Vector4d triangulated_point;
	triangulated_point =
		      design_matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();
	// 齐次向量归一化
	point_3d(0) = triangulated_point(0) / triangulated_point(3);
	point_3d(1) = triangulated_point(1) / triangulated_point(3);
	point_3d(2) = triangulated_point(2) / triangulated_point(3);
}

/**
 * @brief 根据上一帧的位姿通过pnp求解当前帧的位姿
 * 
 * @param[in] R_initial 上一帧的位姿
 * @param[in] P_initial 
 * @param[in] i 	当前帧的索引
 * @param[in] sfm_f 	所有特征点的信息
 * @return true 
 * @return false 
 */

bool GlobalSFM::solveFrameByPnP(Matrix3d &R_initial, Vector3d &P_initial, int i,
								vector<SFMFeature> &sfm_f)
{
	vector<cv::Point2f> pts_2_vector;
	vector<cv::Point3f> pts_3_vector;
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state != true) // 是false就是没有被三角化，pnp是3d到2d求解，因此需要3d点
			continue;
		Vector2d point2d;
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			if (sfm_f[j].observation[k].first == i)
			{
				Vector2d img_pts = sfm_f[j].observation[k].second;  // 这里实际上是去畸变的归一化相机系的x和y坐标
				cv::Point2f pts_2(img_pts(0), img_pts(1));
				pts_2_vector.push_back(pts_2);
				cv::Point3f pts_3(sfm_f[j].position[0], sfm_f[j].position[1], sfm_f[j].position[2]);
				pts_3_vector.push_back(pts_3);
				break;
			}
		}
	}
	if (int(pts_2_vector.size()) < 15)
	{
		printf("unstable features tracking, please slowly move you device!\n");
		if (int(pts_2_vector.size()) < 10)
			return false;
	}
	cv::Mat r, rvec, t, D, tmp_r;
	cv::eigen2cv(R_initial, tmp_r);
	cv::Rodrigues(tmp_r, rvec);
	cv::eigen2cv(P_initial, t);
	/**
	 * notes:
	 *      这里的K矩阵是单位阵的原因是，这里的二维点是归一化相机系坐标而不是像素坐标
	 */
	cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
	bool pnp_succ;
	/**
	 * notes:
	 *      1. 这里使用的是CV_ITERATIVE，也就是将其建模为重投影误差最小化的问题，使用LM进行求解，这里传入的R，t将作为初值进行迭代优化
	 *      2. 这里还可以使用CV_P3P、CV_EPNP
	 */
	pnp_succ = cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1);
	if(!pnp_succ)
	{
		return false;
	}
	cv::Rodrigues(rvec, r);
	//cout << "r " << endl << r << endl;
	MatrixXd R_pnp;
	cv::cv2eigen(r, R_pnp);
	MatrixXd T_pnp;
	cv::cv2eigen(t, T_pnp);
	R_initial = R_pnp;
	P_initial = T_pnp;
	return true;

}

Eigen::Matrix3d Hat(const Eigen::Vector3d &v) {
    Eigen::Matrix3d mat;
    mat << 0, -v(2, 0), v(1, 0),
           v(2, 0), 0, -v(0, 0),
           -v(1, 0), v(0, 0), 0;
    return mat;
}

/**
 * @brief 根据两帧索引和位姿计算对应特征点的三角化位置
 * 
 * @param[in] frame0 
 * @param[in] Pose0 
 * @param[in] frame1 
 * @param[in] Pose1 
 * @param[in] sfm_f 
 */
void GlobalSFM::triangulateTwoFrames(int frame0, Eigen::Matrix<double, 3, 4> &Pose0, 
									 int frame1, Eigen::Matrix<double, 3, 4> &Pose1,
									 vector<SFMFeature> &sfm_f)
{
	assert(frame0 != frame1);
	for (int j = 0; j < feature_num; j++)	// feature_num是特征点总数
	{
		if (sfm_f[j].state == true)	// 已经三角化过了：优化的点在于可以使用多视图来进行优化
			continue;
		bool has_0 = false, has_1 = false;
		Vector2d point0;
		Vector2d point1;
		// 遍历该特征点的观测，看看是不能两帧都能看到：这里可以看出VINS-MONO的框架组织的不如ORB-SLAM，很多地方都需要遍历操作，而不是直接对抽象的对象（如frame）进行操作
		for (int k = 0; k < (int)sfm_f[j].observation.size(); k++)
		{
			if (sfm_f[j].observation[k].first == frame0)
			{
				point0 = sfm_f[j].observation[k].second;	// 取出在该帧的观测
				has_0 = true;
			}
			if (sfm_f[j].observation[k].first == frame1)
			{
				point1 = sfm_f[j].observation[k].second;
				has_1 = true;
			}
		}
		if (has_0 && has_1)	// 如果都能被看到
		{
			Vector3d point_3d;
			// 将这个特征点进行DLT三角化
			triangulatePoint(Pose0, Pose1, point0, point1, point_3d);
			sfm_f[j].state = true;	// 标志位置true
			sfm_f[j].position[0] = point_3d(0);
			sfm_f[j].position[1] = point_3d(1);
			sfm_f[j].position[2] = point_3d(2);
			//cout << "trangulated : " << frame1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;

#if 0
        // 测试：如果认为三维点在第一帧的射线上，只有深度为变量，求解深度值并恢复三维点，与DLT的三维点做比较
        Eigen::Vector3d u1 = Eigen::Vector3d(point0(0, 0), point0(1, 0), 1);
        Eigen::Vector3d u2 = Eigen::Vector3d(point1(0, 0), point1(1, 0), 1);

        Eigen::Matrix3d R1 = Pose0.block<3, 3>(0, 0);
        Eigen::Vector3d t1 = Pose0.block<3, 1>(0, 3);

        Eigen::Matrix3d R2 = Pose1.block<3, 3>(0, 0);
        Eigen::Vector3d t2 = Pose1.block<3, 1>(0, 3);

        Eigen::Vector3d p1 = Hat(u2) * R2 * R1.transpose() * u1;
        Eigen::Vector3d p2 = -Hat(u2) * (t2 - R2 * R1.transpose() * t1);

        double s = (p1.transpose() * p2)(0, 0) / (p1.transpose() * p1)(0, 0);
        std::cout << "s = " << s << std::endl;

        Eigen::Vector3d point = s * R1.transpose() * u1 - R1.transpose() * t1;
        std::cout << "根据其在第一帧下的投影为point0，得到的世界点的坐标为：" << point.transpose() << std::endl;

        std::cout << "DLT三角化的结果为：" << point_3d.transpose() << std::endl;

        Eigen::Vector3d error = point - point_3d;
        if (error.norm() > 0.1)
            std::cout << "此时的误差不可忽略，误差为：" << error.transpose() << ", norm = " << error.norm() << std::endl;

//        Eigen::Vector3d dlt_p1 = R1 * point_3d + t1;
//        Eigen::Vector3d dlt_p2 = R2 * point_3d + t2;
//        std::cout << "dlt_p1 = " << dlt_p1.transpose() << ", dlt_p2 = " << dlt_p2.transpose() << std::endl;
//
//        Eigen::Vector3d s_p1 = R1 * point + t1;
//        Eigen::Vector3d s_p2 = R2 * point + t2;
//        Eigen::Vector3d s_p3 = s * u1;
//        std::cout << "s_p1 = " << s_p1.transpose() << ", s_p2 = " << s_p2.transpose() << ", s_p3 = " << s_p3.transpose() << std::endl;


        std::cout << "\n" << std::endl;
#endif
		}							  
	}
}

// 	 q w_R_cam t w_R_cam
//  c_rotation cam_R_w 
//  c_translation cam_R_w
// relative_q[i][j]  j_q_i
// relative_t[i][j]  j_t_ji  (j < i)
/**
 * @brief 根据已有的枢纽帧和最后一帧的位姿变换，得到各帧位姿和3d点坐标，最后通过ceres进行优化
 * 
 * @param[in] frame_num 滑窗内KF总数 
 * @param[out] q  恢复出来的滑窗中各个姿态
 * @param[out] T  恢复出来的滑窗中各个平移
 * @param[in] l 	枢纽帧的idx
 * @param[in] relative_R 	枢纽帧和最后一帧的旋转
 * @param[in] relative_T 	枢纽帧和最后一帧的平移
 * @param[in] sfm_f 	用来做sfm的特征点集合
 * @param[out] sfm_tracked_points 恢复出来的地图点
 * @return true 
 * @return false 
 */
bool GlobalSFM::construct(int frame_num, Quaterniond* q, Vector3d* T, int l,
			  const Matrix3d relative_R, const Vector3d relative_T,
			  vector<SFMFeature> &sfm_f, map<int, Vector3d> &sfm_tracked_points)
{
    /**
     * sfm_f获得的地图点，仅仅只是用于做一个单目的GBA，并不会将结果直接更新到特征点管理器中的特征点中
     * 后面也只是使用到了位姿，并用于重新多帧三角化，而没有使用这里的地图点信息
     */
	feature_num = sfm_f.size();
	//cout << "set 0 and " << l << " as known " << endl;
	// have relative_r relative_t
	// intial two view
	// 枢纽帧设置为单位帧，也可以理解为世界系
	q[l].w() = 1;
	q[l].x() = 0;
	q[l].y() = 0;
	q[l].z() = 0;
	T[l].setZero();
	// 求得最后一帧的位姿：camera → world
	q[frame_num - 1] = q[l] * Quaterniond(relative_R);
	T[frame_num - 1] = relative_T;
	//cout << "init q_l " << q[l].w() << " " << q[l].vec().transpose() << endl;
	//cout << "init t_l " << T[l].transpose() << endl;

	// 由于纯视觉slam处理都是Tcw,因此下面把Twc转成Tcw
	//rotate to cam frame
	Matrix3d c_Rotation[frame_num];
	Vector3d c_Translation[frame_num];
	Quaterniond c_Quat[frame_num];
	double c_rotation[frame_num][4];
	double c_translation[frame_num][3];
	Eigen::Matrix<double, 3, 4> Pose[frame_num];

	// 将枢纽帧和最后一帧Twc转成Tcw，包括四元数，旋转矩阵，平移向量和增广矩阵
	c_Quat[l] = q[l].inverse();
	c_Rotation[l] = c_Quat[l].toRotationMatrix();
	c_Translation[l] = -1 * (c_Rotation[l] * T[l]);
	Pose[l].block<3, 3>(0, 0) = c_Rotation[l];
	Pose[l].block<3, 1>(0, 3) = c_Translation[l];

	c_Quat[frame_num - 1] = q[frame_num - 1].inverse();
	c_Rotation[frame_num - 1] = c_Quat[frame_num - 1].toRotationMatrix();
	c_Translation[frame_num - 1] = -1 * (c_Rotation[frame_num - 1] * T[frame_num - 1]);
	Pose[frame_num - 1].block<3, 3>(0, 0) = c_Rotation[frame_num - 1];
	Pose[frame_num - 1].block<3, 1>(0, 3) = c_Translation[frame_num - 1];

	// 以上准备工作做好后开始具体实现

	//1: trangulate between l ----- frame_num - 1
	//2: solve pnp l + 1; trangulate l + 1 ------- frame_num - 1; 
	// Step 1 求解枢纽帧到最后一帧之间帧的位姿及对应特征点的三角化处理
	for (int i = l; i < frame_num - 1 ; i++)
	{
		// solve pnp
		if (i > l)
		{
			// 这是依次求解，因此上一帧的位姿是已知量
			Matrix3d R_initial = c_Rotation[i - 1];
			Vector3d P_initial = c_Translation[i - 1];
			if(!solveFrameByPnP(R_initial, P_initial, i, sfm_f))
				return false;
			c_Rotation[i] = R_initial;
			c_Translation[i] = P_initial;
			c_Quat[i] = c_Rotation[i];
			Pose[i].block<3, 3>(0, 0) = c_Rotation[i];
			Pose[i].block<3, 1>(0, 3) = c_Translation[i];
		}
		/**
		 * notes:
		 *      1. 这里按照l~frame-1的顺序进行pnp求解和三角化的一个理由是：前面三角化的点 对后面的帧一定可见，因此pnp求解的匹配数目就会增加，求解的稳定性增强
		 */

		// triangulate point based on the solve pnp result
		// 当前帧和最后一帧进行三角化处理
		triangulateTwoFrames(i, Pose[i], frame_num - 1, Pose[frame_num - 1], sfm_f);
	}
	// Step 2 考虑有些特征点不能被最后一帧看到，因此，fix枢纽帧，遍历枢纽帧到最后一帧进行特征点三角化
	//3: triangulate l-----l+1 l+2 ... frame_num -2
	for (int i = l + 1; i < frame_num - 1; i++)
	    // 对后面求解枢纽帧之前的图像位姿非常重要
		triangulateTwoFrames(l, Pose[l], i, Pose[i], sfm_f);
	/**
	 * notes: 这里并没有三角化枢纽帧到最后一帧之间的这些帧之间的三角化；这些特征点留在了step4中进行三角化
	 */
	// Step 3 处理完枢纽帧到最后一帧，开始处理枢纽帧之前的帧
	//4: solve pnp l-1; triangulate l-1 ----- l
	//             l-2              l-2 ----- l
	for (int i = l - 1; i >= 0; i--)
	{
		//solve pnp
		// 这种情况就是后一帧先求解出来：如果直接求解第0帧的位姿，可能由于没有足够的2D-3D匹配导致位姿求解失败，倒序求解有利于缓解这个问题
		/**
		 * notes: 为什么枢纽帧之前的位姿使用倒序求解，枢纽帧之后的位姿使用顺序求解
		 *      1. 枢纽帧之后顺序求解是因为枢纽帧到最后一帧的位姿已经已知了，那么按照顺序的求解方式，前面帧三角化的地图点会后面图像的求解有利
		 *      2. 枢纽帧之前倒序求解是因为直接求解第0帧的位姿，可能由于没有足够的2D-3D匹配导致位姿求解失败，倒序求解有利于缓解这个问题，也就是后面帧与枢纽帧三角化的地图点有可能被前面帧看到
		 */
		Matrix3d R_initial = c_Rotation[i + 1];
		Vector3d P_initial = c_Translation[i + 1];
		if(!solveFrameByPnP(R_initial, P_initial, i, sfm_f))
			return false;
		c_Rotation[i] = R_initial;
		c_Translation[i] = P_initial;
		c_Quat[i] = c_Rotation[i];
		Pose[i].block<3, 3>(0, 0) = c_Rotation[i];
		Pose[i].block<3, 1>(0, 3) = c_Translation[i];
		//triangulate
		triangulateTwoFrames(i, Pose[i], l, Pose[l], sfm_f);
	}
	// Step 4 得到了所有关键帧的位姿，遍历没有被三角化的特征点，进行三角化
	//5: triangulate all other points
	for (int j = 0; j < feature_num; j++)
	{
		if (sfm_f[j].state == true)
			continue;
		if ((int)sfm_f[j].observation.size() >= 2)	// 只有被两个以上的KF观测到才可以三角化
		{
			Vector2d point0, point1;
			// 取首尾两个KF，尽量保证两KF之间足够位移：为什么不使用所有观测来三角化尼，这样精确度会更高
			int frame_0 = sfm_f[j].observation[0].first;
			point0 = sfm_f[j].observation[0].second;
			int frame_1 = sfm_f[j].observation.back().first;
			point1 = sfm_f[j].observation.back().second;
			Vector3d point_3d;
			triangulatePoint(Pose[frame_0], Pose[frame_1], point0, point1, point_3d);
			sfm_f[j].state = true;
			sfm_f[j].position[0] = point_3d(0);
			sfm_f[j].position[1] = point_3d(1);
			sfm_f[j].position[2] = point_3d(2);
			//cout << "trangulated : " << frame_0 << " " << frame_1 << "  3d point : "  << j << "  " << point_3d.transpose() << endl;
		}		
	}

/*
	for (int i = 0; i < frame_num; i++)
	{
		q[i] = c_Rotation[i].transpose(); 
		cout << "solvePnP  q" << " i " << i <<"  " <<q[i].w() << "  " << q[i].vec().transpose() << endl;
	}
	for (int i = 0; i < frame_num; i++)
	{
		Vector3d t_tmp;
		t_tmp = -1 * (q[i] * c_Translation[i]);
		cout << "solvePnP  t" << " i " << i <<"  " << t_tmp.x() <<"  "<< t_tmp.y() <<"  "<< t_tmp.z() << endl;
	}
*/
#if 0
    std::cout << std::fixed << std::setprecision(15);
    for (int j = 0; j < feature_num; j++)
    {
        if (sfm_f[j].state) {
            int id = sfm_f[j].id;
            if (10 <= id && 30 >= id)
                std::cout << "id = " << id << ", position = " << sfm_f[j].position[0] << " " << sfm_f[j].position[1] << " " << sfm_f[j].position[2] << std::endl;
        }
    }
#endif
	//full BA
	// Step 5 求出了所有的位姿和3d点之后，进行一个视觉slam的global BA
	// 可能需要介绍一下ceres  http://ceres-solver.org/
	ceres::Problem problem;
	// 设置更新方式，不同的表示的更新方式是不同的
	ceres::LocalParameterization* local_parameterization = new ceres::QuaternionParameterization();  // 四元数采用的是右乘更新
	//cout << " begin full BA " << endl;
#if 0
    double temp_position[frame_num][3];
    for (int k = 0; k < feature_num; k++) {
        if (sfm_f[k].state != true)    // 必须是三角化之后的，也就是有对应的世界系坐标
            continue;
        temp_position[k][0] = sfm_f[k].position[0];
        temp_position[k][1] = sfm_f[k].position[1];
        temp_position[k][2] = sfm_f[k].position[2];
    }
#endif
	for (int i = 0; i < frame_num; i++)
	{
		//double array for ceres
		// 这些都是待优化的参数块
		c_translation[i][0] = c_Translation[i].x();
		c_translation[i][1] = c_Translation[i].y();
		c_translation[i][2] = c_Translation[i].z();
		c_rotation[i][0] = c_Quat[i].w();
		c_rotation[i][1] = c_Quat[i].x();
		c_rotation[i][2] = c_Quat[i].y();
		c_rotation[i][3] = c_Quat[i].z();
		// 只有对于特殊的参数块才需要使用AddParameterBlock
		// 旋转是特殊的；后面对最后一帧的平移是固定的SetParameterBlockConstant，因此也需要添加参数块才能取出这一帧，否则会报找不到参数块的错误
		// 需要注意的是，地图点的参数块并没有使用AddParameterBlock来添加，但是这并没有影响
		problem.AddParameterBlock(c_rotation[i], 4, local_parameterization);
		problem.AddParameterBlock(c_translation[i], 3);
		// 由于是单目视觉slam，有七个自由度不可观，因此，fix一些参数块避免在零空间漂移
		// fix设置的世界坐标系第l帧的位姿，同时fix最后一帧的位移用来fix尺度
		if (i == l)
		{
			problem.SetParameterBlockConstant(c_rotation[i]);
		}
		if (i == l || i == frame_num - 1)
		{
			problem.SetParameterBlockConstant(c_translation[i]);
		}
	}

	// 只有视觉重投影构成约束，因此遍历所有的特征点，构建约束
	for (int i = 0; i < feature_num; i++)
	{
		if (sfm_f[i].state != true)	// 必须是三角化之后的，也就是有对应的世界系坐标
			continue;
#if 0
        problem.AddParameterBlock(position[i], 3);
#endif

		// 遍历所有的观测帧，对这些帧建立约束
		for (int j = 0; j < int(sfm_f[i].observation.size()); j++)
		{
			int l = sfm_f[i].observation[j].first;
			ceres::CostFunction* cost_function = ReprojectionError3D::Create(
												sfm_f[i].observation[j].second.x(),
												sfm_f[i].observation[j].second.y());
			// 约束了这一帧位姿和3d地图点
    		problem.AddResidualBlock(cost_function, NULL, c_rotation[l], c_translation[l],
    								sfm_f[i].position);  // 没有使用核函数
#if 0
            problem.AddResidualBlock(cost_function, NULL, c_rotation[l], c_translation[l],
                                     position[i]);  // 没有使用核函数
#endif
		}

	}
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	//options.minimizer_progress_to_stdout = true;
	options.max_solver_time_in_seconds = 0.2;  // 最大求解时间
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	//std::cout << summary.BriefReport() << "\n";
	// 收敛或者损失在阈值范围内
	if (summary.termination_type == ceres::CONVERGENCE || summary.final_cost < 5e-03)
	{
		//cout << "vision only BA converge" << endl;
	}
	else
	{
		//cout << "vision only BA not converge " << endl;
		return false;
	}
	// 优化结束，把double数组的值返回成对应类型的值
	// 同时Tcw -> Twc
	for (int i = 0; i < frame_num; i++)
	{
		q[i].w() = c_rotation[i][0]; 
		q[i].x() = c_rotation[i][1]; 
		q[i].y() = c_rotation[i][2]; 
		q[i].z() = c_rotation[i][3]; 
		q[i] = q[i].inverse();
		//cout << "final  q" << " i " << i <<"  " <<q[i].w() << "  " << q[i].vec().transpose() << endl;
	}
	for (int i = 0; i < frame_num; i++)
	{

		T[i] = -1 * (q[i] * Vector3d(c_translation[i][0], c_translation[i][1], c_translation[i][2]));
		//cout << "final  t" << " i " << i <<"  " << T[i](0) <<"  "<< T[i](1) <<"  "<< T[i](2) << endl;
	}
	// 存储地图点的坐标
	for (int i = 0; i < (int)sfm_f.size(); i++)
	{
		if(sfm_f[i].state)
			sfm_tracked_points[sfm_f[i].id] = Vector3d(sfm_f[i].position[0], sfm_f[i].position[1], sfm_f[i].position[2]);
	}

#if 0
    std::cout << std::fixed << std::setprecision(15);
    for (int j = 0; j < feature_num; j++)
    {
        if (sfm_f[j].state) {
            int id = sfm_f[j].id;
            if (10 <= id && 30 >= id)
                std::cout << "id = " << id << ", position = " << sfm_f[j].position[0] << " " << sfm_f[j].position[1] << " " << sfm_f[j].position[2] << std::endl;
        }
    }
#endif

	return true;

}

