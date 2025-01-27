#include "marginalization_factor.h"

/**
 * @brief 待边缘化的各个残差块计算残差和雅克比矩阵，同时处理核函数的case
 * 
 */
void ResidualBlockInfo::Evaluate()
{
    residuals.resize(cost_function->num_residuals());   // 确定残差的维数

    std::vector<int> block_sizes = cost_function->parameter_block_sizes();  // 确定相关的参数块数目
    raw_jacobians = new double *[block_sizes.size()];   // ceres接口都是double数组，因此这里给雅克比准备数组
    jacobians.resize(block_sizes.size());
    // 这里就是把jacobians每个matrix地址赋给raw_jacobians，然后把raw_jacobians传递给ceres的接口，这样计算结果直接放进了这个matrix
    for (int i = 0; i < static_cast<int>(block_sizes.size()); i++)
    {
        jacobians[i].resize(cost_function->num_residuals(), block_sizes[i]);    // 雅克比矩阵大小 残差×变量
        raw_jacobians[i] = jacobians[i].data();
        //dim += block_sizes[i] == 7 ? 6 : block_sizes[i];
    }
    // 调用各自重载的接口计算残差和雅克比
    cost_function->Evaluate(parameter_blocks.data(), residuals.data(), raw_jacobians);  // 这里实际上结果放在了jacobians

    //std::vector<int> tmp_idx(block_sizes.size());
    //Eigen::MatrixXd tmp(dim, dim);
    //for (int i = 0; i < static_cast<int>(parameter_blocks.size()); i++)
    //{
    //    int size_i = localSize(block_sizes[i]);
    //    Eigen::MatrixXd jacobian_i = jacobians[i].leftCols(size_i);
    //    for (int j = 0, sub_idx = 0; j < static_cast<int>(parameter_blocks.size()); sub_idx += block_sizes[j] == 7 ? 6 : block_sizes[j], j++)
    //    {
    //        int size_j = localSize(block_sizes[j]);
    //        Eigen::MatrixXd jacobian_j = jacobians[j].leftCols(size_j);
    //        tmp_idx[j] = sub_idx;
    //        tmp.block(tmp_idx[i], tmp_idx[j], size_i, size_j) = jacobian_i.transpose() * jacobian_j;
    //    }
    //}
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(tmp);
    //std::cout << saes.eigenvalues() << std::endl;
    //ROS_ASSERT(saes.eigenvalues().minCoeff() >= -1e-6);
    // 如果有核函数，那么就对残差进行相关调整
    if (loss_function)
    {
        double residual_scaling_, alpha_sq_norm_;

        double sq_norm, rho[3];

        sq_norm = residuals.squaredNorm();  // 获得残差的模的平方
        loss_function->Evaluate(sq_norm, rho);  // rho[0]:核函数这个点的值 rho[1]这个点的导数 rho[2]这个点的二阶导数
        //printf("sq_norm: %f, rho[0]: %f, rho[1]: %f, rho[2]: %f\n", sq_norm, rho[0], rho[1], rho[2]);

        double sqrt_rho1_ = sqrt(rho[1]);

        /**
         * notes:
         *      1. 此处的详细论述请见：http://ceres-solver.org/nnls_modeling.html?highlight=loss#_CPPv4N5ceres12LossFunctionE
         *
         * 当residual很大的时候，需要根据robust kernel来降低大的residual的权重
         * rho[1]一般并不会小于1，一般只有在residual比较大的时候开始小于1，这也是为什么ceres-solver文档将其称为outlier region
         * rho[2]一般并不会小于0，一般只有在residual比较大的时候开始小于1，这也是为什么ceres-solver文档将其称为outlier region
         */

        /**
         * Trivial: ρ(s) = s   →   ρ'(s) = 1   →   ρ''(s) = 0
         *
         * Huber: ρ(s) = s                s <= 1   →   ρ'(s) = 1         s <= 1   →   ρ''(s) = 0                s <= 1
         *               2 * sqrt(s) - 1  s > 1                s^(-0.5)  s > 1                -0.5 * s^(-1.5)  s > 1
         *
         * SoftLone: ρ(s) = 2 * ((1 + s)^0.5 - 1)   →   ρ'(s) = (1 + s)^(-0.5)   →   ρ''(s) = -0.5 * (1 + s)^(-1.5)
         *
         * Cauchy: ρ(s) = log(1 + s)   →   ρ'(s) = (1 + s)^(-1)   →   ρ''(s) = -(1 + s)^(-2)
         *
         * Arctan: ρ(s) = arctan(s)   →   ρ'(s) = (1 + s^2)^(-1)   →   ρ''(s) = -2 * s * (1 + s^2)^(-2)
         *
         * Tolerant: ρ(s) = b * log(1 + exp((s - a) / b)) - b * log(1 + exp(-a / b))   →
         *           ρ'(s) = exp((s - a) / b) / (1 + exp((s - a) / b))   →
         *           ρ''(s) = exp((s - a) / b) / (b * (1 + exp((s - a) / b))^2)
         *
         * 从以上的常用核函数看来，一般都是ρ'(s) <= 1, ρ''(s) <= 0
         */

        // https://github.com/ceres-solver/ceres-solver/blob/master/internal/ceres/corrector.cc#L51
        if ((sq_norm == 0.0) || (rho[2] <= 0.0))
        {
            // sq_norm为0没有办法进行else中的操作，不能除0
            // rho[2]<=0是outlier region，直接降低这个residual的权重，降低的幅度是sqrt_rho1_(在outlier的时候，一般小于1)

            // 一般情况下，核函数的二阶导数小于0，因此对于函数及残差都使用同一个缩放因子控制残差的大小
            // 直观上讲，如果残差很大，使用核函数的一阶导数来缩放残差，使得这个残差不会占据主导地位，进而影响优化的方向；
            // 而且对于常用的核函数而言，残差越大，一阶导数越小，在理论上就控制了残差
            residual_scaling_ = sqrt_rho1_;
            alpha_sq_norm_ = 0.0;
        }
        else
        {  // Triggs证明：此种情况下，性能较差，收敛较慢
            // inlier region；一般要求rho[1]与rho[2]同正，那么D一定大于1
            const double D = 1.0 + 2.0 * sq_norm * rho[2] / rho[1];
            // 这里的alpha也就是上述网址中的Theory部分的alpha
            const double alpha = 1.0 - sqrt(D);  // alpha < 0 → 1 - alpha > 1
            // 在inlier region时，sqrt_rho1_一般不小于1，此时降低权重靠1 - alpha > 1
            residual_scaling_ = sqrt_rho1_ / (1 - alpha);
            alpha_sq_norm_ = alpha / sq_norm;
        }
        // 这里就相当于雅克比都乘上sqrt_rho1_，即核函数所在的点的一阶导，基本都是小于1
        for (int i = 0; i < static_cast<int>(parameter_blocks.size()); i++)
        {
            jacobians[i] = sqrt_rho1_ * (jacobians[i] - alpha_sq_norm_ * residuals * (residuals.transpose() * jacobians[i]));
        }

        residuals *= residual_scaling_;
    }
}

MarginalizationInfo::~MarginalizationInfo()
{
    //ROS_WARN("release marginlizationinfo");
    
    for (auto it = parameter_block_data.begin(); it != parameter_block_data.end(); ++it)
        delete[] it->second;

    for (int i = 0; i < (int)factors.size(); i++)
    {

        delete[] factors[i]->raw_jacobians;
        
        delete factors[i]->cost_function;

        delete factors[i];
    }
}

/**
 * @brief 收集各个残差
 * 
 * @param[in] residual_block_info 
 */
void MarginalizationInfo::addResidualBlockInfo(ResidualBlockInfo *residual_block_info)
{
    /**
     * 保存所有的残差块
     * 保存所有的参数块以及其对应的global_size
     * 保存所有的待边缘化的参数块以及其对应的global_size
     */
    factors.emplace_back(residual_block_info);  // 残差块收集起来

    std::vector<double *> &parameter_blocks = residual_block_info->parameter_blocks;    // 这个是和该约束相关的参数块
    std::vector<int> parameter_block_sizes = residual_block_info->cost_function->parameter_block_sizes();   // 各个参数块的大小，注意这里是global size

    for (int i = 0; i < static_cast<int>(residual_block_info->parameter_blocks.size()); i++)
    {
        double *addr = parameter_blocks[i];
        int size = parameter_block_sizes[i];
        // 这里是个map，避免重复添加
        parameter_block_size[reinterpret_cast<long>(addr)] = size;  // 地址->global size
    }
    // 待边缘化的参数块
    for (int i = 0; i < static_cast<int>(residual_block_info->drop_set.size()); i++)
    {
        double *addr = parameter_blocks[residual_block_info->drop_set[i]];
        // 先准备好待边缘化的参数块的map
        parameter_block_idx[reinterpret_cast<long>(addr)] = 0;
    }
}

/**
 * @brief 将各个残差块计算残差和雅克比，同时备份所有相关的参数块内容
 * 
 */
void MarginalizationInfo::preMarginalize()
{
    for (auto it : factors)
    {
        it->Evaluate(); // 调用这个接口计算各个残差块的残差和雅克比矩阵，并根据设定的核函数对残差和jacobian进行缩放

        std::vector<int> block_sizes = it->cost_function->parameter_block_sizes();  // 得到每个残差块的参数块大小
        for (int i = 0; i < static_cast<int>(block_sizes.size()); i++)
        {
            long addr = reinterpret_cast<long>(it->parameter_blocks[i]);    // 得到该参数块的地址
            int size = block_sizes[i];  // 参数块大小
            // 把各个参数块都备份起来，使用map避免重复参数块，之所以备份，是为了后面的状态保留
            // 在IMU预积分的残差构建中使用了视觉位姿的参数块，在视觉重投影的残差中依旧使用了视觉位姿的参数块
            if (parameter_block_data.find(addr) == parameter_block_data.end())  // parameter_block_data中保存了原始的参数块
            {
                double *data = new double[size];
                // 深拷贝
                memcpy(data, it->parameter_blocks[i], sizeof(double) * size);
                parameter_block_data[addr] = data;  // 地址->参数块实际内容的地址
            }
        }
    }
}

// 传入的是global size，返回的是local size
int MarginalizationInfo::localSize(int size) const
{
    return size == 7 ? 6 : size;
}

// 传入的是local size，返回的是global size
int MarginalizationInfo::globalSize(int size) const
{
    return size == 6 ? 7 : size;
}

/**
 * @brief 分线程构造Ax = b
 * 
 * @param[in] threadsstruct 
 * @return void* 
 */
void* ThreadsConstructA(void* threadsstruct)
{
    ThreadsStruct* p = ((ThreadsStruct*)threadsstruct);
    // 遍历这么多分配过来的任务
    for (auto it : p->sub_factors)
    {
        // 遍历参数块
        for (int i = 0; i < static_cast<int>(it->parameter_blocks.size()); i++)
        {
            int idx_i = p->parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[i])]; // 在大矩阵中的id，也就是落座的位置
            int size_i = p->parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[i])];
            // 确保是local size
            if (size_i == 7)
                size_i = 6;
            // 之前pre margin 已经算好了各个残差和雅克比，这里取出来
            Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
            // 和本身以及其他雅克比块构造H矩阵
            // i: 当前参数块， j: 另一个参数块
            for (int j = i; j < static_cast<int>(it->parameter_blocks.size()); j++)
            {
                int idx_j = p->parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[j])];    // 在大矩阵中的id，也就是落座的位置
                int size_j = p->parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[j])];
                if (size_j == 7)
                    size_j = 6;
                Eigen::MatrixXd jacobian_j = it->jacobians[j].leftCols(size_j);
                // 两个雅克比都取出来了
                if (i == j) // 如果是自己，一块田地
                    p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                else    // 和别人，由于对称，两块田地
                {
                    p->A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                    // 第一块田地转置一下
                    p->A.block(idx_j, idx_i, size_j, size_i) = p->A.block(idx_i, idx_j, size_i, size_j).transpose();
                }
            }
            // 然后构建g矩阵
            // xc's todo: 为什么不是-J.t * e：这里没有加负号，所以后面构建残差e的时候也没有加负号；代码最好不要这样写，比较迷惑
            p->b.segment(idx_i, size_i) += jacobian_i.transpose() * it->residuals;
        }
    }
    return threadsstruct;
}

/**
 * @brief 边缘化处理，并将结果转换成残差和雅克比的形式
 * 
 */

/**
 * 边缘化一般是伴随着稀疏求解存在的，也就是对Ax = b，且A是一个分块矩阵，其左上角矩阵和右下角矩阵互为转置，对角线上均为对角矩阵；此时一般利用边缘化的策略来求解方程，从而达到降低运算量的目的
 * 相应的技术还有sherman-morrison-woodbury formula，其可以将一个高维求逆转变为低维求逆
 *
 * 在滑窗优化中，由于窗口不停地滑动，当某一帧滑出窗口的时候，并不能简单的就认为这一帧不再对窗口优化起作用，否则，对每一个窗口优化都毫无约束，窗口优化的结果将并不符合预期
 * 与ORB-SLAM3中类比的话，对局部地图进行优化的时候，会控制局部地图中的地图点的观测（但不在局部地图中的关键帧）的这部分关键帧固定，从而限制局部地图的优化，防止其可以任意偏移和旋转
 *
 * 因此，当某一帧滑出窗口后，依然要保留其对窗口内的帧的影响，从这个含义上将，将与解Ax = b的目的不同，但是二者的做法是接近的
 */
void MarginalizationInfo::marginalize()
{
    int pos = 0;
    // parameter_block_idx key是各个待边缘化参数块地址 value预设都是0
    for (auto &it : parameter_block_idx)
    {
        // 对于待边缘化的参数块而言，第一个的idx为0，第二个待边缘化的idx是第一个待边缘化的参数块的local size，依次递推
        it.second = pos;    // 这就是在所有参数中排序的idx，待边缘化的排在前面
        pos += localSize(parameter_block_size[it.first]);   // 因为要进行求导，因此大小时local size，具体一点就是使用李代数
    }
    // 到这里，待边缘化的参数块已经遍历结束，后面就开始遍历不需要边缘化的参数块了

    m = pos;    // 总共待边缘化的参数块总大小（不是个数）
    // 其他参数块
    for (const auto &it : parameter_block_size)
    {
        // 在addResidualBlockInfo中，parameter_block_idx仅仅保存了待边缘化的参数块，在这里将所有的参数块都添加进去
        if (parameter_block_idx.find(it.first) == parameter_block_idx.end())  // 不需要边缘化的参数块
        {
            parameter_block_idx[it.first] = pos;    // 这样每个参数块的大小都能被正确找到
            pos += localSize(it.second);
        }
    }

    // n表示不需要边缘化的参数块的总计的size
    n = pos - m;    // 其他参数块的总大小

    //ROS_DEBUG("marginalization, pos: %d, m: %d, n: %d, size: %d", pos, m, n, (int)parameter_block_idx.size());

    TicToc t_summing;
    Eigen::MatrixXd A(pos, pos);    // Ax = b预设大小
    Eigen::VectorXd b(pos);
    A.setZero();
    b.setZero();
    /*
    for (auto it : factors)
    {
        for (int i = 0; i < static_cast<int>(it->parameter_blocks.size()); i++)
        {
            int idx_i = parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[i])];
            int size_i = localSize(parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[i])]);
            Eigen::MatrixXd jacobian_i = it->jacobians[i].leftCols(size_i);
            for (int j = i; j < static_cast<int>(it->parameter_blocks.size()); j++)
            {
                int idx_j = parameter_block_idx[reinterpret_cast<long>(it->parameter_blocks[j])];
                int size_j = localSize(parameter_block_size[reinterpret_cast<long>(it->parameter_blocks[j])]);
                Eigen::MatrixXd jacobian_j = it->jacobians[j].leftCols(size_j);
                if (i == j)
                    A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                else
                {
                    A.block(idx_i, idx_j, size_i, size_j) += jacobian_i.transpose() * jacobian_j;
                    A.block(idx_j, idx_i, size_j, size_i) = A.block(idx_i, idx_j, size_i, size_j).transpose();
                }
            }
            b.segment(idx_i, size_i) += jacobian_i.transpose() * it->residuals;
        }
    }
    ROS_INFO("summing up costs %f ms", t_summing.toc());
    */
    //multi thread

    // 往A矩阵和b矩阵中填东西，利用多线程加速
    TicToc t_thread_summing;
    pthread_t tids[NUM_THREADS];
    ThreadsStruct threadsstruct[NUM_THREADS];
    int i = 0;
    for (auto it : factors)
    {
        threadsstruct[i].sub_factors.push_back(it); // 每个线程均匀分配任务
        i++;
        i = i % NUM_THREADS;
    }
    // 每个线程构造一个A矩阵和b矩阵，最后大家加起来
    for (int i = 0; i < NUM_THREADS; i++)
    {
        TicToc zero_matrix;
        // 所以A矩阵和b矩阵大小一样，预设都是0
        threadsstruct[i].A = Eigen::MatrixXd::Zero(pos,pos);
        threadsstruct[i].b = Eigen::VectorXd::Zero(pos);
        // 多线程访问会带来冲突，因此每个线程备份一下要查询的map
        threadsstruct[i].parameter_block_size = parameter_block_size;   // 大小
        threadsstruct[i].parameter_block_idx = parameter_block_idx; // 索引
        // 产生若干线程
        int ret = pthread_create( &tids[i], NULL, ThreadsConstructA ,(void*)&(threadsstruct[i]));
        if (ret != 0)
        {
            ROS_WARN("pthread_create error");
            ROS_BREAK();
        }
    }
    for( int i = NUM_THREADS - 1; i >= 0; i--)  
    {
        // 等待各个线程完成各自的任务
        pthread_join( tids[i], NULL );
        // 把各个子模块拼起来，就是最终的Hx = g的矩阵了 
        A += threadsstruct[i].A;
        b += threadsstruct[i].b;
    }
    //ROS_DEBUG("thread summing up costs %f ms", t_thread_summing.toc());
    //ROS_INFO("A diff %f , b diff %f ", (A - tmp_A).sum(), (b - tmp_b).sum());


    //TODO
    // Amm矩阵的构建是为了保证其正定性
    /**
     * notes:
     *      1. 由于浮点数运算的误差，有可能A不再是一个对称矩阵，非对称矩阵不一定能够相似对角化
     *      2. 利用0.5 * (A + A.transpose())来近似A，理论上A就是一个对称矩阵，其结果就是A
     */
    Eigen::MatrixXd Amm = 0.5 * (A.block(0, 0, m, m) + A.block(0, 0, m, m).transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);   // 特征值分解

    //ROS_ASSERT_MSG(saes.eigenvalues().minCoeff() >= -1e-4, "min eigenvalue %f", saes.eigenvalues().minCoeff());
    // 一个逆矩阵的特征值是原矩阵的倒数，特征向量相同　select类似c++中 ? :运算符
    // 利用特征值取逆来构造其逆矩阵
    // A * P = P * diag()，A.inv = P * diag().inv * P.inv
    /**
     * 当存在为0或者接近0的特征值时，说明Amm奇异或者接近奇异，这个时候对Amm求逆是错误或者不稳定的
     * 详细的原理分析见SLAM14讲边缘化部分的笔记，详细阐述了当特征值为0的时候，可以按照如下的方式计算Amm的逆矩阵，这样做并不会改变问题的解
     */
    Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd((saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() * saes.eigenvectors().transpose();
    //printf("error1: %f\n", (Amm * Amm_inv - Eigen::MatrixXd::Identity(m, m)).sum());

    Eigen::VectorXd bmm = b.segment(0, m);  // 带边缘化的大小
    Eigen::MatrixXd Amr = A.block(0, m, m, n);  // 对应的四块矩阵
    Eigen::MatrixXd Arm = A.block(m, 0, n, m);

    Eigen::MatrixXd Arr = A.block(m, m, n, n);
    Eigen::VectorXd brr = b.segment(m, n); // 剩下的参数
#if 0
    std::cout << "b.row = " << b.rows() << ", b.col = " << b.cols() << std::endl;
#endif
    A = Arr - Arm * Amm_inv * Amr;
    b = brr - Arm * Amm_inv * bmm;  // b的维度发生了变化
#if 0
    std::cout << "b.row = " << b.rows() << ", b.col = " << b.cols() << std::endl;
#endif

    // 这个地方根据Ax = b => JT*J = - JT * e
    // 对A做特征值分解 A = V * S * VT,其中Ｓ是特征值构成的对角矩阵
    // 所以J = S^(1/2) * VT , 这样JT * J = (S^(1/2) * VT)T * S^(1/2) * VT = V * S^(1/2)T *  S^(1/2) * VT = V * S * VT(对角矩阵的转置等于其本身)
    // e = -(JT)-1 * b = - (S^-(1/2) * V^-1) * b = - (S^-(1/2) * VT) * b
    /**
     * notes:
     *      1. 上面为了避免浮点数运算导致A变为非对称矩阵，使用了0.5 * (A + A.T)，这里理论上也应该这样操作，否则也不应该使用SelfAdjointEigenSolver接口
     *      2. 使用SelfAdjointEigenSolver能够加速运算
     */
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(A);
    // 对A矩阵相似对角化后的矩阵求逆
    /**
     * A = p * diag * p.t = J.t * J → J = diag^0.5 * p.t
     * b = -J.t * f(x) → f(x) = -diag^(-0.5) * p.t * b；需要注意的是，前面构建b矩阵的时候省略了负号，因此这里也必须省略；最好不要这样写代码，逻辑耦合的比较强烈
     *
     * 当A奇异或者接近奇异的时候，求逆是错误的或者是不稳定的；
     * 有关奇异情形求逆的理论基础见SLAM14讲边缘化部分的笔记，详细阐述了当特征值为0的时候，可以按照如下的方式计算A的逆矩阵，这样做并不会改变问题的解
     */
    Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
    Eigen::VectorXd S_inv = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));

    // S = S^0.5 * S^0.5
    Eigen::VectorXd S_sqrt = S.cwiseSqrt(); // 这个求得就是 S^(1/2)，不过这里是向量还不是矩阵
    Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();
    // 边缘化为了实现对剩下参数块的约束，为了便于一起优化，就抽象成了残差和雅克比的形式，这样也形成了一种残差约束
    /**
     * 假设参数块的总维度为n，对于维度为2的残差，构成的jacobian为2 * n；对于维度为3的残差，构成的jacobian为3 * n
     * 但是在边缘化之后，剩余的参数块的维度为m，那么我们得到的jacobian为m * m，残差是m * 1
     * 注意：这里的jacobian和残差可以认为是虚构的，不是由物理世界实际存在的约束构建的
     */
    linearized_jacobians = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
    linearized_residuals = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose() * b;
    //std::cout << A << std::endl
    //          << std::endl;
    //std::cout << linearized_jacobians << std::endl;
    //printf("error2: %f %f\n", (linearized_jacobians.transpose() * linearized_jacobians - A).sum(),
    //      (linearized_jacobians.transpose() * linearized_residuals - b).sum());
}

std::vector<double *> MarginalizationInfo::getParameterBlocks(std::unordered_map<long, double *> &addr_shift)
{
    std::vector<double *> keep_block_addr;
    keep_block_size.clear();
    keep_block_idx.clear();
    keep_block_data.clear();

    for (const auto &it : parameter_block_idx)  // 遍历边缘化相关的每个参数块
    {
        if (it.second >= m) // 如果是留下来的，说明后续会对其形成约束
        {
            keep_block_size.push_back(parameter_block_size[it.first]);  // 留下来的参数块大小 global size
            keep_block_idx.push_back(parameter_block_idx[it.first]);    // 留下来的在原向量中的排序
            keep_block_data.push_back(parameter_block_data[it.first]);  // 边缘化前各个参数块的值的备份
            keep_block_addr.push_back(addr_shift[it.first]); // 对应的新地址
        }
    }
    // 留下来的边缘化后的参数块总大小
    sum_block_size = std::accumulate(std::begin(keep_block_size), std::end(keep_block_size), 0);

    return keep_block_addr;
}
/**
 * @brief Construct a new Marginalization Factor:: Marginalization Factor object
 * 边缘化信息结果的构造函数，根据边缘化信息确定参数块总数和大小以及残差维数
 * 
 * @param[in] _marginalization_info 
 */
MarginalizationFactor::MarginalizationFactor(MarginalizationInfo* _marginalization_info):marginalization_info(_marginalization_info)
{
    int cnt = 0;
    for (auto it : marginalization_info->keep_block_size)   // keep_block_size表示上一次边缘化留下来的参数块的大小
    {
        mutable_parameter_block_sizes()->push_back(it); // 调用ceres接口，添加参数块大小信息
        cnt += it;
    }
    //printf("residual size: %d, %d\n", cnt, n);
    set_num_residuals(marginalization_info->n); // 残差维数就是所有剩余状态量的维数和，这里时local size
};

/**
 * @brief 边缘化结果残差和雅克比的计算
 * 
 * @param[in] parameters 
 * @param[in] residuals 
 * @param[in] jacobians 
 * @return true 
 * @return false 
 */
bool MarginalizationFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    //printf("internal addr,%d, %d\n", (int)parameter_block_sizes().size(), num_residuals());
    //for (int i = 0; i < static_cast<int>(keep_block_size.size()); i++)
    //{
    //    //printf("unsigned %x\n", reinterpret_cast<unsigned long>(parameters[i]));
    //    //printf("signed %x\n", reinterpret_cast<long>(parameters[i]));
    //printf("jacobian %x\n", reinterpret_cast<long>(jacobians));
    //printf("residual %x\n", reinterpret_cast<long>(residuals));
    //}
    int n = marginalization_info->n;    // 上一次边缘化保留的残差块的local size的和,也就是残差维数
    int m = marginalization_info->m;    // 上次边缘化的被margin的残差块的local size的总和
    Eigen::VectorXd dx(n);  // 用来存储残差
    // 遍历所有的剩下的有约束的残差块
    for (int i = 0; i < static_cast<int>(marginalization_info->keep_block_size.size()); i++)
    {
        int size = marginalization_info->keep_block_size[i];
        int idx = marginalization_info->keep_block_idx[i] - m;  // idx起点统一到0
        Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(parameters[i], size); // 当前参数块的值，在优化过程中动态变化
        Eigen::VectorXd x0 = Eigen::Map<const Eigen::VectorXd>(marginalization_info->keep_block_data[i], size); // 原始参数块的值
        if (size != 7)
            dx.segment(idx, size) = x - x0; // 不需要local param的直接做差
        else    // 代表位姿的param
        {
            dx.segment<3>(idx + 0) = x.head<3>() - x0.head<3>();    // 位移直接做差
            // 旋转就是李代数做差：q' = q * delta_q → delta_q = q.inv * q'，其中q'对应着x，是一个优化过程中的动态值，q对应着x0
            dx.segment<3>(idx + 3) = 2.0 * Utility::positify(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            /**
             * 确保实部大于0：优化过程都是小步迭代，因此相对旋转对应的轴角很小，因此cos(theta/2)接近1
             *
             * q1 = (s, v).t, q2 = (-1, 0).t, q3 = (-s, -v).t，q1与q3实际上表示同一个旋转
             * q1 * q2 = q3, q1.inv * q1 = (1, 0).t, q1.inv * q3 = q2，尽管q1与q3表示同一个旋转，但是q1.inv * q3并不是单位四元数
             * 以下操作是为了避免出现(-1, 0).t的情形
             */
            if (!((Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).w() >= 0))
            {
                dx.segment<3>(idx + 3) = 2.0 * -Utility::positify(Eigen::Quaterniond(x0(6), x0(3), x0(4), x0(5)).inverse() * Eigen::Quaterniond(x(6), x(3), x(4), x(5))).vec();
            }
        }
    }
    // 更新残差　边缘化后的先验误差 e = e0 + J * dx
    // 个人理解：根据FEJ．雅克比保持不变，但是残差随着优化会变化，因此下面不更新雅克比　只更新残差
    // 关于FEJ的理解可以参考 https://www.zhihu.com/question/52869487
    // 可以参考　https://blog.csdn.net/weixin_41394379/article/details/89975386
#if 0
    std::cout << "linearized_residuals.size = " << marginalization_info->linearized_residuals.size() << std::endl;
    std::cout << "linearized_jacobians.size = " << marginalization_info->linearized_jacobians.size() << std::endl;
    std::cout << "marginalization_info->linearized_jacobians * dx = " << marginalization_info->linearized_jacobians * dx << std::endl;
#endif
#if 0
    std::cout << "m = " << m << ", n = " << n << std::endl;
#endif

    /**
     * 保持jacobian不变，仅仅根据状态量的变动更新残差，也就是说，预积分形成的约束(残差)只会在构成边缘化的时候计算得到的状态量处展开，之后将不会再发生变化
     *
     * 需要注意的是，在这里仅仅保持了边缘化形成的约束对相关状态量的first estimate jacobian，但是这些状态量还可能形成IMU预积分约束和视觉重投影约束，而在这些约束中并没有固定住jacobian
     * 因此，vins-mono实际上并没有严格执行FEJ，对边缘化形成的约束执行FEJ，对其余的约束不做特殊处理；从结果上看，不用FEJ也可以运行的很好
     *
     * 在https://heyijia.blog.csdn.net/article/details/52822104?spm=1001.2014.3001.5502中，engel提供的图像说明，在不同的点处展开，可能会引入人为的伪信息，导致不可观测的状态变的可以观测了
     *
     * 在高翔博士的博客https://zhuanlan.zhihu.com/p/29177540中对零空间和FEJ的描述：
     *          在SLAM中，GN或LM优化过程中，状态变量事实上是存在零空间（nullspace）的。所谓零空间，就是说，如果在这个子空间内改变状态变量的值，不会影响到优化的目标函数取值。
     *      在纯视觉SLAM中，典型的零空间就是场景的绝对尺度，以及相机的绝对位姿。可以想象，如果将相机和地图同时扩大一个尺度，或者同时作一次旋转和平移，整个目标函数应该不发生改变。
     *      零空间的存在也会导致在优化过程中，线性方程的解不唯一。尽管每次迭代的步长都很小，但是时间长了，也容易让整个系统在这个零空间中游荡，令尺度发生漂移。
     *          另一方面，由于边缘化的存在，如果一个点被边缘化，那它的雅可比等矩阵就要被固定于被边缘化之时刻，而其他点的雅可比还会随着时间更新，就使得整个系统中不同变量被线性化的时刻是不同的。
     *      这会导致零空间的降维——本应存在的零空间，由于线性化时刻的不同，反而消失了！而我们的估计也会变得过于乐观。
     *          实际上，优化问题会归结到Hx = b的求解上，如果H有零空间，那么x的变化将不会影响系统的残差，也就是对系统不会有影响，举例来说就是单目视觉的尺度；如果Hy = 0，那么x = x + y不会改变系统
     *
     * 残差更新方式：状态量x的变化引起b的变化，b的变化引起e的变化
     *      b = J.t * e → b' = b + Jb_x * dx
     *      Jb_x = Jb_e * Je_x = J.t * J, b' = J.t * e', b = J.t * e
     *      J.t * e' = J.t * e + J.t * J * dx → e' = e + J * dx
     */
    Eigen::Map<Eigen::VectorXd>(residuals, n) = marginalization_info->linearized_residuals + marginalization_info->linearized_jacobians * dx;
    if (jacobians)
    {
        for (int i = 0; i < static_cast<int>(marginalization_info->keep_block_size.size()); i++)
        {
            if (jacobians[i])
            {
                // 将整个大的jacobian即marginalization_info->linearized_jacobians分配给每一个参数块
                int size = marginalization_info->keep_block_size[i], local_size = marginalization_info->localSize(size);
                int idx = marginalization_info->keep_block_idx[i] - m;
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> jacobian(jacobians[i], n, size);
                jacobian.setZero();
                jacobian.leftCols(local_size) = marginalization_info->linearized_jacobians.middleCols(idx, local_size);
            }
        }
    }
    return true;
}
