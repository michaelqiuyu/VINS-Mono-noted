#include "camodocal/camera_models/Camera.h"
#include "camodocal/camera_models/ScaramuzzaCamera.h"

#include <opencv2/calib3d/calib3d.hpp>

namespace camodocal
{

Camera::Parameters::Parameters(ModelType modelType)
 : m_modelType(modelType)
 , m_imageWidth(0)
 , m_imageHeight(0)
{
    switch (modelType)
    {
    case KANNALA_BRANDT:
        m_nIntrinsics = 8;
        break;
    case PINHOLE:
        m_nIntrinsics = 8;
        break;
    case SCARAMUZZA:
        m_nIntrinsics = SCARAMUZZA_CAMERA_NUM_PARAMS;
        break;
    case MEI:
    default:
        m_nIntrinsics = 9;
    }
}

Camera::Parameters::Parameters(ModelType modelType,
                               const std::string& cameraName,
                               int w, int h)
 : m_modelType(modelType)
 , m_cameraName(cameraName)
 , m_imageWidth(w)
 , m_imageHeight(h)
{
    switch (modelType)
    {
    case KANNALA_BRANDT:
        m_nIntrinsics = 8;
        break;
    case PINHOLE:
        m_nIntrinsics = 8;
        break;
    case SCARAMUZZA:
        m_nIntrinsics = SCARAMUZZA_CAMERA_NUM_PARAMS;
        break;
    case MEI:
    default:
        m_nIntrinsics = 9;
    }
}

/**
 * 返回引用，函数可以作为左值修改对应的变量；
 * 当返回一个引用时，要注意被引用的对象不能超出作用域。所以返回一个对局部变量的引用是不合法的，但是，可以返回一个对静态变量的引用。
 */
Camera::ModelType&
Camera::Parameters::modelType(void)
{
    return m_modelType;
}

std::string&
Camera::Parameters::cameraName(void)
{
    return m_cameraName;
}

int&
Camera::Parameters::imageWidth(void)
{
    return m_imageWidth;
}

int&
Camera::Parameters::imageHeight(void)
{
    return m_imageHeight;
}

Camera::ModelType
Camera::Parameters::modelType(void) const
{
    return m_modelType;
}

const std::string&
Camera::Parameters::cameraName(void) const
{
    return m_cameraName;
}

int
Camera::Parameters::imageWidth(void) const
{
    return m_imageWidth;
}

int
Camera::Parameters::imageHeight(void) const
{
    return m_imageHeight;
}

int
Camera::Parameters::nIntrinsics(void) const
{
    return m_nIntrinsics;
}

cv::Mat&
Camera::mask(void)
{
    return m_mask;
}

const cv::Mat&
Camera::mask(void) const
{
    return m_mask;
}

// 根据2D-3D的匹配关系，调用PnP求解位姿
void
Camera::estimateExtrinsics(const std::vector<cv::Point3f>& objectPoints,
                           const std::vector<cv::Point2f>& imagePoints,
                           cv::Mat& rvec, cv::Mat& tvec) const
{
    std::vector<cv::Point2f> Ms(imagePoints.size());
    for (size_t i = 0; i < Ms.size(); ++i)
    {
        Eigen::Vector3d P;
        // 获取去畸变的相机归一化坐标系
        liftProjective(Eigen::Vector2d(imagePoints.at(i).x, imagePoints.at(i).y), P);

        // 获得相机归一化坐标系坐标
        P /= P(2);

        Ms.at(i).x = P(0);
        Ms.at(i).y = P(1);
    }

    // 在这里已经将图像坐标转化到去畸变后的相机归一化坐标系下了，因此cameraMatrix就是单位阵了
    // assume unit focal length, zero principal point, and zero distortion
    cv::solvePnP(objectPoints, Ms, cv::Mat::eye(3, 3, CV_64F), cv::noArray(), rvec, tvec);
}

// 两个世界点在同一个相机下投影到图像的像素坐标的误差
double
Camera::reprojectionDist(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2) const
{
    Eigen::Vector2d p1, p2;

    // 将三维点投影到图像坐标系
    spaceToPlane(P1, p1);
    spaceToPlane(P2, p2);

    return (p1 - p2).norm();
}

// 标定调用
double
Camera::reprojectionError(const std::vector< std::vector<cv::Point3f> >& objectPoints,
                          const std::vector< std::vector<cv::Point2f> >& imagePoints,
                          const std::vector<cv::Mat>& rvecs,
                          const std::vector<cv::Mat>& tvecs,
                          cv::OutputArray _perViewErrors) const
{
    int imageCount = objectPoints.size();
    size_t pointsSoFar = 0;
    double totalErr = 0.0;

    bool computePerViewErrors = _perViewErrors.needed();
    cv::Mat perViewErrors;
    if (computePerViewErrors)
    {
        _perViewErrors.create(imageCount, 1, CV_64F);
        perViewErrors = _perViewErrors.getMat();
    }

    for (int i = 0; i < imageCount; ++i)
    {
        size_t pointCount = imagePoints.at(i).size();

        pointsSoFar += pointCount;

        std::vector<cv::Point2f> estImagePoints;
        projectPoints(objectPoints.at(i), rvecs.at(i), tvecs.at(i),
                      estImagePoints);

        double err = 0.0;
        for (size_t j = 0; j < imagePoints.at(i).size(); ++j)
        {
            err += cv::norm(imagePoints.at(i).at(j) - estImagePoints.at(j));
        }

        if (computePerViewErrors)
        {
            perViewErrors.at<double>(i) = err / pointCount;
        }

        totalErr += err;
    }

    return totalErr / pointsSoFar;
}

// 计算世界系坐标投影到图像的重投影误差
double
Camera::reprojectionError(const Eigen::Vector3d& P,
                          const Eigen::Quaterniond& camera_q,
                          const Eigen::Vector3d& camera_t,
                          const Eigen::Vector2d& observed_p) const
{
    Eigen::Vector3d P_cam = camera_q.toRotationMatrix() * P + camera_t;

    Eigen::Vector2d p;
    spaceToPlane(P_cam, p);

    return (p - observed_p).norm();
}

// 根据已知的位姿，将世界系坐标投影到图像坐标
void
Camera::projectPoints(const std::vector<cv::Point3f>& objectPoints,
                      const cv::Mat& rvec,
                      const cv::Mat& tvec,
                      std::vector<cv::Point2f>& imagePoints) const
{
    // project 3D object points to the image plane
    imagePoints.reserve(objectPoints.size());

    //double
    cv::Mat R0;
    cv::Rodrigues(rvec, R0);

    Eigen::MatrixXd R(3,3);
    R << R0.at<double>(0,0), R0.at<double>(0,1), R0.at<double>(0,2),
         R0.at<double>(1,0), R0.at<double>(1,1), R0.at<double>(1,2),
         R0.at<double>(2,0), R0.at<double>(2,1), R0.at<double>(2,2);

    Eigen::Vector3d t;
    t << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);

    for (size_t i = 0; i < objectPoints.size(); ++i)
    {
        const cv::Point3f& objectPoint = objectPoints.at(i);

        // Rotate and translate
        Eigen::Vector3d P;
        P << objectPoint.x, objectPoint.y, objectPoint.z;

        P = R * P + t;

        Eigen::Vector2d p;
        spaceToPlane(P, p);

        imagePoints.push_back(cv::Point2f(p(0), p(1)));
    }
}

}
