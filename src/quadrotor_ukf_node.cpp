#include <iostream>

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TransformStamped.h>
#include <tf2_ros/transform_listener.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <tf2_eigen/tf2_eigen.h>
#include <angles/angles.h>

#include <tf2_ros/transform_broadcaster.h>
#include <quadrotor_ukf/quadrotor_ukf.h>
#include <quadrotor_ukf/vio_utils.h>

class QuadrotorUkfNode
{
public:
  QuadrotorUkfNode(std::string ns = "");
  void init();

private:

  void imu_callback(const sensor_msgs::Imu::ConstPtr& msg);
  void slam_callback(const nav_msgs::Odometry::ConstPtr& msg);
  void pose_callback(const geometry_msgs::PoseStamped::ConstPtr& pose);

  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;
  ros::Subscriber subImu_, subSLAM_, subPose_;
  ros::Publisher pubUKF_;
  ros::Publisher pubMavUKF_;
  ros::Publisher pubBias_;

  tf2_ros::Buffer tfBuffer_;
  std::unique_ptr<tf2_ros::TransformListener> tfListener_;
  tf2_ros::TransformBroadcaster tf_broadcaster;

  QuadrotorUKF quadrotorUKF_;


  std::string odom_frame_id_;

  //arma::mat H_C_B = arma::eye<mat>(4,4);//Never use reshape
  Eigen::Matrix<double, 4, 4> H_C_B_;
  Eigen::Matrix<double, 4, 4> H_I_B_;

  int calLimit_;
  int calCnt_;
  Eigen::Matrix<double, 3, 1> average_g_;

  std::string imu_frame_id_, imu_rotated_base_frame_id_, body_frame_id_, body_local_frame_id_;
  bool tf_initialized_;
  geometry_msgs::TransformStamped tf_imu_to_base_;

};

QuadrotorUkfNode::QuadrotorUkfNode(std::string ns): nh_(ns), pnh_("~"), tf_initialized_(false)
{
  tfListener_.reset(new tf2_ros::TransformListener(tfBuffer_));

  H_I_B_.setZero();
  average_g_.setZero();

  H_C_B_(0,0) = 1;
  H_C_B_(0,1) = 0;
  H_C_B_(0,2) = 0;
  H_C_B_(0,3) = 0;
  H_C_B_(1,0) = 0;
  H_C_B_(1,1) = -1;
  H_C_B_(1,2) = 0;
  H_C_B_(1,3) = 0;
  H_C_B_(2,0) = 0;
  H_C_B_(2,1) = 0;
  H_C_B_(2,2) = -1;
  H_C_B_(2,3) = 0.0;
  H_C_B_(3,0) = 0;
  H_C_B_(3,1) = 0;
  H_C_B_(3,2) = 0;
  H_C_B_(3,3) = 1;

  ROS_DEBUG_STREAM("H_C_B\n" << H_C_B_);

  calLimit_ = 100;
  calCnt_   = 0;

  // UKF Parameters and Noise
  double alpha;
  double beta;
  double kappa;
  double stdAcc[3]     = {0,0,0};
  double stdW[3]       = {0,0,0};
  double stdAccBias[3] = {0,0,0};

  pnh_.param("odom", odom_frame_id_, std::string("odom"));
  pnh_.param("imu", imu_frame_id_, std::string("imu"));
  pnh_.param("base_link", body_frame_id_, std::string("base_link"));
  pnh_.param("base_link_frd", body_local_frame_id_, std::string("base_link_frd"));
  pnh_.param("imu_rotated_frame_id", imu_rotated_base_frame_id_, std::string("imu_rotated_base"));
  pnh_.param("alpha", alpha, 0.4);
  pnh_.param("beta" , beta , 2.0);
  pnh_.param("kappa", kappa, 0.0);
  pnh_.param("noise_std/process/acc/x", stdAcc[0], 0.2);
  pnh_.param("noise_std/process/acc/y", stdAcc[1], 0.2);
  pnh_.param("noise_std/process/acc/z", stdAcc[2], 0.2);
  pnh_.param("noise_std/process/w/x", stdW[0], 0.1);
  pnh_.param("noise_std/process/w/y", stdW[1], 0.1);
  pnh_.param("noise_std/process/w/z", stdW[2], 0.1);
  pnh_.param("noise_std/process/acc_bias/x", stdAccBias[0], 0.05);
  pnh_.param("noise_std/process/acc_bias/y", stdAccBias[1], 0.05);
  pnh_.param("noise_std/process/acc_bias/z", stdAccBias[2], 0.05);

  // Fixed process noise
  Eigen::Matrix<double, 9, 9> Rv;
  Rv.setIdentity();// = eye<mat>(9,9);
  Rv(0,0)   = stdAcc[0] * stdAcc[0];
  Rv(1,1)   = stdAcc[1] * stdAcc[1];
  Rv(2,2)   = stdAcc[2] * stdAcc[2];
  Rv(3,3)   = stdW[0] * stdW[0];
  Rv(4,4)   = stdW[1] * stdW[1];
  Rv(5,5)   = stdW[2] * stdW[2];
  Rv(6,6)   = stdAccBias[0] * stdAccBias[0];
  Rv(7,7)   = stdAccBias[1] * stdAccBias[1];
  Rv(8,8)   = stdAccBias[2] * stdAccBias[2];

  // Initialize UKF
  quadrotorUKF_.SetUKFParameters(alpha, beta, kappa);
  quadrotorUKF_.SetImuCovariance(Rv);

  subImu_  = pnh_.subscribe<sensor_msgs::Imu>("imu", 10, boost::bind(&QuadrotorUkfNode::imu_callback, this,_1));
  subSLAM_ = pnh_.subscribe<nav_msgs::Odometry>("odom_slam", 10, boost::bind(&QuadrotorUkfNode::slam_callback, this, _1));
  subPose_ = pnh_.subscribe<geometry_msgs::PoseStamped>("/qvio/pose", 1, boost::bind(&QuadrotorUkfNode::pose_callback, this, _1));

  pubUKF_  = pnh_.advertise<nav_msgs::Odometry>("control_odom", 10);
  pubMavUKF_  = pnh_.advertise<nav_msgs::Odometry>("/mavros/odometry/in", 10);
  pubBias_ = pnh_.advertise<geometry_msgs::Vector3>("imu_bias", 10);
  quadrotorUKF_.PrintxHist();
}

void QuadrotorUkfNode::init()
{
  ros::Rate rate(2.0);
  while (nh_.ok() && !tf_initialized_)
  {
    try
    {
      tf_imu_to_base_ = tfBuffer_.lookupTransform(body_frame_id_, imu_frame_id_, ros::Time(0));

      //Get the rotation matrix and compose Homogenous matrix (without translation)
      Eigen::Affine3d R_I_B = tf2::transformToEigen(tf_imu_to_base_);
      H_I_B_.block(0,0,3,3) = R_I_B.rotation();
      H_I_B_(3,3) = 1;

/*
  H_I_B_(0,0) = 0;
  H_I_B_(0,1) = 1;
  H_I_B_(0,2) = 0;
  H_I_B_(0,3) = 0;

  H_I_B_(1,0) = -1;
  H_I_B_(1,1) = 0;
  H_I_B_(1,2) = 0;
  H_I_B_(1,3) = 0;

  H_I_B_(2,0) = 0;
  H_I_B_(2,1) = 0;
  H_I_B_(2,2) = 1;
  H_I_B_(2,3) = 0.0;

  H_I_B_(3,0) = 0;
  H_I_B_(3,1) = 0;
  H_I_B_(3,2) = 0;
  H_I_B_(3,3) = 1;
  */

      ROS_DEBUG_STREAM("Got imu to imu_rotated_base tf " << tf_imu_to_base_);
      ROS_DEBUG_STREAM("H_I_B\n" << H_I_B_);

      tf_initialized_ = true;
    }
    catch (tf2::TransformException &ex)
    {
      ROS_WARN_THROTTLE(1, "Failed to find transform from [%s] to [%s]",
                        imu_frame_id_.c_str(), imu_rotated_base_frame_id_.c_str());
    }
  }

  tfListener_.reset();
}

void QuadrotorUkfNode::imu_callback(const sensor_msgs::Imu::ConstPtr& msg)
{

  if (isnan(msg->orientation.x) || isnan(msg->orientation.y) || isnan(msg->orientation.z) || isnan(msg->orientation.w))
  {
    ROS_WARN("One or more orientation fields are NaN!");
  }

  // Check if the angular_velocity fields are NaN
  if (isnan(msg->angular_velocity.x) || isnan(msg->angular_velocity.y) || isnan(msg->angular_velocity.z))
  {
    ROS_WARN("One or more angular_velocity fields are NaN!");
  }

  // Check if the linear_acceleration fields are NaN
  if (isnan(msg->linear_acceleration.x) || isnan(msg->linear_acceleration.y) || isnan(msg->linear_acceleration.z))
  {
    ROS_WARN("One or more linear_acceleration fields are NaN!");
  }

  if(!tf_initialized_)
    return;

  //Transform imu into base frame
  geometry_msgs::Vector3 linear_acceleration_rotated;
  geometry_msgs::Vector3 angular_velocity_rotated;

  //ROS kinetic and above
  tf2::doTransform(msg->linear_acceleration, linear_acceleration_rotated, tf_imu_to_base_);
  tf2::doTransform(msg->angular_velocity, angular_velocity_rotated, tf_imu_to_base_);

  //ROS_WARN_STREAM("Orig lin acc " << msg->linear_acceleration << " rotated " << linear_acceleration_rotated);
  //ROS_WARN_STREAM("Orig ang vel " << msg->angular_velocity << " rotated " << angular_velocity_rotated);

  // Assemble control input, and calibration
  Eigen::Matrix<double, 6, 1> u;
  u(0,0) = linear_acceleration_rotated.x;
  u(1,0) = linear_acceleration_rotated.y;
  u(2,0) = linear_acceleration_rotated.z;
  u(3,0) = angular_velocity_rotated.x;
  u(4,0) = angular_velocity_rotated.y;
  u(5,0) = angular_velocity_rotated.z;

  for (int i = 0; i < u.rows(); ++i) {
    if (std::isnan(u(i, 0))) {
      ROS_DEBUG_STREAM("U Matrix Contains NAN!");
    }
  }

  Eigen::Matrix<double, 6, 1> u_old;
  u_old(0,0) = msg->linear_acceleration.x;
  u_old(1,0) = -msg->linear_acceleration.y;
  u_old(2,0) = -msg->linear_acceleration.z;
  u_old(3,0) = msg->angular_velocity.x;
  u_old(4,0) = -msg->angular_velocity.y;
  u_old(5,0) = -msg->angular_velocity.z;

  //ROS_WARN_STREAM("u\n" << u << "\nu old\n " << u_old);

  if (calCnt_ < calLimit_)       // Calibration
  {
    calCnt_++;
    average_g_ += u.block(0,0,3,1);//rows(0,2);
  }
  else if (calCnt_ == calLimit_) // Save gravity vector
  {
    calCnt_++;
    average_g_ /= calLimit_;
    double g = average_g_.norm(); //norm(ag,2);
    quadrotorUKF_.SetGravity(g);
    ROS_DEBUG_STREAM("Set Gravity Calibration!");
  }
  else if (quadrotorUKF_.ProcessUpdate(u, msg->header.stamp))  // Process Update
  {
    ROS_DEBUG_STREAM("Printing xHist from imu_callback...");
    quadrotorUKF_.PrintxHist();
    nav_msgs::Odometry odomUKF;

    // Publish odom
    odomUKF.header.stamp = quadrotorUKF_.GetStateTime();
    odomUKF.header.frame_id = odom_frame_id_;
    odomUKF.child_frame_id = "base_link";
    Eigen::Matrix<double, Eigen::Dynamic, 1> x = quadrotorUKF_.GetState();
    for (int i = 0; i < x.rows(); ++i) {
      if (std::isnan(x(i, 0))) {
        ROS_DEBUG_STREAM("x Matrix Contains NAN!");
        ROS_DEBUG_STREAM("Failure from imu_callback ");
        ros::shutdown();
      }
    }
    odomUKF.pose.pose.position.x = x(0,0);
    odomUKF.pose.pose.position.y = x(1,0);
    odomUKF.pose.pose.position.z = x(2,0);
    ROS_DEBUG_STREAM("Running ypr_to_R from node...");
    Eigen::Matrix<double, 4, 1> q = VIOUtil::MatToQuat(VIOUtil::ypr_to_R(x.block(6,0,3,1)));
    for (int i = 0; i < q.rows(); ++i) {
      if (std::isnan(q(i, 0))) {
        ROS_DEBUG("q Matrix Contains NAN!");
      }
    }
    odomUKF.pose.pose.orientation.w = q(0,0);
    odomUKF.pose.pose.orientation.x = q(1,0);
    odomUKF.pose.pose.orientation.y = q(2,0);
    odomUKF.pose.pose.orientation.z = q(3,0);
    odomUKF.twist.twist.linear.x = x(3,0);
    odomUKF.twist.twist.linear.y = x(4,0);
    odomUKF.twist.twist.linear.z = x(5,0);
    odomUKF.twist.twist.angular.x = u(3,0);
    odomUKF.twist.twist.angular.y = u(4,0);
    odomUKF.twist.twist.angular.z = u(5,0);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  P = quadrotorUKF_.GetStateCovariance();
    for (int j = 0; j < 6; j++)
      for (int i = 0; i < 6; i++)
        odomUKF.pose.covariance[i+j*6] = P((i<3)?i:i+3 , (j<3)?j:j+3);
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        odomUKF.twist.covariance[i+j*6] = P(i+3 , j+3);
    nav_msgs::Odometry odomMavUKF(odomUKF);


    geometry_msgs::TransformStamped static_transformStamped;
    static_transformStamped.header.stamp = ros::Time::now();
    static_transformStamped.header.frame_id = "odom";
    static_transformStamped.child_frame_id = "base_link";

    static_transformStamped.transform.translation.x = odomUKF.pose.pose.position.x;
    static_transformStamped.transform.translation.y = odomUKF.pose.pose.position.y;
    static_transformStamped.transform.translation.z = odomUKF.pose.pose.position.z;

    static_transformStamped.transform.rotation = odomUKF.pose.pose.orientation;

    geometry_msgs::TransformStamped base_link_flat(static_transformStamped); //no orientation
    tf_broadcaster.sendTransform(static_transformStamped);

    base_link_flat.child_frame_id = "base_link_flat";
    base_link_flat.transform.rotation.x = 0;
    base_link_flat.transform.rotation.y = 0;
    base_link_flat.transform.rotation.z = 0;
    base_link_flat.transform.rotation.w = 1;
    //tf_broadcaster.sendTransform(base_link_flat);

    pubUKF_.publish(odomUKF);

    // Publish bias
    geometry_msgs::Vector3 bias;
    bias.x = x(9);
    bias.y = x(10);
    bias.z = x(11);
    pubBias_.publish(bias);
    
  }
}

void QuadrotorUkfNode::slam_callback(const nav_msgs::Odometry::ConstPtr& msg)
{
  if (isnan(msg->pose.pose.position.x) || isnan(msg->pose.pose.position.y) || isnan(msg->pose.pose.position.z) ||
      isnan(msg->pose.pose.orientation.x) || isnan(msg->pose.pose.orientation.y) || isnan(msg->pose.pose.orientation.z) || 
      isnan(msg->pose.pose.orientation.w))
  {
    ROS_WARN("One or more pose fields are NaN!");
  }

  // Check if the twist fields are NaN
  if (isnan(msg->twist.twist.linear.x) || isnan(msg->twist.twist.linear.y) || isnan(msg->twist.twist.linear.z) ||
      isnan(msg->twist.twist.angular.x) || isnan(msg->twist.twist.angular.y) || isnan(msg->twist.twist.angular.z))
  {
    ROS_WARN("One or more twist fields are NaN!");
  }
  if(!tf_initialized_)
    return;

  // Get orientation
  Eigen::Matrix<double, 4, 1> q;
  q(0,0) = msg->pose.pose.orientation.w;
  q(1,0) = msg->pose.pose.orientation.x;
  q(2,0) = msg->pose.pose.orientation.y;
  q(3,0) = msg->pose.pose.orientation.z;
  Eigen::Matrix<double, 3, 1> ypr = VIOUtil::R_to_ypr(VIOUtil::QuatToMat(q));

  // Assemble measurement
  Eigen::Matrix<double, 6, 1> z;
  z(0,0) = msg->pose.pose.position.x;
  z(1,0) = msg->pose.pose.position.y;
  z(2,0) = msg->pose.pose.position.z;
  z(3,0) = ypr(0,0);
  z(4,0) = ypr(1,0);
  z(5,0) = ypr(2,0);

  bool z_nan = false;
  for (int i = 0; i < z.rows(); ++i) {
    if (std::isnan(z(i, 0))) {
      z_nan = true;
    }
  }
  if (z_nan) {
    ROS_DEBUG_STREAM("z Matrix Contains NAN!");
    ROS_DEBUG_STREAM("Shutting Down from within slam_callback");
    ROS_DEBUG_STREAM(z.matrix());
    ros::shutdown();
  }
  // Assemble measurement covariance
  Eigen::Matrix<double, 6, 6> RnSLAM;
  RnSLAM.setZero();
  RnSLAM(0,0) = msg->pose.covariance[0];
  RnSLAM(1,1) = msg->pose.covariance[1+1*6];
  RnSLAM(2,2) = msg->pose.covariance[2+2*6];
  RnSLAM(3,3) = msg->pose.covariance[3+3*6];
  RnSLAM(4,4) = msg->pose.covariance[4+4*6];
  RnSLAM(5,5) = msg->pose.covariance[5+5*6];

  //rotate the measurement for control purpose
  Eigen::Matrix<double, 4, 4> H_C_C0;
  H_C_C0.setIdentity();
  H_C_C0.block(0, 0, 3, 3) = VIOUtil::QuatToMat(q);
  H_C_C0(0,3) = z(0,0);
  H_C_C0(1,3) = z(1,0);
  H_C_C0(2,3) = z(2,0);

  //robot frame
  Eigen::Matrix<double, 4, 4> H_R_R0;
  H_R_R0.setIdentity();
  H_R_R0 = H_I_B_*H_C_C0*H_I_B_.inverse();

  //Set the rotation
  Eigen::Matrix<double, 4, 1> q_R_R0 = VIOUtil::MatToQuat(H_R_R0.block(0, 0, 3, 3));

  // Assemble measurement
  Eigen::Matrix<double, 6, 1> z_new;
  z_new(0,0) = H_R_R0(0,3);
  z_new(1,0) = H_R_R0(1,3);
  z_new(2,0) = H_R_R0(2,3);

  //define the matrix to rotate in the original frame
  Eigen::Matrix<double, 3, 1> ypr_new = VIOUtil::R_to_ypr(VIOUtil::QuatToMat(q_R_R0));
  z_new(3,0) = ypr_new(0,0);
  z_new(4,0) = ypr_new(1,0);
  z_new(5,0) = ypr_new(2,0);

  //rotate the covariance
  Eigen::Matrix<double, 6, 6>  RnSLAM_new;
  RnSLAM_new.setZero();
  RnSLAM_new.block(0,0,3,3) = H_I_B_.block(0, 0, 3, 3)*RnSLAM.block(0,0,3,3)*H_I_B_.block(0, 0, 3, 3).transpose();
  RnSLAM_new.block(3,3,3,3) = H_I_B_.block(0, 0, 3, 3)*RnSLAM.block(3,3,3,3)*H_I_B_.block(0, 0, 3, 3).transpose();

  // Measurement update
  if (quadrotorUKF_.isInitialized())
  {
    quadrotorUKF_.MeasurementUpdateSLAM(z_new, RnSLAM_new, msg->header.stamp);
    ROS_DEBUG_STREAM("Printing xHist from slam_callback...");
    quadrotorUKF_.PrintxHist();
  }
  else
  {
    quadrotorUKF_.SetInitPose(z, msg->header.stamp);
  }

  //-------------------------

  /*
  //Transform odom (in imu_init frame) into base frame
  geometry_msgs::PoseStamped ps_in_imu;
  ps_in_imu.header.frame_id = imu_frame_id_;
  ps_in_imu.pose = msg->pose.pose;
  geometry_msgs::PoseStamped ps_in_imu_rotated;
  tf2::doTransform(ps_in_imu, ps_in_imu_rotated, tf_imu_to_base_);

  tf2::Quaternion q_rot; //Quaternion in rotated base frame
  tf2::convert(ps_in_imu_rotated.pose.orientation , q_rot);
  tf2::Matrix3x3 R_rot(q_rot); //Rotation matrix
  double y, p, r;
  R_rot.getEulerYPR(y, p, r);
  Eigen::Matrix<double, 3, 1> ypr_rot;
  ypr_rot(0,0) = y;
  ypr_rot(1,0) = p;
  ypr_rot(2,0) = r;

  Eigen::Matrix<double, 4, 1> q_rotated_base;
  q_rotated_base(0,0) = ps_in_imu_rotated.pose.orientation.w;
  q_rotated_base(1,0) = ps_in_imu_rotated.pose.orientation.x;
  q_rotated_base(2,0) = ps_in_imu_rotated.pose.orientation.y;
  q_rotated_base(3,0) = ps_in_imu_rotated.pose.orientation.z;
  Eigen::Matrix<double, 3, 1> ypr_rotated = VIOUtil::R_to_ypr(VIOUtil::QuatToMat(q_rotated_base));
  ypr_rotated(2,0) = angles::normalize_angle(ypr_rotated(2,0));

  // Assemble measurement
  Eigen::Matrix<double, 6, 1> z_rotated_base;
  z_rotated_base(0,0) = ps_in_imu_rotated.pose.position.x;
  z_rotated_base(1,0) = ps_in_imu_rotated.pose.position.y;
  z_rotated_base(2,0) = ps_in_imu_rotated.pose.position.z;
  z_rotated_base(3,0) = ypr_rotated(0,0);
  z_rotated_base(4,0) = ypr_rotated(1,0);
  z_rotated_base(5,0) = ypr_rotated(2,0);

  //rotate the covariance
  Eigen::Matrix<double, 6, 6>  RnSLAM_rot;
  RnSLAM_rot.setZero();
  RnSLAM_rot.block(0,0,3,3) = H_I_B_.block(0, 0, 3, 3)*RnSLAM.block(0,0,3,3)*H_I_B_.block(0, 0, 3, 3).transpose();
  RnSLAM_rot.block(3,3,3,3) = H_I_B_.block(0, 0, 3, 3)*RnSLAM.block(3,3,3,3)*H_I_B_.block(0, 0, 3, 3).transpose();

  ROS_ERROR_STREAM("z " << z << "\nz_new\n" << z_new << "\nz_rotated_base\n" << z_rotated_base);
  */

  //ROS_ERROR_STREAM("RnSLAM " << RnSLAM << "\nRnSLAM_new\n" << RnSLAM_new << "\nRnSLAM_rot\n" << RnSLAM_rot);

  //ROS_INFO_STREAM("H_C_B_\n" << H_C_B_ << "\nH_I_B\n" << H_I_B_);
  //ROS_WARN_STREAM("Pose " << ps_in_imu << " rotated " << ps_in_imu_rotated);

/*
  // Measurement update
  if (quadrotorUKF_.isInitialized())
  {
    quadrotorUKF_.MeasurementUpdateSLAM(z_rotated_base, RnSLAM, msg->header.stamp);
  }
  else
  {
    quadrotorUKF_.SetInitPose(z_rotated_base, msg->header.stamp);
  }

*/

  //-------------------------

  /*
  // Get orientation
  Eigen::Matrix<double, 4, 1> q;
  q(0,0) = msg->pose.pose.orientation.w;
  q(1,0) = msg->pose.pose.orientation.x;
  q(2,0) = msg->pose.pose.orientation.y;
  q(3,0) = msg->pose.pose.orientation.z;
  Eigen::Matrix<double, 3, 1> ypr = VIOUtil::R_to_ypr(VIOUtil::QuatToMat(q));

  // Assemble measurement
  Eigen::Matrix<double, 6, 1> z;
  z(0,0) = msg->pose.pose.position.x;
  z(1,0) = msg->pose.pose.position.y;
  z(2,0) = msg->pose.pose.position.z;
  z(3,0) = ypr(0,0);
  z(4,0) = ypr(1,0);
  z(5,0) = ypr(2,0);

  //rotate the measurement for control purpose
  //arma::mat H_C_C0 = arma::eye<mat>(4,4);
  Eigen::Matrix<double, 4, 4> H_C_C0;
  H_C_C0.setIdentity();
  H_C_C0.block(0, 0, 3, 3) = VIOUtil::QuatToMat(q);
  H_C_C0(0,3) = z(0,0);
  H_C_C0(1,3) = z(1,0);
  H_C_C0(2,3) = z(2,0);

  //robot frame
  //mat H_R_R0 = zeros<mat>(4,4);
  Eigen::Matrix<double, 4, 4> H_R_R0;
  H_R_R0.setIdentity();
  H_R_R0 = H_C_B_*H_C_C0*H_C_B_.inverse();

  //Set the rotation
  Eigen::Matrix<double, 4, 1> q_R_R0 = VIOUtil::MatToQuat(H_R_R0.block(0, 0, 3, 3));

  // Assemble measurement
  Eigen::Matrix<double, 6, 1> z_new;
  z_new(0,0) = H_R_R0(0,3);
  z_new(1,0) = H_R_R0(1,3);
  z_new(2,0) = H_R_R0(2,3);

  //define the matrix to rotate in the original frame
  Eigen::Matrix<double, 3, 1> ypr_new = VIOUtil::R_to_ypr(VIOUtil::QuatToMat(q_R_R0));//R_to_ypr(H_R_R0.submat(0, 0, 2, 2));
  z_new(3,0) = ypr_new(0,0);
  z_new(4,0) = ypr_new(1,0);
  z_new(5,0) = ypr_new(2,0);

  //rotate the covariance
  Eigen::Matrix<double, 6, 6>  RnSLAM_new;
  RnSLAM_new.setZero();// = zeros<mat>(6,6);
  //RnSLAM_new.submat(0,0,2,2) = H_C_B.submat(0, 0, 2, 2)*RnSLAM.submat(0,0,2,2)*H_C_B.submat(0, 0, 2, 2).t();
  //RnSLAM_new.submat(3,3,5,5) = H_C_B.submat(0, 0, 2, 2)*RnSLAM.submat(3,3,5,5)*H_C_B.submat(0, 0, 2, 2).t();
  RnSLAM_new.block(0,0,3,3) = H_C_B_.block(0, 0, 3, 3)*RnSLAM.block(0,0,3,3)*H_C_B_.block(0, 0, 3, 3).transpose();
  RnSLAM_new.block(3,3,3,3) = H_C_B_.block(0, 0, 3, 3)*RnSLAM.block(3,3,3,3)*H_C_B_.block(0, 0, 3, 3).transpose();

  //ROS_ERROR_STREAM("z " << z << "\nz_new\n" << z_new << "\nz_rotated_base\n" << z_rotated_base);

  //ROS_ERROR_STREAM("RnSLAM " << RnSLAM << "\nRnSLAM_new\n" << RnSLAM_new << "\nRnSLAM_rot\n" << RnSLAM_rot);

  //ROS_INFO_STREAM("H_C_B_\n" << H_C_B_ << "\nH_I_B\n" << H_I_B_);
  //ROS_WARN_STREAM("Pose " << ps_in_imu << " rotated " << ps_in_imu_rotated);

  Eigen::Matrix<double, 3, 1> ypr_temp;
  ypr_temp(0,0) = 0.0;
  ypr_temp(1,0) = 0.0;
  ypr_temp(2,0) = static_cast<double>(M_PI);

  Eigen::Matrix<double, 4, 1> q_temp = VIOUtil::MatToQuat(VIOUtil::ypr_to_R(ypr_temp));
  ROS_INFO_STREAM("q\n" << q_temp);

  ROS_WARN_STREAM(" ypr\n" << ypr << "\n rot ypr\n" << ypr_rotated << "\nypr_new\n" << ypr_new << "\ncheck\n" << ypr_rot);

  // Measurement update
  if (quadrotorUKF_.isInitialized())
  {
    quadrotorUKF_.MeasurementUpdateSLAM(z_new, RnSLAM_new, msg->header.stamp);
  }
  else
  {
    quadrotorUKF_.SetInitPose(z, msg->header.stamp);
  }*/

}

void QuadrotorUkfNode::pose_callback(const geometry_msgs::PoseStamped::ConstPtr& pose)
{
  geometry_msgs::TransformStamped static_transformStamped;
  static_transformStamped.header.stamp = ros::Time::now();
  static_transformStamped.header.frame_id = "map_ned";
  static_transformStamped.child_frame_id = "base_link_frd";

  static_transformStamped.transform.translation.x = pose->pose.position.x;
  static_transformStamped.transform.translation.y = pose->pose.position.y;
  static_transformStamped.transform.translation.z = pose->pose.position.z;

  static_transformStamped.transform.rotation.x = pose->pose.orientation.x;
  static_transformStamped.transform.rotation.y = pose->pose.orientation.y;
  static_transformStamped.transform.rotation.z = pose->pose.orientation.z;
  static_transformStamped.transform.rotation.w = pose->pose.orientation.w;


  tf_broadcaster.sendTransform(static_transformStamped);
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "quadrotor_ukf");
  ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Debug);

  QuadrotorUkfNode ukf_node;
  ukf_node.init();
  ros::spin();

  return 0;
}
