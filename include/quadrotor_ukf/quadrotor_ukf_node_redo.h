#include <iostream>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <quadrotor_ukf/quadrotor_ukf.h>
#include <quadrotor_ukf/vio_utils.h>

#include <geometry_msgs/PoseStamped.h>
#include <tf2_ros/static_transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>

#define FE_DEBUG(fmt, ...) ROS_DEBUG("[frontier] " fmt, ##__VA_ARGS__)
#define FE_INFO(fmt, ...) ROS_INFO("[frontier] " fmt, ##__VA_ARGS__)
#define FE_WARN(fmt, ...) ROS_WARN("[frontier] " fmt, ##__VA_ARGS__)
#define FE_ERROR(fmt, ...) ROS_ERROR("[frontier] " fmt, ##__VA_ARGS__)
#define FE_FATAL(fmt, ...) ROS_FATAL("[frontier] " fmt, ##__VA_ARGS__)

class quadrotor_ukf_node
{
 public:
  quadrotor_ukf_node();

  void imu_callback(const sensor_msgs::Imu::ConstPtr& msg);
  void slam_callback(const nav_msgs::Odometry::ConstPtr& msg);
  void pose_callback(const geometry_msgs::PoseStamped::ConstPtr& pose_msg);

 private:

  tf2_ros::StaticTransformBroadcaster tf_broadcaster;

  ros::NodeHandle nh, pnh;
  ros::Subscriber subImu, subSlam, pose_sub;
  ros::Publisher pubUKF;

  QuadrotorUKF quadrotorUKF;
  std::string frame_id;

  Eigen::Matrix<double, 4, 4> H_C_B;
  Eigen::Matrix<double, 9, 9> Rv;

  double alpha = 0.0;
  double beta  = 0.0;
  double kappa = 0.0;
  double stdAcc[3]     = {0,0,0};
  double stdW[3]       = {0,0,0};
  double stdAccBias[3] = {0,0,0};

};
