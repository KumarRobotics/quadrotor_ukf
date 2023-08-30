#include <quadrotor_ukf/quadrotor_ukf_node_redo.h>

quadrotor_ukf_node::quadrotor_ukf_node() :
  pnh("~")
{

  subImu = nh.subscribe("imu", 1, &quadrotor_ukf_node::imu_callback, this, ros::TransportHints().tcpNoDelay());
  subSlam = nh.subscribe("odom_slam", 1, &quadrotor_ukf_node::slam_callback, this, ros::TransportHints().tcpNoDelay());
  pose_sub = nh.subscribe("/qvio/pose", 1, &quadrotor_ukf_node::pose_callback, this, ros::TransportHints().tcpNoDelay());

  pubUKF = nh.advertise<nav_msgs::Odometry>("control_odom", 10);

  H_C_B (0,0) = 1;
  H_C_B (0,1) = 0;
  H_C_B (0,2) = 0;
  H_C_B (0,3) = 0;
  H_C_B (1,0) = 0;
  H_C_B (1,1) = -1;
  H_C_B (1,2) = 0;
  H_C_B (1,3) = 0;
  H_C_B (2,0) = 0;
  H_C_B (2,1) = 0;
  H_C_B (2,2) = -1;
  H_C_B (2,3) = 0.0;
  H_C_B (3,0) = 0;
  H_C_B (3,1) = 0;
  H_C_B (3,2) = 0;
  H_C_B (3,3) = 1;

  pnh.param("frame_id", frame_id, std::string("/world"));
  pnh.param("alpha", alpha, 0.4);
  pnh.param("beta" , beta , 2.0);
  pnh.param("kappa", kappa, 0.0);
  pnh.param("noise_std/process/acc/x", stdAcc[0], 0.2);
  pnh.param("noise_std/process/acc/y", stdAcc[1], 0.2);
  pnh.param("noise_std/process/acc/z", stdAcc[2], 0.2);
  pnh.param("noise_std/process/w/x", stdW[0], 0.1);
  pnh.param("noise_std/process/w/y", stdW[1], 0.1);
  pnh.param("noise_std/process/w/z", stdW[2], 0.1);
  pnh.param("noise_std/process/acc_bias/x", stdAccBias[0], 0.05);
  pnh.param("noise_std/process/acc_bias/y", stdAccBias[1], 0.05);
  pnh.param("noise_std/process/acc_bias/z", stdAccBias[2], 0.05);

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
  quadrotorUKF.SetUKFParameters(alpha, beta, kappa);
  quadrotorUKF.SetImuCovariance(Rv);
}


void quadrotor_ukf_node::imu_callback(const sensor_msgs::Imu::ConstPtr& msg)
{
  // Assemble control input, and calibration
  static int calLimit = 100;
  static int calCnt   = 0;
  //static colvec ag = zeros<colvec>(3);
  static Eigen::Matrix<double, 3, 1> ag;
  Eigen::Matrix<double, 6, 1> u;
  u(0,0) = msg->linear_acceleration.x;
  u(1,0) = -msg->linear_acceleration.y;
  u(2,0) = -msg->linear_acceleration.z;
  u(3,0) = msg->angular_velocity.x;
  u(4,0) = -msg->angular_velocity.y;
  u(5,0) = -msg->angular_velocity.z;
  if (calCnt < calLimit)       // Calibration
  {
    calCnt++;
    ag += u.block(0,0,3,1);//rows(0,2);
  }
  else if (calCnt == calLimit) // Save gravity vector
  {
    calCnt++;
    ag /= calLimit;
    double g = ag.norm();//norm(ag,2);
    quadrotorUKF.SetGravity(g);
  }
  else if (quadrotorUKF.ProcessUpdate(u, msg->header.stamp))  // Process Update
  {
    nav_msgs::Odometry odomUKF;
    // Publish odom
    odomUKF.header.stamp = quadrotorUKF.GetStateTime();
    odomUKF.header.frame_id = frame_id;
    Eigen::Matrix<double, Eigen::Dynamic, 1> x = quadrotorUKF.GetState();
    odomUKF.pose.pose.position.x = x(0,0);
    odomUKF.pose.pose.position.y = x(1,0);
    odomUKF.pose.pose.position.z = x(2,0);
    Eigen::Matrix<double, 4, 1> q = VIOUtil::MatToQuat(VIOUtil::ypr_to_R(x.block(6,0,3,1)));
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
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  P = quadrotorUKF.GetStateCovariance();
    for (int j = 0; j < 6; j++)
      for (int i = 0; i < 6; i++)
        odomUKF.pose.covariance[i+j*6] = P((i<3)?i:i+3 , (j<3)?j:j+3);
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
        odomUKF.twist.covariance[i+j*6] = P(i+3 , j+3);
    pubUKF.publish(odomUKF);
/*
    // Publish bias
    geometry_msgs::Vector3 bias;
    bias.x = x(9);
    bias.y = x(10);
    bias.z = x(11);
    pubBias.publish(bias);
*/
  }
}

void quadrotor_ukf_node::slam_callback(const nav_msgs::Odometry::ConstPtr& msg)
{
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
  // Assemble measurement covariance
  Eigen::Matrix<double, 6, 6> RnSLAM;
  RnSLAM.setZero();// = zeros<mat>(6,6);
  //for (int j = 0; j < 3; j++)
    //for (int i = 0; i < 3; i++)
      //RnSLAM(i,j) = msg->pose.covariance[i+j*6];
  RnSLAM(0,0) = msg->pose.covariance[0];
  RnSLAM(1,1) = msg->pose.covariance[1+1*6];
  RnSLAM(2,2) = msg->pose.covariance[2+2*6];
  RnSLAM(3,3) = msg->pose.covariance[3+3*6];
  RnSLAM(4,4) = msg->pose.covariance[4+4*6];
  RnSLAM(5,5) = msg->pose.covariance[5+5*6];

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
  H_R_R0 = H_C_B*H_C_C0*H_C_B.inverse();
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
  RnSLAM_new.block(0,0,3,3) = H_C_B.block(0, 0, 3, 3)*RnSLAM.block(0,0,3,3)*H_C_B.block(0, 0, 3, 3).transpose();
  RnSLAM_new.block(3,3,3,3) = H_C_B.block(0, 0, 3, 3)*RnSLAM.block(3,3,3,3)*H_C_B.block(0, 0, 3, 3).transpose();
  
  // Measurement update
  if (quadrotorUKF.isInitialized())
  {
    quadrotorUKF.MeasurementUpdateSLAM(z_new, RnSLAM_new, msg->header.stamp);
  }
  else
  {
    quadrotorUKF.SetInitPose(z, msg->header.stamp);
  }
}

void quadrotor_ukf_node::pose_callback(const geometry_msgs::PoseStamped::ConstPtr& pose)
{
  geometry_msgs::TransformStamped static_transformStamped;
  static_transformStamped.header.frame_id = "map";
  static_transformStamped.child_frame_id = "qvio";

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
  ros::init(argc, argv, "quadrotor_ukf_node");
  quadrotor_ukf_node quadrotor_ukf_node;
  while(ros::ok())
  {
    ros::spinOnce();
  }
  return 0;
}

