#ifndef QUADROTOR_UKF_H
#define QUADROTOR_UKF_H

#include <iostream>
#include <string.h>
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>
#include <ros/ros.h>
//#include "pose_utils.h"
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <quadrotor_ukf/vio_utils.h>
using namespace std;

#define PI 3.14159265359
#define NUM_INF 999999.9

class QuadrotorUKF
{
  private:

    // State History and Covariance
    std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >  xHist;
    std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >  uHist;
    std::list<ros::Time> xTimeHist;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P;
    Eigen::Matrix<double, 3, 3> meanR;

    // Process Covariance Matrix
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Rv;

    // Instance sigma points
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xa;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Va;

    // Initial process update indicator
    bool initMeasure;
    bool initGravity;

    // Dimemsions
    int stateCnt;
    int procNoiseCnt;
    int measNoiseSLAMCnt;
    int L;

    // Gravity
    double g;

    // UKF Parameters
    double alpha;
    double beta;
    double kappa;
    double lambda;
    double gamma;
    std::vector<Eigen::Matrix<double, 3, 3> > Xa_manifold_in;
    Eigen::Matrix<double, 3, 3>  currR_mean;



    // UKF Weights
    Eigen::Matrix<double, 1, Eigen::Dynamic>  wm;
    Eigen::Matrix<double, 1, Eigen::Dynamic>  wc;

    // Private functions
    void GenerateWeights();
    void GenerateSigmaPoints();
    Eigen::Matrix<double, Eigen::Dynamic, 1> ProcessModel(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, const Eigen::Matrix<double, 6, 1>& u, const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, double dt);
    Eigen::Matrix<double, Eigen::Dynamic, 1> ProcessModelManifold(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, Eigen::Matrix<double, 3, 3>& x_manifold, const Eigen::Matrix<double, 6, 1>& u, const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, double dt);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MeasurementModelSLAM();
    void PropagateAprioriCovariance(const ros::Time time,
                                              std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator& kx, std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator& ku, std::list<ros::Time>::iterator& kt, std::list<Eigen::Matrix<double, 3, 3> >::iterator& kxManHist);
    //void PropagateAposterioriState(std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator kx, std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator ku, std::list<ros::Time>::iterator kt);
    void PropagateAposterioriState(std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator kx, std::list<Eigen::Matrix<double, 3, 3> >::iterator kxManHist, std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator ku, std::list<ros::Time>::iterator kt);

  public:
    std::list<Eigen::Matrix<double, 3, 3> > xManHist;

    QuadrotorUKF();
    ~QuadrotorUKF();

    //manifold option
    bool manifold;
    
    //bool QuadrotorUKF::isInitialized() { return (initMeasure && initGravity); }();
    bool isInitialized();
    Eigen::Matrix<double, Eigen::Dynamic, 1> GetState();
    ros::Time GetStateTime();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> GetStateCovariance();

    void SetGravity(double _g);
    void SetImuCovariance(const Eigen::Matrix<double, Eigen::Dynamic,  Eigen::Dynamic>& _Rv);
    void SetUKFParameters(double _alpha, double _beta, double _kappa);
    void SetInitPose(Eigen::Matrix<double, Eigen::Dynamic, 1> p, ros::Time time);

    bool ProcessUpdate(Eigen::Matrix<double, Eigen::Dynamic, 1> u, ros::Time time);
    bool MeasurementUpdateSLAM(const Eigen::Matrix<double, Eigen::Dynamic, 1>&  z, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&  RnSLAM, ros::Time time);
};

#endif
