#include <quadrotor_ukf/quadrotor_ukf.h>

#include <iostream>
using namespace std;

QuadrotorUKF::QuadrotorUKF()
{
  // Init Dimensions
  stateCnt         = 12;
  procNoiseCnt     = 9;
  measNoiseSLAMCnt = 6;
  L = stateCnt + procNoiseCnt;
  // Init State
  xHist.clear();
  uHist.clear();
  xTimeHist.clear();
  Xa.setZero(stateCnt, 2*L+1);
  Va.setZero(procNoiseCnt, 2*L+1);
  P.setZero(stateCnt,stateCnt);
  P(0,0) = 0.5*0.5;
  P(1,1) = 0.5*0.5;
  P(2,2) = 0.1*0.1;
  P(3,3) = 0.1*0.1;
  P(4,4) = 0.1*0.1;
  P(5,5) = 0.1*0.1;
  P(6,6) = 10*PI/180*10*PI/180;
  P(7,7) = 10*PI/180*10*PI/180;
  P(8,8) = 10*PI/180*10*PI/180;
  P(9,9)   =  0.01*0.01;
  P(10,10) =  0.01*0.01;
  P(11,11) =  0.01*0.01;
  Rv.setIdentity(procNoiseCnt,procNoiseCnt);
  // Init Sigma Points
  alpha = 0.1;
  beta  = 2;
  kappa = 0;
  GenerateWeights();
  // Other Inits
  g = 9.81;
  initMeasure = false;
  initGravity = false;
}

          QuadrotorUKF::~QuadrotorUKF() { }
//bool      QuadrotorUKF::isInitialized() { return (initMeasure && initGravity); }
bool      QuadrotorUKF::isInitialized() { return (initMeasure); }

Eigen::Matrix<double, Eigen::Dynamic, 1>   QuadrotorUKF::GetState()           { return xHist.front(); }
ros::Time QuadrotorUKF::GetStateTime()       { return xTimeHist.front(); }
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>       QuadrotorUKF::GetStateCovariance() { return P;}
void      QuadrotorUKF::SetGravity(double _g) { g = _g; initGravity = true; }
void      QuadrotorUKF::SetImuCovariance(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& _Rv) { Rv = _Rv; }

void QuadrotorUKF::SetUKFParameters(double _alpha, double _beta, double _kappa)
{
  alpha = _alpha;
  beta  = _beta;
  kappa = _kappa;
  GenerateWeights();
}

void QuadrotorUKF::SetInitPose(Eigen::Matrix<double, Eigen::Dynamic, 1> p, ros::Time time)
{
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  x.setZero(stateCnt, 1);
  //x.rows(0,2)  = p.rows(0,2);
  //x.rows(6,8)  = p.rows(3,5);
  x.block<3,1>(0,0) = p.block<3,1>(0,0);
  x.block<3,1>(6,0) = p.block<3,1>(3,0);
  xHist.push_front(x);
  uHist.push_front(Eigen::MatrixXd::Zero(6, 1));
  xManHist.push_front(VIOUtil::expSO3(x.block<3,1>(6,0)));
  xTimeHist.push_front(time);
  initMeasure = true;
  cout<<initMeasure<<endl;
}

bool QuadrotorUKF::ProcessUpdate(Eigen::Matrix<double, Eigen::Dynamic, 1> u, ros::Time time)
{
  if (!initMeasure || !initGravity)
    return false;
  // Just update state, defer covariance update
  double dt = (time-xTimeHist.front()).toSec();
  //Eigen::Matrix<double, Eigen::Dynamic, 1> x = ProcessModel(xHist.front(), u, Eigen::MatrixXd::Zero(procNoiseCnt,1), dt);//ProcessModel(xHist.front(), u, zeros<colvec>(procNoiseCnt), dt);
  Eigen::Matrix<double, 3, 3> xManHisttmp = xManHist.front();
  Eigen::Matrix<double, Eigen::Dynamic, 1> x = ProcessModelManifold(xHist.front(), xManHisttmp, u, Eigen::MatrixXd::Zero(procNoiseCnt,1), dt);
  //Eigen::Matrix<double, Eigen::Dynamic, 1> x = ProcessModel(xHist.front(), u, Eigen::MatrixXd::Zero(procNoiseCnt,1), dt);
  xHist.push_front(x);
  uHist.push_front(u);
  xManHist.push_front(xManHisttmp);
  xTimeHist.push_front(time);

  return true;
}

bool QuadrotorUKF::MeasurementUpdateSLAM(const Eigen::Matrix<double, Eigen::Dynamic, 1>& z, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& RnSLAM, ros::Time time)
{
  // Init
  if (!initMeasure || !initGravity)
   return false;
  
  // A priori covariance
  std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator kx;
  std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator ku;
  std::list<ros::Time>::iterator kt;
  std::list<Eigen::Matrix<double, 3, 3> >::iterator kxManHist;

  PropagateAprioriCovariance(time, kx, ku, kt, kxManHist);

  Eigen::Matrix<double, Eigen::Dynamic, 1> x = *kx;
  Eigen::Matrix<double, 3, 3> x_manifold = *kxManHist;


  // Get Measurement
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H = MeasurementModelSLAM();
  Eigen::Matrix<double, Eigen::Dynamic, 1> za = H * x;
  // Compute Kalman Gain
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S = H * P * H.transpose() + RnSLAM;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K = P * H.transpose() * S.inverse();
  // Innovation
  Eigen::Matrix<double, Eigen::Dynamic, 1> inno = z - za;

  inno.block<3,1>(3,0) = VIOUtil::LogSO3(x_manifold.transpose() * VIOUtil::expSO3(z.block<3,1>(3,0)));//*VIOUtil::expSO3(z.block<3,1>(3,0)));
  // Posteriori Mean
  Eigen::Matrix<double, Eigen::Dynamic, 1> k_inno = K * inno;
  x.block<6,1>(0,0) += k_inno.block<6,1>(0,0);
  x.block<3,1>(6,0) = VIOUtil::LogSO3(x_manifold * VIOUtil::expSO3(k_inno.block<3,1>(6,0)));
  x.block<3,1>(9,0) += k_inno.block<3,1>(9,0);
  x_manifold = x_manifold * VIOUtil::expSO3(k_inno.block<3,1>(6,0));
  Eigen::Matrix<double, 6, 1> k_inno_parallel;
  k_inno_parallel.block<3,1>(0,0) = k_inno.block<3,1>(0,0);
  k_inno_parallel.block<3,1>(3,0) = k_inno.block<3,1>(6,0);
  Eigen::Matrix<double, 6, 6> M = VIOUtil::parallel_transport_trans(k_inno_parallel);

  *kx = x;
  *kxManHist = x_manifold;
  // Posteriori Covariance
  P = P - K * H * P;
  //Parallel transport
  P.block<3,3>(6,6) = M.block<3,3>(3,3) * P.block<3,3>(6,6) * M.block<3,3>(3,3).transpose();
  P.block(6, 0, 3, 6) = M.block<3,3>(3,3) * P.block(6, 0, 3, 6) ;
  P.block(0, 6, 6, 3) = P.block(0, 6, 6, 3) * M.block<3,3>(3,3).transpose();
  P.block(6, 9, 3, 3) = M.block<3,3>(3,3) * P.block(6, 9, 3, 3);
  P.block(9, 6, 3, 3) = P.block(9, 6, 3, 3) * M.block<3,3>(3,3).transpose();
  // Propagate Aposteriori State

  PropagateAposterioriState(kx, kxManHist, ku, kt);


  return true;
}

void QuadrotorUKF::GenerateWeights()
{
  // State + Process noise
  lambda = alpha*alpha*(L+kappa)-L;
  //wm = zeros<rowvec>(2*L+1);
  //wc = zeros<rowvec>(2*L+1);
  wm.setZero(1,2*L+1);
  wc.setZero(1,2*L+1);
  wm(0,0) = lambda / (L+lambda);
  wc(0,0) = lambda / (L+lambda) + (1-alpha*alpha+beta);
  for (int k = 1; k <= 2*L; k++)
  {
    wm(0,k) = 1 / (2 * (L+lambda));
    wc(0,k) = 1 / (2 * (L+lambda));
  }
  gamma = sqrt(L + lambda);
}

void QuadrotorUKF::GenerateSigmaPoints()
{
  // Expand state
  Eigen::Matrix<double, Eigen::Dynamic, 1> x  = xHist.back();
  Eigen::Matrix<double, 3, 3> xman  = xManHist.back();

  Eigen::Matrix<double, Eigen::Dynamic, 1> xaa;
  xaa.setZero(L,1);
  xaa.block(0,0,stateCnt,1) = x;

  //xaa.rows(0,stateCnt-1) = x;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Xaa; 
  Xaa.setZero(L, 2*L+1);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Paa;
  Paa.setZero(L,L);
  //Paa.submat(0, 0, stateCnt-1, stateCnt-1) = P;
  //Paa.submat(stateCnt, stateCnt, L-1, L-1) = Rv;
  Paa.block(0,0,stateCnt,stateCnt) = P;
  Paa.block(stateCnt,stateCnt, L - stateCnt, L - stateCnt) = Rv;
  // Matrix square root
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sqrtPaa = Paa.llt().matrixL();
  cout<<"state:"<<x<<endl;
  cout<<"sqrtPaa:"<<sqrtPaa<<endl;
  // Mean
  //Xaa.col(0) = xaa;
  //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  xaaMat = repmat(xaa,1,L);
  //Xaa.cols(1,L) = xaaMat + gamma * sqrtPaa;
  //Xaa.cols(L+1,L+L) = xaaMat - gamma * sqrtPaa;

  // Create a matrix with columns the increased state 0
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xaaMat;
  xaaMat.setZero(L, L);
  for (int i = 0; i < L; i++)
    xaaMat.col(i) = xaa;

  Xaa.col(0) = xaa;
  Xa_manifold_in.clear();
  Xa_manifold_in.push_back(xman);
  Xaa.block(0,1,L,L) =   xaaMat.block(0,0,L,L) + gamma * sqrtPaa;
  Xaa.block(0,L+1,L,L) = xaaMat.block(0,0,L,L) - gamma * sqrtPaa;
  //redefine the sigma points for the manifold part
  for (int i = 1; i < L+1; i++)
  {
   Xaa.block<3,1>(6,i) =   VIOUtil::LogSO3(xman * VIOUtil::expSO3(gamma * sqrtPaa.block<3,1>(6,i - 1)));
   Xa_manifold_in.push_back(xman * VIOUtil::expSO3(gamma * sqrtPaa.block<3,1>(6,i - 1)));
  }
  for (int i = L+1; i < 2*L+1; i++)
  {
   Xaa.block<3,1>(6,i) =   VIOUtil::LogSO3(xman * VIOUtil::expSO3(-gamma * sqrtPaa.block<3,1>(6,i - L - 1)));
   Xa_manifold_in.push_back(xman * VIOUtil::expSO3(-gamma * sqrtPaa.block<3,1>(6,i - L - 1)));
  }

  //Xa = Xaa.rows(0, stateCnt-1);
  //Va = Xaa.rows(stateCnt, L-1);
  Xa = Xaa.block(0,0,stateCnt, 2*L+1);
  Va = Xaa.block(stateCnt,0,L-stateCnt, 2*L+1);
}

Eigen::Matrix<double, Eigen::Dynamic, 1> QuadrotorUKF::ProcessModel(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, const Eigen::Matrix<double, 6, 1>& u, const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, double dt)
{
  Eigen::Matrix<double, 3, 3> R;
  R = VIOUtil::expSO3(x.block<3,1>(6,0));

  Eigen::Matrix<double, 3, 1> ag;
  ag(0,0) = 0;
  ag(1,0) = 0;
  ag(2,0) = g;
  // Acceleration
  Eigen::Matrix<double, Eigen::Dynamic, 1> a = u.block<3,1>(0,0) + v.block<3,1>(0,0);//u.rows(0,2) + v.rows(0,2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> ddx = R * (a - x.block<3,1>(9,0)) - ag;//R * (a - x.rows(9,11)) - ag;

 // Rotation
  Eigen::Matrix<double, 3, 1> w = u.block<3,1>(3,0) + v.block<3,1>(3,0);//u.rows(3,5) + v.rows(3,5);
  /*Eigen::Matrix<double, 3, 3> dR;//eye<mat>(3,3);
  dR.setIdentity();
  dR(0,1) = -w(2,0) * dt;
  dR(0,2) =  w(1,0) * dt;
  dR(1,0) =  w(2,0) * dt;
  dR(1,2) = -w(0,0) * dt;
  dR(2,0) = -w(1,0) * dt;
  dR(2,1) =  w(0,0) * dt;*/
  Eigen::Matrix<double, 3, 3> dR = VIOUtil::expSO3(w * dt);
  Eigen::Matrix<double, 3, 3> Rt = R * dR;
  // State
  Eigen::Matrix<double, Eigen::Dynamic, 1> xt = x;
  //xt.rows(0,2)  = x.rows(0,2) + x.rows(3,5)*dt + ddx*dt*dt/2;
  //xt.rows(3,5)  =               x.rows(3,5)    + ddx*dt     ;
  //xt.rows(6,8)  = R_to_ypr(Rt);
  //xt.rows(9,11) = x.rows(9,11)  + v.rows(6,8) *dt;
  xt.block<3,1>(0,0)  = x.block<3,1>(0,0) + x.block<3,1>(3,0)*dt + ddx*dt*dt/2;
  xt.block<3,1>(3,0)  =                     x.block<3,1>(3,0)    + ddx*dt     ;
  xt.block<3,1>(6,0)  = VIOUtil::LogSO3(Rt);
  xt.block<3,1>(9,0) = x.block<3,1>(9,0)  + v.block<3,1>(6,0) *dt;
  return xt;
}


Eigen::Matrix<double, Eigen::Dynamic, 1> QuadrotorUKF::ProcessModelManifold(const Eigen::Matrix<double, Eigen::Dynamic, 1>& x, Eigen::Matrix<double, 3, 3>& x_manifold, const Eigen::Matrix<double, 6, 1>& u, const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, double dt)
{
  Eigen::Matrix<double, 3, 3> R;
  R = x_manifold;

  Eigen::Matrix<double, 3, 1> ag;
  ag(0,0) = 0;
  ag(1,0) = 0;
  ag(2,0) = g;
  // Acceleration
  Eigen::Matrix<double, Eigen::Dynamic, 1> a = u.block<3,1>(0,0) + v.block<3,1>(0,0);//u.rows(0,2) + v.rows(0,2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> ddx = R * (a - x.block<3,1>(9,0)) - ag;//R * (a - x.rows(9,11)) - ag;


 // Rotation
  Eigen::Matrix<double, 3, 1> w = u.block<3,1>(3,0) + v.block<3,1>(3,0);//u.rows(3,5) + v.rows(3,5);
  /*Eigen::Matrix<double, 3, 3> dR;//eye<mat>(3,3);
  dR.setIdentity();
  dR(0,1) = -w(2,0) * dt;
  dR(0,2) =  w(1,0) * dt;
  dR(1,0) =  w(2,0) * dt;
  dR(1,2) = -w(0,0) * dt;
  dR(2,0) = -w(1,0) * dt;
  dR(2,1) =  w(0,0) * dt;*/
  Eigen::Matrix<double, 3, 3> dR = VIOUtil::expSO3(w * dt);
  Eigen::Matrix<double, 3, 3> Rt = R * dR;
  // State
  Eigen::Matrix<double, Eigen::Dynamic, 1> xt = x;
  //xt.rows(0,2)  = x.rows(0,2) + x.rows(3,5)*dt + ddx*dt*dt/2;
  //xt.rows(3,5)  =               x.rows(3,5)    + ddx*dt     ;
  //xt.rows(6,8)  = R_to_ypr(Rt);
  //xt.rows(9,11) = x.rows(9,11)  + v.rows(6,8) *dt;
  xt.block<3,1>(0,0)  = x.block<3,1>(0,0) + x.block<3,1>(3,0)*dt + ddx*dt*dt/2;
  xt.block<3,1>(3,0)  =                     x.block<3,1>(3,0)    + ddx*dt     ;
  x_manifold = Rt;
  xt.block<3,1>(6,0)  = VIOUtil::LogSO3(Rt);
  xt.block<3,1>(9,0) = x.block<3,1>(9,0)  + v.block<3,1>(6,0) *dt;
  return xt;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> QuadrotorUKF::MeasurementModelSLAM()
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H;// = zeros<mat>(measNoiseSLAMCnt, stateCnt);
  H.setZero(measNoiseSLAMCnt, stateCnt);
  H(0,0) = 1;
  H(1,1) = 1;
  H(2,2) = 1;
  H(3,6) = 1;
  H(4,7) = 1;
  H(5,8) = 1;
  return H;
}

void QuadrotorUKF::PropagateAprioriCovariance(const ros::Time time,
                                              std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator& kx, std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator& ku, std::list<ros::Time>::iterator& kt, std::list<Eigen::Matrix<double, 3, 3> >::iterator& kxManHist)
{
  // Find aligned state, Time
  double mdt = NUM_INF;
  std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator k1 = xHist.begin();
  std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator k2 = uHist.begin();
  std::list<ros::Time>::iterator k3 = xTimeHist.begin();
  int                       k4 = 0;
  std::list<Eigen::Matrix<double, 3, 3> >::iterator k5 = xManHist.begin();

  for (; k1 != xHist.end(); k1++, k2++, k3++, k4++, k5++)
  {
    double dt = fabs((*k3 - time).toSec());
    if (dt < mdt)
    {
      mdt = dt;
      kx  = k1;
      ku  = k2;
      kt  = k3;
      kxManHist = k5;
    }
    else
    {
      break;
    }
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> cx = *kx;
  Eigen::Matrix<double, 3, 3> cxManHist = *kxManHist;
  ros::Time ct = *kt;
  Eigen::Matrix<double, Eigen::Dynamic, 1> px = xHist.back();
  Eigen::Matrix<double, 3, 3> pxManHist = xManHist.back();
  ros::Time pt = xTimeHist.back();
  double dt = (ct - pt).toSec();
  if (fabs(dt) < 0.001)
  {
    kx = xHist.begin();
    ku = uHist.begin();
    kt = xTimeHist.begin();
    kxManHist = xManHist.begin();
    return;
  }

  // Delete redundant states
  xHist.erase(k1, xHist.end());
  uHist.erase(k2, uHist.end());
  xTimeHist.erase(k3, xTimeHist.end());
  xManHist.erase(k5, xManHist.end());
  // rot, gravity
  Eigen::Matrix<double, 3, 3> pR;
  pR = pxManHist;//VIOUtil::expSO3(px.block<3,1>(6,0));//px.rows(6,8)
  Eigen::Matrix<double, 3, 1> ag;
  ag(0,0) = 0;
  ag(1,0) = 0;
  ag(2,0) = g;

  // Linear Acceleration
  Eigen::Matrix<double, 3, 1> dv = cx.block<3,1>(3,0) - px.block<3,1>(3,0);//rows(3,5);
  Eigen::Matrix<double, 3, 1> a = pR.transpose() * (dv / dt + ag) + px.block<3,1>(9,0);//.rows(9,11);

  // Angular Velocity
  Eigen::Matrix<double, 3, 3> dR;
  dR = pR.transpose() * cxManHist;//VIOUtil::expSO3(cx.block<3,1>(6,0));//cx.rows(6,8)
  Eigen::Matrix<double, 3, 1> w = VIOUtil::LogSO3(dR)/dt;

  //w(0,0) = dR(2,1) / dt;
  //w(1,0) = dR(0,2) / dt;
  //w(2,0) = dR(1,0) / dt;
  // Assemble state and control
  Eigen::Matrix<double, 6, 1> u;
  u.block<3,1>(0,0) = a;
  u.block<3,1>(3,0) = w;

  // Generate sigma points
  GenerateSigmaPoints();
  std::vector<Eigen::Matrix<double, 3, 3> > vec_R;	

  // Mean
  for (int k = 0; k < 2*L+1; k++){
    //Xa.col(k) = ProcessModel(Xa.col(k), u, Va.col(k), dt);
    Xa.col(k) = ProcessModelManifold(Xa.col(k), Xa_manifold_in.at(k), u, Va.col(k), dt);
    //Xa.col(k) = ProcessModel(Xa.col(k), u, Va.col(k), dt);
    vec_R.push_back(Xa_manifold_in.at(k));
  }

   // Now we can get the mean...
  Eigen::Matrix<double, Eigen::Dynamic, 1> xa;
  xa.resize(Xa.rows(),1);
  for (int i = 0; i < 2 * L + 1; i++)
  {
  	xa += wm(0,i) * Xa.col(i);// = sum( repmat(wm,stateCnt,1) % Xa, 1 );
  }
  //compute the mean for the manifold part
  meanR = VIOUtil::MeanofSigmaPointsManifoldSO3(vec_R, wm);//VIOUtil::expSO3(xa.block<3,1>(6,0));//VIOUtil::MeanofSigmaPointsManifoldSO3(vec_R, wm);

  //    Eigen::Matrix<double, 3, 3> meanRtmp = VIOUtil::MeanofSigmaPointsManifoldSO3(vec_R, wm);
  //xa.block<3,1>(6,0)  = VIOUtil::LogSO3(meanR);

  // Covariance
  P.setZero(Xa.rows(), Xa.rows());//.zeros();
  for (int k = 0; k < 2*L+1; k++)
  {
    Eigen::Matrix<double, Eigen::Dynamic, 1> d = Xa.col(k) - xa;
    //redefine the part for the manifold
    d.block<3,1>(6,0) = VIOUtil::LogSO3(Xa_manifold_in.at(k));
    P += wc(0,k) * d * d.transpose();
  }

  return;
}

void QuadrotorUKF::PropagateAposterioriState(std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator kx, std::list<Eigen::Matrix<double, 3, 3> >::iterator kxManHist, std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator ku, std::list<ros::Time>::iterator kt)
{
  for (; kx != xHist.begin(); kx--, ku--, kt--, kxManHist--)
  {
    std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator _kx = kx;
    _kx--;
    std::list<Eigen::Matrix<double, Eigen::Dynamic, 1>>::iterator _ku = ku;
    _ku--;
    std::list<ros::Time>::iterator _kt = kt;
    _kt--;
     std::list<Eigen::Matrix<double, 3, 3> >::iterator _kxManHist  = kxManHist ;
    _kxManHist--;
    //*_kx = ProcessModel(*kx, *_ku, Eigen::MatrixXd::Zero(procNoiseCnt,1), (*_kt - *kt).toSec());

    Eigen::Matrix<double, 3, 3>  xManHisttmp = *kxManHist;
    *_kx = ProcessModelManifold(*kx, xManHisttmp, *_ku,  Eigen::MatrixXd::Zero(procNoiseCnt,1), (*_kt - *kt).toSec());
    *_kxManHist = xManHisttmp;
  }
}

/*void QuadrotorUKF::PropagateAposterioriState(std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator kx, std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator ku, std::list<ros::Time>::iterator kt)
{
  for (; kx != xHist.begin(); kx--, ku--, kt--)
  {
    std::list<Eigen::Matrix<double, Eigen::Dynamic, 1> >::iterator _kx = kx;
    _kx--;
    std::list<Eigen::Matrix<double, Eigen::Dynamic, 1>>::iterator _ku = ku;
    _ku--;
    std::list<ros::Time>::iterator _kt = kt;
    _kt--;
    *_kx = ProcessModel(*kx, *_ku, Eigen::MatrixXd::Zero(procNoiseCnt,1), (*_kt - *kt).toSec());
  }
}*/