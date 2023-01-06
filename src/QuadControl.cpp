#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif


void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{

  float length = L / sqrt(2.f);

  float Ft = collThrustCmd;
  float Fp = momentCmd.x / length;
  float Fq = momentCmd.y / length;
  float Fr = momentCmd.z / kappa;

  float F1 = (Ft + Fp + Fq - Fr) / 4;
  float F2 = (F1 - (Fp - Fr) / 2);
  float F4 = ((Ft - Fp) / 2) - F2;
  float F3 = Ft - F1 - F2 - F4;

  cmd.desiredThrustsN[0] = CONSTRAIN(F1, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[1] = CONSTRAIN(F2, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[2] = CONSTRAIN(F3, minMotorThrust, maxMotorThrust);
  cmd.desiredThrustsN[3] = CONSTRAIN(F4, minMotorThrust, maxMotorThrust);

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  V3F momentCmd;
  V3F I = V3F(Ixx, Iyy, Izz);
  momentCmd = I * kpPQR * (pqrCmd - pqr);

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  if (collThrustCmd > 0.f) {
       float R_11 = R(0, 0);
       float R_12 = R(0, 1);
       float R_13 = R(0, 2);
       float R_21 = R(1, 0);
       float R_22 = R(1, 1);
       float R_23 = R(1, 2);
       float R_33 = R(2, 2);

       float c = collThrustCmd / mass;
       V3F b = V3F(R(0, 2), R(1, 2), 0.f);
       V3F b_c = V3F(accelCmd.x / -c, accelCmd.y / -c, 0.f);
       b_c.constrain(-maxTiltAngle, maxTiltAngle);

       V3F b_error = b_c - b;
       V3F b_c_dot = kpBank * b_error;

       pqrCmd.x = (R_21 * b_c_dot.x - R_11 * b_c_dot.y) / R_33;
       pqrCmd.y = (R_22 * b_c_dot.x - R_12 * b_c_dot.y) / R_33;
       pqrCmd.z = 0.f;
  }
  else {
       pqrCmd.x = 0.f;
       pqrCmd.y = 0.f;
       pqrCmd.z = 0.f;
  }
  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;
  
  float zErr = posZCmd - posZ;
  float zDotCmd = kpPosZ * zErr + velZCmd;
  zDotCmd = CONSTRAIN(zDotCmd, -maxDescentRate, maxDescentRate);
  integratedAltitudeError += zErr * dt;
  float zErrDot = zDotCmd - velZ;
  accelZCmd += KiPosZ * integratedAltitudeError + kpVelZ * zErrDot;
  float BZ = R(2, 2);
  thrust = mass * ((float)CONST_GRAVITY - accelZCmd) / BZ;

  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{

  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  V3F accelCmd = accelCmdFF;

  velCmd = kpPosXY * (posCmd - pos);
  velCmd.constrain(-maxSpeedXY, maxSpeedXY);

  const V3F err_dot = velCmd - vel;
  accelCmd += kpVelXY * err_dot;

  accelCmd.constrain(-maxAccelXY, maxAccelXY);


  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{

  float yawRateCmd=0;
  float yawErr = yawCmd - yaw;
  yawErr = fmodf(yawErr, 2.*F_PI);
  
  if (yawErr > F_PI) {
       yawErr -= 2.f * F_PI;
  }
  else if (yawErr < - F_PI) {
       yawErr += 2.f * F_PI;
  }

  yawRateCmd = this->kpYaw * yawErr;

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
