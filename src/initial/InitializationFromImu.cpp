/* ----------------------------------------------------------------------------
 * Copyright 2017, Massachusetts Institute of Technology,
 * Cambridge, MA 02139
 * All Rights Reserved
 * Authors: Luca Carlone, et al. (see THANKS for the full author list)
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 * @file   InitializationFromImu.cpp
 * @brief  Class to initialize VIO pipeline from IMU measurements only.
 * @author Antoni Rosinol
 */

#include "kimera-vio/initial/InitializationFromImu.h"

#include <cstdlib>
#include <iostream>

#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Unit3.h>

namespace VIO {

namespace {

// JT-ZERO opt-in initialization path. Kimera's generic AlignGravityVectors()
// treats vectors with (1 - dot) < 1e-3 as already aligned. That corresponds
// to roughly 2.56 degrees and is too large for the JT-Zero stationary startup:
// a real ~2.4 degree tilt was being collapsed to identity and absorbed into
// accelerometer bias. Keep the upstream behavior unchanged unless the explicit
// environment switch is enabled.
gtsam::Rot3 AlignGravityVectorsExactForJtZero(
    const gtsam::Vector3& local_gravity_dir,
    const gtsam::Vector3& global_gravity_dir) {
  gtsam::Unit3 local_gravity(local_gravity_dir);
  gtsam::Unit3 global_gravity(global_gravity_dir);

  const double c = local_gravity.dot(global_gravity);
  constexpr double kExactEps = 1e-12;
  if (std::fabs(1.0 - c) < kExactEps) {
    return gtsam::Rot3();
  }

  gtsam::Unit3 cross_product = local_gravity.cross(global_gravity);
  if (std::fabs(1.0 + c) < kExactEps) {
    gtsam::Unit3 perturbed_gravity(
        local_gravity.unitVector() + gtsam::Vector3(1, 2, 3));
    cross_product = local_gravity.cross(perturbed_gravity);
    if (std::isnan(cross_product.unitVector()(0))) {
      perturbed_gravity = gtsam::Unit3(
          local_gravity.unitVector() + gtsam::Vector3(3, 2, 1));
      cross_product = local_gravity.cross(perturbed_gravity);
    }
    return gtsam::Rot3::Expmap(cross_product.unitVector() * M_PI);
  }

  return gtsam::Rot3::AlignPair(
      cross_product, global_gravity, local_gravity);
}

}  // namespace

VioNavState InitializationFromImu::getInitialStateEstimate(
    const ImuAccGyrS& imu_accgyr,
    const gtsam::Vector3& global_gravity,
    const bool& round) {
  LOG(WARNING) << "InitializationFromImu: assumes that the "
                  "vehicle is stationary and upright along some axis,"
                  "and gravity vector is along a single axis!";

  // Compute mean acceleration and angular velocity.
  ImuAccGyr mean_accgyr = computeAverageImuMeasurements(imu_accgyr);

  // Guess initial pose assuming vehicle is stationary (zero acceleration
  // besides negative of gravity), and that gravity is aligned with one IMU
  // axis (vehicle upright).
  gtsam::Pose3 initial_pose_guess =
      guessPoseFromImuMeasurements(mean_accgyr.head(3),  // Mean Acc values.
                                   global_gravity,
                                   round);

  // Zero Velocity assumption!
  gtsam::Vector3 velocity_guess = gtsam::Vector3::Zero();

  // Convert global gravity to local frame of reference.
  // TODO(Toni): this guy should be the same as -1 * mean_acc (aka measured
  // local gravity)...
  gtsam::Vector3 local_gravity =
      initial_pose_guess.rotation().inverse().matrix() * global_gravity;

  // Guess IMU bias. Assumes static vehicle!
  ImuBias imu_bias_guess = guessImuBias(mean_accgyr, local_gravity);

  // JT-ZERO diagnostic: expose the exact stationary IMU packet and the
  // pose/bias seed produced by InitializationFromImu. This is intentionally
  // gated by an environment variable and does not change initialization math.
  if (std::getenv("JTZERO_DIAG_IMU_INIT") != nullptr) {
    constexpr double kRadToDeg = 57.2957795130823208768;

    const gtsam::Vector3 mean_acc = mean_accgyr.head(3);
    const gtsam::Vector3 mean_gyro = mean_accgyr.tail(3);
    const gtsam::Vector3 rpy = initial_pose_guess.rotation().rpy();

    std::cerr
        << "[JT-IMU-INIT]"
        << " initMode="
        << (std::getenv("JTZERO_GRAVITY_ALIGNED_IMU_INIT") != nullptr
                ? "jtzero_exact_gravity"
                : "kimera_default")
        << " samples=" << imu_accgyr.cols()
        << " meanAcc=[" << mean_acc.transpose() << "]"
        << " accNorm=" << mean_acc.norm()
        << " meanGyro=[" << mean_gyro.transpose() << "]"
        << " globalG=[" << global_gravity.transpose() << "]"
        << " localG=[" << local_gravity.transpose() << "]"
        << " initRPYdeg=["
        << rpy.x() * kRadToDeg << " "
        << rpy.y() * kRadToDeg << " "
        << rpy.z() * kRadToDeg << "]"
        << " initBA=[" << imu_bias_guess.accelerometer().transpose() << "]"
        << " initBG=[" << imu_bias_guess.gyroscope().transpose() << "]"
        << std::endl;
  }

  // Return estimated state.
  return VioNavState(initial_pose_guess, velocity_guess, imu_bias_guess);
}

gtsam::Pose3 InitializationFromImu::guessPoseFromImuMeasurements(
    const ImuAcc& mean_acc,
    const gtsam::Vector3& global_gravity,
    const bool& round) {
  // We measure the negative of gravity. Assumes static vehicle.
  gtsam::Vector3 measured_gravity = -1.0 * mean_acc;
  // Align measured gravity with real gravity to figure out our attitude.
  // Assumes gravity aligned along an axis.
  const bool jtzero_gravity_aligned_init =
      std::getenv("JTZERO_GRAVITY_ALIGNED_IMU_INIT") != nullptr;
  gtsam::Rot3 attitude_wrt_gravity =
      jtzero_gravity_aligned_init
          ? AlignGravityVectorsExactForJtZero(measured_gravity, global_gravity)
          : UtilsOpenCV::AlignGravityVectors(
                measured_gravity, global_gravity, round);
  // Absolute translation is unobservable, so return [0, 0, 0].
  return gtsam::Pose3(attitude_wrt_gravity, gtsam::Point3::Zero());
}

ImuBias InitializationFromImu::guessImuBias(
    const ImuAccGyr& mean_accgyr,
    const gtsam::Vector3& local_gravity) {
  // Assumes static vehicle.
  return ImuBias(mean_accgyr.head(3) + local_gravity,  // Acceleration
                 mean_accgyr.tail(3));                 // Gyro
}

}  // namespace VIO
