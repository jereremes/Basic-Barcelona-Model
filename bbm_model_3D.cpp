#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense> // For matrix operations

#include "bbm_model3D.h"

BarcelonaBasicModel::BarcelonaBasicModel(MaterialProperties props) : properties(props) {}

void BarcelonaBasicModel::computeStressStrain(SoilState& state, const Eigen::Matrix3d& d_strain, double d_suction, double dt) {
    // Update preconsolidation pressure using Newton-Raphson iteration
    double tolerance = 1e-6;
    double error = 1.0;
    int max_iter = 10;
    int iter = 0;
    double new_p_c = state.preconsolidation_pressure;
    while (error > tolerance && iter < max_iter) {
        double f = new_p_c - state.preconsolidation_pressure * pow((state.suction + properties.s0) / properties.s0, properties.lambda_s / properties.lambda);
        double df = 1.0 - (properties.lambda_s / properties.lambda) * (state.preconsolidation_pressure * pow((state.suction + properties.s0) / properties.s0, properties.lambda_s / properties.lambda)) / new_p_c;
        double delta = -f / df;
        new_p_c += delta;
        error = std::abs(delta);
        iter++;
    }
    state.preconsolidation_pressure = new_p_c;
    
    // Compute compliance tensor
    double bulk_modulus = (1.0 + properties.kappa) * state.preconsolidation_pressure / properties.kappa;
    double shear_modulus = 3 * (1 - 2 * properties.kappa) / (2 * (1 + properties.kappa)) * bulk_modulus;
    Eigen::Matrix3d compliance_tensor = Eigen::Matrix3d::Identity() / bulk_modulus;
    
    // Compute stress update using implicit formulation
    state.stress += dt * compliance_tensor * d_strain;
    state.suction += dt * d_suction;
}

void BarcelonaBasicModel::computeThermalConduction(SoilState& state, double dt, double thermal_gradient) {
    // Fourier's law: q = -k * dT/dx, implicit update
    double heat_flux = -properties.thermal_conductivity * thermal_gradient;
    double inv_denominator = 1.0 / (1.0 + dt * properties.thermal_conductivity);
    
    // Optimized implicit time integration for temperature update
    state.temperature = (state.temperature + dt * heat_flux / properties.specific_heat) * inv_denominator;
}

void BarcelonaBasicModel::computeHydrology(SoilState& state, double dt) {
    // Nonlinear permeability model: permeability decreases exponentially with suction
    double effective_permeability = properties.permeability * exp(-properties.alpha * state.suction);
    double water_flow = -effective_permeability * state.suction;
    
    // Update water content implicitly
    state.water_content = (state.water_content + dt * water_flow) / (1 + dt * effective_permeability);
}

//int main() {
//    // Define material properties
//    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
//    BarcelonaBasicModel bbm(props);
    
    // Initialize soil state
//    SoilState state;
//    state.stress = Eigen::Matrix3d::Zero();
//    state.strain = Eigen::Matrix3d::Zero();
//    state.suction = 20.0;
//    state.preconsolidation_pressure = 200.0;
//    state.temperature = 298.0;
//    state.water_content = 0.3;
    
    // Time step loop
//    double dt = 0.01;
//    Eigen::Matrix3d d_strain = Eigen::Matrix3d::Identity() * -0.001;
//    for (int step = 0; step < 100; ++step) {
//        bbm.computeStressStrain(state, d_strain, -0.5, dt);
//        bbm.computeThermalConduction(state, dt, -5.0);
//        bbm.computeHydrology(state, dt);
//    }
    
//    std::cout << "Final Mean Stress: " << state.stress.trace() / 3.0 << " kPa\n";
//    std::cout << "Final Temperature: " << state.temperature << " K\n";
//    std::cout << "Final Water Content: " << state.water_content << "\n";
//    return 0;
//}
