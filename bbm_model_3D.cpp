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

    // Debug output
    std::cout << "Stress: \n" << state.stress << "\n";
    std::cout << "Suction: " << state.suction << "\n";
}

void BarcelonaBasicModel::computeThermalConduction(SoilState& state, double dt, double thermal_gradient) {
    // Fourier's law: q = -k * dT/dx, implicit update
    double heat_flux = -properties.thermal_conductivity * thermal_gradient;
    double inv_denominator = 1.0 / (1.0 + dt * properties.thermal_conductivity);
    
    // Optimized implicit time integration for temperature update
    state.temperature = (state.temperature + dt * heat_flux / properties.specific_heat) * inv_denominator;

    // Debug output
    std::cout << "Temperature: " << state.temperature << "\n";
}

void BarcelonaBasicModel::computeHydrology(SoilState& state, double dt) {
    // Nonlinear permeability model: permeability decreases exponentially with suction
    double effective_permeability = properties.permeability * exp(-properties.alpha * state.suction);
    double water_flow = -effective_permeability * state.suction;
    
    // Update water content implicitly
    state.water_content = (state.water_content + dt * water_flow) / (1 + dt * effective_permeability);

    // Debug output
    std::cout << "Water Content: " << state.water_content << "\n";
}
