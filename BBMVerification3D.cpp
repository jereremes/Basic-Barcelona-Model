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

// Function to test an analytical solution for thermal diffusion in 3D
bool testThermalDiffusion3D() {
    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
    BarcelonaBasicModel bbm(props);
    SoilState state;
    state.temperature = 298.0; // Initial temperature
    
    double dt = 0.01;
    double analytical_temperature = state.temperature;
    
    for (int step = 0; step < 100; ++step) {
        bbm.computeThermalConduction(state, dt, -5.0);
        double temperature_change = (props.thermal_conductivity * 5.0 * dt / props.specific_heat);
        analytical_temperature -= temperature_change;
        
        // Detailed debug output
        std::cout << "Step " << step << ": "
                  << "State Temperature = " << state.temperature << ", "
                  << "Analytical Temperature = " << analytical_temperature << ", "
                  << "Temperature Change = " << temperature_change << std::endl;
    }
    
    return std::abs(state.temperature - analytical_temperature) < 1e-3;
}

// Function to test Terzaghi’s consolidation in 3D for hydrology
bool testTerzaghiConsolidation3D() {
    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
    BarcelonaBasicModel bbm(props);
    SoilState state;
    state.water_content = 0.3; // Initial water content
    
    double dt = 0.01;
    double analytical_water_content = state.water_content;
    
    for (int step = 0; step < 100; ++step) {
        bbm.computeHydrology(state, dt);
        analytical_water_content *= exp(-props.permeability * dt);
        
        // Detailed debug output
        std::cout << "Step " << step << ": "
                  << "State Water Content = " << state.water_content << ", "
                  << "Analytical Water Content = " << analytical_water_content << std::endl;
    }
    
    return std::abs(state.water_content - analytical_water_content) < 1e-3;
}

// Function to perform a convergence study on stress-strain behavior in 3D
bool testStressStrainConvergence3D() {
    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
    BarcelonaBasicModel bbm(props);
    SoilState state;
    state.stress = Eigen::Matrix3d::Zero(); // Initial stress
    state.suction = 50.0; // Initial suction
    
    Eigen::Matrix3d d_strain = Eigen::Matrix3d::Identity() * -0.001;
    double dt = 0.01;
    Eigen::Matrix3d previous_stress = state.stress;
    
    for (int step = 0; step < 100; ++step) {
        bbm.computeStressStrain(state, d_strain, -0.5, dt);
        if ((state.stress - previous_stress).norm() < 1e-5) return true;
        previous_stress = state.stress;
        
        // Detailed debug output
        std::cout << "Step " << step << ": "
                  << "State Stress = \n" << state.stress << std::endl;
    }
    
    return false;
}

int main() {
    bool thermal_test = testThermalDiffusion3D();
    bool hydrology_test = testTerzaghiConsolidation3D();
    bool stress_test = testStressStrainConvergence3D();

    std::cout << "Running the 3D tests:" << std::endl;
    
    std::cout << "Thermal Diffusion Test: " << (thermal_test ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Hydrology Test (Terzaghi): " << (hydrology_test ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Stress-Strain Convergence Test: " << (stress_test ? "PASSED" : "FAILED") << std::endl;
    
    return 0;
}