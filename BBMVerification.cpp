#include <iostream>
#include <cmath>
#include "bbm_model.h" // Include the BBM model header file

// Function to test an analytical solution for thermal diffusion
bool testThermalDiffusion() {
    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
    BarcelonaBasicModel bbm(props);
    SoilState state = {100.0, 50.0, 20.0, 200.0, 298.0, 0.3};
    
    double dt = 0.01;
    double analytical_temperature = state.temperature;
    
    for (int step = 0; step < 100; ++step) {
        bbm.computeThermalConduction(state, dt, -5.0);
        analytical_temperature -= (props.thermal_conductivity * 5.0 * dt / props.specific_heat);
    }
    
    return std::abs(state.temperature - analytical_temperature) < 1e-3;
}

// Function to test Terzaghi’s consolidation in 1D for hydrology
bool testTerzaghiConsolidation() {
    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
    BarcelonaBasicModel bbm(props);
    SoilState state = {100.0, 50.0, 20.0, 200.0, 298.0, 0.3};
    
    double dt = 0.01;
    double analytical_water_content = state.water_content;
    
    for (int step = 0; step < 100; ++step) {
        bbm.computeHydrology(state, dt);
        analytical_water_content *= exp(-props.permeability * dt);
    }
    
    return std::abs(state.water_content - analytical_water_content) < 1e-3;
}

// Function to perform a convergence study on stress-strain behavior
bool testStressStrainConvergence() {
    MaterialProperties props = {0.1, 0.02, 0.08, 10.0, 1.2, 1.5, 900, 1e-6, 0.05};
    BarcelonaBasicModel bbm(props);
    SoilState state = {100.0, 50.0, 20.0, 200.0, 298.0, 0.3};
    
    double dt = 0.01;
    double previous_stress = state.mean_stress;
    
    for (int step = 0; step < 100; ++step) {
        bbm.computeStressStrain(state, -0.001, -0.5, dt);
        if (std::abs(state.mean_stress - previous_stress) < 1e-5) return true;
        previous_stress = state.mean_stress;
    }
    
    return false;
}

int main() {
    bool thermal_test = testThermalDiffusion();
    bool hydrology_test = testTerzaghiConsolidation();
    bool stress_test = testStressStrainConvergence();

    std::cout << "Running the tests:" << std::endl;
    
    std::cout << "Thermal Diffusion Test: " << (thermal_test ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Hydrology Test (Terzaghi): " << (hydrology_test ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Stress-Strain Convergence Test: " << (stress_test ? "PASSED" : "FAILED") << std::endl;
    
    return 0;
}
