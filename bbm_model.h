#ifndef BBM_MODEL_H
#define BBM_MODEL_H

#include <iostream>
#include <vector>
#include <cmath>

// Constants and material properties
struct MaterialProperties {
    double lambda;
    double kappa;
    double lambda_s;
    double s0;
    double M;
    double thermal_conductivity;
    double specific_heat;
    double permeability;
    double alpha; // Nonlinear permeability coefficient
};

// State variables
struct SoilState {
    double mean_stress;
    double deviatoric_stress;
    double suction;
    double preconsolidation_pressure;
    double temperature;
    double water_content;
};

// Barcelona Basic Model (BBM) Class
class BarcelonaBasicModel {
public:
    BarcelonaBasicModel(MaterialProperties props);
    
    void computeStressStrain(SoilState& state, double d_eps_v, double d_suction, double dt);
    void computeThermalConduction(SoilState& state, double dt, double thermal_gradient);
    void computeHydrology(SoilState& state, double dt);
    
private:
    MaterialProperties properties;
};

#endif // BBM_MODEL_H
