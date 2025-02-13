#ifndef BBM_MODEL3D_H
#define BBM_MODEL3D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense> // For matrix operations

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
    Eigen::Matrix3d stress;
    Eigen::Matrix3d strain;
    double suction;
    double preconsolidation_pressure;
    double temperature;
    double water_content;
};

// Barcelona Basic Model (BBM) Class
class BarcelonaBasicModel {
    public:
        BarcelonaBasicModel(MaterialProperties props) : properties(props) {}
        
        void computeStressStrain(SoilState& state, const Eigen::Matrix3d& d_strain, double d_suction, double dt);
        void computeThermalConduction(SoilState& state, double dt, double thermal_gradient);
        void computeHydrology(SoilState& state, double dt);
        
    private:
        MaterialProperties properties;
    };

    #endif // BBM_MODEL3D_H