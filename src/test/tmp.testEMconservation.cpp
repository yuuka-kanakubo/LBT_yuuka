#include <iostream>
#include <cmath>
#include <array>
#include "../ParticleInfo.h"       // Assumes Particle struct is defined here
#include "../LBTcl.h"           // Your main LBT class
#include "../LBTConfig.h"       // LBT configuration (mock version or minimal initialization)

bool checkEnergyMomentumConservation(const std::array<double, 4>& p1,
                                     const std::array<double, 4>& p2,
                                     const std::array<double, 4>& p3,
                                     const std::array<double, 4>& p4,
                                     double tol = 1e-6) {
    for (int i = 0; i < 4; ++i) {
        if (std::abs((p1[i] + p2[i]) - (p3[i] + p4[i])) > tol) return false;
    }
    return true;
}

void test_energy_momentum_conservation_collHQ22() {
    Particle p;
    p.P[0] = 10.0;  p.P[1] = 1.0; p.P[2] = 2.0; p.P[3] = 3.0;
    p.Tfrozen = 0.3;
    p.pid = 4;
    p.vcfrozen[1] = 0.0; p.vcfrozen[2] = 0.0; p.vcfrozen[3] = 0.0;

    double qt = 0.0;

    LBTConfig config;  // Ensure tables and parameters are properly initialized in production
    LBTcl lbt(config);
    int channel = 11;

    std::array<double, 4> p_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
    std::array<double, 4> pc_med ={0.,0.,0.,0.}; 
    std::array<double, 4> pc_rec = {0.,0.,0.,0.}; 
    std::array<double, 4> pc_fin ={0.,0.,0.,0.}; 
    lbt.collHQ22(channel, p, pc_rec, pc_med, pc_fin, qt);

    if (!checkEnergyMomentumConservation(p_jet, pc_med, pc_fin, pc_rec)) {
        std::cerr << "[Test Failed] Energy-momentum not conserved in collHQ22" << std::endl;
    } else {
        std::cout << "[Test Passed] Energy-momentum conserved in collHQ22" << std::endl;
    }
}

void test_energy_momentum_conservation_colljet22() {
    Particle p;
    p.P[0] = 20.0; p.P[1] = 5.0; p.P[2] = 0.0; p.P[3] = 18.0;
    p.Tfrozen = 0.3;
    p.pid = 21;
    p.vcfrozen[1] = 0.0; p.vcfrozen[2] = 0.0; p.vcfrozen[3] = 0.0;

    std::array<double, 4> pc_rec, pc_med, pc_fin;
    double qt = 0.0;

    LBTConfig config;
    LBTcl lbt(config);

    int channel = 1;
    lbt.colljet22(channel, p, pc_rec, pc_med, pc_fin, qt);

    std::array<double, 4> p_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
    if (!checkEnergyMomentumConservation(p_jet, pc_med, pc_fin, pc_rec)) {
        std::cerr << "[Test Failed] Energy-momentum not conserved in colljet22" << std::endl;
    } else {
        std::cout << "[Test Passed] Energy-momentum conserved in colljet22" << std::endl;
    }
}

int main() {
    test_energy_momentum_conservation_collHQ22();
    test_energy_momentum_conservation_colljet22();
    // You can add similar functions for collHQ23 or radiation as needed
    return 0;
}

