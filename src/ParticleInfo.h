#ifndef PARTICLE_H
#define PARTICLE_H
#include <iomanip>
#include <iostream>

struct Particle {
        double P[6] = {0.0};        // Four-momentum (E, px, py, pz)
        double V[4] = {0.0};        // Position (t or tau, x, y, z)
        int CAT = 0;                // Category: 0=active, 1=free streaming, 2=recoiled, 4=radiated, 9=negative
        int KATT = 0;               // Flavor code
        double Vfrozen[4] = {0.0};  // Frozen position (for interpolation)
        double Tfrozen = 0.0;       // Temperature when frozen
        double vcfrozen[4] = {0.0}; // Flow velocity when frozen
        double Tint_lrf = 0.0;      // Interaction time in local rest frame
        double radng = 0.0;         // Average radiated gluons
        double WT = 1.0;            // Weight
        double mass = 0.0;          // Mass (0 for light partons, positive for heavy quarks)
        bool isPrimary = false;     // True if parton is among the original initial partons (instead of nj)
        bool isActive = false;     // True if parton exist, and not a past history.
        int mom1 = -1;  // Index of primary mother parton
        int mom2 = -1;  // Index of secondary (e.g. recoil) mother parton
        int index = -1;  // Unique ID for this particle in the vector

	void Print(){
		std::cout 
			<< std::setw(6)  << this->KATT
			<< std::setw(6)  << this->CAT
			<< std::endl
			<< "P  "
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[0]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[1]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[2]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[3]
			<< std::endl
			<< "V  "
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->V[0]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->V[1]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->V[2]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->V[3]
			<< std::endl
			<< "Vfrozen  "
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->Vfrozen[0] 
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->Vfrozen[1] 
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->Vfrozen[2] 
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->Vfrozen[3] 
			<< std::endl
			<< "vcfrozen  "
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->vcfrozen[0] 
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->vcfrozen[1] 
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->vcfrozen[2] 
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->vcfrozen[3] 
			<< std::endl;

	}

};
#endif
