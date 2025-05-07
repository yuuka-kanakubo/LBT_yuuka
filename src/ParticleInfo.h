#ifndef PARTICLE_H
#define PARTICLE_H
#include <iomanip>
#include <iostream>
#include "LBTcl_base.h"

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
        double tot_el_rate = 0.0;     // Total elastic scattering rate (2->2)
        double Xtau_keep = 0.0;     // Xtau calculated by titau before propagation.
        int mom1 = -1;  // Index of primary mother parton
        int mom2 = -1;  // Index of secondary (e.g. recoil) mother parton
        int index = -1;  // Unique ID for this particle in the vector
	double timedilation = 0.0;//time dilation factor between lab frame and fluid rest frame.
	double D2piT = 0.0;//dimention-less diffusion factor
	void get_timedilation(){
		double vMag =  sqrt(pow(vcfrozen[1],2) + pow(vcfrozen[2],2) + pow(vcfrozen[3],2));
		this->timedilation =  (1.0-(P[1]*vcfrozen[1]+P[2]*vcfrozen[2]+P[3]*vcfrozen[3])/P[0])/sqrt(1.0-vMag*vMag);
	};

	void get_D2piT(const double qhat_over_T3){
		if(KATT==21) this->D2piT=8.0*base::pi/(qhat_over_T3/base::CA*base::CF); 
		else this->D2piT=8.0*base::pi/qhat_over_T3;
	}

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
