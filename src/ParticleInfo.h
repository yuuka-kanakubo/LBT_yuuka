#ifndef PARTICLE_H
#define PARTICLE_H
#include <iomanip>
#include <iostream>
#include "LBTcl_base.h"


//Particle createChildParticle(const Particle& parent, const int newpid) {
//    Particle child;
//
//    child.V = parent.V;
//    child.vcfrozen[1] = parent.vcfrozen[1];
//    child.vcfrozen[2] = parent.vcfrozen[2];
//    child.vcfrozen[3] = parent.vcfrozen[3];
//
//    child.Tfrozen = parent.Tfrozen;
//    child.mass = parent.mass;  // update later if needed
//    child.WT = parent.WT;
//    child.isPrimary = false;
//    child.isActive = true;
//
//    child.pid = newpid;
//    child.CAT = 2;  // recoiled
//    child.mom1 = parent.index;
//    child.mom2 = parent.index;
//
//    return child;
//}
//

enum class channel {
      gg2gg, //case 1
      qg2qg, //case 2
      gq2gq, //case 3
      qg2gq, //case 13
      qq2qq, //case 4
      qqbar2qqbar, //case 5
      qqbar2qqbar_flavor_exchange, //case 6
      qqbar2qqbar_flavor_kept, //case 7
      HF_channel11, //case 11
      HF_channel12, //case 12
      qqbar2gg //case 8
};


inline std::string channel_str(const channel cnl) {
	switch (cnl) {
		case channel::gg2gg: // g + g → g + g
			return " channel::gg2gg: // g + g → g + g " ;
		case channel::qg2qg: // q + g → q + g
			return " channel::qg2qg: // q + g → q + g " ;
		case channel::qg2gq: // duplicate of 3, still gluonic → f1
			return " channel::qg2gq: // duplicate of 3, still gluonic → f1 " ;
		case channel::qq2qq: // q + q → q + q
			return " channel::qq2qq: // q + q → q + q " ;
		case channel::qqbar2qqbar: // q + q̄ → q + q̄
			return " channel::qqbar2qqbar: // q + q̄ → q + q̄ " ;
		case channel::qqbar2qqbar_flavor_exchange: // q + q̄ → q + q̄ (flavor-exchange)
			return " channel::qqbar2qqbar_flavor_exchange: // q + q̄ → q + q̄ (flavor-exchange) " ;
		case channel::qqbar2qqbar_flavor_kept: // q + q̄ → q + q̄ (same-flavor)
			return " channel::qqbar2qqbar_flavor_kept: // q + q̄ → q + q̄ (same-flavor) " ;
		case channel::qqbar2gg: // q + q̄ → g + g
			return " channel::qqbar2gg: // q + q̄ → g + g " ;
		case channel::gq2gq: // g + q → g + q 
			return " channel::gq2gq: // g + q → g + q  " ;
		default:
			std::cerr << "ERROR: no other channel available." << std::endl;
			exit(EXIT_FAILURE);
	}
}



struct Particle {
        double P[6] = {0.0};        // Four-momentum (E, px, py, pz, m, pt)
        double V[4] = {0.0};        // Position (t or tau, x, y, z)
        int CAT = -1;                // Category: 0=primary, 1=free streaming, 2=recoiled, 3=medium(negative), 4=radiated
        int pid = 0;               // Flavor code
        double Vfrozen[4] = {0.0};  // Frozen position (for interpolation)
        double Tfrozen = 0.0;       // Temperature when frozen
        double vcfrozen[4] = {0.0}; // Flow velocity when frozen
        double Tint_lrf = 0.0;      // Interaction time in local rest frame
        double radng = 0.0;         // Average radiated gluons
        double WT = 1.0;            // Weight
        double mass = 0.0;          // Mass (0 for light partons, positive for heavy quarks)
        bool isLeading = false;     // Particles which may collide with medium. (instead of nj) or, recoiled, radiated/
        bool isActive = false;     // True if parton exist, and not a past history.
        double tot_el_rate = 0.0;     // Total elastic scattering rate (2->2)
        double Xtau_keep = 0.0;     // Xtau calculated by titau before propagation.
        double max_Ng = 0.0;     // Maximum radiation rate. Used for efficient sampling.
        int parent1 = -1;  // Index of primary parent parton
        int parent2 = -1;  // Index of secondary (e.g. recoil) parent parton
        int kid1 = -1;  // Index of kid parton
        int kid2 = -1;  // Index of secondary (e.g. recoil) kid parton
        int kid3 = -1;  
	double timedilation = 0.0;//time dilation factor between lab frame and fluid rest frame.
	double qhat_over_T3 = 0.0;
	double D2piT = 0.0;//dimention-less diffusion factor

        void copy_thisV_to_Vfrozen(const double Vin[4]){
             for(int i =0; i<4; i++) this->Vfrozen[i] = Vin[i];
        };
        void copy_thisV_to_vcfrozen(const double Vin[4]){
             for(int i =0; i<4; i++) this->vcfrozen[i] = Vin[i];
        };


	void get_timedilation(){
		double vMag =  sqrt(pow(vcfrozen[1],2) + pow(vcfrozen[2],2) + pow(vcfrozen[3],2));
		this->timedilation =  (1.0-(P[1]*vcfrozen[1]+P[2]*vcfrozen[2]+P[3]*vcfrozen[3])/P[0])/sqrt(1.0-vMag*vMag);
	};

	void get_D2piT(const double qhat_over_T3){
		if(pid==21) this->D2piT=8.0*base::pi/(qhat_over_T3/base::CA*base::CF); 
		else this->D2piT=8.0*base::pi/qhat_over_T3;
	}

	void Print(const bool printonlyActive){
		if(printonlyActive==true && this->isActive==false) {
		//std::cout << "------------------ this particle is not actuve: index " << this->index_ << std::endl;
			return;
		}
		if(printonlyActive==true && this->CAT==3) {
			return;
		}
		std::cout << "----" << std::endl;
		std::cout 
			<< "[ index  "
			<< std::setw(6)  << this->index_
			<< "],  [parent1 / parent2  "
			<< std::setw(6)  << this->parent1 << "/" << std::setw(6) << this->parent2
			<< "],  [kid1 / kid2  "
			<< std::setw(6)  << this->kid1 << "/" << std::setw(6) << this->kid2
			<< "],  [pid  "
			<< std::setw(6)  << this->pid
			<< "]   [CAT 1:freestrm, 2:recoiled, 3:medium(negative), 4:radiated"
			<< std::setw(6)  << this->CAT
			<< "] [isActive / isLeading "
			<< std::setw(6)  << this->isActive << "/" << std::setw(6) << this->isLeading
			<< "]"
			<< std::endl
			<< " P  "
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[0]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[1]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[2]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[3]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[4]
			<< std::endl
			<< " V  "
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
	void PrintMed(){
		if(this->CAT!=3) {
			return;
		}
		std::cout << "----" << std::endl;
		std::cout 
			<< "[ index  "
			<< std::setw(6)  << this->index_
			<< "],  [parent1 / parent2  "
			<< std::setw(6)  << this->parent1 << "/" << std::setw(6) << this->parent2
			<< "],  [kid1 / kid2  "
			<< std::setw(6)  << this->kid1 << "/" << std::setw(6) << this->kid2
			<< "],  [pid  "
			<< std::setw(6)  << this->pid
			<< "]   [CAT 1:freestrm, 2:recoiled, 3:medium(negative), 4:radiated"
			<< std::setw(6)  << this->CAT
			<< "] [isActive / isLeading "
			<< std::setw(6)  << this->isActive << "/" << std::setw(6) << this->isLeading
			<< "]"
			<< std::endl
			<< " P  "
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[0]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[1]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[2]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[3]
			<< std::setw(15) << std::setprecision(8) << std::fixed <<  this->P[4]
			<< std::endl
			<< " V  "
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

	void assign_index(){ 
		this->index_ = index_counter; 
                index_counter ++;
	}

        int index() const{return this->index_;}

	private:
	int index_ = -1;  // Unique ID for this particle in the vector

};
#endif
