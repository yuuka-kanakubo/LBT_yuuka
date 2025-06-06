#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "../main_functions.h"           // Your main LBT class
#include "../LBTConfig.h"
#include "../ParticleInfo.h"       // Assumes Particle struct is defined here
#include "../LBTcl.h"           // Your main LBT class
#include "macro_cmn.h"



void test_check_Gamma_el(LBTConfig& config, Particle& p, const bool isQuark){

	//Testing LBTcl::computeRadiationProbability in LBTcl.cpp.
	std::ofstream test_ofs;

	if(isQuark){ 
		p.pid=1;
		test_ofs.open("DATA_Gamma_el_q.dat");
	}else{
		p.pid=21;
		test_ofs.open("DATA_Gamma_el_g.dat");
	}
	for(double E_ini=0.1; E_ini < 200.0; E_ini+=0.5){

		p.P[0] = E_ini;  p.P[1] = E_ini; p.P[2] = 0.; p.P[3] = 0.;

		LBTcl lbt(config);
		lbt.computeScatteringRate(p, E_ini, p.Tfrozen);

		//This is Gamma_el [GeV]
		// ==> converting the unit to fm^-1
		double Gamma_el = p.tot_el_rate/0.1973;
		test_ofs 
			<< std::setw(18) << std::setprecision(10) << E_ini
			<< std::setw(18) << std::setprecision(10) << Gamma_el
			<< std::endl;
	}
	test_ofs.close();


}

void test_check_Gamma_inel(LBTConfig& config, Particle& p, const bool isQuark){

	//Testing LBTcl::computeRadiationProbability in LBTcl.cpp.
	std::ofstream test_ofs;

	if(isQuark){ 
		p.pid=1;
		test_ofs.open("DATA_Gamma_inel_q.dat");
	}else{
		p.pid=21;
		test_ofs.open("DATA_Gamma_inel_g.dat");
	}
	for(double E_ini=0.1; E_ini < 1000.0; E_ini+=0.5){

		p.P[0] = E_ini;  p.P[1] = E_ini; p.P[2] = 0.; p.P[3] = 0.;

		LBTcl lbt(config);
		lbt.computeScatteringRate(p, E_ini, p.Tfrozen);

		//This is Gamma_inel [fm^-1]
		double Gamma_inel = lbt.nHQgluon(p, p.Tfrozen, p.P[0]);
		test_ofs 
			<< std::setw(18) << std::setprecision(10) << E_ini
			<< std::setw(18) << std::setprecision(10) << Gamma_inel
			<< std::endl;
	}
	test_ofs.close();


}



void test_check_qhat(LBTConfig &config, Particle& p, const bool isQuark){

	//Testing LBTcl::computeRadiationProbability in LBTcl.cpp.
	std::ofstream test_ofs;

	if(isQuark){ 
		p.pid=1;
		test_ofs.open("DATA_qhat_q.dat");
	}else{
		p.pid=21;
		test_ofs.open("DATA_qhat_g.dat");
	}
	for(double E_ini=0.1; E_ini < 200.0; E_ini+=0.5){

		p.P[0] = E_ini;  p.P[1] = E_ini; p.P[2] = 0.; p.P[3] = 0.;

		LBTcl lbt(config);
		lbt.computeScatteringRate(p, E_ini, p.Tfrozen);


		double qhat = p.qhat_over_T3 * (p.Tfrozen * p.Tfrozen * p.Tfrozen);
		test_ofs 
			<< std::setw(18) << std::setprecision(10) << E_ini
			<< std::setw(18) << std::setprecision(10) << qhat
			<< std::endl;
	}
	test_ofs.close();



}


void test_check_channel_probability(LBTConfig &config, Particle& p){

	LBTcl lbt(config);
	std::array <int, 2> pid_set = {1, 21};
	for(size_t j=0; j<2; j++){

		int pid = pid_set[j];
		std::array <double, 3> E_set = {10.0, 100.0, 1000.0};

		for(size_t i=0; i<3; i++){
			double E = E_set[i];
			std::ofstream test_ofs;
			if (j==0) test_ofs.open("DATA_channel_"+std::to_string((int)E)+"GeVqJET_vsT.dat");
			else if (j==1) test_ofs.open("DATA_channel_"+std::to_string((int)E)+"GeVgJET_vsT.dat");

			for(double T = 0.30; T<0.80; T+=0.025){

				//Calculating total elastic collision rate
				double RTE;
				{
					p.P[0] = E;  p.P[1] = E; p.P[2] = 0.; p.P[3] = 0.;
					p.Tfrozen = T;
					p.pid = pid;

					lbt.computeScatteringRate(p, E, p.Tfrozen);
					RTE=p.tot_el_rate;
				}


				// Scattering rate components
				double RTEg1, RTEg2, RTEg3;
				double RTEq, RTEq3, RTEq4, RTEq5, RTEq6, RTEq7, RTEq8;
				double RTEHQ11, RTEHQ12;
				lbt.linear(p.pid, E, T,
						RTEg1, RTEg2, RTEg3,
						RTEq, RTEq3, RTEq4, RTEq5, RTEq6, RTEq7, RTEq8,
						RTEHQ11, RTEHQ12);
				if(pid==1){
					double R0 = RTE;
					double R3 = RTEq3, R4 = RTEq4, R5 = RTEq5, R6 = RTEq6, R7 = RTEq7;

					test_ofs 
						<< std::setw(18) << std::setprecision(10) << T
						<< std::setw(18) << std::setprecision(10) << R3 / R0
						<< std::setw(18) << std::setprecision(10) << (R4) / R0
						<< std::setw(18) << std::setprecision(10) << (R5) / R0
						<< std::setw(18) << std::setprecision(10) << (R6) / R0
						<< std::setw(18) << std::setprecision(10) << (R7) / R0
						<< std::setw(18) << std::setprecision(10) << 1-((R3 + R4 + R5 + R6 + R7) / R0)
						<< std::endl;
				}else if(pid==21){
					double R0 = RTE;

					test_ofs 
						<< std::setw(18) << std::setprecision(10) << T
						<< std::setw(18) << std::setprecision(10) << RTEg1 / R0
						<< std::setw(18) << std::setprecision(10) << RTEg2 / R0
						<< std::setw(18) << std::setprecision(10) << 1-((RTEg1 + RTEg2) / R0)
						<< std::endl;
				}

			}

			test_ofs.close();
		}
	}

}







int main() {

	//Settings
	LBTConfig config;
	config.medium.bulkFlag = 2;
	config.jet.fixPosition = 1;
	config.lbtinput.KPamp =0.0;
	config.lbtinput.KPsig = 0.0;
	config.lbtinput.KTamp  = 0.0;
	config.medium.hydro_Tc = 0.165;
	config.lbtinput.KTsig  = 0.0;
	config.physics.fixAlphas  = 0.15;
	config.compute_otherParameter();


	//Read tables
	readTables("../../",config);

	//Analytical results exist.
	const double T_med = 0.3;
	Particle p;
	p.Tfrozen = T_med;
	p.vcfrozen[1] = 0.0; p.vcfrozen[2] = 0.0; p.vcfrozen[3] = 0.0;
	//(t-t_i) time dulation from last interaction
	p.Tint_lrf = 1.0;

	//1. (In)Elastic Scattering  rate
	//=========================
	double isQuark;
	isQuark = false;
	test_check_qhat(config, p, isQuark);
	test_check_Gamma_inel(config, p, isQuark);
	test_check_Gamma_el(config, p, isQuark);
	isQuark = true;
	test_check_qhat(config, p, isQuark);
	test_check_Gamma_inel(config, p, isQuark);
	test_check_Gamma_el(config, p, isQuark);

	//2. Channel
	//=========================
	test_check_channel_probability(config, p);


	std::cout << "All tests completed successfully.\n";
	return 0;
}
