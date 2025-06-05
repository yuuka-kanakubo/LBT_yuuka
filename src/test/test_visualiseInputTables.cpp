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

	double isQuark = false;
	test_check_qhat(config, p, isQuark);
	test_check_Gamma_inel(config, p, isQuark);
	test_check_Gamma_el(config, p, isQuark);
	isQuark = true;
	test_check_qhat(config, p, isQuark);
	test_check_Gamma_inel(config, p, isQuark);
	test_check_Gamma_el(config, p, isQuark);


	std::cout << "All tests completed successfully.\n";
	return 0;
}
