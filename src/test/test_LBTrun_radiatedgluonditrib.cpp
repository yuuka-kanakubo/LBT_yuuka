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



void test_check_Gamma_inel_q(LBTConfig& config){

	//Testing LBTcl::computeRadiationProbability in LBTcl.cpp.
	std::ofstream test_ofs;
	test_ofs.open("DATA_Gamma_inel_q.dat");

	const double T_med = 0.3;
	for(double E_ini=0.1; E_ini < 200.0; E_ini+=0.5){
		Particle p;
		p.P[0] = E_ini;  p.P[1] = E_ini; p.P[2] = 0.; p.P[3] = 0.;
		p.Tfrozen = T_med;
		p.pid = 1;
		p.vcfrozen[1] = 0.0; p.vcfrozen[2] = 0.0; p.vcfrozen[3] = 0.0;
		//(t-t_i) time dulation from last interaction
		p.Tint_lrf = 1.0;


		LBTcl lbt(config);
		lbt.computeScatteringRate(p, E_ini, T_med);

		double KTfactor = 1.0 
			+ config.lbtinput.KTamp * 
			exp(-pow(0.3 - config.medium.hydro_Tc, 2) / (2.0 * config.lbtinput.KTsig * config.lbtinput.KTsig));
		double Gamma_inel_q = lbt.nHQgluon(p, p.Tfrozen, p.P[0]) * KTfactor * config.lbtinput.runKT;
                test_ofs 
			<< std::setw(18) << std::setprecision(10) << E_ini
			<< std::setw(18) << std::setprecision(10) << Gamma_inel_q
			<< std::endl;
	}
	test_ofs.close();


}



void test_check_qhat_q(LBTConfig &config){

	//Testing LBTcl::computeRadiationProbability in LBTcl.cpp.
	std::ofstream test_ofs;
	test_ofs.open("DATA_qhat_q.dat");

	const double T_med = 0.3;
	for(double E_ini=0.1; E_ini < 200.0; E_ini+=0.5){
		Particle p;
		p.P[0] = E_ini;  p.P[1] = E_ini; p.P[2] = 0.; p.P[3] = 0.;
		p.Tfrozen = T_med;
		p.pid = 1;
		p.vcfrozen[1] = 0.0; p.vcfrozen[2] = 0.0; p.vcfrozen[3] = 0.0;
		//(t-t_i) time dulation from last interaction
		p.Tint_lrf = 1.0;

		LBTcl lbt(config);
		lbt.computeScatteringRate(p, E_ini, T_med);


		double qhat = p.qhat_over_T3 * (T_med * T_med * T_med);
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
	config.lbtinput.KPamp =5.0;
	config.lbtinput.KPsig = 5.0;
	config.lbtinput.KTamp  = 0.0;
	config.medium.hydro_Tc = 0.165;
	config.lbtinput.KTsig  = 0.05;
	config.physics.fixAlphas  = 0.15;
	config.compute_otherParameter();


	//Read tables
	readTables("../../",config);

	//Analytical results exist.
	test_check_qhat_q(config);
	test_check_Gamma_inel_q(config);

	std::cout << "All tests completed successfully.\n";
	return 0;
}
