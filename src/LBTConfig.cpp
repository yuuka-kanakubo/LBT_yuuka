#include "LBTConfig.h"
#include "LBTcl_base.h"

void LBTConfig::loadFromFile(const std::string& filename) {
	std::ifstream input(filename);
	if (!input) {
		std::cerr << "Parameter file is missing!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string line;
	while (std::getline(input, line)) {
		if (line.empty() || line[0] == '#') continue;

		char str[1024];
		std::strncpy(str, line.c_str(), sizeof(str));

		char* tokens[5];
		tokens[0] = std::strtok(str, " ");
		int nChar = 0;
		while (tokens[nChar] != nullptr) {
			if (nChar == 3) break;
			nChar++;
			tokens[nChar] = std::strtok(nullptr, " ");
		}

		if (nChar == 3) {
			std::string key(tokens[0]);
			std::string eq(tokens[1]);
			std::string val(tokens[2]);

			if (eq == "=") {
				double fval = std::stod(val);
				int ival = std::stoi(val);

				//Flags (integers)
				if (key == "flag:fixMomentum") jet.fixMomentum = ival;
				else if (key == "flag:fixPosition") jet.fixPosition = ival;
				else if (key == "flag:initHardFlag") jet.initHardFlag = ival;
				else if (key == "flag:flagJetX") jet.flagJetX = ival;
				else if (key == "flag:vacORmed") medium.vacORmed = ival;
				else if (key == "flag:bulkType") medium.bulkFlag = ival;
				else if (key == "flag:Kprimary") physics.Kprimary = ival;
				else if (key == "flag:KINT0") physics.Kinteraction = ival;
				else if (key == "flag:outFormat") output.outFormat = ival;
				else if (key == "flag:heavyOut") output.heavyOut = ival;
				else if (key == "flag:lightOut") output.lightOut = ival;

				// Parameters (double/int)
				else if (key == "para:nEvent") jet.ncall = ival;
				else if (key == "para:nParton") jet.nj = ival;
				else if (key == "para:parID") jet.Kjet = ival;
				else if (key == "para:tau0") clock.tau0 = fval;
				else if (key == "para:tauend") clock.tauend = fval;
				else if (key == "para:dtau") clock.dtau = fval;
				else if (key == "para:temperature") medium.temp0 = fval;
				else if (key == "para:alpha_s") physics.fixAlphas = fval;
				else if (key == "para:hydro_Tc") medium.hydro_Tc = fval; else if (key == "para:pT_min") jet.ipTmin = fval;
				else if (key == "para:pT_max") jet.ipTmax = fval;
				else if (key == "para:eta_cut") jet.eta_cut = fval;
				else if (key == "para:Ecut") physics.Ecut = fval;
				else if (key == "para:ener") jet.ener = fval;
				else if (key == "para:mass") jet.amss = fval;
				else if (key == "para:KPamp") lbtinput.KPamp = fval;
				else if (key == "para:KPsig") lbtinput.KPsig = fval;
				else if (key == "para:KTamp") lbtinput.KTamp = fval;
				else if (key == "para:KTsig") lbtinput.KTsig = fval;
			}
		}
	}

	// Derived values
	if (medium.vacORmed == 0) jet.fixPosition = 1; // position is irrelevant in vacuum

	compute_otherParameter();


	// Optional: echo the loaded parameters for verification
	std::cout << "\n############################################################\n";
	std::cout << "Parameters for this run\n";
	std::cout << "flag:initHardFlag: " << jet.initHardFlag << "\n";
	std::cout << "flag:flagJetX: " << jet.flagJetX << "\n";
	std::cout << "flag:fixMomentum: " << jet.fixMomentum << "\n";
	std::cout << "flag:fixPosition: " << jet.fixPosition << "\n";
	std::cout << "flag:vacORmed: " << medium.vacORmed << "\n";
	std::cout << "flag:bulkType: " << medium.bulkFlag << "\n";
	std::cout << "flag:Kprimary: " << physics.Kprimary << "\n";
	std::cout << "flag:KINT0: " << physics.Kinteraction << "\n";
	std::cout << "flag:outFormat: " << output.outFormat << "\n";
	std::cout << "flag:heavyOut: " << output.heavyOut << "\n";
	std::cout << "flag:lightOut: " << output.lightOut << "\n";
	std::cout << "para:nEvent: " << jet.ncall << "\n";
	std::cout << "para:nParton: " << jet.nj << "\n";
	std::cout << "para:parID: " << jet.Kjet << "\n";
	std::cout << "para:tau0: " << clock.tau0 << "\n";
	std::cout << "para:tauend: " << clock.tauend << "\n";
	std::cout << "para:dtau: " << clock.dtau << "\n";
	std::cout << "para:temperature: " << medium.temp0 << "\n";
	std::cout << "para:alpha_s: " << physics.fixAlphas << "\n";
	std::cout << "para:hydro_Tc: " << medium.hydro_Tc << "\n";
	std::cout << "para:pT_min: " << jet.ipTmin << "\n";
	std::cout << "para:pT_max: " << jet.ipTmax << "\n";
	std::cout << "para:eta_cut: " << jet.eta_cut << "\n";
	std::cout << "para:Ecut: " << physics.Ecut << "\n";
	std::cout << "para:ener: " << jet.ener << "\n";
	std::cout << "para:mass: " << jet.amss << "\n";
	std::cout << "para:KPamp: " << lbtinput.KPamp << "\n";
	std::cout << "para:KPsig: " << lbtinput.KPsig << "\n";
	std::cout << "para:KTamp: " << lbtinput.KTamp << "\n";
	std::cout << "para:KTsig: " << lbtinput.KTsig << "\n";
	std::cout << "############################################################\n" << std::endl;
}

int LBTConfig::checkParameter(int nArg) const {
	int ctErr = 0;
std::cout << "DEBUG: initHardFlag=" << jet.initHardFlag
          << " lightOut=" << output.lightOut
          << " heavyOut=" << output.heavyOut << "\n";


	if (jet.initHardFlag == 1) {
		if (output.lightOut == 1 && output.heavyOut == 1) {
			if (nArg != 5) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file HQ_output light_positive_output light_negative_output\n";
			}
		} else if (output.heavyOut == 1) {
			if (nArg != 3) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file HQ_output\n";
			}
		} else if (output.lightOut == 1) {
			if (nArg != 4) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file light_positive_output light_negative_output\n";
			}
		} else {
			if (nArg != 2) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file\n";
				std::cout << "No output specified, both heavyOut and lightOut are set to 0\n";
			}
		}
	} else if (jet.initHardFlag == 2) {
		if (output.lightOut == 1 && output.heavyOut == 1) {
			if (nArg != 6) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file input_parton_list HQ_output light_positive_output light_negative_output\n";
			}
		} else if (output.heavyOut == 1) {
			if (nArg != 4) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file input_parton_list HQ_output\n";
			}
		} else if (output.lightOut == 1) {
			if (nArg != 5) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file input_parton_list light_positive_output light_negative_output\n";
			}
		} else {
			if (nArg != 2) {
				ctErr++;
				std::cout << "Wrong # of arguments for command line input\n";
				std::cout << "./LBT parameter_file input_parton_list\n";
				std::cout << "No output specified, both heavyOut and lightOut are set to 0\n";
			}
		}
	} else {
		std::cout << "Wrong input for initHardFlag!" << std::endl;
		ctErr++;
	}

	return ctErr;
}

void LBTConfig::compute_otherParameter(){
	lbtinput.KTsig = lbtinput.KTsig * medium.hydro_Tc;
	lbtinput.preKT = physics.fixAlphas / 0.3;
	medium.temp00 = medium.temp0;
	clock.dt = clock.dtau;
	clock.timend = clock.tauend;
	clock.time0 = clock.tau0;
	physics.alphas = alphas0(physics.Kalphas, medium.temp0);
	physics.qhat0 = DebyeMass2(physics.Kqhat0, physics.alphas, medium.temp0);
	lbtinput.runKT = physics.fixAlphas / 0.3;  // running coupling factor
}
