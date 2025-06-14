#include "main_functions.h"
#include "LBTcl.h"

// Declare external functions
void jetInitialize(){
	std::cout << "Currently this is not implemented in the refactored version of LBT." << std::endl;
	std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
	exit(EXIT_FAILURE);
}
void setJetX(int numInitXY){
	std::cout << "Currently this is not implemented in the refactored version of LBT." << std::endl;
	std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
	exit(EXIT_FAILURE);

}



int read_xyMC(LBTConfig &config) {

	int numXY=0;

	std::ifstream fxyMC("hydroProfile/mc_glauber.dat");
	if(!fxyMC.is_open()) {
		std::cout<<"Erro openning date fxyMC!\n";
		exit (EXIT_FAILURE);
	}

	while(true) {
		if(fxyMC.eof()) {
			numXY--;
			break;
		}
		if(numXY>=config.flow.maxMC) break;
		fxyMC >> config.flow.initMCX[numXY] >> config.flow.initMCY[numXY];
		numXY++;
	}

	std::cout << "Number of (x,y) points from MC-Glauber: " << numXY << std::endl;

	fxyMC.close();
	return numXY;

}




void readTables(const std::string w_path, LBTConfig &config) {
    // Read hydro profile depending on bulkFlag
    if (config.medium.bulkFlag == 1) {
        std::string dataFN_in = w_path + "hydroProfile/JetData.dat";
        std::string ctlFN_in = w_path + "hydroProfile/JetCtl.dat";
        std::cout << "ERROR: not refactored yet. " << std::endl;
        exit(EXIT_FAILURE);
    } else if (config.medium.bulkFlag == 2) {
        std::string dataFN_in = w_path + "hydroProfile/bulk.dat";
	char buffer[256];
	strncpy(buffer, dataFN_in.c_str(), sizeof(buffer));
	buffer[sizeof(buffer)-1] = '\0'; // null-terminate

	int len1 = std::strlen(buffer);
	read_ccnu_(buffer, len1);
    }

    if (config.jet.fixPosition != 1) {
        int dummy = read_xyMC(config);
        std::cout << "ERROR: not refactored yet. " << dummy << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read scattering rate table
    std::ifstream f1(w_path + "tables/ratedata");
    int n = 450;
    int it, ie;

    if (!f1.is_open()) {
        std::cerr << "Error opening " << w_path<<  " tables/ratedata!\n";
    } else {
        for (int i = 1; i <= n; ++i) {
            f1 >> it >> ie;
            f1 >> config.tables.qhatG[it][ie]
               >> config.tables.Rg[it][ie]
               >> config.tables.Rg1[it][ie]
               >> config.tables.Rg2[it][ie]
               >> config.tables.Rg3[it][ie]
               >> config.tables.qhatLQ[it][ie]
               >> config.tables.Rq[it][ie]
               >> config.tables.Rq3[it][ie]
               >> config.tables.Rq4[it][ie]
               >> config.tables.Rq5[it][ie]
               >> config.tables.Rq6[it][ie]
               >> config.tables.Rq7[it][ie]
               >> config.tables.Rq8[it][ie];
        }
        f1.close();
    }

    // Read HQ elastic rate data
    std::ifstream f11(w_path + "b-tables/ratedata-HQ");
    if (!f11.is_open()) {
        std::cerr << "Error opening b-tables/ratedata-HQ!\n";
    } else {
        for (int i = 1; i <= n; ++i) {
            f11 >> it >> ie;
            f11 >> config.tables.RHQ[it][ie]
                >> config.tables.RHQ11[it][ie]
                >> config.tables.RHQ12[it][ie]
                >> config.tables.qhatHQ[it][ie];
        }
        f11.close();
    }

    // Read radiation tables for HQ and gluons
    std::ifstream f12(w_path + "b-tables/dNg_over_dt_bD6.dat");
    std::ifstream f13(w_path + "tables/dNg_over_dt_qD6.dat");
    std::ifstream f14(w_path + "tables/dNg_over_dt_gD6.dat");
    if (!f12.is_open() || !f13.is_open() || !f14.is_open()) {
        std::cerr << "Error opening dNg_over_dt tables!\n";
    } else {
        for (int k = 1; k <= config.hqrad.t_gn; ++k) {
		std::string dummy;
		std::streampos pos = f14.tellg();
		std::string headerLine;
		std::getline(f12 >> std::ws, dummy);
		std::getline(f13 >> std::ws, dummy);

		std::getline(f14 >> std::ws, headerLine);
		if (headerLine.empty() || headerLine.find("timestep:") == std::string::npos) {
			std::cerr << "[DEBUG] Unexpected line at k=" << k << ": '" << headerLine << "'\n";
			std::cerr << "[DEBUG] Seek pos before line: " << pos << "\n";
			// Optional: dump the next few lines manually
			for (int i = 0; i < 3; ++i) {
				std::string l;
				std::getline(f14, l);
				std::cerr << "[LOOKAHEAD] " << l << "\n";
			}
			std::exit(EXIT_FAILURE);
		}

            for (int i = 1; i <= config.hqrad.temp_gn; ++i) {
                config.hqrad.dNg_over_dt_c[k + 1][i][0] = 0.0;
                config.hqrad.dNg_over_dt_q[k + 1][i][0] = 0.0;
                config.hqrad.dNg_over_dt_g[k + 1][i][0] = 0.0;
                config.hqrad.max_dNgfnc_c[k + 1][i][0] = 0.0;
                config.hqrad.max_dNgfnc_q[k + 1][i][0] = 0.0;
                config.hqrad.max_dNgfnc_g[k + 1][i][0] = 0.0;

                for (int j = 1; j <= config.hqrad.HQener_gn; ++j) {
                    f12 >> config.hqrad.dNg_over_dt_c[k + 1][i][j] >> config.hqrad.max_dNgfnc_c[k + 1][i][j];
                    f13 >> config.hqrad.dNg_over_dt_q[k + 1][i][j] >> config.hqrad.max_dNgfnc_q[k + 1][i][j];
		    f14 >> config.hqrad.dNg_over_dt_g[k + 1][i][j] >> config.hqrad.max_dNgfnc_g[k + 1][i][j];
		    if( k == 1 && i == config.hqrad.temp_gn && j == config.hqrad.HQener_gn)
			    std::cout << " config.hqrad.dNg_over_dt_g[k + 1][i][j] " << config.hqrad.dNg_over_dt_g[k + 1][i][j]  << std::endl;
		    if( k == config.hqrad.t_gn - 1 && i == config.hqrad.temp_gn && j == config.hqrad.HQener_gn)
			    std::cout << " config.hqrad.dNg_over_dt_g[k + 1][i][j] " << config.hqrad.dNg_over_dt_g[k + 1][i][j]  << std::endl;
		}
            }
	}
        f12.close();
        f13.close();
        f14.close();
    }

    // Zero out time slice 1
    for (int i = 1; i <= config.hqrad.temp_gn; ++i) {
        for (int j = 1; j <= config.hqrad.HQener_gn; ++j) {
            config.hqrad.dNg_over_dt_c[1][i][j] = 0.0;
            config.hqrad.dNg_over_dt_q[1][i][j] = 0.0;
            config.hqrad.dNg_over_dt_g[1][i][j] = 0.0;
            config.hqrad.max_dNgfnc_c[1][i][j] = 0.0;
            config.hqrad.max_dNgfnc_q[1][i][j] = 0.0;
            config.hqrad.max_dNgfnc_g[1][i][j] = 0.0;
        }
    }




     std::ifstream fileB(w_path + "b-tables/distB.dat");
     if(!fileB.is_open()) {
        std::cout << "Erro openning data file distB.dat!" << std::endl;
     } else {
        for(int i=0;i<config.hq22.N_T;i++) {
           for(int j=0;j<config.hq22.N_p1;j++) {
              double dummy_T,dummy_p1;
              fileB>>dummy_T>>dummy_p1;
              if(fabs(config.hq22.min_T+(0.5+i)*config.hq22.bin_T-dummy_T)>1.0e-5 || fabs(config.hq22.min_p1+(0.5+j)*config.hq22.bin_p1-dummy_p1)>1.0e-5) {
                  std::cout << "Erro in reading data file distB.dat!" << std::endl;
                  exit (EXIT_FAILURE);
              }
              fileB>>config.hq22.distFncBM[i][j];
              for(int k=0;k<config.hq22.N_e2;k++) fileB>>config.hq22.distFncB[i][j][k];
              for(int k=0;k<config.hq22.N_e2;k++) fileB>>config.hq22.distMaxB[i][j][k];
           }
        }
     }
     fileB.close();

     std::ifstream fileF(w_path + "b-tables/distF.dat");
     if(!fileF.is_open()) {
        std::cout << "Erro openning data file distF.dat!" << std::endl;
     } else {
        for(int i=0;i<config.hq22.N_T;i++) {
           for(int j=0;j<config.hq22.N_p1;j++) {
              double dummy_T,dummy_p1;
              fileF>>dummy_T>>dummy_p1;
              if(fabs(config.hq22.min_T+(0.5+i)*config.hq22.bin_T-dummy_T)>1.0e-5 || fabs(config.hq22.min_p1+(0.5+j)*config.hq22.bin_p1-dummy_p1)>1.0e-5) {
                  std::cout << "Erro in reading data file distF.dat!" << std::endl;
                  exit (EXIT_FAILURE);
              }
              fileF>>config.hq22.distFncFM[i][j];
              for(int k=0;k<config.hq22.N_e2;k++) fileF>>config.hq22.distFncF[i][j][k];
              for(int k=0;k<config.hq22.N_e2;k++) fileF>>config.hq22.distMaxF[i][j][k];
           }
        }
     }
     fileF.close();


}






bool initialize(LBTConfig& config, int argc, char* argv[], std::time_t& time_start) {
	config.loadFromFile(argv[1]);
	if (config.checkParameter(argc) != 0) {
		std::cerr << "Parameter check failed" << std::endl;
		return false;
	}

	if (config.medium.vacORmed == 1)
		readTables("./", config);


	// Random number initialization
	srand(123); // or: srand((unsigned)time(NULL));
	config.rng.NUM1 = -1 * rand();

	// Print start time
	time_start = std::time(nullptr);
	std::cout << "Program starts at: " << std::asctime(std::localtime(&time_start));

	return true;
}


void open_output(int argc, char* argv[],
		const LBTConfig& config,
		std::ifstream& fpList,
		std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg) {
	const int initHardFlag = config.jet.initHardFlag;
	const int heavyOut = config.output.heavyOut;
	const int lightOut = config.output.lightOut;

	if (initHardFlag == 1) {
		if (argc < 3) {
			std::cerr << "Error: Too few arguments for initHardFlag=1" << std::endl;
			std::exit(EXIT_FAILURE);
		}

		if (heavyOut && lightOut) {
			outHQ.open(argv[2]);
			outLightPos.open(argv[3]);
			outLightNeg.open(argv[4]);
		} else if (heavyOut) {
			outHQ.open(argv[2]);
		} else if (lightOut) {
			outLightPos.open(argv[2]);
			outLightNeg.open(argv[3]);
		}
	} else {
		if (argc < 4) {
			std::cerr << "Error: Too few arguments for initHardFlag=2" << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//Open input parton list!
		fpList.open(argv[2]);
		if (!fpList) {
			std::cerr << "Error: Unable to open input parton list: " << argv[2] << std::endl;
			std::exit(EXIT_FAILURE);
		}

		if (heavyOut && lightOut) {
			outHQ.open(argv[3]);
			outLightPos.open(argv[4]);
			outLightNeg.open(argv[5]);
		} else if (heavyOut) {
			outHQ.open(argv[3]);
		} else if (lightOut) {
			outLightPos.open(argv[3]);
			outLightNeg.open(argv[4]);
		}
	}

	//Headers
	if(heavyOut && outHQ.is_open()){

		outHQ << "# Heavy quark output ----" << std::endl;
		outHQ << "#"
			<< std::setw(5)  << "i"
			<< std::setw(5)  << "index"
			<< std::setw(5)  << "dmy"
			<< std::setw(5)  << "dmy"
			<< std::setw(5)  << "pid"
			<< std::setw(15) << "mass"
			<< std::setw(15) << "e"
			<< std::setw(15) << "px"
			<< std::setw(15) << "py"
			<< std::setw(15) << "pz"
			<< std::setw(15) << "rap"
			<< std::setw(15) << "x"
			<< std::setw(15) << "y"
			<< std::setw(15) << "z"
			<< std::setw(15) << "t"
			<< std::setw(15) << "ft(=t)"
			<< std::setw(5) << "CAT"
			<< std::endl;
	}
	if(lightOut && outLightPos.is_open()){
		outLightPos << "# Positive partons output ----" << std::endl;
		outLightPos << "#"
			<< std::setw(5)  << "i"
			<< std::setw(5)  << "index"
			<< std::setw(5)  << "dmy"
			<< std::setw(5)  << "dmy"
			<< std::setw(5)  << "pid"
			<< std::setw(15) << "mass"
			<< std::setw(15) << "e"
			<< std::setw(15) << "px"
			<< std::setw(15) << "py"
			<< std::setw(15) << "pz"
			<< std::setw(15) << "rap"
			<< std::setw(15) << "x"
			<< std::setw(15) << "y"
			<< std::setw(15) << "z"
			<< std::setw(15) << "t"
			<< std::setw(15) << "ft(=t)"
			<< std::setw(5) << "CAT"
			<< std::endl;
	}
	if(lightOut && outLightNeg.is_open()){

		outLightNeg << "# Negative partons output ----" << std::endl;
		outLightNeg << "#"
			<< std::setw(5)  << "i"
			<< std::setw(5)  << "index"
			<< std::setw(5)  << "dmy"
			<< std::setw(5)  << "dmy"
			<< std::setw(5)  << "pid"
			<< std::setw(15) << "mass"
			<< std::setw(15) << "e"
			<< std::setw(15) << "px"
			<< std::setw(15) << "py"
			<< std::setw(15) << "pz"
			<< std::setw(15) << "rap"
			<< std::setw(15) << "x"
			<< std::setw(15) << "y"
			<< std::setw(15) << "z"
			<< std::setw(15) << "t"
			<< std::setw(15) << "ft(=t)"
			<< std::setw(5) << "CAT"
			<< std::endl;
	}

	// Optional: check if file streams are valid
	if ((heavyOut && !outHQ) ||
			(lightOut && (!outLightPos || !outLightNeg))) {
		std::cerr << "Error: Failed to open one or more output files." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	return;
}


void get_set_ready(std::vector <Particle>& part_event, LBTConfig& config){

	for (auto it = part_event.begin(); it != part_event.end(); ++it) {
		for(int dim=1; dim<=3; dim++){
			it->Vfrozen[dim]=it->Vfrozen[dim]+(it->P[dim]/it->P[0])*it->Vfrozen[0];
			it->V[dim]=it->Vfrozen[dim];
		}
		it->V[0]=-std::log(1.0-ran0(&config.rng.NUM1));
		it->Vfrozen[3]=0.0;
		it->V[3]=0.0;

                // adjust momentum to fit energy and mass
		if(abs(it->pid)==1||abs(it->pid)==2||abs(it->pid)==3||abs(it->pid)==21){
			it->P[4]=0.0;
			it->P[0]=sqrt(it->P[1]*it->P[1]+it->P[2]*it->P[2]+it->P[3]*it->P[3]+it->P[4]*it->P[4]);
			it->P[5]=sqrt(it->P[1]*it->P[1]+it->P[2]*it->P[2]);//transverse momentum
			it->WT=1.0;
		}
	}
	return;
}



void runLBT(std::ifstream& fpList,
		LBTConfig& config,
		std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg
	   ) {

	int& nj = config.jet.nj;
	int& np = config.jet.np;
	const int initHardFlag = config.jet.initHardFlag;
	int numEvent = 0;

	LBTcl lbtcl(config);

	for (int n = 1; ; ++n) {
		std::vector<Particle> partons_event;

		// Stop conditions
		if (initHardFlag == 1 && numEvent >= config.jet.ncall) break;
		if (initHardFlag == 2 && (numEvent >= config.jet.ncall || fpList.eof())) break;

		if (initHardFlag == 1) {
			jetInitialize();
		} else {
			int dummyInt;

			fpList >> dummyInt >> nj;
			std::string line;
			int count = 0;
			if(!std::getline(fpList, line)) break; // consume rest of header

			for (int i = 0; i < nj; ++i) {
				if (!std::getline(fpList, line)) break;

				std::istringstream iss(line);
				Particle p;
				if (!(iss >> dummyInt >> p.pid >> p.P[1] >> p.P[2] >> p.P[3] >> p.P[0]
							>> p.Vfrozen[1] >> p.Vfrozen[2] >> p.Vfrozen[3] >> p.Vfrozen[0])) {
					std::cerr << "[Warning] Skipping malformed line: " << line << "\n";
					continue;
				}

				p.assign_index();
				p.isActive = true;
				p.isLeading = true;
				p.CAT = 0;
				partons_event.push_back(p);
				++count;
			}

			if (nj != count) {
				std::cerr << "ERROR: Input format mismatch in initial parton file." << std::endl;
				std::cerr << "       nj!=count  " << nj << " != " << count << std::endl;
				std::exit(EXIT_FAILURE);
			}
			count = 0;

			// NOTE: particle list is not yet stored in a central manager
			// It can be returned or handled globally
			
		}//End of reading one event info

		np = nj;

		//1. propagation of partons to formation time, 2. energy momentum adjustment.
		get_set_ready(partons_event, config);

		//CHECKING
		std::cout << __FILE__ << "(" << __LINE__ << ")" << "After reading" << std::endl;
		std::cout << "np " << np << std::endl;
		std::cout << "(int) partons_event.size() " << (int) partons_event.size() << std::endl;
		for (auto it = partons_event.begin(); it != partons_event.end(); ++it) {
			it->Print(false);
		}


		// Time evolution (skipping detailed LBT logic)

		lbtcl.set_np_snapshot(np);
		if (config.medium.vacORmed == 1) {
			for (double ti = config.clock.time0 + config.clock.dt;
					ti <= config.clock.timend + base::epsilon;
					ti += config.clock.dt) {
				lbtcl.LBT(partons_event, ti);
			}
                        std::cout << "LBT evolution completed in event# " << n << std::endl;
		}

		++numEvent;

		if (n % config.clock.nprint == 0) {
			std::cout << "n = " << n << ", np = " << np << std::endl;
		}

		// TODO: Output logic and energy accounting can be separated later


                writeout(n, outHQ, outLightPos, outLightNeg, config, partons_event);

		std::vector<Particle>().swap(partons_event);
		if((int) partons_event.size()>0){
			std::cout << "ERROR !!!!" << std::endl;
			exit(1);
		}

	}//Event loop

return;
}


double get_rapidity(const double E, const double pz){
	const double LARGE = 10.0;
	double rap;
	if(E==0.0 && pz==0.0) rap=0.0;
	else if(E==-pz) rap=-LARGE;
	else if(E==pz) rap=LARGE;
	else rap=std::log((E + pz)/(E - pz)) / 2.0;
return rap;
}

void writeout(const int n, std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg, const LBTConfig config,
		const std::vector<Particle> part_event){

	const int heavyOut = config.output.heavyOut;
	const int lightOut = config.output.lightOut;

	int col=0;
	int acol=0;
	double mass=0.0;
	if(heavyOut){
		for(int i = 0; i<(int)part_event.size(); i++){

			const Particle &p = part_event[i];

			if(p.CAT==4 && p.isActive)
				outHQ 
					<< std::setw(5) << i
					<< std::setw(5) << p.index()
					<< std::setw(5) << col//col and acols: dummy
					<< std::setw(5) << acol//col and acols: dummy
					<< std::setw(5) << p.pid
					<< std::setw(15) << std::setprecision(10) << mass
					<< std::setw(15) << std::setprecision(10) << p.P[0]
					<< std::setw(15) << std::setprecision(10) << p.P[1]
					<< std::setw(15) << std::setprecision(10) << p.P[2]
					<< std::setw(15) << std::setprecision(10) << p.P[3]
					<< std::setw(15) << std::setprecision(10) << get_rapidity(p.P[0], p.P[3])
					<< std::setw(15) << std::setprecision(10) << p.V[1]
					<< std::setw(15) << std::setprecision(10) << p.V[2]
					<< std::setw(15) << std::setprecision(10) << p.V[3]
					<< std::setw(15) << std::setprecision(10) << p.V[0]
					<< std::setw(15) << std::setprecision(10) << p.V[0]
					<< std::setw(5) << p.CAT
					<< std::endl;
		}
	}
	if(lightOut){
		for(int i = 0; i<(int)part_event.size(); i++){

			const Particle &p = part_event[i];

			if(p.CAT!=3 && p.isActive)
				outLightPos 
					<< std::setw(5) << i
					<< std::setw(5) << p.index()
					<< std::setw(5) << col//col and acols: dummy
					<< std::setw(5) << acol//col and acols: dummy
					<< std::setw(5) << p.pid
					<< std::setw(15) << std::setprecision(10) << mass
					<< std::setw(15) << std::setprecision(10) << p.P[0]
					<< std::setw(15) << std::setprecision(10) << p.P[1]
					<< std::setw(15) << std::setprecision(10) << p.P[2]
					<< std::setw(15) << std::setprecision(10) << p.P[3]
					<< std::setw(15) << std::setprecision(10) << get_rapidity(p.P[0], p.P[3])
					<< std::setw(15) << std::setprecision(10) << p.V[1]
					<< std::setw(15) << std::setprecision(10) << p.V[2]
					<< std::setw(15) << std::setprecision(10) << p.V[3]
					<< std::setw(15) << std::setprecision(10) << p.V[0]
					<< std::setw(15) << std::setprecision(10) << p.V[0]
					<< std::setw(5) << p.CAT
					<< std::endl;
		}
		for(int i = 0; i<(int)part_event.size(); i++){

			const Particle &p = part_event[i];

			if(p.CAT==3&& p.isActive)
				outLightNeg 
					<< std::setw(5) << i
					<< std::setw(5) << p.index()
					<< std::setw(5) << col//col and acols: dummy
					<< std::setw(5) << acol//col and acols: dummy
					<< std::setw(5) << p.pid
					<< std::setw(15) << std::setprecision(10) << mass
					<< std::setw(15) << std::setprecision(10) << p.P[0]
					<< std::setw(15) << std::setprecision(10) << p.P[1]
					<< std::setw(15) << std::setprecision(10) << p.P[2]
					<< std::setw(15) << std::setprecision(10) << p.P[3]
					<< std::setw(15) << std::setprecision(10) << get_rapidity(p.P[0], p.P[3])
					<< std::setw(15) << std::setprecision(10) << p.V[1]
					<< std::setw(15) << std::setprecision(10) << p.V[2]
					<< std::setw(15) << std::setprecision(10) << p.V[3]
					<< std::setw(15) << std::setprecision(10) << p.V[0]
					<< std::setw(15) << std::setprecision(10) << p.V[0]
					<< std::setw(5) << p.CAT
					<< std::endl;
		}
	}

	if(heavyOut){
		outHQ << "%" << std::setw(5) << n << std::endl;
	}
	if(lightOut){
		outLightPos << "%" << std::setw(5) << n << std::endl;
		outLightNeg << "%" << std::setw(5) << n << std::endl;
	}
return;
}




void finalize(std::ofstream& outHQ,
              std::ofstream& outLightPos,
              std::ofstream& outLightNeg,
              const LBTConfig& config,
              std::time_t time_start) {
    if (config.output.lightOut) {
        outLightPos.close();
        outLightNeg.close();
    }
    if (config.output.heavyOut) {
        outHQ.close();
    }

    std::time_t time_end = std::time(nullptr);
    std::cout << "Program ends at: " << std::asctime(std::localtime(&time_end));

    unsigned duration = std::difftime(time_end, time_start);
    unsigned hours = duration / 3600;
    unsigned minutes = (duration % 3600) / 60;
    unsigned seconds = duration % 60;

    std::cout << "Total run time: "
              << duration << "s (" << hours << "h " << minutes << "m " << seconds << "s)" << std::endl;
    return;
}
