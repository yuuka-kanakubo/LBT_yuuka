#include "main_functions.h"

// Declare external functions
void read_tables(LBTConfig& config){};
void jetClean(){};
void jetInitialize(int numInitXY){};
void setJetX(int numInitXY){};
void LBT(int eventID, double time, LBTConfig& config){};
double alphas0(int Kalphas, double T){return 0.0;};
double DebyeMass2(int Kqhat0, double alphas, double T){return 0.0;};

float ran0(long *idum){
	int j;
	long k;
	//##############
	//static long seed_w=(unsigned) time(NULL);
	static long seed_w=123;

	static long idum2=seed_w;
	//##############
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) { 
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { 
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1; 
	*idum=IA1*(*idum-k*IQ1)-k*IR1; 
	if (*idum < 0) *idum += IM1; 
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV; 
	iy=iv[j]-idum2; 
	iv[j] = *idum; 
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp;
}

bool initialize(LBTConfig& config, int argc, char* argv[], std::time_t& time_start) {
	config.loadFromFile(argv[1]);
	if (config.checkParameter(argc) != 0) {
		std::cerr << "Parameter check failed" << std::endl;
		return false;
	}

	if (config.medium.vacORmed == 1)
		read_tables(config);

	//Set derived quantities
	config.medium.temp00 = config.medium.temp0;
	config.clock.dt = config.clock.dtau;
	config.clock.timend = config.clock.tauend;
	config.clock.time0 = config.clock.tau0;
	config.physics.alphas = alphas0(config.physics.Kalphas, config.medium.temp0);
	config.physics.qhat0 = DebyeMass2(config.physics.Kqhat0, config.physics.alphas, config.medium.temp0);

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
		if(abs(it->KATT)==1||abs(it->KATT)==2||abs(it->KATT)==3||abs(it->KATT)==21){
			it->P[4]=0.0;
			it->P[0]=sqrt(it->P[1]*it->P[1]+it->P[2]*it->P[2]+it->P[3]*it->P[3]+it->P[4]*it->P[4]);
			it->P[5]=sqrt(it->P[1]*it->P[1]+it->P[2]*it->P[2]);//transverse momentum
			it->WT=1.0;
		}
	}
	return;
}



void runLBT(std::ifstream& fpList,
		std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg,
		LBTConfig& config) {

	int& nj = config.jet.nj;
	int& np = config.jet.np;
	const int initHardFlag = config.jet.initHardFlag;
	int numEvent = 0;

	for (int n = 1; ; ++n) {
		std::vector<Particle> partons_event;

		// Stop conditions
		if (initHardFlag == 1 && numEvent >= config.jet.ncall) break;
		if (initHardFlag == 2 && (numEvent >= config.jet.ncall || fpList.eof())) break;

		jetClean();

		if (initHardFlag == 1) {
			jetInitialize(config.counter.numInitXY);
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
				if (!(iss >> dummyInt >> p.KATT >> p.P[1] >> p.P[2] >> p.P[3] >> p.P[0]
							>> p.Vfrozen[1] >> p.Vfrozen[2] >> p.Vfrozen[3] >> p.Vfrozen[0])) {
					std::cerr << "[Warning] Skipping malformed line: " << line << "\n";
					continue;
				}

				p.index = i;
				p.isPrimary = true;
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
		std::cout << "np " << np << std::endl;
		std::cout << "(int) partons_event.size() " << (int) partons_event.size() << std::endl;
		for (auto it = partons_event.begin(); it != partons_event.end(); ++it) {
			it->Print();
		}


		// Time evolution (skipping detailed LBT logic)


		if (config.medium.vacORmed == 1) {
			for (double ti = config.clock.time0 + config.clock.dt;
					ti <= config.clock.timend + LBTConfig::epsilon;
					ti += config.clock.dt) {
				LBT(n, ti, config);
			}
		}

		++numEvent;

		if (n % config.clock.nprint == 0) {
			std::cout << "n = " << n << ", np = " << np << std::endl;
		}

		// TODO: Output logic and energy accounting can be separated later
	}//Event loop
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
