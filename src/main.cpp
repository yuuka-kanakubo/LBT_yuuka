#include "LBTConfig.h"
#include <iostream>
#include <fstream>
#include <ctime>

// Declare external functions
void read_tables(LBTConfig& config);
void jetClean();
void jetInitialize(int numInitXY);
void setJetX(int numInitXY);
void LBT(int eventID, double time, LBTConfig& config);
double alphas0(int Kalphas, double T);
double DebyeMass2(int Kqhat0, double alphas, double T);
double ran0(long* idum);

int main(int argc, char* argv[]) {

	if (argc < 2) {
		std::cerr << "Usage: ./LBT parameter_file [...files]" << std::endl;
		return EXIT_FAILURE;
	}

	LBTConfig config;
	config.loadFromFile(argv[1]);
	if (config.checkParameter(argc) != 0) return EXIT_FAILURE;

	if (config.medium.vacORmed == 1) read_tables(config);

	config.medium.temp00 = config.medium.temp0;
	config.clock.dt = config.clock.dtau;
	config.clock.timend = config.clock.tauend;
	config.clock.time0 = config.clock.tau0;
	config.physics.alphas = alphas0(config.physics.Kalphas, config.medium.temp0);
	config.physics.qhat0 = DebyeMass2(config.physics.Kqhat0, config.physics.alphas, config.medium.temp0);

	srand(123);
	//srand((unsigned)time(NULL));
	config.rng.NUM1 = -1 * rand();

	std::time_t time_start = std::time(nullptr);
	std::cout << "Program starts at: " << std::asctime(std::localtime(&time_start));

	std::ifstream fpList;
	std::ofstream outHQ, positive, negative;
	const int initHardFlag = config.jet.initHardFlag;
	const int heavyOut = config.output.heavyOut;
	const int lightOut = config.output.lightOut;

	if (initHardFlag == 1) {
		if (heavyOut && lightOut) {
			outHQ.open(argv[2]);
			positive.open(argv[3]);
			negative.open(argv[4]);
		} else if (heavyOut) {
			outHQ.open(argv[2]);
		} else if (lightOut) {
			positive.open(argv[2]);
			negative.open(argv[3]);
		}
	} else {
		fpList.open(argv[2]);
		if (heavyOut && lightOut) {
			outHQ.open(argv[3]);
			positive.open(argv[4]);
			negative.open(argv[5]);
		} else if (heavyOut) {
			outHQ.open(argv[3]);
		} else if (lightOut) {
			positive.open(argv[3]);
			negative.open(argv[4]);
		}
	}

	int numEvent = 0;
	for (int n = 1; ; ++n) {
		if (initHardFlag == 1 && numEvent >= config.jet.ncall) break;
		if (initHardFlag == 2 && (numEvent >= config.jet.ncall || fpList.eof())) break;

		jetClean();

		if (initHardFlag == 1) {
			jetInitialize(config.counter.numInitXY);
		} else {
			// Read particle list from fpList and populate config.P[], etc.
			// Not implemented here for brevity
		}

		config.jet.np = config.jet.nj;

		if (config.medium.vacORmed == 1) {
			for (double ti = config.clock.time0 + config.clock.dt;
					ti <= config.clock.timend + LBTConfig::epsilon; 
					ti += config.clock.dt) {
				LBT(n, ti, config);
			}
		}

		// Output logic and energy accounting would go here

		numEvent++;
		if (n % config.clock.nprint == 0) {
			std::cout << "n = " << n << ", np = " << config.jet.np << std::endl;
		}
	}

	if (lightOut) {
		positive.close();
		negative.close();
	}
	if (heavyOut) outHQ.close();

	std::time_t time_end = std::time(nullptr);
	std::cout << "Program ends at: " << std::asctime(std::localtime(&time_end));
	unsigned duration = std::difftime(time_end, time_start);
	std::cout << "Total run time: " << duration << " seconds" << std::endl;

	return 0;
}

