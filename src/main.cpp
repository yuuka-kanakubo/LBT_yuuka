#include <fstream>
#include <iostream>
#include <ctime>
#include "main_functions.h"


int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: ./LBT parameter_file [...files]" << std::endl;
		return EXIT_FAILURE;
	}

	LBTConfig config;

	std::time_t time_start;
	if (!initialize(config, argc, argv, time_start)) return EXIT_FAILURE;
	std::cout << "The end of initialize () " << std::endl;

	std::ifstream fpList;
	std::ofstream outHQ, outLightPos, outLightNeg;
	open_output(argc, argv, config, fpList, outHQ, outLightPos, outLightNeg);

	runLBT(fpList, outHQ, outLightPos, outLightNeg, config);
	finalize(outHQ, outLightPos, outLightNeg, config, time_start);

	std::cout << "SUCCESSFULLY FINISHED LBT RUN! :)" << std::endl;
	return 0;

}
