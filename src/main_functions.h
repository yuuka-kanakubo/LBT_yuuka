#pragma once
#include <fstream>
#include <iostream>
#include <ctime>
#include "LBTConfig.h"
#include "ParticleInfo.h"



		bool initialize(LBTConfig& config, int argc, char* argv[], std::time_t& time_start);
		void open_output(int argc, char* argv[],
				const LBTConfig& config,
				std::ifstream& fpList,
				std::ofstream& outHQ,
				std::ofstream& outLightPos,
				std::ofstream& outLightNeg);
		void runLBT(std::ifstream& fpList,
				std::ofstream& outHQ,
				std::ofstream& outLightPos,
				std::ofstream& outLightNeg,
				LBTConfig& config);
		void finalize(std::ofstream& outHQ,
				std::ofstream& outLightPos,
				std::ofstream& outLightNeg,
				const LBTConfig& config,
				std::time_t time_start);


