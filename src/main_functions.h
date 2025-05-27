#pragma once
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <math.h>
#include "LBTConfig.h"
#include "ParticleInfo.h"
#include "LBTcl.h"

bool initialize(LBTConfig& config, int argc, char* argv[], std::time_t& time_start);
void open_output(int argc, char* argv[],
		const LBTConfig& config,
		std::ifstream& fpList,
		std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg);
void runLBT(std::ifstream& fpList,
		LBTConfig& config,
		std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg);
void writeout(const int n, std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg, const LBTConfig config,
		const std::vector<Particle> part_event);
void finalize(std::ofstream& outHQ,
		std::ofstream& outLightPos,
		std::ofstream& outLightNeg,
		const LBTConfig& config,
		std::time_t time_start);
double get_rapidity(const double E, const double pz);
void readTables(LBTConfig &config);
