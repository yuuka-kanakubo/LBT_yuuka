#pragma once
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <math.h>
#include "LBTConfig.h"
#include "ParticleInfo.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran0(long *idum);
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
