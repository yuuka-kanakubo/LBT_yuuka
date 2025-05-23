#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../LBTConfig.h"
#include "../main_functions.h"


TEST_CASE("initialize works with valid file") {
	LBTConfig config;
	std::time_t start;
	const char* argv[] = { "./LBT", "data/parameters.dat", "data/jet_shower_parton_test.dat", "data/lp-posi_test.dat", "data/lp-nega_test.dat"};
	int argc = 5;

	CHECK(initialize(config, argc, const_cast<char**>(argv), start));
}



TEST_CASE("open_output opens all expected files with initHardFlag=1") {
	LBTConfig config;
	config.jet.initHardFlag = 1;
	config.output.heavyOut = 1;
	config.output.lightOut = 1;

	const char* argv[] = {
		"./LBT",
		"dummy_param.dat",  // argv[1]
		"test_hq.txt",      // argv[2]
		"test_pos.txt",     // argv[3]
		"test_neg.txt"      // argv[4]
	};
	int argc = 5;

	std::ifstream fpList;
	std::ofstream outHQ, outPos, outNeg;

	CHECK_NOTHROW(open_output(argc, const_cast<char**>(argv), config, fpList, outHQ, outPos, outNeg));
	CHECK(outHQ.is_open());
	CHECK(outPos.is_open());
	CHECK(outNeg.is_open());

	outHQ.close(); outPos.close(); outNeg.close();
	std::remove("test_hq.txt");
	std::remove("test_pos.txt");
	std::remove("test_neg.txt");
}
TEST_CASE("open_output with heavyOut only") {
	LBTConfig config;
	config.jet.initHardFlag = 1;
	config.output.heavyOut = 1;
	config.output.lightOut = 0;

	const char* argv[] = {
		"./LBT",
		"dummy_param.dat",
		"test_hq.txt"  // only HQ output
	};
	int argc = 3;

	std::ifstream fpList;
	std::ofstream outHQ, outPos, outNeg;

	open_output(argc, const_cast<char**>(argv), config, fpList, outHQ, outPos, outNeg);

	CHECK(outHQ.is_open());
	CHECK(!outPos.is_open());
	CHECK(!outNeg.is_open());

	outHQ.close();
	std::remove("test_hq.txt");
}

TEST_CASE("open_output with lightOut only") {
	LBTConfig config;
	config.jet.initHardFlag = 1;
	config.output.heavyOut = 0;
	config.output.lightOut = 1;

	const char* argv[] = {
		"./LBT",
		"dummy_param.dat",
		"test_pos.txt",
		"test_neg.txt"
	};
	int argc = 4;

	std::ifstream fpList;
	std::ofstream outHQ, outPos, outNeg;

	open_output(argc, const_cast<char**>(argv), config, fpList, outHQ, outPos, outNeg);

	CHECK(!outHQ.is_open());
	CHECK(outPos.is_open());
	CHECK(outNeg.is_open());

	outPos.close(); outNeg.close();
	std::remove("test_pos.txt");
	std::remove("test_neg.txt");
}

TEST_CASE("open_output with initHardFlag == 2 and both outputs") {
	LBTConfig config;
	config.jet.initHardFlag = 2;
	config.output.heavyOut = 1;
	config.output.lightOut = 1;

	// Create dummy input list file to satisfy fpList.open
	std::ofstream dummyInput("input_list.txt");
	dummyInput << "0 0\n";  // header line
	dummyInput.close();

	const char* argv[] = {
		"./LBT",
		"dummy_param.dat",
		"input_list.txt",
		"test_hq.txt",
		"test_pos.txt",
		"test_neg.txt"
	};
	int argc = 6;

	std::ifstream fpList;
	std::ofstream outHQ, outPos, outNeg;

	open_output(argc, const_cast<char**>(argv), config, fpList, outHQ, outPos, outNeg);

	CHECK(fpList.is_open());
	CHECK(outHQ.is_open());
	CHECK(outPos.is_open());
	CHECK(outNeg.is_open());

	fpList.close(); outHQ.close(); outPos.close(); outNeg.close();
	std::remove("input_list.txt");
	std::remove("test_hq.txt");
	std::remove("test_pos.txt");
	std::remove("test_neg.txt");
}


