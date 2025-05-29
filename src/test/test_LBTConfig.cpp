#include "../LBTConfig.h"
#include "macro_cmn.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>

void test_checkParameter() {
    std::cout << "Running test_checkParameter()...\n";

    // Save current cout buffer
    std::streambuf* originalCout = std::cout.rdbuf();
    // Redirect cout to /dev/null or nul
#ifdef _WIN32
    std::ofstream nullStream("nul");
#else
    std::ofstream nullStream("/dev/null");
#endif
    std::cout.rdbuf(nullStream.rdbuf());

    LBTConfig config;

    config.jet.initHardFlag = 1;
    config.output.heavyOut = 1;
    config.output.lightOut = 1;
    ASSERT_EQUAL(config.checkParameter(5), 0);
    ASSERT_EQUAL(config.checkParameter(4), 1);

    config.jet.initHardFlag = 2;
    config.output.heavyOut = 1;
    config.output.lightOut = 0;
    ASSERT_EQUAL(config.checkParameter(4), 0);
    ASSERT_EQUAL(config.checkParameter(3), 1);

    config.output.heavyOut = 0;
    config.output.lightOut = 1;
    ASSERT_EQUAL(config.checkParameter(5), 0);
    ASSERT_EQUAL(config.checkParameter(4), 1);

    config.output.lightOut = 0;
    ASSERT_EQUAL(config.checkParameter(2), 0);
    ASSERT_EQUAL(config.checkParameter(1), 1);

    config.jet.initHardFlag = 999;
    ASSERT_EQUAL(config.checkParameter(2), 1);

    std::cout.rdbuf(originalCout);
    std::cout << "All checkParameter() tests passed.\n\n";
}

void test_loadFromFile() {
    std::cout << "Running test_loadFromFile()...\n";
    const char* testFile = "test_params.dat";
    std::ofstream ofs(testFile);
    ofs << "flag:initHardFlag = 1\n";
    ofs << "flag:fixPosition = 1\n";
    ofs << "flag:vacORmed = 0\n";
    ofs << "para:alpha_s = 0.27\n";
    ofs << "para:KTsig = 0.05\n";
    ofs << "para:hydro_Tc = 0.165\n";
    ofs.close();

    // Save current cout buffer
    std::streambuf* originalCout = std::cout.rdbuf();
    // Redirect cout to /dev/null or nul
#ifdef _WIN32
    std::ofstream nullStream("nul");
#else
    std::ofstream nullStream("/dev/null");
#endif
    std::cout.rdbuf(nullStream.rdbuf());

    LBTConfig config;
    config.loadFromFile(testFile);
    std::cout.rdbuf(originalCout);

    ASSERT_EQUAL(config.jet.initHardFlag, 1);
    ASSERT_EQUAL(config.jet.fixPosition, 1);
    ASSERT_EQUAL(config.medium.vacORmed, 0);
    ASSERT_EQUAL(config.physics.fixAlphas, 0.27);
    ASSERT_EQUAL(config.medium.hydro_Tc, 0.165);

    double expectedKTsig = 0.05 * config.medium.hydro_Tc;
    double expectedPreKT = config.physics.fixAlphas / 0.3;
    ASSERT_EQUAL(std::round(config.lbtinput.KTsig * 1e6), std::round(expectedKTsig * 1e6));
    ASSERT_EQUAL(std::round(config.lbtinput.preKT * 1e6), std::round(expectedPreKT * 1e6));

    std::remove(testFile);
    std::cout << "loadFromFile() test passed.\n\n";
}

int main() {
    test_checkParameter();
    test_loadFromFile();
    std::cout << "All tests completed successfully.\n";
    return 0;
}

