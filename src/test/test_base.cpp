#include "../LBTConfig.h"
#include "../LBTcl_base.h"
#include "macro_cmn.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>




void test_DebyeMass2(){

	LBTConfig config;
	config.physics.Kqhat0 = 2;
	config.physics.Kalphas = 1;//Constant alphas 0.3
        double alphas = alphas0(config.physics.Kalphas, 0.30);
	double tol = 1e10;
	ASSERT_EQUAL(std::round(DebyeMass2(config.physics.Kqhat0, alphas, 0.30) * tol), std::round(0.5089380098815464*tol));


}


void test_alpha_s(){

	LBTConfig config;
	config.physics.Kalphas = 1;//Constant alphas 0.3

        ASSERT_EQUAL(alphas0(config.physics.Kalphas, 0.30), 0.30);

}




int main() {
    test_alpha_s();//running coupling constant
    test_DebyeMass2();
    std::cout << "All tests completed successfully.\n";
    return 0;
}

