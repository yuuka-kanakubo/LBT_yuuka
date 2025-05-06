#ifndef LBT_CONFIG_H
#define LBT_CONFIG_H
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <string>
#include <string>
#include <array>
#include <vector>


struct MediumParameters {
    int vacORmed = 1;
    int bulkFlag = 0;
    double temp0 = 0.25;
    double hydro_Tc = 0.165;
    double temp00 = 0.0;  // derived
};

struct LBTInputParameters {
    double KPamp = 0.0;
    double KPsig = 5.0;
    double KTamp = 0.0;
    double KTsig = 0.05 * 0.165; // initialized with hydro_Tc, should update
    double preKT = 1.0;//0.3/0.3
    double scaleAK = 2.0;
    double KPfactor = 0.0, KTfactor = 0.0, Kfactor = 0.0, runKT = 0.0; // derived
};

struct JetInitializationParameters {
    int initHardFlag = 2;
    int fixMomentum = 0;
    int fixPosition = 0;
    int flagJetX = 0;
    int Kjet = 21;
    double ipTmin = 0.0;
    double ipTmax = 800.0;
    double eta_cut = 0.5;
    double ener = 50.0;
    double amss = 0.0;
    int nj = 1000;
    int ncall = 1;
    int np = 0; // derived
};

struct ClockParameters {
    int nprint = 100;
    int tauswitch = 0;
    double tau0 = 0.0;
    double dtau = 0.1;
    double tauend = 4.0;
    double time0 = 0.0;
    double dt = 0.0;
    double timend = 0.0;
};

struct PhysicsSwitches {
    int Kprimary = 1;
    int KINT0 = 1;
    int Kqhat0 = 2;
    int Kalphas = 1;
    double Ecut = 0.0;
    double fixAlphas = 0.3;
    int KINT = 0;
    double alphas = 0.0;
    double qhat0 = 0.0;
    double qhat00 = 0.0;
};

struct RNGSettings {
    long NUM1 = -33;
};

struct RadiationSwitches {
    int icl22 = 0;
    int icl23 = 0;
    int iclrad = 0;
    int isp = 0;
};

struct Counters {
    int n_coll22 = 0;
    int n_coll23 = 0;
    int ng_coll23 = 0;
    int ng_nrad = 0;
    int n_radiation = 0;
    int ng_radiation = 0;
    int n_gluon = 0;
    int n_sp1 = 0;
    int n_sp2 = 0;
    int loopN = 1000;
    int counth100 = 0;
};

struct OutputSettings {
    int lightOut = 1;
    int heavyOut = 0;
    int outFormat = 1;
    double cutOut = 1e-6; // same as epsilon
};

struct ScatteringTables {
    double Rg[60][20]{};
    double Rg1[60][20]{};
    double Rg2[60][20]{};
    double Rg3[60][20]{};
    double Rq[60][20]{};
    double Rq3[60][20]{};
    double Rq4[60][20]{};
    double Rq5[60][20]{};
    double Rq6[60][20]{};
    double Rq7[60][20]{};
    double Rq8[60][20]{};
    double qhatLQ[60][20]{};
    double qhatG[60][20]{};
    double RHQ[60][20]{};
    double RHQ11[60][20]{};
    double RHQ12[60][20]{};
    double qhatHQ[60][20]{};
    double qhat_over_T3 = 0.0;
    double D2piT = 0.0;
};

struct HeavyQuarkRadiation {
    static constexpr int HQener_gn = 500;
    static constexpr int t_gn = 75;
    static constexpr int temp_gn = 100;

    std::vector<std::vector<std::vector<double>>> dNg_over_dt_c;
    std::vector<std::vector<std::vector<double>>> dNg_over_dt_q;
    std::vector<std::vector<std::vector<double>>> dNg_over_dt_g;

    std::vector<std::vector<std::vector<double>>> max_dNgfnc_c;
    std::vector<std::vector<std::vector<double>>> max_dNgfnc_q;
    std::vector<std::vector<std::vector<double>>> max_dNgfnc_g;

    double HQener_max = 1000.0;
    double t_max = 15.0;
    double temp_max = 0.65;
    double temp_min = 0.15;

    double delta_tg;
    double delta_temp;
    double delta_HQener;

    HeavyQuarkRadiation() {

	    delta_tg = t_max / t_gn;
	    delta_temp = (temp_max - temp_min) / temp_gn;
	    delta_HQener = HQener_max / HQener_gn;

	    dNg_over_dt_c.resize(t_gn + 2, std::vector<std::vector<double>>(temp_gn + 1, std::vector<double>(HQener_gn + 1, 0.0)));
	    dNg_over_dt_q.resize(t_gn + 2, std::vector<std::vector<double>>(temp_gn + 1, std::vector<double>(HQener_gn + 1, 0.0)));
	    dNg_over_dt_g.resize(t_gn + 2, std::vector<std::vector<double>>(temp_gn + 1, std::vector<double>(HQener_gn + 1, 0.0)));

	    max_dNgfnc_c.resize(t_gn + 2, std::vector<std::vector<double>>(temp_gn + 1, std::vector<double>(HQener_gn + 1, 0.0)));
	    max_dNgfnc_q.resize(t_gn + 2, std::vector<std::vector<double>>(temp_gn + 1, std::vector<double>(HQener_gn + 1, 0.0)));
	    max_dNgfnc_g.resize(t_gn + 2, std::vector<std::vector<double>>(temp_gn + 1, std::vector<double>(HQener_gn + 1, 0.0)));

    }
};

struct FlowAndPositionVectors {
    static constexpr int maxMC = 2000000;
    std::vector<double> initMCX = std::vector<double>(maxMC, 0.0);
    std::vector<double> initMCY = std::vector<double>(maxMC, 0.0);

    double Vtemp[4]{};
    double xwm[3]{};
    double wkt2m = 0.0;
    double vf[4]{};
    double vfcar[4]{};
    double vp[4]{};
    double vc0[4]{};
    double vc[4]{};
    double pc0[4]{};
    double pc2[4]{};
    double pc3[4]{};
    double pc4[4]{};
    double p0[4]{};
    double p2[4]{};
    double p3[4]{};
    double p4[4]{};
    double pc00[4]{};
    double pc30[4]{};
    double pc01[4]{};
    double pb[4]{};
    double Pj0[4]{};
    double PGm[4]{};
    double vfrozen[4]{};
    double vfrozen0[4]{};
    double vcfrzn[4]{};
    double vcfrzn0[4]{};
    double Vx = 0.0, Vy = 0.0, Veta = 0.0, Xtau = 0.0;
    double Ecmcut = 2.0;
};

struct HQ22Distributions {
    static constexpr int N_p1 = 500;
    static constexpr int N_T = 60;
    static constexpr int N_e2 = 75;

    std::vector<std::vector<std::vector<double>>> distFncB;
    std::vector<std::vector<std::vector<double>>> distFncF;
    std::vector<std::vector<std::vector<double>>> distMaxB;
    std::vector<std::vector<std::vector<double>>> distMaxF;

    std::vector<std::vector<double>> distFncBM;
    std::vector<std::vector<double>> distFncFM;

    double min_p1 = 0.0;
    double max_p1 = 1000.0;
    double bin_p1;

    double min_T = 0.1;
    double max_T = 0.7;
    double bin_T;

    double min_e2 = 0.0;
    double max_e2 = 15.0;
    double bin_e2;

    HQ22Distributions() {
	    bin_p1 = (max_p1 - min_p1) / N_p1;
	    bin_T = (max_T - min_T) / N_T;
	    bin_e2 = (max_e2 - min_e2) / N_e2;

	    distFncB.resize(N_T, std::vector<std::vector<double>>(N_p1, std::vector<double>(N_e2, 0.0)));
	    distFncF.resize(N_T, std::vector<std::vector<double>>(N_p1, std::vector<double>(N_e2, 0.0)));
	    distMaxB.resize(N_T, std::vector<std::vector<double>>(N_p1, std::vector<double>(N_e2, 0.0)));
	    distMaxF.resize(N_T, std::vector<std::vector<double>>(N_p1, std::vector<double>(N_e2, 0.0)));

	    distFncBM.resize(N_T, std::vector<double>(N_p1, 0.0));
	    distFncFM.resize(N_T, std::vector<double>(N_p1, 0.0));
    }


};

struct GluonEmissionTable {
    double dng0[101][101]{};
};




class LBTConfig {

public:
    MediumParameters medium;
    LBTInputParameters lbtinput;
    JetInitializationParameters jet;
    ClockParameters clock;
    PhysicsSwitches physics;
    RNGSettings rng;
    RadiationSwitches radiation;
    Counters counter;
    OutputSettings output;
    ScatteringTables tables;
    HeavyQuarkRadiation hqrad;
    FlowAndPositionVectors flow;
    HQ22Distributions hq22;
    GluonEmissionTable emission;

    static constexpr double pi = 3.1415926;
    static constexpr double epsilon = 1e-6;
    static constexpr double CA = 3.0;
    static constexpr double CF = 4.0 / 3.0;
    static constexpr double sctr = 0.1973;  // fm to GeV^-1


    void loadFromFile(const std::string& filename);
    int checkParameter(int nArg) const;



};

#endif // LBT_CONFIG_H

