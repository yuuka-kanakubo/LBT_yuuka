///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define variables and functions for LBT
//
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef common_H
#define common_H

class common{

public:
//...constant
//const double    pi=3.1415926;
//const double    epsilon=1e-6;
//const double    CA=3.0; 
//const double    CF=4.0/3.0;
//const double    sctr=0.1973;	     //fm to GeV^-1

//...fixed parameters (cannot change using input parameter file)
//int lightOut=1;     // 1 -> write out light parton information
//int heavyOut=0;     // 1 -> write out heavy quark information
//int outFormat=1;    // different output formats: 1 - my preferred, 2 - for JETSCAPE patch code
//onst double cutOut=epsilon;  // neglect light partons in output files for E < cutOut;

//...paramters for the bulk, can be changed from input parameter file
//int       vacORmed=1; // 0-vacuum; 1-medium
//int       bulkFlag=0; // 0-static medium; 1-OSU hydro; 2-CCNU hydro 
//double    temp0=0.25; // medium temperature                 
//double    hydro_Tc=0.165;
////...derived quantities
//double    temp00; 

//...input parameters for LBT
//double    KPamp=0.0;
//double    KPsig=5.0;
//double    KTamp=0.0;
//double    KTsig=0.05*hydro_Tc;
//double    preKT=0.3/0.3;
//double    scaleAK=2.0;
////...derived quantities
//double    KPfactor,KTfactor,Kfactor,runKT;

//...initialization parameters for jet partons
//int       initHardFlag=2;  //1 initialize by the code itself; 2 initialize by reading in particle list
//int       fixMomentum=0;
//int       fixPosition=0;
//int       flagJetX=0;      // 0: do nothing; 1: keep momentum but reset jet position within LBT
//int       Kjet=21;         //initial flavor of the jet parton
//double    ipTmin=0.0;
//double    ipTmax=800.0;
//double    eta_cut=0.5;
//double    ener=50.0;       //initial energy of the jet parton/heavy quark
//double    amss=0.0;        //initial mass of the jet parton/heavy quark 
//int       nj=1000;         //number of initial partons per event      
//int       ncall=1;         // number of events
////...derived quantities
//int       np;              //number of partons at this time step ti				

//...clock information
//int       nprint=100;	// print out infomation per # of events	
//int       tauswitch=0;  // 0 for t-z coordinate, 1 (not included in this version) for tau-eta 
//double    tau0=0.0; // initial time		
//double    dtau=0.1; // time interval		
//double    tauend=4.0; // time when program ends		
////...derived quantities
//double    time0;
//double    dt;        //dtime when tauswitch is turned off (0) / dtau when tauswitch is turned on (1)
//double    timend;    //end time or tau RENAME
//
//...switches and cuts
//int       Kprimary=1;  // 0 keep all partons, 1 keep leading jet parton only (switch off other partons)
//int       KINT0=1;     // 0 no radiation, 1 elastic + inelastic
//int       Kqhat0=2;    //Debye screening mass switch RENAME
//int       Kalphas=1;   //alphas switch		
//double    Ecut=0.0;    //energy cut of the recoiled partons
//double    fixAlphas=0.3; 
////...derived quantities
//int       KINT;                     //radiation switch
//double    alphas;
//double    qhat0;                    //Debye mass RENAME
//double    qhat00;		

		
//...input with current machine time
//...random number seed (any negative integer)	  
//  long  NUM1=-33;                  
//long  NUM1;                  

//    scattering rate
//double Rg[60][20];         //total gluon scattering rate as functions of initial energy and temperature 
//double Rg1[60][20];        //gg-gg              CT1
//double Rg2[60][20];        //gg-qqbar           CT2
//double Rg3[60][20];        //gq-qg              CT3
//double Rq[60][20];         //total gluon scattering rate as functions of initial energy and temperature
//double Rq3[60][20];        //qg-qg              CT13
//double Rq4[60][20];        //qiqj-qiqj          CT4
//double Rq5[60][20];        //qiqi-qiqi          CT5
//double Rq6[60][20];        //qiqibar-qjqjbar    CT6
//double Rq7[60][20];        //qiqibar-qiqibar    CT7
//double Rq8[60][20];        //qqbar-gg           CT8	  
//double qhatLQ[60][20];
//double qhatG[60][20];
//
//double RHQ[60][20];        //total scattering rate for heavy quark
//double RHQ11[60][20];      //Qq->Qq
//double RHQ12[60][20];      //Qg->Qg
//double qhatHQ[60][20];     //qhat of heavy quark
//double qhat_over_T3;        //qhat/T^3 for heavy quark as fnc of (T,p)
//double D2piT;

// for heavy quark radiation table
//const int HQener_gn=500;
//const int t_gn=75;
//const int temp_gn=100;
//
//double dNg_over_dt_c[t_gn+2][temp_gn+1][HQener_gn+1]={0.0};
//double dNg_over_dt_q[t_gn+2][temp_gn+1][HQener_gn+1]={0.0};
//double dNg_over_dt_g[t_gn+2][temp_gn+1][HQener_gn+1]={0.0};
//double max_dNgfnc_c[t_gn+2][temp_gn+1][HQener_gn+1]={0.0};
//double max_dNgfnc_q[t_gn+2][temp_gn+1][HQener_gn+1]={0.0};
//double max_dNgfnc_g[t_gn+2][temp_gn+1][HQener_gn+1]={0.0};
//
//const double HQener_max=1000.0;
//const double t_max=15.0;
//const double temp_max=0.65;
//const double temp_min=0.15;
//double delta_tg=t_max/t_gn;
//double delta_temp=(temp_max-temp_min)/temp_gn;
//double delta_HQener=HQener_max/HQener_gn;

// for MC initialization of jet partons
//const int maxMC=2000000;
//double initMCX[maxMC],initMCY[maxMC];
//		
//...radiation block		
//int    icl22;                   
//int    icl23;                    //the numerical switch in colljet23
//int    iclrad;                   //the numerical switch in radiation 	  
//int    isp;                      //the splitting function switch
//
//...global variable qhat	
//int counth100=0;
  
/////double qhat;	                 //transport parameter

//double dng0[101][101]={0.0};	 //table of dn/dkperp2/dx 	  

//double Vtemp[4]={0.0};

const int dimParList=500000;

//...time system
double tirad[dimParList]={0.0};
double tiscatter[dimParList]={0.0};
double tiform[dimParList]={0.0};	         //pythia undone 
double Tint_lrf[dimParList]={0.0};          //for heavy quark
//double eGluon=0.0;
//double nGluon=0.0;

//....radiated gluon
double radng[dimParList]={0.0};	  
//double xwm[3]={0.0};
//double wkt2m;
	  
//double vf[4]={0.0};              //flow velocity in tau-eta coordinate
//double vfcar[4]={0.0};           //flow velocity in t-z coordinate
//
//double vp[4]={0.0};              //position of particle		
//double vc0[4]={0.0};             //flow velocity     
//	
////...dimensions in subroutine colljet and twcoll
//double vc[4]={0.0};              
//double pc0[4]={0.0};
//double pc2[4]={0.0};
//double pc3[4]={0.0};
//double pc4[4]={0.0};		
//double p0[4]={0.0};
//double p2[4]={0.0};
//double p3[4]={0.0};		
//double p4[4]={0.0};

//double pc00[4]={0.0};
//double pc30[4]={0.0};	
//		
//double pc01[4]={0.0};
//double pb[4]={0.0};
//		
//double Pj0[4]={0.0};
		
double V[4][dimParList]={0.0};     //parton position 
double P[6][dimParList]={0.0};     //parton 4-momentum
double V0[4][dimParList]={0.0};    //negative parton position
double P0[6][dimParList]={0.0};    //negative parton 4-momentum

double Prad[4][dimParList]={0.0};		

double WT[dimParList]={0};
double WT0[dimParList]={0};
		
int NR[dimParList]={0};            //scattering rank 
int KATT1[dimParList]={0};         //parton flavor  
int KATT10[dimParList]={0};        //negative parton flavor
//

//double PGm[4]={0.0};
double tjp[dimParList]={0.0};
double Vfrozen[4][dimParList]={0.0};     //parton final 4 coordinate
double Vfrozen0[4][dimParList]={0.0};    //negative parton final 4 coordinate 
//double Ecmcut=2.0;                   //energy cut for free streaming   
double Tfrozen[dimParList]={0.0};
double Tfrozen0[dimParList]={0.0};
double vcfrozen[4][dimParList]={0.0};
double vcfrozen0[4][dimParList]={0.0};

//........................................................................................................NCUT
double VV[4][dimParList]={0.0};
double VV0[4][dimParList]={0.0};
double PP[4][dimParList]={0.0};
double PP0[4][dimParList]={0.0};
int CAT[dimParList]={0};
int CAT0[dimParList]={0};
		
int ncut;
int ncut0;

//........................................................................................................NCUTEND

//int n_coll22=0;
//int n_coll23=0;
//int ng_coll23=0;
//int ng_nrad=0;
//int n_radiation=0;
//int ng_radiation=0;
//int n_gluon=0;
//int n_sp1=0;
//int n_sp2=0;
//
// Variables for HQ 2->2
//const int N_p1=500;
//const int N_T=60;
//const int N_e2=75;
//double distFncB[N_T][N_p1][N_e2],distFncF[N_T][N_p1][N_e2],distMaxB[N_T][N_p1][N_e2],distMaxF[N_T][N_p1][N_e2];
//double distFncBM[N_T][N_p1],distFncFM[N_T][N_p1];
//double min_p1=0.0;
//double max_p1=1000.0;
//double bin_p1=(max_p1-min_p1)/N_p1;
//double min_T=0.1;
//double max_T=0.7;
//double bin_T=(max_T-min_T)/N_T;
//double min_e2=0.0;
//double max_e2=15.0;
//double bin_e2=(max_e2-min_e2)/N_e2;
 
//int loopN=10000;
//int loopN=1000;
//
//int numInitXY=0; 
 
//...functions	  

////void LBT(int &n, double &ti);

void trans(double v[4],double p[4]);
void transback(double v[4],double p[4]);
  
float ran0(long *idum);
   
double alphas0(int &Kalphas,double temp0);
double DebyeMass2(int &Kqhat0,double alphas,double temp0);
	  
void lam(int KATT0,double &RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2);	  

void flavor(int &CT,int &KATT0,int &KATT2,int &KATT3,double RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2);
void linear(int KATT,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2,double &RTEg,double &RTEg1,double &RTEg2,double &RTEg3,double &RTEq,double &RTEq3,double &RTEq4,double &RTEq5,double &RTEq6,double &RTEq7,double &RTEq8,double &RTEHQ,double &RTEHQ11,double &RTEHQ12,double &qhatTP);
void twflavor(int &CT,int &KATT0,int &KATT2,double E,double T);
void twlinear(int KATT,double E,double T,double &RTEg1,double &RTEg2,double &RTEq6,double &RTEq7,double &RTEq8);
void colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
void twcoll(int CT,double qhat0ud,double v0[4],double p0[4],double p2[4]);
	  
void titau(double ti,double vf[4],double vp[4],double p0[4],double &Vx,double &Vy,double &Veta,double &Xtau);


void rotate(double px,double py,double pz,double p[4],int icc);
	  
int KPoisson(double alambda);	  

void radiationHQ(int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int &ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt);  
void collHQ22(int CT,double temp,double qhat,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
double Mqc2qc(double s, double t, double M);
double Mgc2gc(double s, double t, double M);

void collHQ23(int parID, double temp_med, double qhat0ud, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double qt, int &ic, double Tdiff, double HQenergy, double max_Ng, double xLow, double xInt);
double dNg_over_dxdydt(int parID, double x0g, double y0g, double HQenergy, double HQmass, double temp_med, double Tdiff);
double tau_f(double x0g, double y0g, double HQenergy, double HQmass);
double splittingP(int parID, double z0g);
double lambdas(double kTFnc);
double nflavor(double kTFnc);
double alphasHQ(double kTFnc, double tempFnc);
double nHQgluon(int parID,double dtLRF,double &time_gluon,double &temp_med,double &HQenergy,double &max_Ng);

void read_xyMC(int &numXY);
void jetInitialize(int numXY);
void setJetX(int numXY);
void read_tables();
void jetClean();
void setParameter(string fileName);
int  checkParameter(int nArg);

extern "C" {
    void read_ccnu_(char *dataFN_in, int len1);
    void hydroinfoccnu_(double *Ct, double *Cx, double *Cy, double *Cz, double *Ctemp, double *Cvx, double *Cvy, double *Cvz, int *Cflag);

    void sethydrofilesez_(int *dataID_in, char *dataFN_in, int *ctlID_in, char *ctlFN_in, int *bufferSize, int len1, int len2);
    void readhydroinfoshanshan_(double *t, double *x, double *y, double *z, double *e, double *s, double *temp, double *vx, double *vy, double *vz, int *flag);
}


};
#endif
