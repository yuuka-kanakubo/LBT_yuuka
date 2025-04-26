#ifndef LBTcl_H
#define LBTcl_H


class LBTcl{

  int CT;                        //collision type 1 2 3 13 4 5 6 7 8
  int KATTC0;                    //flavor code 1/d 2/u 3/s -1/dbar -2/ubar -3/sbar 21/gluon 
  int KATT2;
  int KATT3;
  int KATT30;
		
  double RTE;                    //scattering rate (energy temperature)
  double E;                      //parton energy
  double PLen;                   //parton momentum
  double T;                      //local temperature		

  double T1;                     //index for difference method
  double T2;
  double E1;
  double E2;
  int iT1;
  int iT2;
  int iE1;
  int iE2;	
		
  int nrad;
  int idlead;
  int idlead1;
  int idlead2;
  double Xtau;                   //....................main & titau	
  double Vx;
  double Vy;		
  double Veta;

  double tcar;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar;		
  double ycar;
  double zcar;					

  double tcar0;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar0;		
  double ycar0;
  double zcar0;


  double rans;

  int KATTx;
  int KATTy;

  double Ejp;
  double Elab;
		
  double qt;
		
  double px0;
  double py0;		
  double pz0;

  double Ncoll22;                    //...................average number of elastic scattering for particle i in the current step  
  int Nscatter;                   //...................number of elastic scattering for particle i in the current step


  int free=0;
  int free0=0;
  double fraction=0.0;
  double fraction0=0.0;
  double vc0b[4]={0.0};             //flow velocity     
  double pMag,vMag,flowFactor;

  double kt2;
  
  double vp0[4]={0.0};
  double p0temp[4]={0.0}; 

  double p0temp1=0.0; 
  double p0temp2=0.0; 
		
  double pcx[4]={0.0};
  double pcy[4]={0.0};		

  double pcx1[4]={0.0};
  double pcy1[4]={0.0};		

  // for heavy quark
  double dt_lrf,maxFncHQ,Tdiff,lim_low,lim_high,lim_int;
  double probCol,probRad,probTot;


public:
  LBTcl(){};
  ~LBTcl(){};

};
#endif
