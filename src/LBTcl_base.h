#pragma once

#include <cmath>

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

inline int index_counter =0;

namespace base{

    static constexpr double pi = 3.141592653589793;
    static constexpr double epsilon = 1e-6;
    static constexpr double tol = 1e-10;
    static constexpr double CA = 3.0;
    static constexpr double CF = 4.0 / 3.0;
    static constexpr double GEVFM = 0.19732698;
 
}

inline double get_magnitude(const std::array<double, 4>& p_) {
    return std::sqrt(p_[1]*p_[1] + p_[2]*p_[2] + p_[3]*p_[3]);
}


inline double get_energy(const std::array<double, 4>& p) {
    return p[0];
}

inline double get_mass(const std::array<double, 4>& p) {
    double E2 = p[0] * p[0];
    double p2 = p[1]*p[1] + p[2]*p[2] + p[3]*p[3];
    double m2 = E2 - p2;
    return (m2 > 1e-12) ? std::sqrt(m2) : 0.0;
}



inline double tau_f(const double x0g, const double y0g, const double energy, const double mass) {
	return 2.0*energy*x0g*(1.0-x0g)/((x0g*y0g*energy)*(x0g*y0g*energy)+x0g*x0g*mass*mass);
}

inline double splittingP(const int parID, const double z0g) {

      if(parID==21) return 2.0*pow(1.0-z0g+pow(z0g,2),3)/z0g/(1.0-z0g);
      else return (1.0-z0g)*(2.0-2.0*z0g+z0g*z0g)/z0g;
}


inline double dNg_over_dxdydt(int parID, double x0g, double y0g, double energy, double mass, double temp_med, double alpha_s, double qhat_over_T3, double Tdiff) {

	double tauFnc = tau_f(x0g,y0g,energy,mass);
	double qhatFnc = qhat_over_T3*pow(temp_med,3);   // no longer have CF factor, taken out in splittingP too.

	return 4.0/base::pi*base::CA*alpha_s*splittingP(parID,x0g)*qhatFnc*pow(y0g,5)*pow(sin(Tdiff/2.0/tauFnc/base::GEVFM),2)*pow((energy*energy/(y0g*y0g*energy*energy+mass*mass)),4)/x0g/x0g/energy/energy/base::GEVFM;

	//TODO: Maybe in the future update on alphas for HQ?
	//return 4.0/base::pi*base::CA*alphasHQ(x0g*y0g*energy,temp_med)*splittingP(parID,x0g)*qhatFnc*pow(y0g,5)*pow(sin(Tdiff/2.0/tauFnc/base::GEVFM),2)*pow((energy*energy/(y0g*y0g*energy*energy+mass*mass)),4)/x0g/x0g/energy/energy/base::GEVFM;

}



inline std::array <double, 4> get_centerofmass(const std::array<double, 4> p0_, const std::array<double, 4> p2_){

	std::array <double, 4> vc_;
	vc_[0]=0.;
	vc_[1]=(p0_[1]+p2_[1])/(p0_[0]+p2_[0]);
	vc_[2]=(p0_[2]+p2_[2])/(p0_[0]+p2_[0]);
	vc_[3]=(p0_[3]+p2_[3])/(p0_[0]+p2_[0]);

	return vc_;
}




inline void rotate(
    double px, double py, double pz,               // target vector to align with z-axis
    std::array<double, 4>& vec,                    // vector to be rotated (in-place)
    int direction                                  // +1 = align with z-axis, -1 = inverse
) {
    // Input: vec = (E, x, y, z)
    // Output: vec[1,2,3] rotated so that target (px,py,pz) aligns with z-axis
    // If direction == 1: rotate (px,py,pz) -> (0,0,E)
    // If direction == -1: inverse rotation

    double x = vec[1];
    double y = vec[2];
    double z = vec[3];

    const double p_mag = std::sqrt(px * px + py * py + pz * pz);    // |p|
    const double pT = std::sqrt(px * px + py * py);                 // transverse part of p

    // Define rotation angles around z and y axes
    const double cos_phi = (pT == 0.0) ? 1.0 : px / pT;
    const double sin_phi = (pT == 0.0) ? 0.0 : py / pT;

    const double cos_theta = pz / p_mag;
    const double sin_theta = pT / p_mag;

    // Rotated vector components
    double x_rot, y_rot, z_rot;

    if (direction == 1) {
        // Rotate (x,y,z) into z-axis
        x_rot =  x * cos_theta * cos_phi + y * cos_theta * sin_phi - z * sin_theta;
        y_rot = -x * sin_phi + y * cos_phi;
        z_rot =  x * sin_theta * cos_phi + y * sin_theta * sin_phi + z * cos_theta;
    } else {
        // Inverse rotation: z-axis -> original direction (px, py, pz)
        x_rot = x * cos_phi * cos_theta - y * sin_phi + z * cos_phi * sin_theta;
        y_rot = x * sin_phi * cos_theta + y * cos_phi + z * sin_phi * sin_theta;
        z_rot = -x * sin_theta + z * cos_theta;
    }

    vec[1] = x_rot;
    vec[2] = y_rot;
    vec[3] = z_rot;

    // Optionally recompute energy (vec[0]) from massless approximation:
    // vec[0] = std::sqrt(vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3]);
}





inline float ran0(long *idum){
	int j;
	long k;
	//##############
	//static long seed_w=(unsigned) time(NULL);
	static long seed_w=123;

	static long idum2=seed_w;
	//##############
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) { 
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { 
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1; 
	*idum=IA1*(*idum-k*IQ1)-k*IR1; 
	if (*idum < 0) *idum += IM1; 
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV; 
	iy=iv[j]-idum2; 
	iv[j] = *idum; 
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp;
}



inline float CLCG(long *idum){
	int j;
	long k;
	//##############
	//static long seed_w=(unsigned) time(NULL);
	static long seed_w=123;

	static long idum2=seed_w;
	//##############
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) { 
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { 
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1; 
	*idum=IA1*(*idum-k*IQ1)-k*IR1; 
	if (*idum < 0) *idum += IM1; 
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV; 
	iy=iv[j]-idum2; 
	iv[j] = *idum; 
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp;
}





inline double alphas0(const int Kalphas, const double T){
       double X;
       if(Kalphas==1)
       {
               X=0.3;
       }else{
               X=0.3;
               std::cout << "Currently Kalphas!=1 is not implemented. temperature " << T << " is not used."<< std::endl;
               std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
               exit(EXIT_FAILURE);
       }
       return X;
}

// functions for gluon radiation from heavy quark
//double alphasHQ(const double kTFnc, const double tempFnc) {
//
//	//      double kTEff,error_para,resultFnc;
//	//
//	//      error_para=1.0;
//	//
//	//      if(kTFnc<pi*tempFnc*error_para) {
//	//         kTEff=pi*tempFnc*error_para;
//	//      } else {
//	//         kTEff=kTFnc;
//	//      }
//	//
//	//      resultFnc=4.0*pi/(11.0-2.0*nflavor(kTEff)/3.0)/2.0/log(kTEff/lambdas(kTEff));
//	//
//	//      return(resultFnc);
//	///return(alphas);
//}




inline double DebyeMass2(const int Kqhat0, const double alphas, const double T){

       double qhat;
       if(Kqhat0==1){
               qhat=4.0*base::pi*alphas*pow(T,2);
               return qhat;
       }else if(Kqhat0==2){
               qhat=(3.0/2.0)*4.0*base::pi*alphas*pow(T,2);
               return qhat;
       }else if(Kqhat0==3){
               qhat=1.0;
               return qhat;
       }else{
               std::cout << "Currently this is not implemented in the refactored version of LBT." << std::endl;
               std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
               exit(EXIT_FAILURE);
       }
}



inline  double Mqc2qc(const double s, const double t, const double M) {

      double m2m=M*M;
      double u=2.0*m2m-s-t;

      double MM=64.0/9.0*(pow((m2m-u),2)+pow((s-m2m),2)+2.0*m2m*t)/t/t;

      return MM;

  }

inline  double Mgc2gc(const double s, const double t, const double M) {

      double m2m=M*M;
      double u=2.0*m2m-s-t;
      double MM = 32.0*(s-m2m)*(m2m-u)/t/t;
      MM+=64.0/9.0*((s-m2m)*(m2m-u)+2.0*m2m*(s+m2m))/pow((s-m2m),2);
      MM+=64.0/9.0*((s-m2m)*(m2m-u)+2.0*m2m*(u+m2m))/pow((u-m2m),2);
      MM+=16.0/9.0*m2m*(4.0*m2m-t)/((s-m2m)*(m2m-u));
      MM+=16.0*((s-m2m)*(m2m-u)+m2m*(s-u))/(t*(s-m2m));
      MM+=16.0*((s-m2m)*(m2m-u)-m2m*(s-u))/(t*(u-m2m));

      return MM;

  }

