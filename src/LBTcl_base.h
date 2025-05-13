#pragma once
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

    static constexpr double pi = 3.1415926;
    static constexpr double epsilon = 1e-6;
    static constexpr double CA = 3.0;
    static constexpr double CF = 4.0 / 3.0;
    static constexpr double sctr = 0.1973;  // fm to GeV^-1
    static constexpr double GEVFM = 0.1970;
 
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





inline double alphas0(const int Kalphas, const double T){
       double X;
       if(Kalphas==1)
       {
               X=0.3;
       }else{
               std::cout << "Currently Kalphas!=1 is not implemented." << std::endl;
               std::cout << __FILE__ << "(" << __LINE__ << ")" << std::endl;
               exit(EXIT_FAILURE);
       }
       return X;
}

inline double DebyeMass2(const int Kqhat0, const double alphas, const double T){

       double Y;
       if(Kqhat0==1){
               Y=4.0*base::pi*alphas*pow(T,2);
               return Y;
       }else if(Kqhat0==2){
               Y=(3.0/2.0)*4.0*base::pi*alphas*pow(T,2);
               return Y;
       }else if(Kqhat0==3){
               Y=1.0;
               return Y;
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

