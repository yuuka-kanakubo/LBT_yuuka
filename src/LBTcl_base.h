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


namespace base{

    static constexpr double pi = 3.1415926;
    static constexpr double epsilon = 1e-6;
    static constexpr double CA = 3.0;
    static constexpr double CF = 4.0 / 3.0;
    static constexpr double sctr = 0.1973;  // fm to GeV^-1
    static constexpr double GEVFM = 0.1970;

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
