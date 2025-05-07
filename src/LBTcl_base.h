#pragma once

namespace base{

    static constexpr double pi = 3.1415926;
    static constexpr double epsilon = 1e-6;
    static constexpr double CA = 3.0;
    static constexpr double CF = 4.0 / 3.0;
    static constexpr double sctr = 0.1973;  // fm to GeV^-1
    static constexpr double GEVFM = 0.1970;

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
