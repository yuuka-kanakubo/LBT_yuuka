#ifndef LBTCL_H
#define LBTCL_H
#include <vector>
#include <cmath>
#include "ParticleInfo.h"
#include "LBTConfig.h"
#include "LBTcl_base.h"


class LBTcl{

	private:

		LBTConfig& config;
		double computeScatteringRate(Particle &p, const double PLen, const double T);
		double computeRadiationProbability(Particle &p, double T, double E);
		void handleElasticCollision(Particle &p, std::vector<Particle> &particles);
		void handleRadiation(Particle &p, std::vector<Particle> &particles, int &icl23, int &iclrad);
		void propagateParticle(Particle &p, double ti, int &free, double &fraction);
		double computeCollisionProbability(
				const Particle &p,
				double qhat,
				double pLen,
				double T,
				double fraction
				);



		void titau(double ti,
				const std::array<double, 4> &vf,
				const std::array<double, 4> &vp,
				const std::array<double, 4> &p0,
				double &Vx, double &Vy, double &Veta, double &Xtau) {

			double gamma = 1.0 / sqrt(1.0 - (vf[1] * vf[1] + vf[2] * vf[2]));

			double mt = sqrt(p0[1]*p0[1] + p0[2]*p0[2]);
			double Yp = 0.5 * log((p0[0] + p0[3]) / (p0[0] - p0[3]));

			double etas = vp[3];
			double etaf = atanh(vf[3]) + etas;

			double pper = mt;  // same as sqrt(p0[1]^2 + p0[2]^2)
			double vper = sqrt(vf[1]*vf[1] + vf[2]*vf[2]);
			double pvper = p0[1]*vf[1] + p0[2]*vf[2];

			Vx   = p0[1] / (pper * cosh(Yp - etas));
			Vy   = p0[2] / (pper * cosh(Yp - etas));
			Veta = tanh(Yp - etas) / ti;

			Xtau = (gamma * mt * cosh(Yp - etaf) - pvper * vper) / (mt * cosh(Yp - etas));
			return;
		}

		void trans(const std::array<double, 4> &v, std::array<double, 4> &p) {
			double vv = sqrt(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
			if (vv < 1e-12) return;

			double ga = 1.0 / sqrt(1.0 - vv * vv);
			double ppar = p[1]*v[1] + p[2]*v[2] + p[3]*v[3];
			double gavv = (ppar * ga / (1.0 + ga) - p[0]) * ga;

			p[0] = ga * (p[0] - ppar);
			p[1] += v[1] * gavv;
			p[2] += v[2] * gavv;
			p[3] += v[3] * gavv;
		}

		void transback(const std::array<double, 4> &v, std::array<double, 4> &p) {
			double vv = sqrt(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
			if (vv < 1e-12) return;

			double ga = 1.0 / sqrt(1.0 - vv * vv);
			double ppar = p[1]*v[1] + p[2]*v[2] + p[3]*v[3];
			double gavv = (-ppar * ga / (1.0 + ga) - p[0]) * ga;

			p[0] = ga * (p[0] + ppar);
			p[1] -= v[1] * gavv;
			p[2] -= v[2] * gavv;
			p[3] -= v[3] * gavv;
		}


		void lam2(const int flavour,double &RTE1, double &RTE2,const double T,double &T1,double &T2,int &iT1,int &iT2,int &iE1,int &iE2){   
			// --- Step 2: Interpolate in Temperature first ---
			if (flavour == 21) {
				// Gluon
				RTE1 = (config.tables.qhatG[iT2][iE1] - config.tables.qhatG[iT1][iE1]) * (T - T1) / (T2 - T1) + config.tables.qhatG[iT1][iE1];
				RTE2 = (config.tables.qhatG[iT2][iE2] - config.tables.qhatG[iT1][iE2]) * (T - T1) / (T2 - T1) + config.tables.qhatG[iT1][iE2];
			} else if (abs(flavour) == 4) {
				// Heavy quark
				RTE1 = (config.tables.qhatHQ[iT2][iE1] - config.tables.qhatHQ[iT1][iE1]) * (T - T1) / (T2 - T1) + config.tables.qhatHQ[iT1][iE1];
				RTE2 = (config.tables.qhatHQ[iT2][iE2] - config.tables.qhatHQ[iT1][iE2]) * (T - T1) / (T2 - T1) + config.tables.qhatHQ[iT1][iE2];
			} else {
				// Light quark
				RTE1 = (config.tables.qhatLQ[iT2][iE1] - config.tables.qhatLQ[iT1][iE1]) * (T - T1) / (T2 - T1) + config.tables.qhatLQ[iT1][iE1];
				RTE2 = (config.tables.qhatLQ[iT2][iE2] - config.tables.qhatLQ[iT1][iE2]) * (T - T1) / (T2 - T1) + config.tables.qhatLQ[iT1][iE2];
			}
			return;
		}

		void lam(const int flavour,double &RTE,const double E,const double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2){   

			double dtemp=0.02;
			iT1=(int)((T-0.1)/dtemp);
			iT2=iT1+1;
			iE1=(int)(log(E)+2);
			if(iE1<1) iE1=1;
			iE2=iE1+1;

			T1=0.12+(iT1-1)*0.02;
			T2=T1+dtemp;
			E1=exp(iE1-2.0);
			E2=exp(iE2-2.0);

			if(flavour==21) {	
				double RTE1=(config.tables.Rg[iT2][iE1]-config.tables.Rg[iT1][iE1])*(T-T1)/(T2-T1)+config.tables.Rg[iT1][iE1];
				double RTE2=(config.tables.Rg[iT2][iE2]-config.tables.Rg[iT1][iE2])*(T-T1)/(T2-T1)+config.tables.Rg[iT1][iE2];
				RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
			} else if (flavour==4||flavour==-4) { // add heavy quark channel
				double RTE1=(config.tables.RHQ[iT2][iE1]-config.tables.RHQ[iT1][iE1])*(T-T1)/(T2-T1)+config.tables.RHQ[iT1][iE1];
				double RTE2=(config.tables.RHQ[iT2][iE2]-config.tables.RHQ[iT1][iE2])*(T-T1)/(T2-T1)+config.tables.RHQ[iT1][iE2];
				RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
			} else {    
				double RTE1=(config.tables.Rq[iT2][iE1]-config.tables.Rq[iT1][iE1])*(T-T1)/(T2-T1)+config.tables.Rq[iT1][iE1];
				double RTE2=(config.tables.Rq[iT2][iE2]-config.tables.Rq[iT1][iE2])*(T-T1)/(T2-T1)+config.tables.Rq[iT1][iE2];
				RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
			}

			std::cout << "RTE: " << RTE << std::endl;
			return;
		}







		double nHQgluon(const Particle &p, const double dtLRF,
				const double temp_med_,const double HQenergy_,double &max_Ng){
			// gluon radiation probability for heavy quark       

			int flavour = p.KATT;
			double time_gluon = p.Tint_lrf;
			double temp_med = temp_med_;
			double HQenergy = HQenergy_;
			if(time_gluon>config.hqrad.t_max) {
				std::cout << "accumulated time exceeds t_max" << std::endl;
				std::cout << time_gluon << "    " << temp_med << "    " << HQenergy << std::endl;
				time_gluon=config.hqrad.t_max;
			}

			if(temp_med>config.hqrad.temp_max) {
				std::cout << "temperature exceeds temp_max -- extrapolation is used" << std::endl;
				std::cout << time_gluon << "    " << temp_med << "    " << HQenergy << std::endl;
				//     temp_med=temp_max;
			}

			if(HQenergy>config.hqrad.HQener_max) {
				std::cout << "HQenergy exceeds HQener_max -- extrapolation is used" << std::endl;
				std::cout << time_gluon << "    " << temp_med << "    " << HQenergy << std::endl;
				//     HQenergy=HQener_max;
			}

			if(temp_med<config.hqrad.temp_min) {
				std::cout << "temperature drops below temp_min" << std::endl;
				std::cout << time_gluon << "    " << temp_med << "    " << HQenergy << std::endl;
				temp_med=config.hqrad.temp_min;
			}

			int time_num=(int)(time_gluon/config.hqrad.delta_tg+0.5)+1;
			//  temp_num=(int)((temp_med-temp_min)/delta_temp+0.5);
			//  HQenergy_num=(int)(HQenergy/delta_HQener+0.5); // use linear interpolation instead of finding nearest point for E and T dimensions
			int temp_num=(int)((temp_med-config.hqrad.temp_min)/config.hqrad.delta_temp);
			int HQenergy_num=(int)(HQenergy/config.hqrad.delta_HQener); // normal interpolation

std::cout << time_num << std::endl;
std::cout << temp_num << std::endl;
std::cout << HQenergy_num << std::endl;

			if(HQenergy_num >= config.hqrad.HQener_gn) HQenergy_num=config.hqrad.HQener_gn-1; // automatically become extrapolation
			if(temp_num >= config.hqrad.temp_gn) temp_num=config.hqrad.temp_gn-1;


			double rate_T1E1,rate_T1E2,rate_T2E1,rate_T2E2,max_T1E1,max_T1E2,max_T2E1,max_T2E2;
			if(flavour==21) {
				rate_T1E1 = config.hqrad.dNg_over_dt_g[time_num][temp_num][HQenergy_num];
				rate_T1E2 = config.hqrad.dNg_over_dt_g[time_num][temp_num][HQenergy_num+1];
				rate_T2E1 = config.hqrad.dNg_over_dt_g[time_num][temp_num+1][HQenergy_num];
				rate_T2E2 = config.hqrad.dNg_over_dt_g[time_num][temp_num+1][HQenergy_num+1];
				max_T1E1 = config.hqrad.max_dNgfnc_g[time_num][temp_num][HQenergy_num];
				max_T1E2 = config.hqrad.max_dNgfnc_g[time_num][temp_num][HQenergy_num+1];
				max_T2E1 = config.hqrad.max_dNgfnc_g[time_num][temp_num+1][HQenergy_num];
				max_T2E2 = config.hqrad.max_dNgfnc_g[time_num][temp_num+1][HQenergy_num+1];
			} else if(abs(flavour)==4) {
				rate_T1E1 = config.hqrad.dNg_over_dt_c[time_num][temp_num][HQenergy_num];
				rate_T1E2 = config.hqrad.dNg_over_dt_c[time_num][temp_num][HQenergy_num+1];
				rate_T2E1 = config.hqrad.dNg_over_dt_c[time_num][temp_num+1][HQenergy_num];
				rate_T2E2 = config.hqrad.dNg_over_dt_c[time_num][temp_num+1][HQenergy_num+1];
				max_T1E1 = config.hqrad.max_dNgfnc_c[time_num][temp_num][HQenergy_num];
				max_T1E2 = config.hqrad.max_dNgfnc_c[time_num][temp_num][HQenergy_num+1];
				max_T2E1 = config.hqrad.max_dNgfnc_c[time_num][temp_num+1][HQenergy_num];
				max_T2E2 = config.hqrad.max_dNgfnc_c[time_num][temp_num+1][HQenergy_num+1];
			} else {
				rate_T1E1 = config.hqrad.dNg_over_dt_q[time_num][temp_num][HQenergy_num];
				rate_T1E2 = config.hqrad.dNg_over_dt_q[time_num][temp_num][HQenergy_num+1];
				rate_T2E1 = config.hqrad.dNg_over_dt_q[time_num][temp_num+1][HQenergy_num];
				rate_T2E2 = config.hqrad.dNg_over_dt_q[time_num][temp_num+1][HQenergy_num+1];
				max_T1E1 = config.hqrad.max_dNgfnc_q[time_num][temp_num][HQenergy_num];
				max_T1E2 = config.hqrad.max_dNgfnc_q[time_num][temp_num][HQenergy_num+1];
				max_T2E1 = config.hqrad.max_dNgfnc_q[time_num][temp_num+1][HQenergy_num];
				max_T2E2 = config.hqrad.max_dNgfnc_q[time_num][temp_num+1][HQenergy_num+1];
			} 

			double rate_EGrid_low = rate_T1E1+(temp_med-config.hqrad.temp_min-temp_num*config.hqrad.delta_temp)/config.hqrad.delta_temp*(rate_T2E1-rate_T1E1);
			double rate_EGrid_high = rate_T1E2+(temp_med-config.hqrad.temp_min-temp_num*config.hqrad.delta_temp)/config.hqrad.delta_temp*(rate_T2E2-rate_T1E2);
			double max_EGrid_low = max_T1E1+(temp_med-config.hqrad.temp_min-temp_num*config.hqrad.delta_temp)/config.hqrad.delta_temp*(max_T2E1-max_T1E1);
			double max_EGrid_high = max_T1E2+(temp_med-config.hqrad.temp_min-temp_num*config.hqrad.delta_temp)/config.hqrad.delta_temp*(max_T2E2-max_T1E2);

			double delta_Ng;
			delta_Ng = rate_EGrid_low+(HQenergy-HQenergy_num*config.hqrad.delta_HQener)/config.hqrad.delta_HQener*(rate_EGrid_high-rate_EGrid_low);
			max_Ng = max_EGrid_low+(HQenergy-HQenergy_num*config.hqrad.delta_HQener)/config.hqrad.delta_HQener*(max_EGrid_high-max_EGrid_low);


std::cout << " D2piT " << p.D2piT << std::endl;
std::cout << " dtLRF " << dtLRF << std::endl;
std::cout << " rate_T2E1 " << rate_T2E1 << std::endl;


			delta_Ng*=6.0/p.D2piT*dtLRF;
			max_Ng*=6.0/p.D2piT;

			//  if(delta_Ng>1) {
			//     std::cout << "Warning: Ng greater than 1   " << time_gluon << "  " << delta_Ng << std::endl;
			//  }

			return delta_Ng;

		}	  












	public:


		void LBT(std::vector<Particle> &particles, double ti);

		LBTcl(LBTConfig& config_in):config(config_in){};
		~LBTcl(){};

};

extern "C" {
	void read_ccnu_(char *dataFN_in, int len1);
	void hydroinfoccnu_(double *Ct, double *Cx, double *Cy, double *Cz, double *Ctemp, double *Cvx, double *Cvy, double *Cvz, int *Cflag);

	void sethydrofilesez_(int *dataID_in, char *dataFN_in, int *ctlID_in, char *ctlFN_in, int *bufferSize, int len1, int len2);
	void readhydroinfoshanshan_(double *t, double *x, double *y, double *z, double *e, double *s, double *temp, double *vx, double *vy, double *vz, int *flag);
}
#endif
