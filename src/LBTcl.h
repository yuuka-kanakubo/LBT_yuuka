#ifndef LBTCL_H
#define LBTCL_H
#include <vector>
#include <cmath>
#include <algorithm>
#include "ParticleInfo.h"
#include "LBTConfig.h"
#include "LBTcl_base.h"


class LBTcl{

	private:

		int np_snapshot;
		bool belowCutOff(const double Eloc, const Particle &p);

		LBTConfig& config;
		double computeScatteringRate(Particle &p, const double PLen, const double T);
		double computeRadiationProbability(Particle &p, double T, double E);
		double handleElasticCollision(Particle &p, const double PLen, std::vector<Particle> &particles_current, 
				std::vector<Particle> &particles);
		int handleRadiation(Particle &p, const double qt, std::vector<Particle> &particles_current, 
				std::vector<Particle> &particles);
		void propagateParticle(Particle &p, double ti, int &free, double &fraction);
		double computeCollisionProbability(
				Particle &p,
				double qhat,
				double pLen,
				double T,
				double fraction
				);
		void FinalTouch(Particle &p, std::vector<Particle> & part_current);




		int KPoisson(const double alambda){
			//....Generate numbers according to Poisson distribution
			//    P(lambda)=lambda**k exp(-lambda)/k!
			//    input: average number of radiated gluons lambda
			//    output: number of radiated gluons in this event


			double target=exp(-alambda);
			double p=ran0(&config.rng.NUM1);

			double KKPoisson=0;
			while(p>target)
			{
				p=p*ran0(&config.rng.NUM1);
				KKPoisson++;
			}		
			return KKPoisson;
		}




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

			return;
		}







		double nHQgluon(Particle &p, const double dtLRF,
				const double temp_med_,const double HQenergy_){
			// gluon radiation probability for heavy quark       
			int flavour = p.pid;
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
			double max_Ng = max_EGrid_low+(HQenergy-HQenergy_num*config.hqrad.delta_HQener)/config.hqrad.delta_HQener*(max_EGrid_high-max_EGrid_low);



			delta_Ng*=6.0/p.D2piT*dtLRF;
			max_Ng*=6.0/p.D2piT;
			p.max_Ng = max_Ng;
			p.Tint_lrf = time_gluon;


			return delta_Ng;

		}	  



		void linear(
				const int pid, const double E, const double T,
				double& RTEg1, double& RTEg2, double& RTEg3,
				double& RTEq, double& RTEq3, double& RTEq4, double& RTEq5,
				double& RTEq6, double& RTEq7, double& RTEq8,
				double& RTEHQ11, double& RTEHQ12
			   ) {
			// Temperature and energy grid binning
			double dtemp = 0.02;
			int iT1 = static_cast<int>((T - 0.1) / dtemp);
			int iT2 = iT1 + 1;
			int iE1 = static_cast<int>(std::log(E) + 2.0);
			int iE2 = iE1 + 1;

			double T1 = 0.12 + (iT1 - 1) * dtemp;
			double T2 = T1 + dtemp;
			double E1 = std::exp(iE1 - 2.0);
			double E2 = std::exp(iE2 - 2.0);

			// Bilinear interpolation per flavor channel
			if (pid == 21) {
				auto interp = [&](auto& table, int j) {
					return (table[iT2][j] - table[iT1][j]) * (T - T1) / (T2 - T1) + table[iT1][j];
				};
				RTEg1 = (interp(config.tables.Rg1, iE2) - interp(config.tables.Rg1, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rg1, iE1);
				RTEg2 = (interp(config.tables.Rg2, iE2) - interp(config.tables.Rg2, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rg2, iE1);
				RTEg3 = (interp(config.tables.Rg3, iE2) - interp(config.tables.Rg3, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rg3, iE1);
			}
			else if (std::abs(pid) == 4) {
				auto interp = [&](auto& table, int j) {
					return (table[iT2][j] - table[iT1][j]) * (T - T1) / (T2 - T1) + table[iT1][j];
				};
				RTEHQ11 = (interp(config.tables.RHQ11, iE2) - interp(config.tables.RHQ11, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.RHQ11, iE1);
				RTEHQ12 = (interp(config.tables.RHQ12, iE2) - interp(config.tables.RHQ12, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.RHQ12, iE1);

				//        qhatTP  = (interp(config.tables.qhatHQ, iE2) - interp(config.tables.qhatHQ, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.qhatHQ, iE1);
			}
			else {
				auto interp = [&](auto& table, int j) {
					return (table[iT2][j] - table[iT1][j]) * (T - T1) / (T2 - T1) + table[iT1][j];
				};
				RTEq3 = (interp(config.tables.Rq3, iE2) - interp(config.tables.Rq3, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rq3, iE1);
				RTEq4 = (interp(config.tables.Rq4, iE2) - interp(config.tables.Rq4, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rq4, iE1);
				RTEq5 = (interp(config.tables.Rq5, iE2) - interp(config.tables.Rq5, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rq5, iE1);
				RTEq6 = (interp(config.tables.Rq6, iE2) - interp(config.tables.Rq6, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rq6, iE1);
				RTEq7 = (interp(config.tables.Rq7, iE2) - interp(config.tables.Rq7, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rq7, iE1);
				RTEq8 = (interp(config.tables.Rq8, iE2) - interp(config.tables.Rq8, iE1)) * (E - E1) / (E2 - E1) + interp(config.tables.Rq8, iE1);

				RTEq = RTEq3 + RTEq4 + RTEq5 + RTEq6 + RTEq7 + RTEq8;
			}
		}




		void flavor(
				const int initialPid, const double RTE, const double E, const double T,
				int& channel, int& pid2, int& pid3
			   ) {
			static const std::array<int, 6> valencePids = {1, 2, 3, -1, -2, -3};

			// Scattering rate components
			double RTEg1, RTEg2, RTEg3;
			double RTEq, RTEq3, RTEq4, RTEq5, RTEq6, RTEq7, RTEq8;
			double RTEHQ11, RTEHQ12;

			linear(initialPid, E, T, 
					RTEg1, RTEg2, RTEg3,
					RTEq, RTEq3, RTEq4, RTEq5, RTEq6, RTEq7, RTEq8,
					RTEHQ11, RTEHQ12);

			double a = ran0(&config.rng.NUM1);
			int b = 0;


			if (initialPid == 21) {
				double R0 = RTE;
				if (a <= RTEg1 / R0) {
					channel = 1; pid2 = pid3 = 21;
				} else if (a <= (RTEg1 + RTEg2) / R0) {
					channel = 2;
					do {
						b = static_cast<int>(ran0(&config.rng.NUM1) * 6);
					} while (b >= 6);
					pid2 = valencePids[b];
					pid3 = 21;
				} else {
					channel = 3;
					do {
						b = static_cast<int>(ran0(&config.rng.NUM1) * 6);
					} while (b >= 6);
					pid2 = pid3 = valencePids[b];
				}
			}
			else if (abs(initialPid) == 4) {
				double R0 = RTE;
				if (a <= RTEHQ11 / R0) {
					channel = 11;
					b = static_cast<int>(ran0(&config.rng.NUM1) * 6);
					if (b >= 6) b = 5;
					pid2 = pid3 = valencePids[b];
				} else {
					channel = 12;
					pid2 = pid3 = 21;
				}
			}
			else {
				double R0 = RTE;
				double R3 = RTEq3, R4 = RTEq4, R5 = RTEq5, R6 = RTEq6, R7 = RTEq7, R8 = RTEq8;

				if (a <= R3 / R0) {
					channel = 13;
					pid2 = pid3 = 21;
				} else if (a <= (R3 + R4) / R0) {
					channel = 4;
					do {
						b = static_cast<int>(ran0(&config.rng.NUM1) * 6);
						if (b >= 6) b = 5;
						pid2 = pid3 = valencePids[b];
					} while (pid2 == initialPid);
				} else if (a <= (R3 + R4 + R5) / R0) {
					channel = 5;
					pid2 = pid3 = initialPid;
				} else if (a <= (R3 + R4 + R5 + R6) / R0) {
					channel = 6;
					pid3 = -initialPid;
					do {
						b = static_cast<int>(ran0(&config.rng.NUM1) * 3);
						if (b >= 3) b = 2;
						pid2 = -initialPid / abs(initialPid) * valencePids[b];
					} while (abs(pid2) == abs(pid3));
				} else if (a <= (R3 + R4 + R5 + R6 + R7) / R0) {
					channel = 7;
					pid2 = pid3 = -initialPid;
				} else {
					channel = 8;
					pid3 = -initialPid;
					pid2 = 21;
				}
			}
		}


		bool sampleThermalParton(
				const int channel,
				const std::array<double, 4>& pc_jet,
				const double T,
				double& maxWeight,
				double& e2_out
				) {

			// Compute HQ mass
			double HQmass2 = pc_jet[0]*pc_jet[0] - pc_jet[1]*pc_jet[1] - pc_jet[2]*pc_jet[2] - pc_jet[3]*pc_jet[3];
			double HQmass = (HQmass2 > 1e-12) ? std::sqrt(HQmass2) : 0.0;

			// Get momentum magnitude and temperature bin indices
			double P = std::sqrt(pc_jet[1]*pc_jet[1] + pc_jet[2]*pc_jet[2] + pc_jet[3]*pc_jet[3]);
			int idx_P = std::clamp(static_cast<int>((P - config.hq22.min_p1) / config.hq22.bin_p1), 0, config.hq22.N_p1 - 1);
			int idx_T = std::clamp(static_cast<int>((T - config.hq22.min_T) / config.hq22.bin_T), 0, config.hq22.N_T - 1);

			double fBmax = config.hq22.distFncBM[idx_T][idx_P];
			double fFmax = config.hq22.distFncFM[idx_T][idx_P];

			for (int attempt = 0; attempt < 1e6; ++attempt) {
				double xw = config.hq22.max_e2 * ran0(&config.rng.NUM1);
				int idx_e2 = std::min(static_cast<int>((xw - config.hq22.min_e2) / config.hq22.bin_e2), config.hq22.N_e2 - 1);
				double fval = 0.0;

				if (channel == 11) {
					fval = config.hq22.distFncF[idx_T][idx_P][idx_e2] / fFmax;
					maxWeight =config.hq22. distMaxF[idx_T][idx_P][idx_e2];
				} else if (channel == 12) {
					fval = config.hq22.distFncB[idx_T][idx_P][idx_e2] / fBmax;
					maxWeight = config.hq22.distMaxB[idx_T][idx_P][idx_e2];
				}

				if (ran0(&config.rng.NUM1) < fval) {
					e2_out = xw * T;
					return true;
				}
			}

			e2_out = 0.0;
			return false;
		}

		bool sampleCollisionAngles(
				const int channel,
				const std::array <double, 4> &pc_jet,
				const double T,
				const double e2,
				const double qhat0ud,
				const double maxWeight,
				double& e4,
				double& theta2,
				double& theta4,
				double& phi24
				) {

			double ff = 0.0;
			double HQmass2 = pc_jet[0]*pc_jet[0] - pc_jet[1]*pc_jet[1] - pc_jet[2]*pc_jet[2] - pc_jet[3]*pc_jet[3];
			double HQmass = (HQmass2 > 1e-12) ? std::sqrt(HQmass2) : 0.0;

			// Get momentum magnitude and temperature bin indices
			double P = std::sqrt(pc_jet[1]*pc_jet[1] + pc_jet[2]*pc_jet[2] + pc_jet[3]*pc_jet[3]);
			double E1 = std::sqrt(P * P + HQmass2);  // correct heavy quark energy

			for (int i = 0; i < 1e6; ++i) {
				theta2 = base::pi * ran0(&config.rng.NUM1);
				theta4 = base::pi * ran0(&config.rng.NUM1);
				phi24  = 2.0 * base::pi * ran0(&config.rng.NUM1);

				double cos24 = std::sin(theta2) * std::sin(theta4) * std::cos(phi24) + std::cos(theta2) * std::cos(theta4);
				double down = E1 - P * std::cos(theta4) + e2 - e2 * cos24;

				e4 = (E1 * e2 - P * e2 * std::cos(theta2)) / down;

				double sigFactor = std::sin(theta2) * std::sin(theta4) * e2 * e4 / down;

				// Mandelstam variables in CM frame
				double s = 2.0 * E1 * e2 + HQmass * HQmass - 2.0 * P * e2 * std::cos(theta2);
				double t = -2.0 * e2 * e4 * (1.0 - cos24);
				double u = 2.0 * HQmass * HQmass - s - t;

				// kinematic cutoffs
				if (s <= 2.0 * qhat0ud || t >= -qhat0ud || u >= -qhat0ud) continue;

				double msq = 0.0;
				if (channel == 11) {
					ff = 1.0 / (std::exp(e2 / T) + 1.0) * (1.0 - 1.0 / (std::exp(e4 / T) + 1.0));
					msq = Mqc2qc(s, t, HQmass);
				} else if (channel == 12) {
					ff = 1.0 / (std::exp(e2 / T) - 1.0) * (1.0 + 1.0 / (std::exp(e4 / T) - 1.0));
					msq = Mgc2gc(s, t, HQmass);
				}

				sigFactor *= ff;
				double rank = ran0(&config.rng.NUM1);
				if (rank <= (msq / maxWeight) * sigFactor) {
					return true;
				}
			}

			return false;

		}



		double Mgg2gg_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) * (
					3.0 - tmin * (ss - tmin) / (ss * ss)
					+ (ss - tmin) * ss / (tmin * tmin)
					+ tmin * ss / ((ss - tmin) * (ss - tmin))
					);
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				(3.0 - tt * uu / (ss * ss)
				 + uu * ss / (tt * tt)
				 + tt * ss / (uu * uu)) / mmax;
			return msq;
		}





		double Mqg2qg_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) * (
					(4.0 / 9.0) * (tmin * tmin + (ss - tmin) * (ss - tmin)) / (tmin * (ss - tmin))
					- (tmin * tmin + (ss - tmin) * (ss - tmin)) / (ss * ss)
					);
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((4.0 / 9.0) * (tt * tt + uu * uu) / (tt * uu)
				 - (tt * tt + uu * uu) / (ss * ss)) / (mmax + 4.0);
			return msq;
		}





		double Mgq2gq_approx(double ss, double tt, double uu, double tmin, double tmax, double p0E, double p2E) {
			double mmax_a = 4.0 / (ss * ss) * (
					(ss * ss + (ss - tmin) * (ss - tmin)) / (tmin * tmin)
					+ (4.0 / 9.0) * (ss * ss + (ss - tmin) * (ss - tmin)) / (ss * (ss - tmin))
					);
			double mmax_b = 4.0 / (ss * ss) * (
					(ss * ss + (ss - tmax) * (ss - tmax)) / (tmax * tmax)
					+ (4.0 / 9.0) * (ss * ss + (ss - tmax) * (ss - tmax)) / (ss * (ss - tmax))
					);
			double mmax = std::max(mmax_a, mmax_b);

			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((ss * ss + uu * uu) / (tt * tt)
				 + (4.0 / 9.0) * (ss * ss + uu * uu) / (ss * uu)) / mmax;
			return msq;
		}




		double Mqq2qq_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) *
				((4.0 / 9.0) * (ss * ss + (ss - tmin) * (ss - tmin)) / (tmin * tmin));
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((4.0 / 9.0) * (ss * ss + uu * uu) / (tt * tt)) / mmax;
			return msq;
		}




		double Mqqbar2qqbar_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) * (
					(4.0 / 9.0) * (ss * ss + (ss - tmin) * (ss - tmin)) / (tmin * tmin)
					+ (ss * ss + tmin * tmin) / ((ss - tmin) * (ss - tmin))
					- (2.0 / 3.0) * ss * ss / (tmin * (ss - tmin))
					);
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((4.0 / 9.0) * ((ss * ss + uu * uu) / (tt * tt)
					+ (ss * ss + tt * tt) / (uu * uu))
				 - (2.0 / 3.0) * ss * ss / (tt * uu)) / mmax;
			return msq;
		}




		double Mqqbar2qqbar_diff_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) * (
					(4.0 / 9.0) * (std::pow(ss, 2) + std::pow((ss - tmin), 2)) / std::pow(tmin, 2)
					);
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((4.0 / 9.0) * (std::pow(ss, 2) + std::pow(uu, 2)) / std::pow(tt, 2)) / mmax;
			return msq;
		}


		double Mqqbar2qqbar_same_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) * (
					(4.0 / 9.0) * (std::pow(ss, 2) + std::pow((ss - tmin), 2)) / std::pow(tmin, 2)
					+ (std::pow(ss, 2) + std::pow(tmin, 2)) / std::pow((ss - tmin), 2)
					- (2.0 / 3.0) * std::pow(ss, 2) / (tmin * (ss - tmin))
					);
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((4.0 / 9.0) * ((std::pow(ss, 2) + std::pow(uu, 2)) / std::pow(tt, 2)
					+ (std::pow(ss, 2) + std::pow(tt, 2)) / std::pow(uu, 2))
				 - (2.0 / 3.0) * std::pow(ss, 2) / (tt * uu)) / mmax;
			return msq;
		}



		double Mqqbar2gg_approx(double ss, double tt, double uu, double tmin, double p0E, double p2E) {
			double mmax = 4.0 / (ss * ss) * (
					(32.0 / 27.0) * (std::pow(ss, 2) + std::pow((ss - tmin), 2)) / (tmin * (ss - tmin))
					- (8.0 / 3.0) * (std::pow(ss, 2) + std::pow((ss - tmin), 2)) / std::pow(ss, 2)
					);
			double msq = std::pow(1.0 / (2.0 * p0E * p2E), 2) *
				((32.0 / 27.0) * (std::pow(tt, 2) + std::pow(uu, 2)) / (tt * uu)
				 - (8.0 / 3.0) * (std::pow(tt, 2) + std::pow(uu, 2)) / std::pow(ss, 2)) / mmax;
			return msq;
		}




		double computeMatrixElement(
				const int channel,
				const double ss, const double tt, const double uu,
				const double tmin, const double tmax,
				const double p0E, const double p2E
				) {
			switch (channel) {
				case 1: return Mgg2gg_approx(ss, tt, uu, tmin, p0E, p2E);
				case 2: return Mqg2qg_approx(ss, tt, uu, tmin, p0E, p2E);
				case 3:
				case 13:
					return Mgq2gq_approx(ss, tt, uu, tmin, tmax, p0E, p2E);
				case 4: return Mqq2qq_approx(ss, tt, uu, tmin, p0E, p2E);
				case 5: return Mqqbar2qqbar_approx(ss, tt, uu, tmin, p0E, p2E);
				case 6: return Mqqbar2qqbar_diff_approx(ss, tt, uu, tmin, p0E, p2E);
				case 7: return Mqqbar2qqbar_same_approx(ss, tt, uu, tmin, p0E, p2E);
				case 8: return Mqqbar2gg_approx(ss, tt, uu, tmin, p0E, p2E);
				default:
					std::cerr << "Unsupported channel = " << channel << " in computeMatrixElement()\n";
					exit(EXIT_FAILURE);
			}
		}



		double getFinalStateStatFactor(const int channel, const double f1, const double f2) {
			switch (channel) {
				// Final state includes gluons → use Bose statistics
				case 1: // g + g → g + g
				case 2: // q + g → q + g
				case 13: // duplicate of 3, still gluonic → f1
					return f1;

					// Final state includes only quarks/antiquarks → use Fermi statistics
				case 4: // q + q → q + q
				case 5: // q + q̄ → q + q̄
				case 6: // q + q̄ → q + q̄ (flavor-exchange)
				case 7: // q + q̄ → q + q̄ (same-flavor)
				case 8: // q + q̄ → g + g
				case 3: // g + q → g + q //TODO: why f2?
					return f2;

				default:
					std::cerr << "Warning: channel=" << channel << " not recognized for ff assignment.\n";
					return 0.0;
			}
		}




		void colljet22(
				const int channel,
				const Particle &p,
				std::array<double, 4> &pc_rec,
				std::array<double, 4> &pc_med,
				std::array<double, 4> &pc_fin,
				double& qt                               // output: transverse momentum transfer
			      ) {

			//In original fnc.h,
			//p2[4] (p_rec) = Incoming thermal parton → sampled from the thermal distribution (e.g. Bose-Einstein)
			//p3[4] (p_med) = Recoil parton from the medium → final state of the thermal parton
			//p0[4] (p and p_fin) = colljet22 takes p0 and save it to p4, then modify p0 (final state p_jet)
			//*p4 in original fnc.h is just p0 (incoming jet).

			std::array<double, 4> pc_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
			std::array<double, 4> v_fluid = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};


			double alpha_s = alphas0(config.physics.Kalphas, p.Tfrozen);  // Assuming alphas0() computes coupling
			double qhat0ud = DebyeMass2(config.physics.Kqhat0, alpha_s, p.Tfrozen);  // qhat_0: Calculated by  \mu_D^2 = 4\pi \alpha_s T^2

			trans(v_fluid, pc_jet);

			double t;
			for (int nloop = 0; nloop < 1e6; ++nloop) {

				//Reasonable momentum sampling for pc_med (thermal parton);
				//Sample until one gets reasonable t.
				double s, f1, f2;
				for (int mloop = 0; mloop < 1e6; ++mloop) {
					double EoverT = 15.0 * ran0(&config.rng.NUM1);//Dimensionless energy x = E / T, max x=15
					double phi24  = 2.0 * base::pi * ran0(&config.rng.NUM1);//Azimuthal angle φ ∈ [0, 2π]
					double cos_ = 1.0 - 2.0 * ran0(&config.rng.NUM1);// cos(θ) ∈ [-1, 1]
					double sin_ = sqrt(1.0 - cos_ * cos_);        // sin(θ) from cos(θ)

					// Construct 4-momentum of the thermal parton
					pc_med[0] = EoverT * p.Tfrozen;                     // Energy E = x·T
					pc_med[1] = pc_med[0] * sin_ * cos(phi24);     // px = E·sin(θ)·cos(φ)
					pc_med[2] = pc_med[0] * sin_ * sin(phi24);     // py = E·sin(θ)·sin(φ)
					pc_med[3] = pc_med[0] * cos_;                  // pz = E·cos(θ)


					//Normalised distribution function
					f1 = pow(EoverT, 3) / (exp(EoverT) - 1) / 1.4215;   // Bose-Einstein weight (gluon)
					f2 = pow(EoverT, 3) / (exp(EoverT) + 1) / 1.2845;   // Fermi-Dirac weight (quark)

					//Mandelstam s = (p0 + p2)**2  = 2 p0 /dot p2 
					//(massless is assumed! TODO check mass somewhere) 
					s = 2.0 * (pc_jet[0]*pc_med[0] - pc_jet[1]*pc_med[1] - pc_jet[2]*pc_med[2] - pc_jet[3]*pc_med[3]);

					//Randomly sample t from [0, s]
					double r_ = ran0(&config.rng.NUM1);
					t = r_ * s;

					if ((t >= qhat0ud) && (t <= (s - qhat0ud))){
						//Accepted!
						break;
					}


				}

				double tmin = qhat0ud;
				double tmid = s / 2.0;
				double tmax = s - qhat0ud;


				//TODO ??? double u = -s - t  
				double u = s - t;

				// Matrix element
				double msq = computeMatrixElement(channel, s, t, u, tmin, tmax, pc_jet[0], pc_med[0]);
				double ff = getFinalStateStatFactor(channel, f1, f2);

				double accept = msq * ff;
				if (ran0(&config.rng.NUM1) <= accept) {
					break;
				}

			}

			//Passing...just to initialize pc_rec and pc_fin
			pc_rec = pc_med;

			//Calculate cm velocity in 
			std::array <double, 4> v_cm =  get_centerofmass(pc_jet, pc_rec);
			trans(v_cm, pc_jet);
			trans(v_cm, pc_rec);

			//Copy initial info
			pc_fin = pc_jet;
			LongitudinalMomentumTransfer(t, pc_jet, pc_rec, pc_fin);

			transback(v_cm, pc_jet);
			transback(v_cm, pc_rec);
			transback(v_cm, pc_fin);

			//     calculate qt in the rest frame of medium
			rotate(pc_jet[1], pc_jet[2], pc_jet[3], pc_fin, 1);
			double qT = sqrt(pc_fin[1] * pc_fin[1] + pc_fin[2] * pc_fin[2]);
			qt = qT;
			rotate(pc_jet[1], pc_jet[2], pc_jet[3], pc_fin, -1);


			// Transform all back to lab frame
			transback(v_fluid, pc_rec);
			transback(v_fluid, pc_med);
			transback(v_fluid, pc_fin);
			transback(v_fluid, pc_jet);
			//std::cout << "pc_rec " <<  pc_rec[0] << "  " << pc_rec[1] << "  " << pc_rec[2] << "  " << pc_rec[3] << std::endl;
			//std::cout << "pc_med " <<  pc_med[0] << "  " << pc_med[1] << "  " << pc_med[2] << "  " << pc_med[3] << std::endl;



			return;
		}



		void LongitudinalMomentumTransfer(const double t, const std::array<double, 4> pc_jet, std::array<double, 4>& pc_rec, std::array<double, 4>& pc_fin){


			double pcm = pc_rec[0];
			double ran_p_=2.0*base::pi*ran0(&config.rng.NUM1);

			double q_T=sqrt(pcm*pcm-(t/2.0/pcm-pcm)*(t/2.0/pcm-pcm));
			double qx=q_T*cos(ran_p_);
			double qy=q_T*sin(ran_p_);

			double q_L=t/2.0/pcm;

			// Compute velocity components of particle 2 in the lab frame
			const double E2 = pc_rec[0];
			const double px2 = pc_rec[1];
			const double py2 = pc_rec[2];
			const double pz2 = pc_rec[3];

			const double transverseVelocity = std::sqrt(px2 * px2 + py2 * py2) / E2;
			const double vx = px2 / E2;
			const double vy = py2 / E2;
			const double vz = pz2 / E2;

			// Subtract longitudinal momentum transfer
			pc_rec[1] -= q_L * vx;
			pc_rec[2] -= q_L * vy;

			// Add transverse momentum transfer, if transverse velocity is non-zero
			if (transverseVelocity != 0.0) {
				pc_rec[1] += (vz * vx * qy + vy * qx) / transverseVelocity;
				pc_rec[2] += (vz * vy * qy - vx * qx) / transverseVelocity;
			}

			// Final z-component momentum update (labelled as s2 in original code)
			pc_rec[3] -= q_L * vz + transverseVelocity * qy;


			//Final state jet particle
			pc_fin[1] = -pc_rec[1]; 
			pc_fin[2] = -pc_rec[2]; 
			pc_fin[3] = -pc_rec[3]; 

		}



		void collHQ22(
				int channel,
				const Particle& p, // incoming particle
				Particle& p_rec,  // output: recoiled thermal parton
				Particle& p_med,  // output: initial thermal medium parton
				Particle& p_fin, 
				double &qt         // output: transverse momentum transfer
			     ) {

			//Only here I am using p.P.
			//Momentum should be put back at the end of this function.
			std::array<double, 4> pc_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
			std::array<double, 4> v_fluid = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
			std::array<double, 4> pc_rec = {0., 0., 0., 0.};
			std::array<double, 4> pc_med = {0., 0., 0., 0.};
			std::array<double, 4> pc_fin = {0.,0.,0.,0.};

			// Transform heavy quark to fluid rest frame
			trans(v_fluid, pc_jet);


			//Keep info
			pc_fin = pc_jet;
			double mass2 = pc_jet[0]*pc_jet[0] 
				- pc_jet[1]*pc_jet[1] 
				- pc_jet[2]*pc_jet[2] 
				- pc_jet[3]*pc_jet[3];
			double mass = (mass2 > 1e-12) ? std::sqrt(mass2) : 0.0;

			double maxWeight, e2;
			bool sampleOK = sampleThermalParton(channel, pc_jet, p.Tfrozen, maxWeight, e2);
			if (!sampleOK) {
				qt = 0.0;
				pc_rec.fill(0.0);
				pc_med.fill(0.0);
				transback(v_fluid, pc_jet);
				//transback(v0, pc_init);
				return;
			}

			// Sample angles
			double qhat0ud, e4, theta2, theta4, phi24;
			bool anglesOK = sampleCollisionAngles(
					channel, pc_jet, p.Tfrozen, e2, qhat0ud, maxWeight,
					e4, theta2, theta4, phi24
					);
			if (!anglesOK) {
				qt = 0.0;
				pc_rec.fill(0.0);
				pc_med.fill(0.0);
				transback(v_fluid, pc_jet);
				//transback(v0, pc_init);
				return;
			}


			// Final kinematics
			pc_med[0] = e2;
			pc_med[1] = e2 * std::sin(theta2);
			pc_med[2] = 0.0;
			pc_med[3] = e2 * std::cos(theta2);

			pc_rec[0] = e4;
			pc_rec[1] = e4 * std::sin(theta4) * std::cos(phi24);
			pc_rec[2] = e4 * std::sin(theta4) * std::sin(phi24);
			pc_rec[3] = e4 * std::cos(theta4);

			// Back-rotate from jet axis
			rotate(pc_jet[1], pc_jet[2], pc_jet[3], pc_rec, -1);
			rotate(pc_jet[1], pc_jet[2], pc_jet[3], pc_med, -1);

			// Update HQ momentum via conservation
			for (int j = 1; j <= 3; ++j)
				pc_fin[j] = pc_jet[j] + pc_med[j] - pc_rec[j];
			pc_fin[0] = std::sqrt(pc_fin[1]*pc_fin[1] 
					+ pc_fin[2]*pc_fin[2] 
					+ pc_fin[3]*pc_fin[3] 
					+ mass*mass);

			// Transverse momentum transfer
			rotate(pc_jet[1], pc_jet[2], pc_jet[3], pc_fin, 1);
			qt = std::sqrt(pc_fin[1]*pc_fin[1] + pc_fin[2]*pc_fin[2]);
			rotate(pc_jet[1], pc_jet[2], pc_jet[3], pc_fin, -1);

			// Back-transform all to lab frame
			transback(v_fluid, pc_rec);
			transback(v_fluid, pc_med);
			transback(v_fluid, pc_jet);
			transback(v_fluid, pc_fin);

			for(int i=0; i<4; i++){
				p_rec.V[i] = pc_rec[i];
				p_med.V[i] = pc_med[i];
			}


			std::cout << "pc_med(p3) " <<  pc_med[0] << "  " << pc_med[1] << "  " << pc_med[2] << "  " << pc_med[3] << std::endl;
			std::cout << "pc_rec(p2) " <<  pc_rec[0] << "  " << pc_rec[1] << "  " << pc_rec[2] << "  " << pc_rec[3] << std::endl;
			std::cout << "pc_fin(p0) " <<  pc_fin[0] << "  " << pc_fin[1] << "  " << pc_fin[2] << "  " << pc_fin[3] << std::endl;

exit(1);
		}

		bool collHQ23(const Particle &p, const double qt, 
				const std::array<double, 4> pc_med, 
				std::array<double, 4> &pc_rec, 
				std::array<double, 4> &pc_rad,
				std::array<double, 4> &pc_fin
			     ){

			std::array<double, 4> pc_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
			std::array<double, 4> v_flow = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};


			//Copying pc_med to pc_rec(p2: initial thermal momentum, will be final thermal momentum)
			pc_rec = pc_med;
			pc_fin = pc_jet;
			trans(v_flow, pc_fin);
			trans(v_flow, pc_rec);

			//At rest frame check if energy is too small.
			double alpha_s = alphas0(config.physics.Kalphas, p.Tfrozen);  // Assuming alphas0() computes coupling
			double qhat0 = DebyeMass2(config.physics.Kqhat0, alpha_s, p.Tfrozen);  // qhat_0: Calculated by  \mu_D^2 = 4\pi \alpha_s T^2
			double Eloc = pc_fin[0];//Eloc == HQenergy (energy of the jet particle at fluid rest frame)
			double mass2 = pc_fin[0]*pc_fin[0]  
				- pc_fin[1]*pc_fin[1]  
				- pc_fin[2]*pc_fin[2]  
				- pc_fin[3]*pc_fin[3];
			double mass = (mass2>1e-12)? 
				sqrt(mass2):0.0;
			if(Eloc<2.0*sqrt(qhat0)){
				return false;//no radiation!
			}

			//Rotate pc_rec(med)  with pc_jet
			std::array <double, 4> zDir = {pc_fin[0], pc_fin[1], pc_fin[2], pc_fin[3]};
			rotate(zDir[1], zDir[2], zDir[3], pc_rec, 1);



			// --- Phase space limitation ---
			double lim_low = sqrt(6.0 * base::pi * alpha_s) * p.Tfrozen / Eloc;
			double lim_high = (abs(p.pid) == 4) ? 1.0 : (1.0 - lim_low);  // heavy quarks allow full range
			double lim_int = lim_high - lim_low;

			int nloopOut=0;
			bool DONE = false;
			//Sampling!
			do {

				double randomX, randomY;
				do {
					randomX=lim_low+lim_int*ran0(&config.rng.NUM1);
					randomY=ran0(&config.rng.NUM1);
				} while(tau_f(randomX,randomY,Eloc,mass)<1.0/base::pi/p.Tfrozen);

				int count = 0;
				while(p.max_Ng*ran0(&config.rng.NUM1)>dNg_over_dxdydt(p.pid,randomX,randomY,Eloc,mass,p.Tfrozen, alpha_s, p.qhat_over_T3, p.Tint_lrf)) {
					count++;
					if(count>1e+5) {
						std::cout << "Give up loop at point 1 ..." << std::endl;
						return false;
					}

					do {
						randomX=lim_low+lim_int*ran0(&config.rng.NUM1);
						randomY=ran0(&config.rng.NUM1);
					} while(tau_f(randomX,randomY,Eloc,mass)<1.0/base::pi/p.Tfrozen);
				}


				if(p.pid==21&&randomX>0.5) randomX=1.0-randomX;
				double theta_gluon=2.0*base::pi*ran0(&config.rng.NUM1);
				double kperp_gluon=randomX*randomY*Eloc;
				pc_rad[1]=kperp_gluon*cos(theta_gluon);
				pc_rad[2]=kperp_gluon*sin(theta_gluon);
				pc_rad[3]=randomX*Eloc*sqrt(1.0-randomY*randomY);
				pc_rad[0]=sqrt(pc_rad[1]*pc_rad[1]+pc_rad[2]*pc_rad[2]+pc_rad[3]*pc_rad[3]);


				if(pc_rad[0]>(Eloc-mass)) {
					pc_rad[1]=0.0;
					pc_rad[2]=0.0;
					pc_rad[3]=0.0;
					pc_rad[0]=0.0;
					nloopOut++;
					continue;
				}


				//Solve energy conservation
				double sqx,sqy,sqzA,sq0A,sqzB,sq0B,sqz,sq0,sqtheta;
				int nloop1=0;
				int nloop2=0;
				bool doneA = false;
				bool doneB = false;
				bool doneAB = false;


				double sE1=pc_fin[0];
				double sp1x=0.0;
				double sp1y=0.0;
				double sp1z=sqrt(pc_fin[1]*pc_fin[1]+pc_fin[2]*pc_fin[2]+pc_fin[3]*pc_fin[3]);
				double sE2=pc_rec[0];
				double sp2x=pc_rec[1];
				double sp2y=pc_rec[2];
				double sp2z=pc_rec[3];
				double sk0=pc_rad[0];
				double skx=pc_rad[1];
				double sky=pc_rad[2];
				double skz=pc_rad[3];

				do {
					double sqtheta=2.0*base::pi*ran0(&config.rng.NUM1);
					sqx=qt*cos(sqtheta);
					sqy=qt*sin(sqtheta);
					double sAA=(sE1+sE2-sk0)/(sp1z+sp2z-skz);
					double sBB=(pow(sE2,2)-pow((sE1-sk0),2)+pow((sp1x-sqx-skx),2)
							+pow((sp1y-sqy-sky),2)+pow((sp1z-skz),2)+pow(mass,2)
							-pow((sp2x+sqx),2)-pow((sp2y+sqy),2)-pow(sp2z,2))/2.0/(sp1z+sp2z-skz);
					double aaa=sAA*sAA-1.0;
					double bbb=2.0*(sAA*sp2z+sAA*sBB-sE2);
					double ccc=pow((sp2x+sqx),2)+pow((sp2y+sqy),2)+sp2z*sp2z+2.0*sp2z*sBB+sBB*sBB-sE2*sE2;
					double abc2=bbb*bbb-4.0*aaa*ccc;

					if(abc2<0.0) {
						nloop1++;
						continue;
					} else {
						nloop1=0;
					}



					double abc=sqrt(abc2);
					sq0A=(-bbb+abc)/2.0/aaa;
					sq0B=(-bbb-abc)/2.0/aaa;
					sqzA=sAA*sq0A+sBB;
					sqzB=sAA*sq0B+sBB;

					if(sq0A*sq0A-sqx*sqx-sqy*sqy-sqzA*sqzA<0 && sE1-sq0A-sk0>mass && sE2+sq0A>0) {
						doneA=true;

					} else {
						doneA=false;
					}
					if(sq0B*sq0B-sqx*sqx-sqy*sqy-sqzB*sqzB<0 && sE1-sq0B-sk0>mass && sE2+sq0B>0) {
						doneB=true;
					} else {
						doneB=false;
					}

					if(!doneA && !doneB) {
						nloop2++;
						continue;
					} else if(doneA && doneB) {
						if(abs(sq0A)<abs(sq0B)) {
							sq0=sq0A;
							sqz=sqzA;
						} else {
							sq0=sq0B;
							sqz=sqzB;
						}
					} else if(doneA) {
						sq0=sq0A;
						sqz=sqzA;
					} else {
						sq0=sq0B;
						sqz=sqzB;
					}

					doneAB=true;

				} while (!doneAB && nloop1<config.counter.loopN && nloop2<config.counter.loopN);


				if(!doneAB) {
					nloopOut++;
					continue;
				} else {
					pc_fin[0]=sE1-sq0-sk0;
					pc_fin[1]=sp1x-sqx-skx;
					pc_fin[2]=sp1y-sqy-sky;
					pc_fin[3]=sp1z-sqz-skz;

					pc_rec[0]=sE2+sq0;
					pc_rec[1]=sp2x+sqx;
					pc_rec[2]=sp2y+sqy;
					pc_rec[3]=sp2z+sqz;


					rotate(zDir[1],zDir[2],zDir[3],pc_fin,-1); // rotate p0 into global system
					rotate(zDir[1],zDir[2],zDir[3],pc_rec,-1); // rotate p2 into global system
					rotate(zDir[1],zDir[2],zDir[3],pc_rad,-1); // rotate p4 into global system                                                                                                                                                                                  transback(v0,p0);
					transback(v_flow,pc_fin);
					transback(v_flow,pc_rec);
					transback(v_flow,pc_rad);

					//TODO: debug: check on-shell condition

					DONE=true;
				}

			} while(!DONE && nloopOut<config.counter.loopN);

			//std::cout << "pc_rec " <<  pc_rec[0] << "  " << pc_rec[1] << "  " << pc_rec[2] << "  " << pc_rec[3] << std::endl;
			//std::cout << "pc_fin " <<  pc_fin[0] << "  " << pc_fin[1] << "  " << pc_fin[2] << "  " << pc_fin[3] << std::endl;
			//std::cout << "pc_rad " <<  pc_rad[0] << "  " << pc_rad[1] << "  " << pc_rad[2] << "  " << pc_rad[3] << std::endl;

			if(!DONE) return false;
			else return true;

		}




		//int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int &ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt
		bool radiationHQ(const Particle &p,
				std::array<double, 4> pc_rad0, 
				std::array<double, 4> pc_fin0, 
				std::array<double, 4> &pc_rad1,
				std::array<double, 4> &pc_rad2,
				std::array<double, 4> &pc_fin12
			     ){

			std::array<double, 4> pc_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
			std::array<double, 4> v_flow = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};

			//    work in the rest frame of medium
			//    return the 4-momentum of final states in 1->3 radiation
			//    input: P2(4)-momentum of radiated gluon from 2->3 (p_rad0 = p_rad in previous radiation)
			//           P3(4)-momentum of daughter parton from 2->3 (p_fin0 = p_fin in previous radiation)
			//           Pj0(4)-inital momentum of jet before 2->3 (p_jet incoming leading parton in 2->2)
			//           v0-local velocity of medium
			//    output:P2(4)-momentum of 1st radiated gluon (pc_rad1)
			//           P3(4)-momentum of daughter parton (pc_fin12)
			//           P4(4)-momentum of 2nd radiated gluon(pc_rad2)
			//           i=1: no radiation; 0:successful radiation

			std::cout << "pc_jet " <<  pc_jet[0] << "  " << pc_jet[1] << "  " << pc_jet[2] << "  " << pc_jet[3] << std::endl;
			std::cout << "pc_fin0 " <<  pc_fin0[0] << "  " << pc_fin0[1] << "  " << pc_fin0[2] << "  " << pc_fin0[3] << std::endl;
			std::cout << "pc_rad0 " <<  pc_rad0[0] << "  " << pc_rad0[1] << "  " << pc_rad0[2] << "  " << pc_rad0[3] << std::endl;
exit(1);

			double mass2=(pc_fin0[0]*pc_fin0[0]
					-pc_fin0[1]*pc_fin0[1]
					-pc_fin0[2]*pc_fin0[2]
					-pc_fin0[3]*pc_fin0[3]);
			double mass=(p.pid==4 && mass2>1e-12)? sqrt(mass2): 0.0;

                        if (p.pid!=4){
				if(fabs(pc_fin0[0]-sqrt(pc_fin0[1]*pc_fin0[1]+pc_fin0[2]*pc_fin0[2]+pc_fin0[3]*pc_fin0[3]))>base::epsilon){
					std::cout << "ERROR!  EE - PP != 0 :  " << pc_fin0[0]*pc_fin0[0] 
						<< "   " <<  pc_fin0[1]*pc_fin0[1]+pc_fin0[2]*pc_fin0[2]+pc_fin0[3]*pc_fin0[3] 
						<< std::endl;
					exit(EXIT_FAILURE);
				}
			}


			// transform to local comoving frame of the fluid
			trans(v_flow,pc_rad0);
			trans(v_flow,pc_fin0);
			trans(v_flow,pc_jet);


			//At rest frame check if energy is too small.
			double alpha_s = alphas0(config.physics.Kalphas, p.Tfrozen);  // Assuming alphas0() computes coupling
			double qhat0 = DebyeMass2(config.physics.Kqhat0, alpha_s, p.Tfrozen);  // qhat_0: Calculated by  \mu_D^2 = 4\pi \alpha_s T^2
			double Eloc = pc_fin0[0];//Eloc == HQenergy (energy of the jet particle at fluid rest frame)



			double lim_low = sqrt(6.0 * base::pi * alpha_s) * p.Tfrozen / Eloc;
			double lim_int=(pc_fin0[0]-mass)/Eloc-lim_low;

			if(lim_int<=0.0) return false;

			double px0=pc_rad0[1]+pc_fin0[1];
			double py0=pc_rad0[2]+pc_fin0[2];
			double pz0=pc_rad0[3]+pc_fin0[3];

			// rotate to the frame in which jet moves along z-axis
			rotate(px0,py0,pz0,pc_fin0,1);
			rotate(px0,py0,pz0,pc_rad0,1);
			rotate(px0,py0,pz0,pc_jet,1);


			int nloopOut=0;
			bool DONE=false;
			//Sampling!
			do {
				double randomX, randomY;
				do {
					randomX=lim_low+lim_int*ran0(&config.rng.NUM1);
					randomY=ran0(&config.rng.NUM1);
				} while(tau_f(randomX,randomY,Eloc,mass)<1.0/base::pi/p.Tfrozen);   

				int count=0;
				while(p.max_Ng*ran0(&config.rng.NUM1)>dNg_over_dxdydt(p.pid,randomX,randomY,Eloc,mass,p.Tfrozen, alpha_s, p.qhat_over_T3, p.Tint_lrf)) {
					count++;
					if(count>1e+5) {
						std::cout << "Give up loop at point 1 ..." << std::endl;
						return false;
					}

					do {
						randomX=lim_low+lim_int*ran0(&config.rng.NUM1);
						randomY=ran0(&config.rng.NUM1);
					} while(tau_f(randomX,randomY,Eloc,mass)<1.0/base::pi/p.Tfrozen);
				}

				if(p.pid==21&&randomX>0.5) randomX=1.0-randomX;
				double theta_gluon=2.0*base::pi*ran0(&config.rng.NUM1);
				double kperp_gluon=randomX*randomY*Eloc;
				pc_rad2[1]=kperp_gluon*cos(theta_gluon);
				pc_rad2[2]=kperp_gluon*sin(theta_gluon);
				pc_rad2[3]=randomX*Eloc*sqrt(1.0-randomY*randomY);
				pc_rad2[0]=sqrt(pc_rad2[1]*pc_rad2[1]+pc_rad2[2]*pc_rad2[2]+pc_rad2[3]*pc_rad2[3]);

				rotate(pc_jet[1],pc_jet[2],pc_jet[3],pc_rad2,-1);

				if(pc_rad2[0]>(pc_fin0[0]-mass)) {
					// which should be impossible due to reset of lim_int above
					std::cout << "Something is weired ... pc_rad2[0] > (pc_fin0[0]-mass) " << pc_rad2[0] << " < " << (pc_fin0[0]-mass) << std::endl;
					exit(EXIT_FAILURE);
				}

				// solve energy-momentum conservation
				// p0 is re-constructed off-shell parton, correponding to P5
				// p1 is the final heavy quark from 2->3, corresponding to P3 (input). P3 (output) is the final heavy quark after 1->3.
				// k1 is the first gluon, corresponding to P2
				// k2 is the second gluon, corresponding to P4
				// assume k10 and p10 unchanged and modify their other components while p0 and k2 are fixed


				double sk1z,sk1p,sk1x,sk1y,sktheta,stheta12;
				int nloop1=0;
				int nloop2=0;
				bool doneA = false;
				bool doneB = false;
				bool doneAB=false;

				std::array <double, 4> Psum = {
					pc_rad0[0]+pc_fin0[0],
					pc_rad0[1]+pc_fin0[1],
					pc_rad0[2]+pc_fin0[2],
					pc_rad0[3]+pc_fin0[3],
				};
				double sp0z=sqrt(Psum[1]*Psum[1]+Psum[2]*Psum[2]+Psum[3]*Psum[3]);
				double sk10=pc_rad0[0];
				double sk1zOld=pc_rad0[3];
				double sp10=pc_fin0[0];
				double sk20=pc_rad2[0];
				double sk2x=pc_rad2[1];
				double sk2y=pc_rad2[2];
				double sk2z=pc_rad2[3];
				double sk2p=sqrt(sk2x*sk2x+sk2y*sk2y);

				double sk1z1,sk1z2,sk1p1,sk1p2;

				double sAA=sk10*sk10+sk2p*sk2p+sp0z*sp0z+sk2z*sk2z-2.0*sp0z*sk2z+mass*mass-(sp10-sk20)*(sp10-sk20);
				double cos12Min2=(sAA*sAA/4.0/sk10/sk10-(sp0z-sk2z)*(sp0z-sk2z))/sk2p/sk2p;

				if(cos12Min2>1.0) {
					nloopOut++;
					continue;
				}

				do {
					stheta12=2.0*base::pi*ran0(&config.rng.NUM1); // theta between k1 and k2
					double aaa=4.0*((sp0z-sk2z)*(sp0z-sk2z)+sk2p*sk2p*cos(stheta12)*cos(stheta12));
					double bbb=-4.0*sAA*(sp0z-sk2z);
					double ccc=sAA*sAA-4.0*sk10*sk10*sk2p*sk2p*cos(stheta12)*cos(stheta12);

					double abc2=bbb*bbb-4.0*aaa*ccc;

					if(abc2<0.0) {
						nloop1++;
						continue;
					} else {
						nloop1=0;
					}

					double abc=sqrt(abc2);
					sk1z1=(-bbb+abc)/2.0/aaa;
					sk1z2=(-bbb-abc)/2.0/aaa;
					sk1p1=sqrt(sk10*sk10-sk1z1*sk1z1);
					sk1p2=sqrt(sk10*sk10-sk1z2*sk1z2);

					// Since we have squared both sides of the original equation during solving, the solutions may not satisfy the original equation. Double check is needed.


					// require time-like of p1 and k10>k1z;
					if(2.0*sk1p1*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z1+sAA<0.000001 && sk10>abs(sk1z1) && sp10*sp10-(sp0z-sk1z1)*(sp0z-sk1z1)-(sk10*sk10-sk1z1*sk1z1)>mass*mass) {
						doneA=true; 
					} else {
						doneA=false;
					}
					if(2.0*sk1p2*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z2+sAA<0.000001 && sk10>abs(sk1z2) && sp10*sp10-(sp0z-sk1z2)*(sp0z-sk1z2)-(sk10*sk10-sk1z2*sk1z2)>mass*mass) {
						doneB=true; 
					} else {
						doneB=false;
					}

					// select appropriate solution
					if(!doneA && !doneB) {
						//          cout << "Solutions fail ..." << endl;
						nloop2++;
						continue;
					} else if(doneA && doneB) {
						//          cout << "Both solutions work!" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
						if(abs(sk1z1-sk1zOld)<abs(sk1z2-sk1zOld)) {
							sk1z=sk1z1;
						} else {
							sk1z=sk1z2;
						}

					} else if(doneA) {
						//          cout << "pass A ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
						sk1z=sk1z1;
					} else {
						//          cout << "pass B ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
						sk1z=sk1z2;
					}

					doneAB=true;   

					sk1p=sqrt(sk10*sk10-sk1z*sk1z);
					//           cout << "check solution: " << pow(sp10-sk20,2)-pow(sk1p,2)-pow(sk2p,2)-2.0*sk1p*sk2p*cos(stheta12)-pow(sp0z-sk1z-sk2z,2)-pow(mass,2) <<"  "<< aaa*sk1z*sk1z+bbb*sk1z+ccc <<"  "<<2.0*sk1p*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z+sAA<<"  "<<2.0*sk1p*sk2p*cos(stheta12)+2.0*(sp0z-sk2z)*sk1z-sAA << endl;

				} while (!doneAB && nloop1<config.counter.loopN && nloop2<config.counter.loopN);

				if(!doneAB) { 
					nloopOut++;
					continue;
				} else {
					sktheta=atan2(sk2y,sk2x); // theta for k2
					sktheta=sktheta+stheta12; // theta for k1
					sk1p=sqrt(sk10*sk10-sk1z*sk1z);
					sk1x=sk1p*cos(sktheta);
					sk1y=sk1p*sin(sktheta);

					pc_rad1[0]=sk10;
					pc_rad1[1]=sk1x;
					pc_rad1[2]=sk1y;
					pc_rad1[3]=sk1z;

					pc_fin12[0]=sp10-sk20;
					pc_fin12[1]=-sk1x-sk2x;
					pc_fin12[2]=-sk1y-sk2y;
					pc_fin12[3]=sp0z-sk1z-sk2z;

					// rotate back to the local coordinate
					rotate(px0,py0,pz0,pc_rad1,-1);
					rotate(px0,py0,pz0,pc_fin12,-1);
					rotate(px0,py0,pz0,pc_rad2,-1);
					rotate(px0,py0,pz0,pc_jet,-1);

					// boost back to the global frame 
					transback(v_flow,pc_rad1);
					transback(v_flow,pc_fin12);
					transback(v_flow,pc_rad2);
					transback(v_flow,pc_jet);

					// debug: check on-shell condition of P2, pc_fin12 and P4
					if(abs(pc_rad1[0]*pc_rad1[0]
								-pc_rad1[1]*pc_rad1[1]
								-pc_rad1[2]*pc_rad1[2]
								-pc_rad1[3]*pc_rad1[3]) >0.000001 || 
							abs(pc_fin12[0]*pc_fin12[0]
								-pc_fin12[1]*pc_fin12[1]
								-pc_fin12[2]*pc_fin12[2]
								-pc_fin12[3]*pc_fin12[3]
								-mass*mass)>0.000001 
							|| abs(pc_rad2[0]*pc_rad2[0]
								-pc_rad2[1]*pc_rad2[1]
								-pc_rad2[2]*pc_rad2[2]
								-pc_rad2[3]*pc_rad2[3])>0.000001) {
						std::cout << "Wrong solution -- not on shell" << "  " << sk10 << "  " << sk1x << "  " << sk1y << "  " << sk1z << "  " << sk20 << "  " << sk2x << "  " << sk2y << "  " << sk2z << "  " << stheta12 << "  " << sp10 << "  " << sp10-sk20 << "  " << -sk1x-sk2x << "  " << -sk1y-sk2y << "  " << sp0z << "  " << sp0z-sk1z-sk2z << "  " <<mass<< "  "<<pow(sp10-sk20,2)-pow(sk1x+sk2x,2)-pow(sk1y+sk2y,2)-pow(sp0z-sk1z-sk2z,2)-pow(mass,2)<<std::endl;
						std::cout << abs(pc_rad1[0]*pc_rad1[0]-pc_rad1[1]*pc_rad1[1]-pc_rad1[2]*pc_rad1[2]-pc_rad1[3]*pc_rad1[3]) <<"  "<< abs(pc_fin12[0]*pc_fin12[0]-pc_fin12[1]*pc_fin12[1]-pc_fin12[2]*pc_fin12[2]-pc_fin12[3]*pc_fin12[3]-mass*mass) <<"  "<< abs(pc_rad2[0]*pc_rad2[0]-pc_rad2[1]*pc_rad2[1]-pc_rad2[2]*pc_rad2[2]-pc_rad2[3]*pc_rad2[3]) <<std::endl;
					}

					DONE=true;

					// SC: check energy-momentum conservation
					for(int i=0; i<=3; i++) {
						if(0 && abs(pc_rad1[i]+pc_fin12[i]-pc_rad1[i]-pc_fin12[i]-pc_rad2[i])>0.000001) {
							std::cout << "Warning: Violation of E.M. conservation!  " << i << " " << abs(pc_rad1[i]+pc_fin12[i]-pc_rad1[i]-pc_fin12[i]-pc_rad2[i]) << std::endl;
						}
					}

					// SC: check on-shell condition
					double shell2=abs(pc_rad1[0]*pc_rad1[0]-pc_rad1[1]*pc_rad1[1]-pc_rad1[2]*pc_rad1[2]-pc_rad1[3]*pc_rad1[3]);
					double shell3=abs(pc_fin12[0]*pc_fin12[0]-pc_fin12[1]*pc_fin12[1]-pc_fin12[2]*pc_fin12[2]-pc_fin12[3]*pc_fin12[3]-mass*mass);
					double shell4=abs(pc_rad2[0]*pc_rad2[0]-pc_rad2[1]*pc_rad2[1]-pc_rad2[2]*pc_rad2[2]-pc_rad2[3]*pc_rad2[3]);
					if(shell2>0.000001 || shell3>0.000001 || shell4>0.000001) {
						std::cout << "Warning: Violation of on-shell: " << shell2 << "  " << shell3 << "  " << shell4 << std::endl;
					}

				}

			} while(!DONE && nloopOut<config.counter.loopN);

			if(!DONE) return false;
			else return true;

		}











	public:


		void LBT(std::vector<Particle> &particles, double ti);
		void set_np_snapshot(const int np_snapshot_in){this->np_snapshot = np_snapshot_in;}

		LBTcl(LBTConfig& config_in):config(config_in){

		};
		~LBTcl(){};

};

extern "C" {
	void read_ccnu_(char *dataFN_in, int len1);
	void hydroinfoccnu_(double *Ct, double *Cx, double *Cy, double *Cz, double *Ctemp, double *Cvx, double *Cvy, double *Cvz, int *Cflag);

	void sethydrofilesez_(int *dataID_in, char *dataFN_in, int *ctlID_in, char *ctlFN_in, int *bufferSize, int len1, int len2);
	void readhydroinfoshanshan_(double *t, double *x, double *y, double *z, double *e, double *s, double *temp, double *vx, double *vy, double *vz, int *flag);
}
#endif
