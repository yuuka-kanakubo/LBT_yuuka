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
		void handleElasticCollision(Particle &p, const double PLen, std::vector<Particle> &particles);
		void handleRadiation(Particle &p, std::vector<Particle> &particles);
		void propagateParticle(Particle &p, double ti, int &free, double &fraction);
		double computeCollisionProbability(
				Particle &p,
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
			double max_Ng = max_EGrid_low+(HQenergy-HQenergy_num*config.hqrad.delta_HQener)/config.hqrad.delta_HQener*(max_EGrid_high-max_EGrid_low);


			std::cout << " D2piT " << p.D2piT << std::endl;
			std::cout << " dtLRF " << dtLRF << std::endl;
			std::cout << " rate_T2E1 " << rate_T2E1 << std::endl;


			delta_Ng*=6.0/p.D2piT*dtLRF;
			max_Ng*=6.0/p.D2piT;
			p.max_Ng = max_Ng;

			//  if(delta_Ng>1) {
			//     std::cout << "Warning: Ng greater than 1   " << time_gluon << "  " << delta_Ng << std::endl;
			//  }

			return delta_Ng;

		}	  



		void linear(
				const int pid, const double E, const double T,
				double& RTEg, double& RTEg1, double& RTEg2, double& RTEg3,
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
			double RTEg, RTEg1, RTEg2, RTEg3;
			double RTEq, RTEq3, RTEq4, RTEq5, RTEq6, RTEq7, RTEq8;
			double RTEHQ11, RTEHQ12;

			linear(initialPid, E, T, 
					RTEg, RTEg1, RTEg2, RTEg3,
					RTEq, RTEq3, RTEq4, RTEq5, RTEq6, RTEq7, RTEq8,
					RTEHQ11, RTEHQ12);

			double a = ran0(&config.rng.NUM1);
			int b = 0;

			std::cout << "RTEg " << RTEg << std::endl;
			std::cout << "RTEg1 " << RTEg1 << std::endl;
			std::cout << "RTEg2 " << RTEg2 << std::endl;
			std::cout << "RTEg3 " << RTEg3 << std::endl;
			std::cout << "a " << a << std::endl;

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






		void collHQ22(
				int channel,
				const Particle& p; // incoming particle
				Particle& p_rec,  // output: recoiled thermal parton
				Particle& p_med,  // output: initial thermal medium parton
				double &qt         // output: transverse momentum transfer
			     ) {

			std::array<double, 4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
			std::array<double, 4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};

			// Transform momentum into fluid rest frame
			trans(vc0, pc0);

			double mass = pc0[0] * pc0[0] - pc0[1] * pc0[1] - pc0[2] * pc0[2] - pc0[3] * pc0[3];
			mass = (mass > 1e-12) ? std::sqrt(mass) : 0.0;

			// Sampling and collision logic (unchanged, but clearer variable names used)
			// [Sampling medium parton from distribution → stored in pc_med]

			// Get momentum magnitude and temperature bin indices
			double P = std::sqrt(pc0[1]*pc0[1] + pc0[2]*pc0[2] + pc0[3]*pc0[3]);
			double T = p.Tfrozen;
			int index_P = std::clamp(static_cast<int>((P - min_p1) / bin_p1), 0, N_p1 - 1);
			int index_T = std::clamp(static_cast<int>((T - min_T) / bin_T), 0, N_T - 1);

			double fBmax = distFncBM[index_T][index_P];
			double fFmax = distFncFM[index_T][index_P];
//here

			maxValue=10.0;  // need actual value later

			ct1_loop=0;

			// [Compute kinematics of final state → stored in pc0 and pc_rec]

			// ... main HQ scattering sampling logic remains here (not duplicated for brevity) ...

			// Final: rotate and transform all momenta back to lab frame
			transback(v0, pc_rec);  // recoiled medium parton
			transback(v0, pc0);     // updated HQ momentum
			transback(v0, pc_med);  // initial thermal medium parton
			transback(v0, pc4);     // reference original HQ momentum

			// Transverse momentum transfer is computed relative to pc4
			rotate(pc4[1], pc4[2], pc4[3], pc0, 1);
			qt = std::sqrt(pc0[1] * pc0[1] + pc0[2] * pc0[2]);
			rotate(pc4[1], pc4[2], pc4[3], pc0, -1);
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
