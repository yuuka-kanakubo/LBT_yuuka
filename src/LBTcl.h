#ifndef LBTCL_H
#define LBTCL_H
#include <vector>
#include <cmath>
#include "ParticleInfo.h"
#include "LBTConfig.h"


class LBTcl{

	private:

		LBTConfig& config;
		double computeScatteringRate(int flavor, double PLen, double T);
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
