#include "LBTcl.h"
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the key function of LBT
// Update elastic and inelastic evolution of the particle list for one time step
//
///////////////////////////////////////////////////////////////////////////////////////////////////





double LBTcl::computeScatteringRate(int flavor, double PLen, double T) {
//	// Lookup tables qhatG/qhatLQ/qhatHQ depending on flavor
//	// Interpolate based on PLen, T
//	// Return RTE (rate)
//
//	// --- Step 1: Find interpolation points ---
//	double T1, T2, E1, E2;
//	int iT1, iT2, iE1, iE2;
//
//	// Assuming lam() is available (sets T1, T2, E1, E2 and index iT1, iT2, iE1, iE2)
//	lam(flavor, /*out*/ T1, T2, E1, E2, iT1, iT2, iE1, iE2);
//
//	double RTE1, RTE2;
//
//	// --- Step 2: Interpolate in Temperature first ---
//	if (flavor == 21) {
//		// Gluon
//		RTE1 = (qhatG[iT2][iE1] - qhatG[iT1][iE1]) * (T - T1) / (T2 - T1) + qhatG[iT1][iE1];
//		RTE2 = (qhatG[iT2][iE2] - qhatG[iT1][iE2]) * (T - T1) / (T2 - T1) + qhatG[iT1][iE2];
//	} else if (abs(flavor) == 4) {
//		// Heavy quark
//		RTE1 = (qhatHQ[iT2][iE1] - qhatHQ[iT1][iE1]) * (T - T1) / (T2 - T1) + qhatHQ[iT1][iE1];
//		RTE2 = (qhatHQ[iT2][iE2] - qhatHQ[iT1][iE2]) * (T - T1) / (T2 - T1) + qhatHQ[iT1][iE2];
//	} else {
//		// Light quark
//		RTE1 = (qhatLQ[iT2][iE1] - qhatLQ[iT1][iE1]) * (T - T1) / (T2 - T1) + qhatLQ[iT1][iE1];
//		RTE2 = (qhatLQ[iT2][iE2] - qhatLQ[iT1][iE2]) * (T - T1) / (T2 - T1) + qhatLQ[iT1][iE2];
//	}
//
//	// --- Step 3: Interpolate in Energy (PLen) ---
//	double qhat_over_T3 = (RTE2 - RTE1) * (PLen - E1) / (E2 - E1) + RTE1;
//
//	// --- Step 4: Apply K-factors ---
//	double runKT = fixAlphas / 0.3;  // running coupling factor
//	double KPfactor = 1.0 + KPamp * exp(-PLen * PLen / (2.0 * KPsig * KPsig));
//	double KTfactor = 1.0 + KTamp * exp(-pow(T - hydro_Tc, 2) / (2.0 * KTsig * KTsig));
//
//	double Kfactor = KPfactor * KTfactor * KTfactor * runKT * preKT;  // full correction
//
//	qhat_over_T3 *= Kfactor;  // Apply correction
//
//	// --- Step 5: Multiply by T^3 to get real qhat ---
//	return qhat_over_T3 * pow(T, 3);  // Final scattering rate (momentum broadening)
return 0.;
}


double LBTcl::computeRadiationProbability(Particle &p, double T, double E) {
//	//Calculate the probability that a radiation event happens for the current particle over time dt_lrf.
//	//Use radng[i] from the original code, which tracks average number of gluons radiated.
//	//Need to use heavy quark or light parton radiation tables, e.g., nHQgluon().
//	// Return probability
//
//	// --- Inputs: particle p, local temperature T, energy E ---
//
//	// Calculate running alphas
//	double alpha_s = alphas0(Kalphas, T);  // Assuming alphas0() computes coupling
//	double runKT = fixAlphas / 0.3;         // running coupling scaling
//
//	// Calculate K-factors
//	double KTfactor = 1.0 + KTamp * exp(-pow(T - hydro_Tc, 2) / (2.0 * KTsig * KTsig));
//
//	// Local interaction time
//	double dt_lrf = dt * sqrt(1.0 - pow(p.vcfrozen[1], 2) - pow(p.vcfrozen[2], 2) - pow(p.vcfrozen[3], 2));
//	p.Tint_lrf += dt_lrf;
//
//	// --- Update expected number of gluons radiated ---
//	if (p.KATT == 21) {
//		// Gluon: need 1/2 suppression
//		p.radng += nHQgluon(p.KATT, dt_lrf, p.Tint_lrf, T, E, maxFncHQ) * KTfactor * runKT / 2.0;
//	} else {
//		// Quarks
//		p.radng += nHQgluon(p.KATT, dt_lrf, p.Tint_lrf, T, E, maxFncHQ) * KTfactor * runKT;
//	}
//
//	// --- Phase space limitation ---
//	double lim_low = sqrt(6.0 * M_PI * alpha_s) * T / E;
//	double lim_high = (abs(p.KATT) == 4) ? 1.0 : (1.0 - lim_low);  // heavy quarks allow full range
//	double lim_int = lim_high - lim_low;
//
//	// --- Final radiation probability ---
//	if (lim_int > 0.0) {
//		return 1.0 - exp(-p.radng);
//	} else {
//		return 0.0;
//	}
//
return 0.;
}


double LBTcl::computeCollisionProbability(
    const Particle &p,
    double qhat,
    double pLen,
    double T,
    double fraction
) {
//    double runKT = 1.0;
//    if (KrunKT == 1) {
//        runKT = log(1.0 + pLen / (T + 1e-12));  // avoid log(0)
//    }
//
//    double KTfactor = 1.0;
//    if (abs(p.KATT) == 4)      KTfactor = 1.0;
//    else if (abs(p.KATT) == 5) KTfactor = 1.0;
//    else if (p.KATT == 21)     KTfactor = 9.0 / 4.0;
//
//    double probCol = fraction * qhat * runKT * KTfactor * (M_PI / 4.0) / (pLen * T + 1e-12);
//    return probCol;
return 0.;
}


void LBTcl::handleElasticCollision(Particle &p, std::vector<Particle> &particles) {
//    int parentIndex = p.index;
//
//    // Step 1: Prepare momentum and flow for transformation
//    double pc0[4] = {p.P[0], p.P[1], p.P[2], p.P[3]};
//    double vc0[4] = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
//    trans(vc0, pc0);
//
//    // Step 2: Determine channel and output momenta
//    int CT, KATT2, KATT3;
//    double temp0 = p.Tfrozen;
//    double PLen = sqrt(pc0[1]*pc0[1] + pc0[2]*pc0[2] + pc0[3]*pc0[3]);
//
//    flavor(CT, p.KATT, KATT2, KATT3, 0.0, PLen, temp0);
//
//    double pc2[4] = {0.0}; // scattered leading parton
//    double pc3[4] = {0.0}; // thermal recoil
//    double pc4[4] = {0.0}; // unused
//    double qt = 0.0;
//
//    if (CT == 11 || CT == 12) {
//        collHQ22(CT, temp0, qhat0, vc0, pc0, pc2, pc3, pc4, qt);
//    } else {
//        colljet22(CT, temp0, qhat0, vc0, pc0, pc2, pc3, pc4, qt);
//    }
//
//    transback(vc0, pc0);
//    transback(vc0, pc2);
//    transback(vc0, pc3);
//
//    // Step 3: Mark parent as inactive and update its momentum
//    for (int j = 0; j < 4; ++j) {
//        p.P[j] = pc0[j];
//    }
//    p.KATT = KATT2;
//    p.isActive = false;
//
//    // Step 4: Create new leading parton (scattered)
//    Particle scattered;
//    for (int j = 0; j < 4; ++j) {
//        scattered.P[j] = pc2[j];
//        scattered.V[j] = p.V[j];
//    }
//    scattered.KATT = KATT2;
//    scattered.CAT = 0;
//    scattered.Tfrozen = temp0;
//    scattered.vcfrozen[1] = p.vcfrozen[1];
//    scattered.vcfrozen[2] = p.vcfrozen[2];
//    scattered.vcfrozen[3] = p.vcfrozen[3];
//    scattered.WT = p.WT;
//    scattered.mass = 0.0;
//    scattered.isPrimary = false;
//    scattered.isActive = true;
//    scattered.mom1 = parentIndex;
//    scattered.mom2 = parentIndex;  // use same if ghost not tracked
//
//    scattered.index = particles.size();  // assign current index
//    particles.push_back(scattered);
//
//    // Step 5: Create thermal recoil parton
//    Particle recoil;
//    for (int j = 0; j < 4; ++j) {
//        recoil.P[j] = pc3[j];
//        recoil.V[j] = p.V[j];
//    }
//    recoil.KATT = KATT3;
//    recoil.CAT = 2;  // recoil
//    recoil.Tfrozen = temp0;
//    recoil.vcfrozen[1] = p.vcfrozen[1];
//    recoil.vcfrozen[2] = p.vcfrozen[2];
//    recoil.vcfrozen[3] = p.vcfrozen[3];
//    recoil.WT = p.WT;
//    recoil.mass = 0.0;
//    recoil.isPrimary = false;
//    recoil.isActive = true;
//    recoil.mom1 = parentIndex;
//    recoil.mom2 = parentIndex;
//
//    recoil.index = particles.size();  // assign current index
//    particles.push_back(recoil);
}




void LBTcl::handleRadiation(Particle &p, std::vector<Particle> &particles, int &icl23, int &iclrad) {
//    int parentIndex = p.index;
//
//    // Step 1: Prepare momentum and flow
//    double pc0[4] = {p.P[0], p.P[1], p.P[2], p.P[3]};
//    double vc0[4] = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
//    trans(vc0, pc0);
//    double Ejp = pc0[0];
//    transback(vc0, pc0);
//
//    if (Ejp <= 2.0 * sqrt(qhat0)) return;
//
//    double pc2[4] = {0.0};  // Updated radiator
//    double pc3[4] = {0.0};  // Ghost (unused)
//    double pc4[4] = {0.0};  // Radiated gluon
//    double qt = 0.0;
//
//    double T = p.Tfrozen;
//    double Tdiff = p.Tint_lrf;
//    double alpha_s = alphas0(Kalphas, T);
//    double lim_low = sqrt(6.0 * M_PI * alpha_s) * T / Ejp;
//    double lim_high = (abs(p.KATT) == 4) ? 1.0 : (1.0 - lim_low);
//    double lim_int = lim_high - lim_low;
//
//    // Step 2: Call radiation kernel
//    collHQ23(p.KATT, T, qhat0, vc0, pc0, pc2, pc3, pc4, qt,
//		    icl23, Tdiff, Ejp, maxFncHQ, lim_low, lim_int, iclrad);
//    //colljet23(T, qhat0, vc0, pc0, pc2, pc3, pc4, qt, icl23, Tdiff, Ejp, iclrad);
//
//    // Step 3: Process radiation if successful
//    if (icl23 != 1 && iclrad != 1) {
//        // Deactivate parent and update momentum
//        for (int j = 0; j < 4; ++j) {
//            p.P[j] = pc0[j];
//        }
//        p.isActive = false;
//
//        // Add radiated gluon
//        Particle gluon;
//        for (int j = 0; j < 4; ++j) {
//            gluon.P[j] = pc4[j];
//            gluon.V[j] = p.V[j];
//        }
//        gluon.KATT = 21;
//        gluon.CAT = 4;
//        gluon.Tfrozen = T;
//        gluon.vcfrozen[1] = p.vcfrozen[1];
//        gluon.vcfrozen[2] = p.vcfrozen[2];
//        gluon.vcfrozen[3] = p.vcfrozen[3];
//        gluon.WT = p.WT;
//        gluon.mass = 0.0;
//        gluon.isPrimary = false;
//        gluon.isActive = true;
//        gluon.mom1 = parentIndex;
//        gluon.mom2 = parentIndex;
//        gluon.index = particles.size();
//        particles.push_back(gluon);
//
//        // Step 4: Handle multiple gluons (Poisson) for HQ
//        int nrad = KPoisson(p.radng);
//        while (--nrad > 0) {
//            double pc2_more[4] = {0.0};
//            double pc4_more[4] = {0.0};
//            double pb[4] = {0.0};
//
//            if (abs(p.KATT) == 4) {
//                radiationHQ(p.KATT, qhat0, vc0, pc4, pc2_more, pc4_more, pb,
//                            iclrad, Tdiff, Ejp, maxFncHQ, T, lim_low, lim_int);
//            } else {
//                radiation(qhat0, vc0, pc4, pc2_more, pc4_more, pb,
//                          iclrad, Tdiff, Ejp);
//            }
//
//            if (iclrad != 1) {
//                Particle extraGluon;
//                for (int j = 0; j < 4; ++j) {
//                    extraGluon.P[j] = pc4_more[j];
//                    extraGluon.V[j] = p.V[j];
//                }
//                extraGluon.KATT = 21;
//                extraGluon.CAT = 4;
//                extraGluon.Tfrozen = T;
//                extraGluon.vcfrozen[1] = p.vcfrozen[1];
//                extraGluon.vcfrozen[2] = p.vcfrozen[2];
//                extraGluon.vcfrozen[3] = p.vcfrozen[3];
//                extraGluon.WT = p.WT;
//                extraGluon.mass = 0.0;
//                extraGluon.isPrimary = false;
//                extraGluon.isActive = true;
//                extraGluon.mom1 = parentIndex;
//                extraGluon.mom2 = parentIndex;
//                extraGluon.index = particles.size();
//                particles.push_back(extraGluon);
//            } else {
//                break;  // Stop emitting if emission fails
//            }
//        }
//    }
}



void LBTcl::propagateParticle(Particle &p, double ti, int &free, double &fraction) {

	if (config.clock.tauswitch == 0) {
		// --- Propagation in t-z coordinates ---

		// Update spatial position
		for (int j = 1; j <= 3; ++j) {
			p.V[j] += config.clock.dt * p.P[j] / p.P[0];
		}

		// Cartesian coordinates
		double tcar = ti;
		double xcar = p.V[1];
		double ycar = p.V[2];
		double zcar = p.V[3];

		// Query hydro fields
		double VX, VY, VZ, ed, sd;
		double &temp = (p.CAT!=5)? config.medium.temp0: config.medium.temp00; 
		int hydro_ctl = 0;

		if (config.medium.bulkFlag == 1) {
			readhydroinfoshanshan_(&tcar,&xcar,&ycar,&zcar,&ed,&sd,&temp,&VX,&VY,&VZ,&hydro_ctl);
		} else if (config.medium.bulkFlag == 2) {
			hydroinfoccnu_(&tcar, &xcar, &ycar, &zcar, &temp, &VX, &VY, &VZ, &hydro_ctl);
		} else {
			VX = VY = VZ = 0.0;
			hydro_ctl = 0;
		}

		// If medium is alive
		if (hydro_ctl == 0 && temp >= config.medium.hydro_Tc) {
			p.Vfrozen[0] = ti;
			for (int j = 1; j <= 3; ++j) {
				p.Vfrozen[j] = p.V[j];
			}
			p.Tfrozen = temp;
			p.vcfrozen[1] = VX;
			p.vcfrozen[2] = VY;
			p.vcfrozen[3] = VZ;
			fraction = 1.0;
			free = 0;
		} else {
			// Free-streaming
			fraction = 0.0;
			free = 1;
		}

	} else {
		// --- Propagation in tau-eta coordinates ---

		std::array<double, 4> vp = {0.0, p.V[1], p.V[2], p.V[3]};
		std::array<double, 4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
		std::array<double, 4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};

		double Vx, Vy, Veta, Xtau;
		this->titau(ti, vc0, vp, pc0, Vx, Vy, Veta, Xtau);

		// Update spatial position
		p.V[1] += config.clock.dt * Vx;
		p.V[2] += config.clock.dt * Vy;
		p.V[3] += config.clock.dt * Veta;

		double tcar = ti * std::cosh(p.V[3]);
		double zcar = ti * std::sinh(p.V[3]);
		double xcar = p.V[1];
		double ycar = p.V[2];

		// Query hydro fields
		double VX, VY, VZ, ed, sd;
		double &temp = (p.CAT!=5)? config.medium.temp0: config.medium.temp00; 
		int hydro_ctl = 0;

		if (config.medium.bulkFlag == 1) {
			readhydroinfoshanshan_(&tcar,&xcar,&ycar,&zcar,&ed,&sd,&temp,&VX,&VY,&VZ,&hydro_ctl);
		} else if (config.medium.bulkFlag == 2) {
			hydroinfoccnu_(&tcar, &xcar, &ycar, &zcar, &temp, &VX, &VY, &VZ, &hydro_ctl);
		} else {
			VX = VY = VZ = 0.0;
			hydro_ctl = 0;
		}

		if (hydro_ctl == 0 && temp >= config.medium.hydro_Tc) {
			p.Vfrozen[0] = ti;
			for (int j = 1; j <= 3; ++j) {
				p.Vfrozen[j] = p.V[j];
			}
			p.Tfrozen = temp;
			p.vcfrozen[1] = VX;
			p.vcfrozen[2] = VY;
			p.vcfrozen[3] = VZ;
			fraction = 1.0;
			free = 0;
		} else {
			fraction = 0.0;
			free = 1;
		}
	}//tauswitch
return;
}






void LBTcl::LBT(std::vector<Particle> &part_event, double ti) {
    int np_snapshot = part_event.size();  // snapshot of part_event at start of this step "np0"

    // Assign indices to existing part_event for ancestry tracking
    for (int i = 0; i < (int)part_event.size(); ++i) {
        part_event[i].index = i;
    }

    // Loop over all active part_event at this step
    for (int i = 0; i < (int)part_event.size(); ++i) {
        Particle &p = part_event[i];

std::cout << __FILE__ << "(" << __LINE__ << ")" << "i, Vfrozen " << i << " " << p.Vfrozen[0] << std::endl;
        // Skip frozen or inactive part_event
        if (!p.isActive || p.Vfrozen[0] >= ti) continue;

        int free = 0;
        double fraction = 0.0;

std::cout << __FILE__ << "(" << __LINE__ << ")" << "Before propagation " << std::endl;
p.Print();

        // Propagate parton
	this->propagateParticle(p, ti, free, fraction);

std::cout << __FILE__ << "(" << __LINE__ << ")" << "After propagation " << std::endl;
p.Print();

	if (p.CAT != 1) {

            if (free == 0) {
                // Boost momentum into local fluid frame
                std::array<double,4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
                std::array<double,4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
                this->trans(vc0, pc0);
                double E = pc0[0];
                double PLen = sqrt(pc0[1]*pc0[1] + pc0[2]*pc0[2] + pc0[3]*pc0[3]);
                this->transback(vc0, pc0);
//
//                // Query medium
//                double T = p.Tfrozen;
//                double qhat = computeScatteringRate(p.KATT, PLen, T);
//
//                // Compute probabilities
//                double probCol = computeCollisionProbability(p, qhat, PLen, T, fraction);
//                double probRad = computeRadiationProbability(p, T, E);
//                double probTot = probCol + probRad;
//
//                // Sample scattering
//                if (ran0(&NUM1) < probTot) {
//                    if (ran0(&NUM1) < probRad / probTot) {
//                        handleRadiation(p, part_event, icl23, iclrad);
//                    } else {
//                        handleElasticCollision(p, part_event);
//                    }
//                }
//
//                // Reset radiation tracker
//                p.radng = 0.0;
            }//free ==0
        }//for positive partons
    }//particle loop
return;
}


