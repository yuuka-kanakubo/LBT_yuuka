#include "LBTcl.h"
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the key function of LBT
// Update elastic and inelastic evolution of the particle list for one time step
//
///////////////////////////////////////////////////////////////////////////////////////////////////




double LBTcl::computeScatteringRate(Particle &p, const double PLen_in, const double T_in) {
	// Lookup tables qhatG/qhatLQ/qhatHQ depending on flavor
	// Interpolate based on PLen(E), T
	// Return RTE (rate)

	// --- Step 1: Find interpolation points ---
	double T1, T2, E1, E2;
	int iT1, iT2, iE1, iE2;

	// Assuming lam() is available (sets T1, T2, E1, E2 and index iT1, iT2, iE1, iE2)
        int flavor = p.pid;
	lam(flavor, p.tot_el_rate, PLen_in, T_in, T1, T2, E1, E2, iT1, iT2, iE1, iE2);
	//Assuming lam returned T1, T2, E1, E2, iT1, iT2, iE1, iE2.
	//They will be used in lam2 to get RTE1 and RTE2.

	// --- Step 2: Interpolate in temperature (T_in) ---
	double RTE1, RTE2;
	lam2(flavor, RTE1, RTE2, T_in, T1, T2, iT1, iT2, iE1, iE2);

	// --- Step 3: Interpolate in Energy (PLen) ---
	double qhat_over_T3 = (RTE2 - RTE1) * (PLen_in - E1) / (E2 - E1) + RTE1;

	// --- Step 4: Apply K-factors ---
	double KPfactor = 1.0 + config.lbtinput.KPamp * exp(-PLen_in * PLen_in / (2.0 * config.lbtinput.KPsig * config.lbtinput.KPsig));
	double KTfactor = 1.0 + config.lbtinput.KTamp * exp(-pow(T_in - config.medium.hydro_Tc, 2) / (2.0 * config.lbtinput.KTsig * config.lbtinput.KTsig));

	double Kfactor = KPfactor * KTfactor * KTfactor * config.lbtinput.runKT * config.lbtinput.preKT;  // full correction

	qhat_over_T3 *= Kfactor;  // Apply correction

	p.get_D2piT(qhat_over_T3);

	// --- Step 5: Multiply by T^3 to get real qhat ---
	return qhat_over_T3 * pow(T_in, 3);  // Final scattering rate (momentum broadening)
}


double LBTcl::computeRadiationProbability(Particle &p, const double T, const double E) {
	//Calculate the probability that a radiation event happens for the current particle over time dt_lrf.
	//Use radng[i] from the original code, which tracks average number of gluons radiated.
	//Need to use heavy quark or light parton radiation tables, e.g., nHQgluon().
	// Return probability

	// --- Inputs: particle p, local temperature T, energy E ---

	// Calculate running alphas
	double alpha_s = alphas0(config.physics.Kalphas, T);  // Assuming alphas0() computes coupling

	// Calculate K-factors
	double KTfactor = 1.0 + config.lbtinput.KTamp * exp(-pow(T - config.medium.hydro_Tc, 2) / (2.0 * config.lbtinput.KTsig * config.lbtinput.KTsig));

	// Local interaction time
	//double dt_lrf = config.clock.dt * sqrt(1.0 - pow(p.vcfrozen[1], 2) - pow(p.vcfrozen[2], 2) - pow(p.vcfrozen[3], 2));
	double dt_lrf = config.clock.dt*p.timedilation;
	p.Tint_lrf += dt_lrf;

	// --- Update expected number of gluons radiated ---
	if (p.pid == 21) {
		// Gluon: need 1/2 suppression
		p.radng += nHQgluon(p, dt_lrf, T, E) * KTfactor * config.lbtinput.runKT / 2.0;
	} else {
		// Quarks
		p.radng += nHQgluon(p, dt_lrf, T, E) * KTfactor * config.lbtinput.runKT;
	}

	// --- Phase space limitation ---
	double lim_low = sqrt(6.0 * M_PI * alpha_s) * T / E;
	double lim_high = (abs(p.pid) == 4) ? 1.0 : (1.0 - lim_low);  // heavy quarks allow full range
	double lim_int = lim_high - lim_low;

	// --- Final radiation probability ---
	if (lim_int > 0.0) {
		return 1.0 - exp(-p.radng);
	} else {
		return 0.0;
	}

}


double LBTcl::computeCollisionProbability(
    Particle &p,
    const double PLen_in,
    const double T_in,
    const double fraction,
    const double probRad
) {

	double dt_lrf = config.clock.dt*p.timedilation;
	double probCol = (config.clock.tauswitch==0)? fraction*dt_lrf*p.tot_el_rate/base::GEVFM: dt_lrf*p.tot_el_rate*p.Xtau_keep/base::GEVFM;

        //Modify vertices to take into account time-ordering collisions
	if(config.clock.tauswitch==0) p.V[0]=p.V[0]-fraction*config.clock.dt*p.tot_el_rate/base::GEVFM*sqrt( pow(p.P[1],2)+pow(p.P[2],2)+pow(p.P[3],2) )/p.P[0];
	else p.V[0]=p.V[0]-config.clock.dt*p.tot_el_rate*p.Xtau_keep/base::GEVFM*sqrt(pow(p.P[1],2)+pow(p.P[2],2)+pow(p.P[3],2))/p.P[0];
	//TODO: Is this config.clock.dt -> config.clock.dt * p.timedilation?

	double KPfactor = 1.0 + config.lbtinput.KPamp * exp(-PLen_in * PLen_in / (2.0 * config.lbtinput.KPsig * config.lbtinput.KPsig));
	double KTfactor = 1.0 + config.lbtinput.KTamp * exp(-pow(T_in - config.medium.hydro_Tc, 2) / (2.0 * config.lbtinput.KTsig * config.lbtinput.KTsig));

	probCol*=KPfactor*KTfactor*config.lbtinput.runKT;
	probCol=(1.0-exp(-probCol))*(1.0-probRad);

    return probCol;
}


void LBTcl::handleElasticCollision(Particle &p, const double PLenloc, std::vector<Particle> &particles) {
    int parentIndex = p.index;

    // Step 1: Prepare momentum and flow for transformation
    std::array <double, 4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
    std::array <double, 4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};

    // Step 2: Determine channel and output momenta
    // Decide flavor outcome
    int channel, pid_med, pid_rec;
    flavor(p.pid, p.tot_el_rate, PLenloc, p.Tfrozen, channel, pid_med, pid_rec);
std::cout << "pid med " << pid_med <<  std::endl;
std::cout << "pid rec " << pid_rec <<  std::endl;
std::cout << "channel " << channel <<  std::endl;
exit(1);

    // Deactivate incoming parton
    p.isActive = false;


    // Create new outgoing particles
    Particle p_med;//incoming medium parton
    Particle p_rec;

    p_med.pid = pid_med;
    p_rec.pid = pid_rec;

    p_med.CAT = 2;  // recoiled
    p_rec.CAT = 2;

    p_med.mom1 = p.index;
    p_rec.mom1 = p.index;

    p_med.mom2 = p.index;
    p_rec.mom2 = p.index;


    // Optionally: log the event
    std::cout << "Elastic collision (channel=" << channel << ") at t=" << p.V[0]
              << ": pid=" << p.pid << " → pid_med=" << pid_med << ", pid_rec=" << pid_rec << std::endl;


    if (channel == 11 || channel == 12) {
	    collHQ22(channel, p, p_rec, p_med, qt);
    } else {
	    colljet22(CT, p, p_rec, p_med, qt);
    }
//
//    transback(vc0, pc0);
//    transback(vc0, pc2);
//    transback(vc0, pc3);
//
//    // Step 3: Mark parent as inactive and update its momentum
//    for (int j = 0; j < 4; ++j) {
//        p.P[j] = pc0[j];
//    }
//    p.pid = pid2;
//    p.isActive = false;
//
//    // Step 4: Create new leading parton (scattered)
//    Particle scattered;
//    for (int j = 0; j < 4; ++j) {
//        scattered.P[j] = pc2[j];
//        scattered.V[j] = p.V[j];
//    }
//    scattered.pid = pid2;
//    scattered.CAT = 0;
//    scattered.Tfrozen = p.Tfrozen;
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
//    recoil.pid = pid3;
//    recoil.CAT = 2;  // recoil
//    recoil.Tfrozen = p.Tfrozen;
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



    out3.index = particles.size();
    particles.push_back(out3);

    out4.index = particles.size();
    particles.push_back(out4);

}




void LBTcl::handleRadiation(Particle &p, std::vector<Particle> &particles) {
    int parentIndex = p.index;

    // Step 1: Prepare momentum and flow
    std::array<double, 4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
    std::array<double, 4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
    trans(vc0, pc0);
    double Eloc = pc0[0];
    transback(vc0, pc0);

    double alpha_s = alphas0(config.physics.Kalphas, p.Tfrozen);  // Assuming alphas0() computes coupling
    double qhat0 = DebyeMass2(config.physics.Kalphas, alpha_s, p.Tfrozen);  // qhat_0: Calculated by  \mu_D^2 = 4\pi \alpha_s T^2

    if (Eloc <= 2.0 * sqrt(qhat0)) return;

    double lim_low = sqrt(6.0 * M_PI * alpha_s) * p.Tfrozen / Eloc;
    double lim_high = (abs(p.pid) == 4) ? 1.0 : (1.0 - lim_low);
    double lim_int = lim_high - lim_low;

    // Step 2: Call radiation kernel
    double qt = 0.0;
    double pc2[4] = {0.0};  // Updated radiator
    double pc3[4] = {0.0};  // Ghost (unused)
    double pc4[4] = {0.0};  // Radiated gluon

//TODO
//    collHQ23(p.pid, p.Tfrozen, qhat0, vc0, pc0, pc2, pc3, pc4, qt,
//		    icl23, p.Tint_lrf, Eloc, p.max_Ng, lim_low, lim_int, iclrad);
//    //colljet23(T, qhat0, vc0, pc0, pc2, pc3, pc4, locqt, icl23, p.Tint_lrf, Ejp, iclrad);

    // Step 3: Process radiation if successful
 //   if (icl23 != 1 && iclrad != 1) {
        // Deactivate parent and update momentum
        for (int j = 0; j < 4; ++j) {
            p.P[j] = pc0[j];
        }
        p.isActive = false;

        // Add radiated gluon
        Particle gluon;
        for (int j = 0; j < 4; ++j) {
            gluon.P[j] = pc4[j];
            gluon.V[j] = p.V[j];
        }
        gluon.pid = 21;
        gluon.CAT = 4;
        gluon.Tfrozen = p.Tfrozen;
        gluon.vcfrozen[1] = p.vcfrozen[1];
        gluon.vcfrozen[2] = p.vcfrozen[2];
        gluon.vcfrozen[3] = p.vcfrozen[3];
        gluon.WT = p.WT;
        gluon.mass = 0.0;
        gluon.isPrimary = false;
        gluon.isActive = true;
        gluon.mom1 = parentIndex;
        gluon.mom2 = parentIndex;
        gluon.index = particles.size();
        particles.push_back(gluon);
//
//        // Step 4: Handle multiple gluons (Poisson) for HQ
//        int nrad = KPoisson(p.radng);
//        while (--nrad > 0) {
//            double pc2_more[4] = {0.0};
//            double pc4_more[4] = {0.0};
//            double pb[4] = {0.0};
//
//            if (abs(p.pid) == 4) {
//                radiationHQ(p.pid, qhat0, vc0, pc4, pc2_more, pc4_more, pb,
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
//                extraGluon.pid = 21;
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
  //  }
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

		//Caclutate some variables used to caluculate interaction rate, radiation rate etc.
		p.get_timedilation();


	} else {
		// --- Propagation in tau-eta coordinates ---

		std::array<double, 4> vp = {0.0, p.V[1], p.V[2], p.V[3]};
		std::array<double, 4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
		std::array<double, 4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};

		double Vx, Vy, Veta, Xtau;
		this->titau(ti, vc0, vp, pc0, Vx, Vy, Veta, Xtau);
		p.Xtau_keep = Xtau;//archived for later use.

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
std::cout << "=============== " << i << std::endl;
std::cout << ":):):):) Particle ....... " << i << "   at time " << ti << std::endl;

        // Skip frozen or inactive part_event
        if (!p.isActive || p.Vfrozen[0] >= ti) continue;

        int free = 0;
        double fraction = 0.0;


        // Propagate parton
	this->propagateParticle(p, ti, free, fraction);


	if (p.CAT != 1 && free == 0) {
                // Boost momentum into local fluid frame
                std::array<double,4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
                std::array<double,4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
                this->trans(vc0, pc0);
                double Eloc = pc0[0];
                double PLenloc = sqrt(pc0[1]*pc0[1] + pc0[2]*pc0[2] + pc0[3]*pc0[3]);
                this->transback(vc0, pc0);


                // Query medium
                double T = p.Tfrozen;
                double qhat = computeScatteringRate(p, PLenloc, T);
		std::cout << "qhat " << qhat << std::endl;

                // Compute probabilities
                double probRad = computeRadiationProbability(p, T, Eloc);
                double probCol = computeCollisionProbability(p, PLenloc, T, fraction, probRad);
		std::cout << "probRad " << probRad << std::endl;
		std::cout << "probCol " << probCol << std::endl;
                double probTot = probCol + probRad;


                // Sample scattering
		//TODO: should be like, ran0 < probCol, ran0 < probRad in parallel?
		//      For now, I am just trying to reproduce the results.
		if (ran0(&config.rng.NUM1) < probTot) {
			std::cout << __FILE__ << "(" << __LINE__ << ")" << "Calling handleElasticCollision. " << std::endl;
			std::cout << "ti " << ti << "  " << i << std::endl;
			handleElasticCollision(p, PLenloc, part_event);
			if (ran0(&config.rng.NUM1) < probRad / probTot) {
				std::cout << __FILE__ << "(" << __LINE__ << ")" << "Calling handleRadiation. " << std::endl;
				std::cout << "ti " << ti << "  " << i << std::endl;
				exit(1);
				handleRadiation(p, part_event);
				//handleRadiation(p, part_event, icl23, iclrad);
			}
		}
		// Reset radiation tracker
                p.radng = 0.0;
            }//positive and free ==0(in medium)
    }//particle loop
    return;
}


