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
	double lim_low = sqrt(6.0 * base::pi * alpha_s) * T / E;
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


double LBTcl::handleElasticCollision(Particle &p, const double PLenloc, std::vector<Particle> &part_event) {
    int parentIndex = p.index();

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

    // Deactivate incoming parton
    p.isActive = false;


    // Create new outgoing particles
    Particle p_med = p;//incoming medium parton
    Particle p_rec = p;

    p_med.pid = pid_med;
    p_rec.pid = pid_rec;

    Particle p_fin = p;

    //index
    p_med.assign_index();
    p_rec.assign_index();
    p_fin.assign_index();


    // Optionally: log the event
    std::cout << "Elastic collision (channel=" << channel << ") at t=" << p.V[0]
              << ": pid=" << p.pid << " → pid_med=" << pid_med << ", pid_rec=" << pid_rec << std::endl;

    double qt;
    if (channel == 11 || channel == 12) {
	    collHQ22(channel, p, p_rec, p_med, p_fin, qt);
            std::cout << "ERROR: refactoring is not completed yet. " << std::endl;
	    exit(EXIT_FAILURE);
    } else {
	    colljet22(channel, p, p_rec, p_med, p_fin, qt);
    }

    //(De)Activate particles after collision
    p.isActive = false;

    p_fin.isPrimary = true;//lost energy, still primary.
    p_med.isPrimary = false;
    p_rec.isPrimary = false;
    p_fin.isActive = true;
    p_med.isActive = true;
    p_rec.isActive = true;

    std::cout << "p_rec " << std::endl;
    p_rec.Print();
    std::cout << "p_med " <<  std::endl;
    p_med.Print();
    std::cout << "p_fin " <<  std::endl;
    p_fin.Print();

    p_med.CAT = 3;//negative
    p_rec.CAT = 2;//recoiled 

    //p_med.parent1 = -1;//I am just commenting out since by default it will be initialised with -1
    //p_med.parent2 = -1;//N/A since it's medium parton

    p_rec.parent1 = p.index();
    p_rec.parent2 = p_med.index();


    p_med.kid1 = p_rec.index();
    p_med.kid2 = p_fin.index();
    p.kid1 = p_rec.index();
    p.kid2 = p_fin.index();

   //Vfrozen, Tfrozen = Vertices and Temperature of interaction point?
    p_med.copy_thisV_to_Vfrozen(p.V);
    p_med.Tfrozen = p.Tfrozen;
    p_med.copy_thisV_to_vcfrozen(p.vcfrozen);

    p_rec.copy_thisV_to_Vfrozen(p.V);
    p_rec.Tfrozen = p.Tfrozen;
    p_rec.copy_thisV_to_vcfrozen(p.vcfrozen);

    part_event.push_back(p_fin);
    part_event.push_back(p_med);
    part_event.push_back(p_rec);

    return qt;
}




void LBTcl::handleRadiation(Particle &p, std::vector<Particle> &part_event) {

	//In original LBT...
	//pc01 is pc4 (initial jet momentum in 2->2)
	//pc2: "recoiled" particle. will be overwritten by pc3 anyways.
	//pc3: scattering "medium" particle. ==> apparently reused from the handleCollision
	//pc4: radiated gluon momentum will be calculated in this function.
	std::array<double, 4> pc_rad = {0., 0., 0., 0.};
	std::array<double, 4> pc_rec = {0., 0., 0., 0.};
	std::array<double, 4> pc_fin = {0., 0., 0., 0.};

	//Looking for the missing partner!
	//Then, what I need to do here is to find the "scattering medium particle" in part_event.
	//I am tracing pc4 (the incoming jet). So, I know daughters of pc4.
	//Once I know the daughters of pc4, I can find their moms indices. 
	//The one has medium tag (CAT = "medium") is going to be pc3. 

	//Here I hope kid1 is allocated in part_event[i_kid1]
	//Is this scattered medium parton?
	std::array<double, 4> pc_med = {0., 0., 0., 0.};
	if(part_event[p.kid1].CAT == 3){
		for(int i = 1; i<=3 ; i++) pc_med[i] =  part_event[p.kid1].P[i];
	}else if(part_event[p.kid2].CAT == 3){
		for(int i = 1; i<=3 ; i++) pc_med[i] =  part_event[p.kid1].P[i];
	}else{
		std::cout << "ERROR! Couldn't find the missing medium partner... " << std::endl; 
		std::cout << "p.kid1 " << p.kid1 << std::endl; 
		std::cout << "p.kid2 " << p.kid2 << std::endl; 
		std::cout << "part_event[p.kid1].CAT " << part_event[p.kid1].CAT << std::endl; 
		std::cout << "part_event[p.kid2].CAT " << part_event[p.kid2].CAT << std::endl; 
		exit(EXIT_FAILURE);
	}


// Step 1: Prepare momentum and flow
    std::array<double, 4> pc_jet = {p.P[0], p.P[1], p.P[2], p.P[3]};
    std::array<double, 4> v_flow = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
			std::cout << "pc_jet(p4) " <<  pc_jet[0] << "  " << pc_jet[1] << "  " << pc_jet[2] << "  " << pc_jet[3] << std::endl;
			std::cout << "pc_med(p2) " <<  pc_med[0] << "  " << pc_med[1] << "  " << pc_med[2] << "  " << pc_med[3] << std::endl;
			//std::cout << "pc_rec(p3) " <<  pc_rec[0] << "  " << pc_rec[1] << "  " << pc_rec[2] << "  " << pc_rec[3] << std::endl;
			//std::cout << "pc_fin(p0) " <<  pc_fin[0] << "  " << pc_fin[1] << "  " << pc_fin[2] << "  " << pc_fin[3] << std::endl;
exit(1);
    trans(v_flow, pc_jet);
    double Eloc = pc_jet[0];
    transback(v_flow, pc_jet);

    double alpha_s = alphas0(config.physics.Kalphas, p.Tfrozen);  // Assuming alphas0() computes coupling
    double qhat0 = DebyeMass2(config.physics.Kqhat0, alpha_s, p.Tfrozen);  // qhat_0: Calculated by  \mu_D^2 = 4\pi \alpha_s T^2

    if (Eloc <= 2.0 * sqrt(qhat0)) return;

    // Step 2: Call radiation kernel


//TODO
//    collHQ23(p.pid, p.Tfrozen, qhat0, vc0, pc0, pc2, pc3, pc4, qt,
//		    icl23, p.Tint_lrf, Eloc, p.max_Ng, lim_low, lim_int, iclrad);
//    //colljet23(T, qhat0, vc0, pc0, pc2, pc3, pc4, locqt, icl23, p.Tint_lrf, Ejp, iclrad);

    // Step 3: Process radiation if successful
 //   if (icl23 != 1 && iclrad != 1) {
        // Deactivate parent and update momentum
        for (int j = 0; j < 4; ++j) {
            p.P[j] = pc_jet[j];
        }
        p.isActive = false;

        // Add radiated gluon
        Particle gluon;
        for (int j = 0; j < 4; ++j) {
            gluon.P[j] = pc_jet[j];
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
        gluon.parent1 = p.index();
        //gluon.parent2 = p.index();
        gluon.assign_index();
        part_event.push_back(gluon);
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

/*
		if (ran0(&config.rng.NUM1) < probTot) {

			if (ran0(&config.rng.NUM1) < probRad / probTot) {
                                          do Radiation
			  }else{
					  do Collision
			  } 
	        }


		if (ran0(&config.rng.NUM1) < probTot) {

					  do Collision

			if (ran0(&config.rng.NUM1) < probRad / probTot) {
                                          do Radiation
			  } 
	        }


*/


                // Sample scattering
		//TODO: should be like, ran0 < probCol, ran0 < probRad in parallel?
		//      For now, I am just trying to reproduce the results.
		if (ran0(&config.rng.NUM1) < probTot) {
			std::cout << __FILE__ << "(" << __LINE__ << ")" << "Calling handleElasticCollision. " << std::endl;
			double qt = handleElasticCollision(p, PLenloc, part_event);

		//CHECKING
		std::cout << __FILE__ << "(" << __LINE__ << ")" << "After handleElastic" << std::endl;
		for (auto it = part_event.begin(); it != part_event.end(); ++it) {
			it->Print();
		}
exit(1);


                        if(config.physics.Kinteraction==0 || (qt<base::epsilon)) continue;
			if (ran0(&config.rng.NUM1) < probRad / probTot) {
				std::cout << __FILE__ << "(" << __LINE__ << ")" << "Calling handleRadiation. " << std::endl;
				std::cout << "ti " << ti << "  " << i << std::endl;
				handleRadiation(p, part_event);
				exit(1);
			}
		}
		// Reset radiation tracker
                p.radng = 0.0;
            }//positive and free ==0(in medium)
    }//particle loop
    return;
}


