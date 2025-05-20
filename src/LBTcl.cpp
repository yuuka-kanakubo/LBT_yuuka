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
        p.qhat_over_T3 = qhat_over_T3;

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


double LBTcl::handleElasticCollision(Particle &p, const double PLenloc, std::vector<Particle> &part_event, std::vector<Particle> &part_current) {
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
    Particle p_fin = p;

    //index
    p_med.assign_index();
    p_rec.assign_index();
    p_fin.assign_index();

    // Optionally: log the event
    std::cout << "Elastic collision (channel=" << channel << ") at t=" << p.V[0]
              << ": pid=" << p.pid << " → pid_med=" << pid_med << ", pid_rec=" << pid_rec << std::endl;

    double qt;
    std::array<double, 4> pc_rec = {0.,0.,0.,0.};// output: final medium parton momentum
    std::array<double, 4> pc_med = {0.,0.,0.,0.};// output: initial medium parton
    std::array<double, 4> pc_fin = {0.,0.,0.,0.};
    if (channel == 11 || channel == 12) {
	    collHQ22(channel, p, p_rec, p_med, p_fin, qt);
            std::cout << "ERROR: refactoring is not completed yet. " << std::endl;
	    exit(EXIT_FAILURE);
    } else {
	    colljet22(channel, p, pc_rec, pc_med, pc_fin, qt);
    }

    //I want to select parton with larger energy
    //Then, put it to the original location of incoming jet.
    //IN otiginal LBT, they are exchanging pc0 and pc2 and flavor0 flavor2.
    if(pc_fin[0]>pc_rec[0]){
	    for(int i=0; i<=3; i++){
		    p_fin.P[i] = pc_fin[i];
		    p_rec.P[i] = pc_rec[i];
		    p_med.P[i] = pc_med[i];
	    }
	    p_fin.pid = p.pid;
	    p_rec.pid = pid_rec;
	    p_med.pid = pid_med;

    }else{
	    for(int i=0; i<=3; i++){
		    p_fin.P[i] = pc_rec[i];
		    p_rec.P[i] = pc_fin[i];
		    p_med.P[i] = pc_med[i];
	    }
	    p_fin.pid = pid_rec;
	    p_rec.pid = p.pid;
	    p_med.pid = pid_med;

    }


    //(De)Activate particles after collision
    p.isActive = false;
    p_fin.isActive = true;
    p_med.isActive = true;
    p_rec.isActive = true;

    p_fin.isPrimary = true;//lost energy, still primary.
    p_med.isPrimary = false;
    p_rec.isPrimary = false;

    p_med.CAT = 3;//negative
    p_rec.CAT = 2;//recoiled 

    //p_med.parent1 = -1;//I am just commenting out since by default it will be initialised with -1
    //p_med.parent2 = -1;//N/A since it's medium parton

    p_rec.parent1 = p.index();
    p_rec.parent2 = p_med.index();
    p_fin.parent1 = p.index();
    p_fin.parent2 = p_med.index();

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


    part_current.push_back(p_fin);
    part_current.push_back(p_rec);
    part_current.push_back(p_med);

    return qt;
}




void LBTcl::handleRadiation(Particle &p, const double qt, std::vector<Particle> &part_event, std::vector<Particle> &part_current) {

	// Step 1: Prepare momentum and flow

	//In original LBT...
	//pc01 is pc4 (initial jet momentum in 2->2), i.e., &p (even before the momentum exchange in 2->2)
	//pc2: "recoiled" particle. will be overwritten by pc3 anyways.
	//pc3: scattering "medium" particle. ==> apparently reused from the handleCollision
	//pc4: radiated gluon momentum will be calculated in this function.
	std::array<double, 4> pc_rad = {0., 0., 0., 0.};
	std::array<double, 4> pc_rec = {0., 0., 0., 0.};
	std::array<double, 4> pc_fin = {0., 0., 0., 0.};

	//Looking for the missing partner!
	//What I need to do here is to find the "scattering medium particle" in part_event.
	//I am tracing all the family tree of p So, I know kids of pc4.
	//The one has medium tag (CAT = "medium") is going to be pc3. 

	//Is this scattered medium parton?
	std::array<double, 4> pc_med = {0., 0., 0., 0.};
	int partner=-1;
	Particle p_med;
	for (auto it = part_current.rbegin(); it != part_current.rend(); ++it) {
		if(it->CAT==3 && (it->kid1==p.kid1 || it->kid1==p.kid2)){
			partner = it->index();
			p_med = *it;
			std::cout << "FOUND partner! -->  " << partner << std::endl;  
		}
	}
	if(partner<0){
		std::cout << "ERROR! Couldn't find the missing medium partner... " << std::endl; 
		std::cout << "p.kid1 " << p.kid1 << std::endl; 
		std::cout << "p.kid2 " << p.kid2 << std::endl; 
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i<=3 ; i++) pc_med[i] =  p_med.P[i];


	//Find recoiled parton in 2->2, i.e. kid of p (leading parton in 2->2)
	int rec_kid=-1;
	int fin_kid=-1;
	Particle *p_rec_ptr = nullptr;
	Particle *p_fin_ptr = nullptr;
	bool FOUND1 = false;
	bool FOUND2 = false;
	for (auto it = part_current.rbegin(); it != part_current.rend(); ++it) {
		if((it->index()==p.kid1 || it->index()==p.kid2)){
			if(it->CAT==2){
				rec_kid = it->index();
				p_rec_ptr = &(*it);
				std::cout << "FOUND kid(rec)! -->  " << rec_kid << std::endl;  
				FOUND1 = true;
			}else if(it->CAT==0){
				fin_kid = it->index();
				p_fin_ptr = &(*it);
				std::cout << "FOUND kid(fin)! -->  " << fin_kid << std::endl;  
				FOUND2 = true;
			}
			if(FOUND1 && FOUND2) break;
		}
	}
	if(rec_kid<0 || fin_kid<0){
		std::cout << "ERROR! Couldn't find the missing kids... " << std::endl; 
		std::cout << "p.kid1 " << p.kid1 << std::endl; 
		std::cout << "p.kid2 " << p.kid2 << std::endl; 
		exit(EXIT_FAILURE);
	}


	// Step 2: Call radiation kernel
	//returns pc_rec, pc_rad, pc_fin
	bool success23 = collHQ23(p, qt, pc_med, pc_rec, pc_rad, pc_fin);
	//TODO: why?    //colljet23(T, qhat0, vc0, pc0, pc2, pc3, pc4, locqt, icl23, p.Tint_lrf, Ejp, iclrad);


	//Step 2.5: Additinal radiation or not?
       int nrad = KPoisson(p.radng);



	// Step 3: Process radiation if successful
	//   if (icl23 != 1 && iclrad != 1) {

        // Add radiated gluon: pc_rad
        Particle gluon;
        gluon.assign_index();
        for (int j = 0; j < 4; ++j) {
            gluon.P[j] = pc_rad[j];
            gluon.V[j] = p.V[j];
            gluon.Vfrozen[j] = p.V[j];
        }
        gluon.V[0] = -log(1.0-ran0(&config.rng.NUM1));
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
        gluon.parent2 = partner;//medium partner of p (leading parton)
        part_current.push_back(gluon);


	//Overwriting recoiled parton info in 2->2
	if (p_rec_ptr) {
		Particle& p_rec = *p_rec_ptr;  // Create a reference when needed
		for (int j = 0; j < 4; ++j) {
			p_rec.P[j] = pc_rec[j];
		}
	}else{
		std::cout << "ERROR: p_rec_ptr = nullptr. " << std::endl;
	}

	//Overwriting final state parton info in 2->2
	if (p_fin_ptr) {
		Particle& p_fin = *p_fin_ptr;  // Create a reference when needed
		for (int j = 0; j < 4; ++j) {
			p_fin.P[j] = pc_fin[j];
		}
	}else{
		std::cout << "ERROR: p_fin_ptr = nullptr. " << std::endl;
	}


	//Modify family tree and diactivate parent
	p.kid1 = gluon.index();
	p.kid2 = fin_kid;
	p.kid3 = rec_kid;
	p.isActive = false;



			//CHECKING
			std::cout << __FILE__ << "(" << __LINE__ << ")" << "After handleElastic" << std::endl;
			for (auto it = part_event.begin(); it != part_event.end(); ++it) {
				it->Print(true);
			}
			std::cout << "=======event part / current part =======" << std::endl;
			for (auto it = part_current.begin(); it != part_current.end(); ++it) {
				it->Print(true);
			}
			//CHECKING


        // Step 4: Handle multiple gluons (Poisson) for HQ
        while (--nrad > 0) {
//            double pc2_more[4] = {0.0};
//            double pc4_more[4] = {0.0};
//            double pb[4] = {0.0};
//
		std::array <double, 4> pc_rad1 = {0.,0.,0.,0.};
		std::array <double, 4> pc_rad2 = {0.,0.,0.,0.};
		std::array <double, 4> pc_fin12 = {0.,0.,0.,0.};
		bool successAdditional23 = radiationHQ(p, pc_rad, pc_fin, pc_rad1, pc_rad2, pc_fin12);
                //radiation(qhat0, vc0, pc4, pc2_more, pc4_more, pb,
                //          iclrad, Tdiff, Ejp);
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
        }
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


bool LBTcl::belowCutOff(const double Eloc, const Particle &p){

		double alpha_s = alphas0(config.physics.Kalphas, p.Tfrozen);  // Assuming alphas0() computes coupling
		double qhat0 = DebyeMass2(config.physics.Kqhat0, alpha_s, p.Tfrozen);  // qhat_0: Calculated by  \mu_D^2 = 4\pi \alpha_s T^2
		if(Eloc<sqrt(qhat0)) {
		      std::cout << __LINE__ << "Eloc < sqrt(qhat0)" << Eloc << "<" << sqrt(qhat0) << std::endl;
                      return true;
		  }
		if(Eloc<config.flow.Ecmcut) {
			std::cout << __LINE__ << " pc0[0]<Ecmcut " << Eloc <<  "<" << config.flow.Ecmcut<< std::endl;
			return true;
		}
		return false;
}



void LBTcl::LBT(std::vector<Particle> &part_event, double ti) {

    archive_np_snapshot((int)part_event.size());  // snapshot of part_event at start of this step "np0"

    // Loop over all active part_event at this step
    for (int i = 0; i < np_snapshot; ++i) {
        Particle &p = part_event[i];

	std::cout << "===============  i " << i << std::endl;
	std::cout << ".......             CAT     " << p.CAT << std::endl;
	std::cout << ".......             index " << p.index() << std::endl;
        // Skip frozen or inactive part_event
	//Only Vfrozen<ti AND active particle (Active AND (primary OR recoiled OR medium)) particles may join the loop
        //medium particle just propagates.
	if ((!p.isActive || (!p.isPrimary && p.CAT!=2 && p.CAT!=3)) || p.Vfrozen[0] >= ti) continue;
        if (p.P[0]<base::epsilon){
		std::cout << "  now this is going to be TERMINATED " << std::endl;
		for(int k=0;k<=3;k++) p.V[i]=0.0;
		p.CAT=1;
	}

	std::cout << ":):):):) Particle ....... " << i << "   at time " << ti << std::endl;
	std::cout << ".......             P     " << p.P[0] << std::endl;
	std::cout << ".......             index " << p.index() << std::endl;

        int free = 0;
        double fraction = 0.0;


        //Propagate parton.
	//In the future this should be just this->propagateParticle(p, ti, free, fraction);
	//Following is for the consistency between this code and original LBT.
	//=========================================
	if(p.CAT!=3 && p.CAT!=1){
		this->propagateParticle(p, ti, free, fraction);
		if((p.parent1>0 || p.parent2>0) && p.CAT==2){
			//find medium parent of this recoiled one
                        int medpart=-1;
			for(int j=0; j<(int)part_event.size(); j++){
				if (part_event[j].index() == p.parent1 || part_event[j].index() == p.parent2){
                                          if(part_event[j].isActive && part_event[j].CAT==3){
						  this->propagateParticle(part_event[j], ti, free, fraction);
						  medpart = part_event[j].index();
						  std::cout << "CHEK NEGATIVE PARTON PROPAGATION ========" << std::endl;
						  part_event[j].Print(true);
					  }
				}
			}
                        if (medpart<0){
				std::cout << "ERROR: Missing partner " << std::endl;
				exit(EXIT_FAILURE);
			}else{
				std::cout << "Missing partner FOUND --> " << medpart << std::endl;
			}
		}
	}
	if(p.CAT==3){
		std::cout << "THIS IS MEDIUM PARTON -- just propagation. " << std::endl;
		free=0;
		continue;
	}



	if (p.CAT != 1 && free == 0) {

                // Boost momentum into local fluid frame
                std::array<double,4> pc0 = {p.P[0], p.P[1], p.P[2], p.P[3]};
                std::array<double,4> vc0 = {0.0, p.vcfrozen[1], p.vcfrozen[2], p.vcfrozen[3]};
                this->trans(vc0, pc0);
                double Eloc = pc0[0];
		//See if below cut off 
		if(belowCutOff(Eloc,p)) continue;
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
		double throwdice = ran0(&config.rng.NUM1);
		std::cout << "--- throw dice!---  " << throwdice << std::endl;
		if (throwdice < probTot) {


			std::vector<Particle> part_current;
			int n_newparticle=0;


			std::cout << __FILE__ << "(" << __LINE__ << ")" << "Calling handleElasticCollision. " << std::endl;
			double qt = handleElasticCollision(p, PLenloc, part_event, part_current);
			n_newparticle++;//increment for recoiled parton. Medium parton is not counted.

			//CHECKING
			std::cout << __FILE__ << "(" << __LINE__ << ")" << "After handleElastic" << std::endl;
			for (auto it = part_event.begin(); it != part_event.end(); ++it) {
				it->Print(true);
			}
			std::cout << "=======event part / current part =======" << std::endl;
			for (auto it = part_current.begin(); it != part_current.end(); ++it) {
				it->Print(true);
			}
			//CHECKING

			std::cout << "--- throw dice for radiation??---  qt " << qt << std::endl;
			if(config.physics.Kinteraction==0 || (qt<base::epsilon)) continue;
			double throwdice_rad = ran0(&config.rng.NUM1);
			std::cout << "--- throw dice for radiation!---  " << throwdice_rad << std::endl;
			if (throwdice_rad < probRad / probTot) {
				std::cout << __FILE__ << "(" << __LINE__ << ")" << "Calling handleRadiation. " << std::endl;
				handleRadiation(p, qt, part_event, part_current);
				exit(1);
			}

			//Set V[0] (=t) for newly created particles.
			for (auto it = part_current.begin(); it != part_current.end(); ++it) {
                                    it->V[0]=-log(1.0-ran0(&config.rng.NUM1));
			}


			if(n_newparticle==0){
				std::cout << "ERROR n_newparticle = " << n_newparticle << std::endl;
				exit(EXIT_FAILURE);
			}
			part_event.insert(part_event.begin() + p.index(),  part_current.begin(), part_current.begin() + n_newparticle);
			part_event.insert(part_event.end(), part_current.begin() + n_newparticle, part_current.end());
			std::vector<Particle>().swap(part_current);
			np_snapshot += n_newparticle;


		}
		// Reset radiation tracker
		p.radng = 0.0;


            }//positive and free ==0(in medium)
    }//particle loop
    return;
}


