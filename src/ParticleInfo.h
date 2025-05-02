#ifndef PARTICLE_H
#define PARTICLE_H

struct Particle {
        double P[4] = {0.0};        // Four-momentum (E, px, py, pz)
        double V[4] = {0.0};        // Position (t or tau, x, y, z)
        double P0[4] = {0.0};       // Negative particle momentum
        double V0[4] = {0.0};       // Negative particle position
        int CAT = 0;                // Category: 0=active, 1=free streaming, 2=recoiled, 4=radiated
        int CAT0 = 0;               // Negative particle category
        int KATT = 0;               // Flavor code
        int KATT0 = 0;              // Negative flavor code
        double Vfrozen[4] = {0.0};  // Frozen position (for interpolation)
        double Vfrozen0[4] = {0.0}; // Frozen negative
        double Tfrozen = 0.0;       // Temperature when frozen
        double Tfrozen0 = 0.0;      // Temperature when frozen (negative)
        double vcfrozen[4] = {0.0}; // Flow velocity when frozen
        double vcfrozen0[4] = {0.0};
        double Tint_lrf = 0.0;      // Interaction time in local rest frame
        double radng = 0.0;         // Average radiated gluons
        double WT = 1.0;            // Weight
        double WT0 = 1.0;           // Negative particle weight
        double mass = 0.0;          // Mass (0 for light partons, positive for heavy quarks)
        bool isPrimary = false;     // True if parton is among the original initial partons (instead of nj)
        int mom1 = -1;  // Index of primary mother parton
        int mom2 = -1;  // Index of secondary (e.g. recoil) mother parton
        int index = -1;  // Unique ID for this particle in the vector
};
#endif
