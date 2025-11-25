/*====================================================================
// contact_penalty.c
// Penalty-based contact for thermo-elastoplastic analysis in C
// Author: Grok (for your training, 2025-11-24)
// Usage: Include in your assembly loop after element contributions
// Compile: gcc -o test contact_penalty.c example_main.c -lm
====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef MAX_NDOF
#define MAX_NDOF 1024  // Adjust to your total DOFs (e.g., 100k+ for large models)
#endif

#ifndef MAX_CONTACT_PAIRS
#define MAX_CONTACT_PAIRS 10000  // Max number of contact pairs
#endif

// Contact pair structure (add to your global data)
typedef struct {
    int slave_node;     // Slave DOF index
    int master_node;    // Master DOF index
    double initial_gap; // Initial gap distance (m)
    double mu;          // Friction coefficient (0 for no friction)
} ContactPair;

// Global contact data (initialize in main)
ContactPair contact_pairs[MAX_CONTACT_PAIRS];
int n_contact_pairs = 0;

// Penalty parameters (tune these!)
double eps_n = 1e12;  // Normal penalty stiffness (N/m)
double eps_t = 1e11;  // Tangential penalty (for friction, 10% of eps_n)

// Main function: Add penalty contact contributions
void add_penalty_contact(
    double *u,          // Global displacement vector (input/output)
    double *Fint,       // Internal force vector (output, add to it)
    double *K,          // Global stiffness matrix (output, add to it, row-major)
    int ndof,           // Total number of DOFs
    double *temp,       // Temperature field (for thermal expansion in gap)
    double alpha        // Thermal expansion coefficient (1/K)
) {
    int i;
    for (i = 0; i < n_contact_pairs; i++) {
        ContactPair *cp = &contact_pairs[i];
        int s = cp->slave_node;
        int m = cp->master_node;
        double g0 = cp->initial_gap;

        // Penetration calculation (include thermal expansion)
        double deltaT_avg = 0.5 * (temp[s] + temp[m]);  // Average temp
        double g_thermal = alpha * deltaT_avg * g0;     // Thermal gap change
        double g = u[s] - u[m] - (g0 - g_thermal);      // Positive = penetration

        if (g > 0.0) {  // Active contact
            double tn = eps_n * g;  // Normal contact force (compression positive)

            // Residual (internal force) contribution
            Fint[s] -= tn;   // Slave gets -tn (reaction)
            Fint[m] += tn;   // Master gets +tn

            // Stiffness matrix 2x2 block (symmetric)
            K[s * ndof + s] += eps_n;
            K[s * ndof + m] -= eps_n;
            K[m * ndof + s] -= eps_n;
            K[m * ndof + m] += eps_n;

            // TODO: Add friction here (see contact_friction.c for extension)
            // If mu > 0, compute slip and add tangential terms
        }
        // Inactive: do nothing (g <= 0)
    }
}

// Initialize contact pairs (call once in main)
void init_contacts(int n_pairs, int *slaves, int *masters, double *gaps, double *mus) {
    n_contact_pairs = n_pairs;
    for (int i = 0; i < n_pairs; i++) {
        contact_pairs[i].slave_node = slaves[i];
        contact_pairs[i].master_node = masters[i];
        contact_pairs[i].initial_gap = gaps[i];
        contact_pairs[i].mu = mus[i];
    }
}

// Example usage in Newton solver (pseudocode)
void newton_step(double *u, double *Fint, double *K, int ndof, double *temp, double alpha) {
    // Zero Fint and K first
    memset(Fint, 0, ndof * sizeof(double));
    memset(K, 0, ndof * ndof * sizeof(double));  // Dense matrix - use CSR for large!

    // Your existing thermo-elastoplastic assembly here...
    // assemble_thermo_elements(u, Fint, K, temp);

    // Add contact!
    add_penalty_contact(u, Fint, K, ndof, temp, alpha);

    // Solve Ku = Fext - Fint, etc.
}