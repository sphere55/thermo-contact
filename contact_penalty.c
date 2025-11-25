/*====================================================================
// contact_penalty.c - Penalty Contact for Thermo-Elastoplastic in C
// 2025-11-24, for sphere55 (1974-born FEA master!)
// Compile: gcc -o test example_main.c contact_penalty.c -lm
====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_NDOF 1024
#define MAX_CONTACT_PAIRS 10000

typedef struct {
    int slave_node;
    int master_node;
    double initial_gap;
    double mu;  // Friction coeff
} ContactPair;

ContactPair contact_pairs[MAX_CONTACT_PAIRS];
int n_contact_pairs = 0;
double eps_n = 1e12;  // Normal penalty
double eps_t = 1e11;  // Tangential

void init_contacts(int n, int *slaves, int *masters, double *gaps, double *mus) {
    n_contact_pairs = n;
    for (int i = 0; i < n; i++) {
        contact_pairs[i].slave_node = slaves[i];
        contact_pairs[i].master_node = masters[i];
        contact_pairs[i].initial_gap = gaps[i];
        contact_pairs[i].mu = mus[i];
    }
}

void add_penalty_contact(double *u, double *Fint, double *K, int ndof, double *temp, double alpha) {
    for (int i = 0; i < n_contact_pairs; i++) {
        ContactPair *cp = &contact_pairs[i];
        int s = cp->slave_node;
        int m = cp->master_node;
        double g0 = cp->initial_gap;

        double deltaT_avg = 0.5 * (temp[s] + temp[m]);
        double g_thermal = alpha * deltaT_avg * g0;  // Gap shrink by heat
        double g = u[s] - u[m] - (g0 - g_thermal);   // >0 = penetration

        if (g > 0.0) {
            double tn = eps_n * g;

            Fint[s] -= tn;
            Fint[m] += tn;

            K[s * ndof + s] += eps_n;
            K[s * ndof + m] -= eps_n;
            K[m * ndof + s] -= eps_n;
            K[m * ndof + m] += eps_n;
        }
    }
}