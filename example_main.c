/*====================================================================
// example_main.c - 1D Two Bars + Thermal + Contact Test
// Run: gcc -o test example_main.c contact_penalty.c -lm && ./test
// Expected: Reaction force ~1000 N if compressed beyond gap
====================================================================*/

#include <stdio.h>
#include <math.h>
#include "contact_penalty.c"  // Include (or compile together)

int main() {
    int ndof = 2;  // u[0]=master, u[1]=slave
    double u[2] = {0.0, -0.0015};  // Initial: slight penetration
    double Fint[2] = {0.0, 0.0};
    double K[4] = {0.0};  // 2x2 matrix (row-major)
    double temp[2] = {150.0, 20.0};  // Heat on master bar
    double alpha = 1.2e-5;  // Thermal exp coeff

    // Setup one contact pair: slave=1, master=0, gap=1mm, no friction
    int slaves[1] = {1};
    int masters[1] = {0};
    double gaps[1] = {0.001};
    double mus[1] = {0.0};
    init_contacts(1, slaves, masters, gaps, mus);

    // Simple assembly: add structural stiffness (e.g., two springs k=1e6)
    double k_struct = 1e6;
    K[0] += k_struct; K[0 + 1] -= k_struct / 2;  // Simplified
    K[1] += k_struct; K[1 + 0] -= k_struct / 2;

    // Add contact
    add_penalty_contact(u, Fint, K, ndof, temp, alpha);

    // Simple solve: du = K^{-1} * (-Fint)  (demo only)
    double det = K[0]*K[3] - K[1]*K[2];
    if (fabs(det) > 1e-6) {
        double du0 = (K[3] * (-Fint[0]) - K[1] * (-Fint[1])) / det;
        double du1 = (K[2] * (-Fint[0]) - K[0] * (-Fint[1])) / det;
        u[0] += du0;
        u[1] += du1;
    }

    // Results
    double g_final = u[1] - u[0] - 0.001 + alpha * 0.5*(temp[0]+temp[1])*0.001;
    double reaction = eps_n * fmax(g_final, 0.0);
    printf("=== Thermo-Contact Test Results ===\n");
    printf("Final penetration (with thermal): %.2e m\n", g_final);
    printf("Contact reaction force: %.2f N\n", reaction);
    printf("Stiffness K (after contact):\n");
    printf("  %.2e   %.2e\n", K[0], K[1]);
    printf("  %.2e   %.2e\n", K[2], K[3]);

    // Success check: If ~1000N, contact works!
    if (reaction > 500) {
        printf("SUCCESS! Contact activated. Ready for your full code.\n");
    }
    return 0;
}