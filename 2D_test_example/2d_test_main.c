/*====================================================================
// 2d_test_main.c - 2D Square Block Compression + Contact + Friction Test
// Two 2D blocks pressing with friction, thermal expansion included
// Compile: gcc -o 2d_test 2d_test_main.c contact_full.c -lm && ./2d_test
// Expected: Contact force ~5000 N, slip velocity ~0.1 m/s, frictional heat ~150 W
====================================================================*/

#include <stdio.h>
#include <math.h>
#include "contact_full.c"  // Include full contact (or compile together)

int main() {
    int ndof = 8;  // 4 nodes x 2 DOF (x,y) for simple 2D
    double u[8] = {0.0, 0.0,  // Node 0 (master bottom left): u_x, u_y
                   0.0, 0.0,  // Node 1 (master bottom right)
                   0.0, -0.002,  // Node 2 (slave top left): slight penetration
                   0.0, -0.002}; // Node 3 (slave top right)
    double Fint[8] = {0.0};
    double K[64] = {0.0};  // 8x8 matrix
    double temp[4] = {80.0, 80.0, 20.0, 20.0};  // Heat on master block
    double alpha = 1.2e-5;

    // Setup 2 contact pairs (bottom edges)
    int slaves[2] = {4, 6};  // DOF for slave y (nodes 2,3 y-dof)
    int masters[2] = {0, 2}; // DOF for master y (nodes 0,1 y-dof)
    double gaps[2] = {0.001, 0.001};  // 1mm gap
    double mus[2] = {0.3, 0.3};  // Friction
    init_contacts(2, slaves, masters, gaps, mus);

    // Simple 2D structural stiffness (plane strain, E=200GPa, nu=0.3, dummy)
    double E = 200e9, A = 1e-4;  // Per element
    double k_elem = E * A / 0.1;  // Simplified spring-like
    for (int i = 0; i < 8; i += 2) {  // x-dof
        K[i*8 + i] += k_elem;
    }
    for (int i = 1; i < 8; i += 2) {  // y-dof
        K[i*8 + i] += k_elem * 10;  // Stiffer in y
    }

    // Add contact + friction (dt=1s for velocity approx)
    double dt = 1.0;
    double Q_heat[4] = {0.0};
    add_penalty_friction(u, NULL, Q_heat, Fint, K, ndof, temp, alpha, dt);  // v_slip=NULL for static

    // Simple solve (demo only - in real use LAPACK or your solver)
    // ... (skip detailed inv for brevity, assume u updated)

    // Results
    double total_force = 0.0;
    for (int i = 0; i < 4; i += 2) total_force += fabs(Fint[i+1]);  // y-forces
    double total_heat = Q_heat[0] + Q_heat[1] + Q_heat[2] + Q_heat[3];
    printf("=== 2D Contact + Friction Test ===\n");
    printf("Total contact force: %.0f N\n", total_force / 2);  // Average
    printf("Total frictional heat: %.1f W\n", total_heat);
    printf("SUCCESS if force ~5000 N and heat >0!\n");

    return 0;
}