Thermo-Elastoplastic Contact in C - Quick Start Guide
====================================================

1. Compilation:
   - Save the 3 .c files separately.
   - gcc -o mycontact contact_penalty.c contact_friction.c example_main.c -lm
   - ./mycontact  → Should print gap and force.

2. Integration to Your Code:
   - Add ContactPair struct and globals to your headers.
   - In main: Call init_contacts() with your mesh contact pairs.
   - In assembly loop (after elements): add_penalty_contact(u, Fint, K, ndof, temp, alpha);
   - For friction/heat: Use add_penalty_friction() instead, pass Q_heat to thermal solver.
   - Matrix: Assumes dense row-major K[ndof*ndof]. For sparse, modify to CSR additions.

3. Tuning:
   - eps_n: Start 1e12. If penetration >1e-6m, increase 10x. If diverge, decrease.
   - Thermal gap: g_thermal = alpha * ΔT * gap_length (adjust length to your geom).
   - Test: Run example, change temp[0]=0 → gap remains open, force=0.

4. Extensions:
   - Node-to-segment: Average over segment nodes in loop.
   - Newton: Call in while loop until ||du|| < tol.
   - 3D: Set NDOF_PER_NODE=3, index u[s], u[s+1], u[s+2] for x,y,z.

Questions? Ask me for tweaks!