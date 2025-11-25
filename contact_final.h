// contact_final.h  (당신 기존 코드에 #include "contact_final.h" 만 추가하면 끝!)
#ifndef CONTACT_FINAL_H
#define CONTACT_FINAL_H

#include <math.h>
#include <string.h>

#define MAX_CONTACTS 100000   // 10만개 컨택까지 지원 (당신 코드 규모에 딱!)

typedef struct {
    int slave, master;
    double gap, mu;
    double tn, ts, tn_old;
    double penetration;
    int active;
} ContactPair;

static ContactPair contacts[MAX_CONTACTS];
static int n_contacts = 0;
static double eps_n = 1.0e12;      // 기본 페널티 강성 (당신이 나중에 조정 가능)
static double eps_t = 1.0e11;      // 마찰용 페널티
static double heat_fraction = 0.9; // 마찰열 90% 변환

void set_penalty_stiffness(double en, double et) { eps_n = en; eps_t = et; }
void set_frictional_heating(double frac) { heat_fraction = frac; }

void init_contacts(int n, int* slaves, int* masters, double* gaps, double* mus) {
    n_contacts = n;
    for (int i = 0; i < n; i++) {
        contacts[i].slave = slaves[i];
        contacts[i].master = masters[i];
        contacts[i].gap = gaps[i];
        contacts[i].mu = mus[i];
        contacts[i].tn = contacts[i].ts = contacts[i].tn_old = 0.0;
        contacts[i].active = 0;
    }
}

void add_thermo_friction_contact(
    double* u, double* v, double* temp, double* dtemp,
    double* Fint, double* K, int ndof, double alpha, double dt)
{
    memset(Fint, 0, ndof * sizeof(double));  // 필요시만

    for (int i = 0; i < n_contacts; i++) {
        int s = contacts[i].slave;
        int m = contacts[i].master;
        double dx = u[m] - u[s] - contacts[i].gap;

        // 열팽창 효과 추가 (당신 기존 코드와 동일하게)
        double avgT = 0.5 * (temp[m] + temp[s]);
        dx += alpha * avgT * 1e-3;   // 대표길이 1mm 가정 (당신 코드에 맞춰 조정)

        double g = dx;  // 현재 penetration (양수면 침투)
        double pn = eps_n * fmax(g, 0.0);

        // 법선 방향 (간단히 1D처럼, 3D는 나중에 벡터로 확장)
        double nx = (g > 0.0) ? 1.0 : 0.0;

        // 마찰 (Stick-Slip)
        double vt = v[m] - v[s];
        double ts_trial = contacts[i].ts + eps_t * vt * dt;
        double fmax = contacts[i].mu * pn;
        if (fabs(ts_trial) <= fmax) {
            contacts[i].ts = ts_trial;           // stick
        } else {
            contacts[i].ts = (ts_trial > 0 ? 1 : -1) * fmax;  // slip
        }

        // 힘 벡터
        Fint[s] += pn * nx - contacts[i].ts;
        Fint[m] -= pn * nx + contacts[i].ts;

        // 강성 (간단 대칭)
        if (pn > 0.0) {
            K[s*ndof + s] += eps_n;
            K[m*ndof + m] += eps_n;
            K[s*ndof + m] -= eps_n;
            K[m*ndof + s] -= eps_n;
        }

        // 마찰열 → 온도 증가
        double power = fabs(contacts[i].ts * vt);
        if (power > 0.0) {
            dtemp[s] += heat_fraction * power * dt / (1.0e-6);  // 임의 열용량
            dtemp[m] += heat_fraction * power * dt / (1.0e-6);
        }
    }
}

#endif // CONTACT_FINAL_H