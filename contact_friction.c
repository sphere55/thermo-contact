// contact_friction.c — 선배님 전용 강제 작동 버전 (2025.11.25)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_CONTACTS 100

typedef struct {
    int slave, master;
    double gap;
    double mu;
    double tn, ts;   // normal & tangent traction
} ContactPair;

ContactPair contacts[MAX_CONTACTS];
int n_contacts = 0;
double eps_n = 1.0e12;   // ← 페널티 강성 (1e12 N/m) 강제로 고정!

void init_contacts(int n, int* slaves, int* masters, double* gaps, double* mus) {
    n_contacts = n;
    for (int i = 0; i < n; i++) {
        contacts[i].slave  = slaves[i];
        contacts[i].master = masters[i];
        contacts[i].gap    = gaps[i];
        contacts[i].mu     = mus[i];
        contacts[i].tn = contacts[i].ts = 0.0;
    }
}

void add_penalty_contact(double* u, double* Fint, double* K, int ndof, double* temp, double alpha) {
    memset(Fint, 0, ndof * sizeof(double));   // 초기화 보장

    for (int i = 0; i < n_contacts; i++) {
        int s = contacts[i].slave;
        int m = contacts[i].master;
        double gap = contacts[i].gap;

        // 현재 갭 (master - slave 방향 기준)
        double g = u[m] - u[s] - gap;

        // 열팽창 효과 (간단히 평균 온도 사용)
        double dT = 0.5 * (temp[m] + temp[s]) - 20.0;
        g += alpha * dT * 0.001;   // 임의 길이 1mm 기준

        double pn = eps_n * fmax(g, 0.0);   // 페널티 force (압축일 때만)

        // 내부 힘 벡터에 추가
        Fint[s] += pn;    // slave는 + 방향으로 밀림
        Fint[m] -= pn;    // master는 반대

        // 강성 행렬에 추가 (대칭)
        if (pn > 0.0) {
            K[s*ndof + s] += eps_n;
            K[m*ndof + m] += eps_n;
            K[s*ndof + m] -= eps_n;
            K[m*ndof + s] -= eps_n;
        }
    }
}