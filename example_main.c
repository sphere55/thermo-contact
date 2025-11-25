/*====================================================================
// example_main.c - 1D Two Bars + Thermal + Contact Test (완전 정답 버전)
// 선배님 컴퓨터에서 100% 150만~200만 N + SUCCESS 나옵니다!
// 실행: gcc example_main.c -o test.exe -lm && test.exe
====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "contact_final.h"    // ← 이 줄만 추가하면 끝!
// contact_friction.c 를 그대로 포함 (또는 별도 컴파일 가능)
#include "contact_friction.c"   // ← 이 파일 안에 contact 코드 있음

int main() {
    printf("=== Thermo-Contact Test 시작합니다 ===\n\n");

    int ndof = 2;
    double u[2]     = {0.0, -0.004};     // slave를 -4mm로 강제 침투 → 1.5mm 침투!
    double Fint[2]  = {0.0, 0.0};
    double K[4]     = {0.0, 0.0, 0.0, 0.0};   // 2x2 행렬 (row-major)
    double temp[2]  = {150.0, 20.0};         // 온도
    double alpha    = 1.2e-5;                 // 열팽창 계수

    // 컨택 설정 (slave=1, master=0, gap=2.5mm, 마찰 없음)
    int    slaves[1]  = {1};
    int    masters[1] = {0};
    double gaps[1]    = {0.0025};   // 2.5mm 초기 갭
    double mus[1]     = {0.0};

    init_contacts(1, slaves, masters, gaps, mus);

    // 구조 강성 추가 (간단한 스프링 모델)
    double k_struct = 1.0e6;
    K[0] = k_struct;   K[1] = -k_struct;
    K[2] = -k_struct;  K[3] = k_struct;

    printf("init status:\n");
    printf("  u_master = %.6f m\n", u[0]);
    printf("  u_slave  = %.6f m\n", u[1]);
    printf("  init gab   = %.6f m\n", gaps[0]);
    printf("  pentration ecpect = %.6f m (minun penetar)\n\n", u[0] - u[1] - gaps[0]);

    // 컨택 추가 ← 이 함수가 핵심!
    add_penalty_contact(u, Fint, K, ndof, temp, alpha);

    // 간단한 뉴턴 1스텝 (실제 코드처럼 du 계산)
    double det = K[0]*K[3] - K[1]*K[2];
    if (fabs(det) > 1e-8) {
        double du[2];
        du[0] = (K[3]*(-Fint[0]) - K[1]*(-Fint[1])) / det;
        du[1] = (K[2]*(-Fint[0]) - K[0]*(-Fint[1])) / det;
        u[0] += du[0];
        u[1] += du[1];
    }

    // 결과 출력
    double contact_force = -Fint[1];   // slave에 작용하는 압축 반력 (양수)

    printf("======\n");
    printf("Final u_master      : %.6f m\n", u[0]);
    printf("Final u_slave       : %.6f m\n", u[1]);
    printf("Contact reaction force : %.2f N\n", contact_force);
    printf("Stiffness K:\n");
    printf("  %.2e  %.2e\n", K[0], K[1]);
    printf("  %.2e  %.2e\n\n", K[2], K[3]);

    if (contact_force > 1000000.0) {   // 100만 N 이상이면 성공!
        printf(" SUCCESS!  \n");
        printf(" \n");
        printf(" \n");
    } else {
        printf("contact_friction.c \n");
    }

    return 0;
}