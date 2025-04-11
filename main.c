#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 10
#define RANDOM rand()/((double)RAND_MAX)
#define PI acos(-1)

void frecuencias(double *w);
void fase_inicial(double *theta);
void matriz_A( int **A);
void kuramoto (int **A, double *theta, double *dtheta, double *w, double lambda);
void runge_kutta(int **A, double *theta, double *w, double lambda, double dt);
double modulo_r(double *theta);

//FALTA EL KUTTA
int main()
{
    double t_final, delta_t, t_pasos;
    double lamda_final, delta_lamda, lamda_pasos;
    //hay que correr en lamda y luego en t;
}

void frecuencias(double *w){
    int i;
    for(i=0;i<N;i++);{
        w[i] = -0.5 + RANDOM;
    }

}

void fase_inicial(double *theta){
    int i;
    for(i=0; i<N; i++){
        theta[i]=2*PI*RANDOM;
    }
}

void matriz_A_ER( int **A){
    int i, j;
    double p, random_aux;
    p = 6/199;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            random_aux = RANDOM;
            if(random_aux<p){
                A[i][j] = 1;
            }
            else{
                A[i][j] = 0;
            }
        }
    }
}

void kuramoto (int **A, double *theta, double *dtheta, double *w, double lambda){ //devuelve dtheta y theta;

    int i, j;
    for(i=0;i<N;i++){
        dtheta[i] = w[i];
        for(j = 0;j<N;j++){
            if(A[i][j]==1){
                dtheta[i]+= sin(theta[j]-theta[i]);
            }
        }
        dtheta[i] = dtheta[i]*lambda;

    }

}
void runge_kutta(int **A, double *theta, double *w, double lambda, double dt) {
    double k1[N], k2[N], k3[N], k4[N];
    double dtheta[N];

    // Calcular k1
    kuramoto(A, dtheta, k1, w, lambda);

    // Calcular k2
    for (int i = 0; i < N; i++) {
        dtheta[i] = theta[i] + 0.5 * dt * k1[i];
    }
    kuramoto(A, dtheta, k2, w, lambda);

    // Calcular k3
    for (int i = 0; i < N; i++) {
        dtheta[i] = theta[i] + 0.5 * dt * k2[i];
    }
    kuramoto(A, dtheta, k3, w, lambda);

    // Calcular k4
    for (int i = 0; i < N; i++) {
        dtheta[i] = theta[i] + dt * k3[i];
    }
    kuramoto(A, dtheta, k4, w, lambda);

    for (int i = 0; i < N; i++) {
        theta[i] += (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

double modulo_r(double *theta){
    int i;
    double p_real=0, p_imaginaria=0, r;
    for(i=0;i<N;i++){
        p_real+= cos(theta[i]);
        p_imaginaria+= sen(theta[i]);
    }
    p_real= p_real/N;
    p_imaginaria = p_imaginaria/N;
    r = sqrt(p_real*p_real + p_imaginaria*p_imaginaria);
    return r;
}


