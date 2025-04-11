#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 10
#define RANDOM rand()/((double)RAND_MAX+1)
#define PI acos(-1)

void frecuencias(double *w);
void fase_inicial(double *theta);
void matriz_A_ER( int A[N][N]);
void matriz_A_BA( int A[N][N]);

void kuramoto(int A[N][N], double *theta, double *dtheta, double *w, double lambda);
void runge_kutta(int A[N][N], double *theta, double *w, double lambda, double dt);
double modulo_r(double *theta);
void printear_file(char *lugar, double variable_x, double *variable_y);

//FALTA EL KUTTA
int main()
{
    char lugar="results/r.txt";
    double w[N], theta[N], dtheta[N];
    double r;
    double ri[N]
    int A[N][N];
    double t_final, t_inicial, delta_t;
    double lambda_final, lambda_inicial, delta_lambda;
    int i, j;
    int pasos_t = (t_final-t_inicial)/delta_t;
    int pasos_lambda = (lambda_final-lambda_inicial)/delta_lambda;
    matriz_A_ER(A);
    frecuencias(w);
    //¿QUÉ THETAS INICIALES COGEMOS?
    for(i = 0;i<delta_lambda;i++){
        for(j = 0;j<delta_t; j++){
            runge_kutta(theta, w, lambda_inicial+delta_lambda*i, t_inicial + delta_t*j);
            r = modulo_r(theta);
            kuramoto(A, theta, dtheta, w, lambda_inicial+delta_lambda*i);

        }
        // ANTES TENEMOS QUE VER CUANDO TERMALIZA  printear_file("results/r.txt", lambda_inicial+delta_lambda*i, r[i]);
        ri[i]=r;
        printear_file(lugar,lambda_inicial+delta_lambda*i,ri[i])

    }
    //hay que correr en lambda y luego en t;
    
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

void matriz_A_ER( int A[N][N]){
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

void matriz_A_BA( int A[N][N]){
    int i, j;
    double p, random_aux;
    int nodo_random;
    //empezamos con 7 links todos conectados entre ellos
    for(i=0;i<7;i++){
        for(j=0;j<7;j++){
            if(i==j){
                A[i][i] = 0;
            }
            else{
                A[j][i] = 1;
            }
        }
    }
    for(i=0;i<N;i++){//Todos los elementos de la diagonal son 0
        A[i][i]=0;
    }
    for(i=7;i<N;i++){//hacemos el resto de nodos
        for(j=0;j<6;j++){//escogemos 6 nodos aleatorios entre los ya existentes
            nodo_random=(int)(i*RANDOM);
            if(A[nodo_random,i]==1){
                j--;
            }
            else{
                A[nodo_random, i]= 1;
                A[i, nodo_random]= 1;
            }
        }
    }
    
}

void kuramoto (int A[N][N], double *theta, double *dtheta, double *w, double lambda){ //devuelve dtheta y theta;

    int i, j;
    for(i=0;i<N;i++){
        dtheta[i] = w[i];
        for(j = 0;j<N;j++){
            if(A[i][j]==1){
                dtheta[i]+= lambda*sin(theta[j]-theta[i]);
            }
        };

    }

}
void runge_kutta(int A[N][N], double *theta, double *w, double lambda, double dt) {
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

void printear_file(char *lugar, double variable_x, double *variable_y){
    int i;
    FILE *f;
    f=fopen(lugar,"a");
    for(i=0;i<N;i++){
        f = fprintf(lugar, "%lf %lf \n",variable_x, variable_y[i]);
    }
    fclose(f);
}


