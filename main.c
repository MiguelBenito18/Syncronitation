#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 10
#define RANDOM rand()/((double)RAND_MAX+1)
#define PI acos(-1)


int conectividad(int A[N][N], int nodo){
    int i;
    int conect = 0;
    for(i=0;i<N;i++){
        conect+=A[nodo][i];
    }
    return conect;
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

void matriz_A_ER( int A[N][N], int k){
    int i, j;
    double p, random_aux;
    p = k/199;
    for(i=0;i<N;i++){//Los elementos de la diagonal son ceros
        A[i][i]=0;
    }    
    for(i=0;i<(N-1);i++){
        for(j=(i+1);j<N;j++){
            random_aux = RANDOM;
            if(random_aux<p){
                A[i][j] = 1;
                A[j][i] = 1;
            }
            else{
                A[i][j] = 0;
                A[j][i] = 0;
            }
        }
    }
}

void matriz_A_BA( int A[N][N], int k){//funciona solo para k par
    int i, j,m;
    double p, random_aux;
    int nodo_random;
    
    for(i=0;i<N;i++){//Empezamos con todo ceros
        for(j=0;j<N;j++){
            A[i][j]=0;
        }
    }
    //empezamos con k+1 nodos todos conectados entre ellos
    for(i=0;i<(k+1);i++){
        for(j=0;j<(k+1);j++){
            if(i!=j){
                A[i][j]=1;
            }
        }
    }
    //hacemos el resto de nodos
    for(i=(k+1);i<N;i++){
        for(j=0;j<(k/2);j++){//hacemos k/2 conexiones por cada nodo nuevo creado
            random_aux=RANDOM;
            m=0
            while(random_aux>0){//método para elejir el nodo segun la probabilidad Barabasi Albert
                random_aux=random_aux-conectividad(int A[N][N],m)/(k*i/2+j);
                m++;
            }
            if(A[m-1][i]==1){//Nos aseguramos de no volver a escoger el mismo nodo
                j--;
            } else{
                A[m-1][i]=1;
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
    kuramoto(A, theta, k1, w, lambda);

    // Calcular k2
    for (int i = 0; i < N; i++) {
        theta[i] = theta[i] + 0.5 * dt * k1[i];
    }
    kuramoto(A, theta, k2, w, lambda);

    // Calcular k3
    for (int i = 0; i < N; i++) {
        theta[i] = theta[i] + 0.5 * dt * k2[i];
    }
    kuramoto(A, theta, k3, w, lambda);

    // Calcular k4
    for (int i = 0; i < N; i++) {
        theta[i] = theta[i] + dt * k3[i];
    }
    kuramoto(A, theta, k4, w, lambda);

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

void printear_file(char *lugar, double *variable_x, double *variable_y){
    int i;
    FILE *f;
    f=fopen(lugar,"a");
    for(i=0;i<N;i++){
        f = fprintf(lugar, "%lf %lf \n",variable_x[i], variable_y[i]);
    }
    fclose(f);
}

//FALTA EL KUTTA
int main()
{
    char *lugar="results/r.txt";
    double w[N], theta[N], dtheta[N], fase_comienzo[N];
    double r;
    double ri[N]
    int A[N][N];
    double t_final=10; double t_inicial=0; double delta_t=0.01;//Por probar xd
    double lambda_final=1.6; lambda_inicial=0.4; delta_lambda=0.02;//Más o menos como en el artículo
    int i, j;
    int pasos_t = (t_final-t_inicial)/delta_t;
    int pasos_lambda = (lambda_final-lambda_inicial)/delta_lambda;
    matriz_A_ER(A);
    frecuencias(w);
    /*
    //¿QUÉ THETAS INICIALES COGEMOS?
    fase_inicial(theta);
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

*/
    
    //hay que correr en lambda y luego en t
    //EVOLUCIÓN TEMPORAL Y CÁLCULO DE r PARA CADA LAMBDA
    fase_inicial(fase_comienzo);

    for(i=0;i<pasos_lambda;i++){
        lambdas[i]=lambda_inicial+i*delta_lambda;//guardamos las lambdas para las que calculamos r
        for(int k = 0; k < N; k++) {//fijamos las fases a las fases de comienzo
            theta[k] = fase_comienzo[k];
        }

        for(j=0;j<pasos_t;j++){//evolucion temporal en pasos_t pasos de delta_t
            runge_kutta(A, theta, w,(lambda_inicial+i*delta_lambda), delta_t);
        }
        r[i]=modulo_r(theta);//cálculo de r (asumo que ya hemos llegado al equilibrio)
    }

    //Pasar los datos a un fichero
    FILE *f = fopen(lugar, "w");
    if (f == NULL) {
        printf("Error: No se pudo crear el archivo %s\n", lugar);
        return;
    }

    for(int i = 0; i < pasos_lambda; i++) {
        fprintf(f, "%lf %lf\n", lambdas[i], r[i]);
    }

    fclose(f);
    
}
