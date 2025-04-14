#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 200
#define RANDOM rand()/((double)RAND_MAX)
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
    for(i=0;i<N;i++){
        w[i] = -0.5 + RANDOM;
    }

}

void fase_inicial(double *theta){
    int i;
    for(i=0; i<N; i++){
        theta[i]=-PI * 2*PI*RANDOM;
    }
}

void matriz_A_ER( int A[N][N], int k){
    int i, j;
    double p, random_aux;
    p = ((double)k)/(N-1);
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
    for(i=0;i<N;i++){//Los elementos de la diagonal son ceros
        A[i][i]=0;
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
            m=0;
            while(random_aux>0){//método para elejir el nodo segun la probabilidad Barabasi Albert
                random_aux=random_aux-conectividad(A,m)/(k*i/2+j);
                m++;
            }
            if(A[m-1][i]==1){//Nos aseguramos de no volver a escoger el mismo nodo
                j--;
            } else{
                A[m-1][i]=1;
            }
        }
    }
    for(i=0;i<N;i++){//Los elementos de la diagonal son ceros
        A[i][i]=0;
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
    double temp[N], theta_orig[N];

    for (int i = 0; i < N; i++) {
        theta_orig[i] = theta[i];
    }

    kuramoto(A, theta_orig, k1, w, lambda);

    for (int i = 0; i < N; i++) {
        temp[i] = theta_orig[i] + 0.5 * dt * k1[i];
    }
    kuramoto(A, temp, k2, w, lambda);

    for (int i = 0; i < N; i++) {
        temp[i] = theta_orig[i] + 0.5 * dt * k2[i];
    }
    kuramoto(A, temp, k3, w, lambda);

    for (int i = 0; i < N; i++) {
        temp[i] = theta_orig[i] + dt * k3[i];
    }
    kuramoto(A, temp, k4, w, lambda);

    for (int i = 0; i < N; i++) {
        theta[i] = theta_orig[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

double modulo_r(double *theta){
    int i;
    double p_real=0, p_imaginaria=0, r;
    for(i=0;i<N;i++){
        p_real+= cos(theta[i]);
        p_imaginaria+= sin(theta[i]);
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
        fprintf(lugar, "%lf %lf \n",variable_x[i], variable_y[i]);
    }
    fclose(f);
}


void genera_histograma(double *x, int x_data, double *histo, int n_histo, double *min, double *delta){
    int i,indice;
    double norma;
    double  histo_min= x[0], histo_max=x[0];
    for(i = 0;i<x_data;i++){
        if(histo_min>x[i]){
            histo_min = x[i];
        }
        if(histo_max<x[i]){
            histo_max = x[i];
        }
    }
    double longitud = (histo_max-histo_min)/n_histo;
    norma = x_data*longitud;

    for(i=0;i<n_histo;i++){
        histo[i]=0;
    }
    for(i=0;i<x_data;i++){
        indice = (x[i]-histo_min)/longitud;
        if(indice == n_histo){
            indice = n_histo-1;
        }
        histo[indice]++;
    }
    for(i=0;i<x_data;i++){
        histo[i] = histo[i]/norma;
    }
    *min = histo_min;
    *delta = longitud;




}
void distribucion_grado_parte_0_2(int A[N][N],int k){
    /*CON N 200 SE VE BIEN, PARA N GRANDES NO CORRE*/
    matriz_A_ER(A, k);
    FILE *f1;
    f1 = fopen("conexiones_cada_nodo_ER.txt", "w");
    double conexiones_vector[N];
    double histo[10];
    int conexiones_aux = 0;
    for(int i=0;i<N;i++){
        for(int j = 0;j<N;j++){
            conexiones_aux+=A[i][j];
            }
        conexiones_vector[i] = ((double)conexiones_aux);
        fprintf(f1, "%d\n", conexiones_aux);
        conexiones_aux = 0;

        }
    fclose(f1);
    double histo_min, delta_histo;
    genera_histograma(conexiones_vector,N, histo, 10, &histo_min, &delta_histo);
    FILE *f;
    f = fopen("histograma_ER.txt", "w");
    for(int i=0;i<10;i++){
        fprintf(f, "%lf %lf\n",histo_min+delta_histo*i,histo[i]);
    }
    fclose(f);

}

int main()
{
    FILE *f;
    f = fopen("resultados_r.txt", "w");
    double w[N], theta[N], dtheta[N], fase_comienzo[N];
    double r;
    int A[N][N];
    double t_final=100000, t_inicial=0, delta_t=1;//Por probar xd
    double lambda_final=2, lambda_inicial=0.5 , delta_lambda=0.2;//Más o menos como en el artículo
    int i, j;
    int pasos_t = (t_final-t_inicial)/delta_t;
    int pasos_lambda = (lambda_final-lambda_inicial)/delta_lambda;
    matriz_A_ER(A,6);
    frecuencias(w);
    int t_thermal = 10000;
    //¿QUÉ THETAS INICIALES COGEMOS?
    fase_inicial(theta);
    for(i = 0;i<pasos_lambda;i++){
        for(j = 0;j<pasos_t; j++){
            runge_kutta(A, theta, w, lambda_inicial+delta_lambda*i, t_inicial + delta_t*j);
            if (j >= t_thermal) { //MEDIMOS A PARTIR DE AQUI, QUE SUPONGO QUE LLEGA A ESTADO ESTACIONARIO PERO NO ES ASÍ
            r += modulo_r(theta);
            }

        }
        r /= (t_final - t_thermal); //HACEMOS LA MEDIA DE LAS MEDIDAS DENTRO DEL IF
        printf("%lf ", r);
        fprintf(f, "%lf %lf \n",lambda_inicial+delta_lambda*i, r );
        r = 0;

    }


    return 0;

}
