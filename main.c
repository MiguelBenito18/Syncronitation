#include <stdio.h>
#include <iostream>
#include <math.h>

#define N 10000 //Numero de nodos de la red
#define L 9999 //Numero de intervalos del histograma
#define Medidas 100
#define PI acos(-1)
#define redes

using namespace std;

double genNumRandom (double a, double b){
    double ran;
    ran=a+(b-a)*rand()/((double)RAND_MAX+1);
    return ran;
}
void histograma(int *data, double *H, int Ndata, int Nintervalos, double *d, int *m, int *M){
    double delta=0.0;
    int minimo, maximo;
    double norm;
    int i,indice;
    minimo=data[0];
    maximo=minimo;
    //calculo del valor maximo y minimo
    for (i=1; i<Ndata;i++){
        if (data[i]>maximo){
            maximo=data[i];
        }
        if (data[i]<minimo){
            minimo=data[i];
        }
    }
    *m=minimo;
    *M=maximo;
    //calculo de la delta
    delta=(double)(maximo-minimo)/Nintervalos;
    *d=delta;
    //Queremos leer los k conexiones posibles que van desde 0 hasta N-1
    delta=1.0;
    //creacion de histograma
    for (i=0;i<Nintervalos;i++){
        H[i]=0;
    }
    for (i=0;i<Ndata;i++){
        indice=(data[i])/delta;
        if (indice==Nintervalos){
            indice--;
        }
        H[indice]++;
    }

    //Normalizar Histograma
    norm=1.0/(delta*Ndata);
    for (i=0;i<Nintervalos;i++){
        H[i]*=norm;
    }

}
double promedio(double *x, int d){
    int i;
    double med=0.0;
    for (i=0;i<d;i++){
        med+=x[i];
    }
    med/=d;
    return med;
}
#ifdef Kuramoto
void matrizAKuramoto(int *A[N]){
    int i,j;
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            if (i==j){
                A[i][j]=0;
            }
            else{
                A[i][j]=1;
            }
        }
    }
}
#endif // Kuramoto
#ifdef redes
void cuentaVecinos(int *A[N], int *k, int i, int j){
    if (A[i][j]==1){
        k[i]++;
        k[j]++;
    }
}
/* Redes a mano, emplearemos las de Python
double prob(int *k,double alpha, int d,int j){
    int i;
    double p,sum;
    sum=0.0;
    for (i=0;i<(d-1);i++){
        sum+=k[i]+alpha;
    }
    p=(k[j]+alpha)/sum;
    return p;
}
void matrizAredes(int* A[N], double alpha, int conectividad){
    int i,j,cuenta;
    int k[N];
    double p,num;
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            A[i][j]=0;
        }
        k[i]=0;
    }
    if (alpha==1.0){
        p=(double) conectividad/(N-1.0);
        for (i=0;i<N;i++){
            for (j=0;j<(i+1);j++){
                if (j==i){
                    A[i][j]=0;
                }
                else{
                    num=genNumRandom(0,1);
                    if (num<p){
                        A[i][j]=1;
                    }
                    else{
                        A[i][j]=0;
                    }
                }
                A[j][i]=A[i][j];
            }
        }
    }
    else{
        for (i=0;i<N;i++){
            cuenta=0;
            for (j=0;j<(i+1);j++){
                if (j==i){
                    A[i][j]=0;
                }
                else{
                    if (i<5){
                        A[i][j]=1;
                    }
                    else{
                        if (cuenta<conectividad/2){
                            p=prob(k,alpha,i+1,j);
                            num=genNumRandom(0.0,1.0);
                            if (num<=p){
                                A[i][j]=1;
                                cuenta++;
                            }
                            else{
                                A[i][j]=0;
                            }
                        }
                        else{
                            A[i][j]=0;
                        }
                    }
                }
                A[j][i]=A[i][j];
                cuentaVecinos(A,k,i,j);
                printf("%d",A[i][j]);
            }
            printf("\n%d\n",i+1);
        }
    }
}
*/
//Redes las generamos en Python, cada linea del fichero de texto es una arista
void nombre(char *name, double alpha){
    if (alpha==1.0){
        strcpy(name,"red_ER.txt");
    }
    else{
        strcpy(name,"red_BA.txt");
    }
}
void matrizAredes(int *A[N], double alpha, int *k){
    int i,j;
    char name[L];
    nombre(name,alpha);
    FILE *f;
    f=fopen(name,"r");
    if (f==NULL){
        printf("No se ha podido abrir el fichero\n");
    }
    //Inicializamos la red a 0
    for (i=0; i<N;i++){
        for (j=0;j<(i+1);j++){
            A[i][j]=0;
            A[j][i]=0;
        }
        k[i]=0;
    }
    while (fscanf(f,"%d %d",&i,&j)==2){
        A[i][j]=1;
        A[j][i]=1;
        k[i]++;
        k[j]++;
    }
    fclose(f);
}
#endif // redes
void kuramoto (int *A[N], double *theta, double *dtheta, double *w, double lambda){ //devuelve dtheta y theta;
    int i, j;
    for(i=0;i<N;i++){
        dtheta[i] = w[i];
        for(j = 0;j<N;j++){
            if(A[i][j]==1){
                dtheta[i]+= lambda*sin(theta[j]-theta[i]);
            }
        }
    }
}
void igualaVector(double *v1, double *v2, int d){
    for (int i=0;i<d;i++){
        v1[i]=v2[i];
    }
}
double moduloImaginario(double preal, double pimaginaria){
    double modulo;
    modulo=sqrt(preal*preal+pimaginaria*pimaginaria);
    return modulo;
}

void RK4 (double h,double *w, double lambda, double *theta, int *A[N]){
    int i;
    double k1[N], k2[N], k3[N], k4[N];
    double thetacopia[N];
    //Calculo k1
    kuramoto(A,theta,k1,w,lambda);
    //Calculo k2
    for (i=0;i<N;i++){
        thetacopia[i]=theta[i]+h/2.0*k1[i];
    }
    kuramoto(A,thetacopia,k2,w,lambda);
    //Calculo k3
    for (i=0;i<N;i++){
        thetacopia[i]=theta[i]+h/2.0*k2[i];
    }
    kuramoto(A,thetacopia,k3,w,lambda);
    //Calculo k4
    for (i=0;i<N;i++){
        thetacopia[i]=theta[i]+h*k3[i];
    }
    kuramoto(A,thetacopia,k4,w,lambda);
    //Calculo de la nueva theta
    for (i=0;i<N;i++){
        theta[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }

}

double calculoR(double *theta, int d){
    double rreal=0.0;
    double rim=0.0;
    double r;
    int i;
    for (i=0;i<d;i++){
        rreal+=cos(theta[i]);
        rim+=sin(theta[i]);
    }
    rreal/=(double)d;
    rim/=(double)d;
    r=moduloImaginario(rreal,rim);
    return r;
}
//FUNCIONES FICHERO DE REPRESENTACION LAMBDA FRENTE A R
void abrirFichero(FILE *f){
    f=fopen("results.txt","w");
    fclose(f);
}
void escribirFichero(FILE *f, double lambda, double r){
    f=fopen("results.txt","a");
    fprintf(f,"%lf\t%lf\n",lambda,r);
    fclose(f);
}
//FICHERO PARA LA W_EFF
void abrirFichero_weff(FILE *f){
    f=fopen("w_eff.txt","w");
    fclose(f);
}
void escribirFichero_weff(FILE *f, double lambda, double *w_eff){
    f=fopen("w_eff.txt","a");
    fprintf(f,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",lambda,w_eff[0],w_eff[1],w_eff[2],w_eff[3],w_eff[4],w_eff[5],w_eff[6],w_eff[7],w_eff[8],w_eff[9]);
    fclose(f);
}
//FUNCIONES FICHERO DE REPRESENTACION HISTOGRAMA P(k)
void abrirHistograma (FILE*f){
    f=fopen("pk.txt","w");
    fclose(f);
}
void escribirHistograma(FILE*f, double H, int i){
    f=fopen("pk.txt","a");
    fprintf(f,"%d\t%lf\n",i,H);
    fclose(f);
}
void kuramoto_star(double *theta, double *dtheta, double *w, double lambda){ //devuelve dtheta y theta;
    int i, j;
    double omega = (10+10)/(11.0);
    double suma_h=0;
    for(i = 1;i<11;i++){
        suma_h+=sin(theta[i]-theta[0]);
        dtheta[i] = w[i]- omega + lambda*sin(theta[i]-theta[0]);
    }
    dtheta[0] = w[0]-omega+lambda*suma_h;

}
void RK4_star (double h,double *w, double lambda, double *theta){
    int i;
    double k1[11], k2[11], k3[11], k4[11];
    double thetacopia[11];
    //Calculo k1
    kuramoto_star(theta,k1,w,lambda);
    //Calculo k2
    for (i=0;i<11;i++){
        thetacopia[i]=theta[i]+h/2.0*k1[i];
    }
    kuramoto_star(thetacopia,k2,w,lambda);
    //Calculo k3
    for (i=0;i<11;i++){
        thetacopia[i]=theta[i]+h/2.0*k2[i];
    }
    kuramoto_star(thetacopia,k3,w,lambda);
    //Calculo k4
    for (i=0;i<11;i++){
        thetacopia[i]=theta[i]+h*k3[i];
    }
    kuramoto_star(thetacopia,k4,w,lambda);
    //Calculo de la nueva theta
    for (i=0;i<11;i++){
        theta[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }

}

void red_estrella(){
    int i, j;
    FILE *f;
    f = fopen("estrella.txt", "wt");
    double theta[11];
    double  w[11];
    int medidas, deltalambda, deltat, lambdapasos, tpasos;
    double lambda, lambdai, lambdaf, ti, tf, h;
    //inicializamos wi y thetas;
    for(i = 0;i<11;i++){
        theta[i] = genNumRandom(-PI, PI);
        w[i] = 1;
    }
    w[0] = 10;
    int sentido = 0;
    lambdai = 0, lambdaf = 2.6, ti = 0, tf = 10, deltalambda = 0.1, h = 0.01;
    lambdapasos = (int)((lambdaf-lambdai)/deltalambda);
    tpasos = (int)((tf-ti)/h);
    double r[2*lambdapasos], r_aux = 0;
    lambda = lambdai;
    for(sentido = 0;sentido<2;sentido++){
        for(i = 0;i<lambdapasos;i++){
            for(j = 0;j<tpasos;j++){
                RK4_star(h,w, lambda, theta);
                if(j>=tpasos/3){
                    r_aux += calculoR(theta, 11);
                }
            }
            r[sentido*lambdapasos+ i] = r_aux/(tpasos-tpasos/3);
            printf("%lf %lf \n", lambda, r[sentido*lambdapasos+ i]);
            fprintf(f, "%lf %lf \n", 0.0, 0.0);
            r_aux = 0;
            lambda+=deltalambda;
        }
        lambda = lambdaf;
        deltalambda = -deltalambda;

    }
    fclose(f);



}
//Nuevas ideas
int mayorConectividad(int *k){
    int i,m;
    m=0;
    for (i=1;i<N;i++){
        if (k[i]>=k[m]){
            m=i;
        }
    }
    return m;
}
int main(){
    FILE *f;
    FILE *g;
    double alpha=0.0;
    double theta[N], w[N],thetacopia[N],w_eff[N],thetapunto[N]; //fase y frecuencia angular
    double h;//paso de tiempo RK
    int i,j,k,t,s,sum;//para los arrays
    //double erre[Medidas]; //Array de r para poder hacer promedio Descartable
    double r;//Calculo de la r en promedio
    double lambdai,lambdaf,deltalambda,lambda;//para simular las graficas de r frente a lambda
    double tfinal, tinicial;//Simulacion temporal
    int pasoslambda,pasost; //pasos dados para lambda mas pasos de estacionamiento al variar lambda
    /*
    //Variables para el calculo de r
    double rsum;
    int medidasestacionarias;
    */
    //Variables para contruccion de histograma
    int m,M;
    double d,H[N];
    //Matriz de adyacencia y array de vecinos
    int ki[N];
    int** A = new int*[N];       // Crear un array de punteros (filas)
    for (i = 0; i < N; i++) {
        A[i] = new int[N];       // Crear cada fila
    }
    //Calculo de la Matriz de adyacencia para la red
    #ifdef Kuramoto
    matrizAKuramoto(A);
    #endif // Kuramoto
    #ifdef redes
    matrizAredes(A,alpha,ki);
    #endif // redes
    //He quitado la parte del estado inicial aquí
    /*
    Escribir P(k)
    histograma(ki,H,N,L,&d,&m,&M);
    abrirHistograma(f);
    for (i=0;i<L;i++){
        escribirHistograma(f,H[i],i+1);
        printf("%d\t%lf\n",i+1,H[i]);
    }
    */
    for (j=0;j<N;j++){
        theta[j]=genNumRandom(-PI,PI);
        w[j]=ki[j];
        //w[j]=genNumRandom(-0.5,0.5);
    }

    //ALGORITMO PARA REPRESENTAR LAMBDA EN FUNCION DE R
    //ABRIR FICHERO
    abrirFichero(f);
    abrirFichero_weff(g);
    //PASOS DEL ARRAY
    tinicial=0.0;
    tfinal=100;
    h=0.01;
    lambdai=0.0;
    lambdaf=2.0;
    deltalambda=0.02;
    pasoslambda=(int)(lambdaf-lambdai)/deltalambda;
    pasost=(int)(tfinal-tinicial)/h;
    lambda=lambdai;
    //Luego quito los dos printf, es para asegurarme que funcionan
    printf("%lf\n",h);
    printf("%d",pasost);
    double w_promedio_k[10]; //la calculamos solo para 10 conectividades diferentes
    int conectividades_weff[]={1,2,3,4,5,6,7,8,9,10};
    for(int i=0;i<10;i++){
        w_promedio_k[i]=0;
    }
    r=0;
    sum=0;
    for (j=0;j<2;j++){
        for (i=0;i<pasoslambda;i++){
            r=0;
            for(int nodo=0;nodo<N;nodo++){
                w_eff[nodo]=0;
            }
            for (t=0;t<pasost;t++){
                RK4(h,w,lambda,theta,A);
                //Calculo de w_eff
                if(t>(pasost/2)){//empezamos a integrar cuando el sistema ya est� evolucionado
                    kuramoto(A, theta, thetapunto, w, lambda);//Para calcular w_eff necesitamos thetapunto
                    for(int nodo=0;nodo<N;nodo++){
                        w_eff[nodo]+=h*thetapunto[nodo];
                    }
                }
            }
            //calculo de r como promedio
            for(t=0;t<100;t++){
                RK4(h,w,lambda,theta,A);
                r+=calculoR(theta,N);
            }
            r=r/100;//promediamos 100 r's en el equilibrio
            escribirFichero(f,lambda,r);

            //calculo w_eff para una lambda concreta
            for(int nodo=0;nodo<N;nodo++){
                w_eff[nodo]=2*w_eff[nodo]/tfinal;//weff/T donde T=tfinal/2
            }

            for(int k=0;k<9;k++){
                w_promedio_k[k]=0;
                int suma=0;
                for(int nodo=0;nodo<N;nodo++){
                    if(ki[nodo]==conectividades_weff[k]){
                        suma++;
                        w_promedio_k[k]+=w_eff[nodo];
                    }
                }
                if(suma!=0){
                    w_promedio_k[k]=w_promedio_k[k]/suma;
                }
            }
            if(j==0){
                escribirFichero_weff(g,lambda, w_promedio_k);
            }


            lambda+=deltalambda;
        }


        deltalambda=-deltalambda;
        lambda+=deltalambda;
    }


    /*
    for (j=0;j<2;j++){
        for (i=0;i<pasoslambda;i++){
            r=0;
            for(int nodo=0;nodo<N;nodo++){
                w_eff[nodo]=0;
            }
            for (t=0;t<pasost;t++){
                //Estabilizamos el sistema
                RK4(h,w,lambda,theta,A);
                //Calculo de w_eff
                kuramoto(A, theta, thetapunto, w, lambda);//Para calcular w_eff necesitamos thetapunto
                for(int nodo=0;nodo<N;nodo++){
                    w_eff[nodo]+=h*thetapunto[nodo];
                }
            }
            //calculo de r como promedio
            for(t=0;t<100;t++){
                RK4(h,w,lambda,theta,A);
                r+=calculoR(theta,N);
            }
            r=r/100;//promediamos 100 r's en el equilibrio
            escribirFichero(f,lambda,r);

            //calculo w_eff para una lambda concreta
            for(int nodo=0;nodo<N;nodo++){
                w_eff[nodo]=w_eff[nodo]/tfinal;
            }

            for(int k=0;k<10;k++){
                int suma=0;
                for(int nodo=0;nodo<N;nodo++){
                    if(ki[nodo]==conectividades_weff[k]){
                        suma++;
                        w_promedio_k[k]+=w_eff[nodo];
                    }
                }
                if(suma!=0){
                    w_promedio_k[k]=w_promedio_k[k]/suma;
                }
            }
            if(j==0){
                escribirFichero_weff(g,lambda, w_promedio_k);
            }


            lambda+=deltalambda;
        }


        deltalambda=-deltalambda;
        lambda+=deltalambda;
    }
    */
    for (i=0;i<N;i++){
        delete[] A[i];
    }
    return 0;
}
