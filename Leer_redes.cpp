#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>

#define N 10000
#define L 9999

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
    delta=1.0;
    *d=delta;
    //creacion de histograma
    for (i=0;i<Nintervalos;i++){
        H[i]=0;
    }
    for (i=0;i<Ndata;i++){
        indice=(data[i]-minimo)/delta;
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

void escribeFichero(FILE *f, int i, double H){
    f=fopen("P(k)_BA_Bueno.txt", "a");
    fprintf(f,"%d\t%lf\n", i, H);
    fclose(f);
}

using namespace std;

int main(){

FILE *in, *out;

int i, j;

int k[N];  //sabemos que hay N datos
double H[L];

double delta;
int minimo, maximo;


for(i=0; i<N; i++){
    k[i]=0;
}

in=fopen("red_BA.txt", "r");

if(in==NULL){
    printf("Error al abrir el fichero de entrada");
    return 1;
  }else{
      //printf("El fichero de entrada se ha abierto correctamente");

      while(fscanf(in, "%d  %d", &i, &j)==2){
        k[i]++;
        k[j]++;
      }
      fclose(in);

      /*
      double suma=0.0;
      for(i=0; i<N; i++){
        //printf("%d\t %d\n", i, k[i]);
        suma+=k[i];
      }
      suma=suma/N;
      printf("La media es: %lf", suma);
     */

      histograma(k, H, N, L, &delta, &minimo, &maximo);
      /*
      double norma=0.0;
      for(i=0; i<L; i++){
        norma+=H[i];
      }
      printf("La norma es: %lf", norma);
      */
      out=fopen("P(k)_BA_Bueno.txt", "w");
      if(out==NULL){
         printf("Error al abrir el fichero de salida");
         return 1;
      }else{
          //printf("El fichero de salida se ha abierto correctamente");
          fclose(out);

         for(j=0; j<L; j++){
            escribeFichero(out, j+1, H[j]);
         }
         fclose(out);
      }
    }
return 0;
}
