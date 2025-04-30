#include <stdio.h>
#include <iostream>
#include <math.h>

#define N 10000 //Numero de nodos de la red
#define PI acos(-1)

using namespace std;

void abreFicheroEscritura(FILE *g){
g=fopen("pk_ER_Normalizado.txt", "w");
fclose(g);
}

void escribeFichero(FILE *g, int i, double H){
    g=fopen("pk_ER_Normalizado.txt", "a");
    fprintf(g,"%d\t%lf\n", i, H);
    fclose(g);
}

int main(){

int m, j;
int num[N];
double H[N];
double suma=0.0;
double norma=0.0;
FILE *f, *g;

f=fopen("pk_ER10000.txt", "r");

m=0;

while(fscanf(f, "%d %lf", &num[m], &H[m]) == 2) {
        //printf("Entero: %d, Decimal: %lf\n", num[j], H[j]);
        m++;
    }
fclose(f);


for(j=0; j<m+1; j++){
    suma+=H[j];
}

abreFicheroEscritura(g);

for (j=0; j<m+1; j++){
    H[j]=H[j]/suma;
    escribeFichero(g, num[j], H[j]);
    norma+=H[j];   //la norma tiene que salir 1
    //printf("%lf", H[j]);
}
 
//printf("%lf    %lf", suma, norma); para comprobar que sale bien

}
