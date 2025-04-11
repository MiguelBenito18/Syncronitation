#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define NIP_SEED 000000
#define N_SHAPES 17890
#define MAX_STRING_LENGTH 256
#define CITY_CENTER ((Point){.x = 676452.6, .y = 4613323})
#define MAX_DISTANCE_CITY 3000
#define GAUSSIAN_FACTOR 0.00005
#define MAX_STEPS 100000

typedef struct Point {
    double x;
    double y;
} Point;


