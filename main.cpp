//
//  main.cpp
//  nematiccrystals
//
//  Created by Oleksii Furman on 17/05/2019.
//  Copyright © 2019 Oleksii Furman. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>
#include <fstream>
#include <map>

#define size 20
int lattice[size][size];

using namespace std;

double calculate_energy(int trial, int i, int j);

int modulo(int a,int b) {
    int c = (int)a % (int)b;
    if (c < 0) {
        return c+b;
    } else {
        return c;
    }
}

void initialize(mt19937& eng) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            lattice[i][j] = 0;
        }
    }
}

void cy_ising_step(mt19937& eng, uniform_real_distribution<double>& dist, float T) {
    int i, j;
    double U1, U0, dU;
    for (i=0; i < size; i++) {
        for (j=0; j<size; j++) {
            int fiNew = lattice[i][j] + round((rand()%1000/1000.0-0.5)*10);
            if (fiNew >= 90) {
                fiNew -= 180;
            } else if (fiNew <= -90){
                fiNew += 180;
            }
            
            U1 = calculate_energy(fiNew, i, j);
            
            U0 = calculate_energy(lattice[i][j], i, j);
            
            dU = U1 - U0;
            
            if (dU <= 0) {
                lattice[i][j] = fiNew;
            } else if (exp(-dU/T) > rand()%1000/1000.0) {
                lattice[i][j] = fiNew;
            }
        }
    }
}

double calculate_Qxx(){
    double Qxx = 0.0;
    for(int i = 0; i<size; i++)
        for(int j = 0; j<size; j++) Qxx += 2 * pow(cos(lattice[i][j]*M_PI/180), 2) - 1;
    return Qxx/pow(size, 2);
}
double calculate_Qxy(){
    double Qxy = 0.0;
    for(int i = 0; i<size; i++)
        for(int j = 0; j<size; j++) Qxy += 2 * cos(lattice[i][j]*M_PI/180) * sin(lattice[i][j]*M_PI/180);
    return Qxy/pow(size, 2);
}

double calculate_S(){
    double Qxx = calculate_Qxx();
    double Qxy = calculate_Qxy();
    return sqrt(pow(Qxx, 2) + pow(Qxy, 2));
}

double calculate_energy(int trial, int i, int j){
    return  -1.5*(pow(cos(M_PI*(trial-lattice[modulo((i - 1),size)][j])/180), 2)
                  +pow(cos(M_PI*(trial-lattice[modulo((i + 1),size)][j])/180), 2)
                  +pow(cos(M_PI*(trial-lattice[i][modulo((j - 1), size)])/180), 2)
                  +pow(cos(M_PI*(trial-lattice[i][modulo((j + 1), size)])/180), 2))+2;
}

void parallel_temperature(unsigned long steps, unsigned long K0, double T) {
    random_device rd;
    mt19937 eng(rd());
    uniform_real_distribution<double> dist(0,1);
    
    initialize(eng);
    
    int factor = 1000;
    double S = 0;
    int i;
    for (i = 0; i < steps; i++) {
        cy_ising_step(eng, dist, T);
        if (i < K0) {
            continue;
        }
        if (modulo(i, factor) == 0) {
            S += calculate_S();
        }
    }
    //TEMPERATURE ITERATOR
    S = S/((steps-K0)/factor);
    cout << S << "," << T << "\n";
}

int main(int argc, const char * argv[]) {
//    string arg = argv[1];
//    double temperature = atof(arg.c_str());
    srand(time(NULL));
    parallel_temperature(230000, 30000, 0.8);
    return 0;
}
