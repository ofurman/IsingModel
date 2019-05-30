//
//  main.cpp
//  isingmodel
//
//  Created by Oleksii Furman on 16/05/2019.
//  Copyright © 2019 Oleksii Furman. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

using namespace std;

int modulo(int a,int b);
int delta_energy(int spin,int nnspins_sum,int J=1);
void calculate_magnetisation_and_energy(int size, int** lattice, double &magn, double &energy);
float calculate_susceptibility(float  mean_magn, float  mean_magn_squared, int N, float T = 1.7);
float calculate_heat_cap(float mean_energy, float mean_energy_squared, int N, float T = 1.7);
float calculate_bc(float mean_M_squared, float mean_M_power_four);

int** initialize(int size);

void _cy_ising_update(int size, int** &field, int i, int j, float T=1, int J=1);
void cy_ising_step(int size, int** &lattice, float T);
void mcs_simulation(int size, unsigned long steps, unsigned long K0);

int main(int argc, const char * argv[]) {
    mcs_simulation(10, 230000, 30000);
    mcs_simulation(20, 230000, 30000);
    mcs_simulation(50, 230000, 30000);
    return 0;
}

int modulo(int a,int b) {
    int c = (int)a % (int)b;
    if (c < 0) {
        return c+b;
    } else {
        return c;
    }
}

int delta_energy(int spin,int nnspins_sum,int J){
    return 2*J*spin*nnspins_sum;
}

void calculate_magnetisation_and_energy(int size, int** lattice, double &magn, double &energy) {
    int S, nnspins_sum;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            magn += lattice[i][j];
            S = lattice[i][j];
            nnspins_sum = lattice[modulo((i - 1),size)][j] + lattice[modulo((i + 1), size)][j] +
                          lattice[i][modulo((j - 1), size)] + lattice[i][modulo((j + 1), size)];
            energy += -nnspins_sum*S;
        }
    }
    energy = energy/(2*size*size);
    
}

float calculate_susceptibility(float  mean_magn, float  mean_magn_squared, int N, float T) {
    float susceptibility = (1.0/(N*T))*(mean_magn_squared-pow(mean_magn, 2));
    return susceptibility;
}

float calculate_heat_cap(float  mean_energy, float mean_energy_squared, int N, float T) {
    float heat_cap = N/(float)(pow(T, 2))*(mean_energy_squared-pow(mean_energy, 2));
    return heat_cap;
}

float calculate_bc(float mean_M_squared, float mean_M_power_four) {
    float bc = 1-(mean_M_power_four/(3*(pow(mean_M_squared, 2))));
    return bc;
}

int** initialize(int size) {
    mt19937 gen(234);
    uniform_int_distribution<int> brandom(0, 1);
    int choose[2] = {-1,1};
    int** lattice = new int*[size];
    for(int i = 0; i < size; ++i)
        lattice[i] = new int[size];
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            lattice[i][j] = choose[brandom(gen)];
        }
    }
    
    return lattice;
}

void _cy_ising_update(int size, int**& field, int i, int j, float T, int J) {
    int nnspins_sum = field[modulo((i - 1),size)][j] + field[modulo((i + 1), size)][j] +
                      field[i][modulo((j - 1), size)] + field[i][modulo((j + 1), size)];
    int dU = delta_energy(field[i][j], nnspins_sum);
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    if (dU <= 0) {
        field[i][j] *= -1;
    } else if (exp(-dU/T) > r) {
        field[i][j] *= -1;
    }
}

void cy_ising_step(int size, int**& lattice, float T) {
    int n_offset, m_offset, n, m;
    for (n_offset = 0; n_offset < 2; n_offset++) {
        for (m_offset = 0; m_offset < 2; m_offset++) {
            for (n = n_offset; n < size; n=n+2) {
                for (m = m_offset; m < size; m=m+2) {
                    _cy_ising_update(size, lattice, n, m, T);
                }
            }
        }
    }
}

void mcs_simulation(int size, unsigned long steps, unsigned long K0) {
    float temp[34] = {1.2,1.35, 1.5, 1.65, 1.79, 1.94, 2.        , 2.03571429, 2.07142857, 2.10714286, 2.14285714,
        2.17857143, 2.21428571, 2.25      , 2.269, 2.28571429, 2.32142857,
        2.35714286, 2.39285714, 2.42857143, 2.46428571, 2.5, 2.55, 2.67, 2.78, 2.89, 3., 3.12, 3.26, 3.41, 3.55, 3.7, 3.85, 4.};
    int** lattice = initialize(size);
    ofstream datafile;
    datafile.open("isingmodeldata" + to_string(size) + ".csv");
    datafile << "magnet" << "," << "energy" << "," << "susc" << "," << "heat_cap" << "," << "bc" << "," << "temp" << "\n";
    for (int t = 0; t < 34; t++) {
        float T = temp[t];
        int factor = 100;
        float magnetization[(steps-K0)/factor];
        float energy[(steps-K0)/factor];
        float susceptibility[(steps-K0)/factor];
        float heat_cap[(steps-K0)/factor];
        float bc[(steps-K0)/factor];
        for (int i = 0; i < steps; i++) {
            cy_ising_step(size, lattice, T);
            if (i < K0) {
                continue;
            }
            if (modulo(i, factor) == 0) {
                double m = 0;
                double en = 0;
                int N = pow(size, 2);
                float mean_magnetization = 0;
                float mean_squared_magnetization = 0;
                float mean_energy = 0;
                float mean_squared_energy = 0;
                float mean_M_power_four = 0;
                calculate_magnetisation_and_energy(size, lattice, m, en);
                magnetization[(i-K0)/factor] = m;
                energy[(i-K0)/factor] = en;
                
                for (int n = 0; n <= (i-K0)/factor;  n++) {
                    mean_magnetization += abs(magnetization[n]);
                    mean_energy += energy[n];
                    mean_squared_magnetization += pow(magnetization[n], 2);
                    mean_squared_energy += pow(energy[n], 2);
                    mean_M_power_four += pow(magnetization[n], 4);
                }
                mean_magnetization = mean_magnetization/(((i-K0)/factor)+1);
                mean_energy = mean_energy/(((i-K0)/factor)+1);
                mean_squared_magnetization = mean_squared_magnetization/(((i-K0)/factor)+1);
                mean_squared_energy = mean_squared_energy/(((i-K0)/factor)+1);
                mean_M_power_four = mean_M_power_four/(((i-K0)/factor)+1);
                susceptibility[(i-K0)/factor] = calculate_susceptibility(mean_magnetization, mean_squared_magnetization, N, T);
                heat_cap[(i-K0)/factor] = calculate_heat_cap(mean_energy, mean_squared_energy, N, T);
                bc[(i-K0)/factor] = calculate_bc(mean_squared_magnetization, mean_M_power_four);
            }
        }
        
        //TEMPERATURE ITERATOR
        float mean_magnetization_for_temp = 0;
        float mean_energy_for_temp = 0;
        float mean_susceptibility_for_temp = 0;
        float mean_heat_cap_for_temp = 0;
        float mean_bc = 0;
        for (int i = 0 ; i < (steps-K0)/factor; i ++) {
            mean_magnetization_for_temp += abs(magnetization[i]);
            mean_energy_for_temp += energy[i];
            mean_susceptibility_for_temp += susceptibility[i];
            mean_heat_cap_for_temp += heat_cap[i];
            mean_bc += bc[i];
        }
        mean_magnetization_for_temp = mean_magnetization_for_temp/(((steps-K0)/factor)*(size*size));
        mean_energy_for_temp = mean_energy_for_temp/((steps-K0)/factor);
        mean_susceptibility_for_temp = mean_susceptibility_for_temp/((steps-K0)/factor);
        mean_heat_cap_for_temp = mean_heat_cap_for_temp/((steps-K0)/factor);
        mean_bc = mean_bc/((steps-K0)/factor);
        // Write data to file
        datafile << mean_magnetization_for_temp << "," << mean_energy_for_temp << "," <<mean_susceptibility_for_temp << "," << mean_heat_cap_for_temp << "," << mean_bc << "," << T << "\n";
    }
    datafile.close();
}

