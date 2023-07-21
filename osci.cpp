#include <cstdio>
#include <cstdlib>
#include <cmath>
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// Runtime constants
const static double Epsilon = 1e-10; // Defines the precision of energy calculations
const static int N_of_Divisions = 1000;
const static int N_max = 5; // Number of calculated Eigenstates

// Physical parameters using dimensionless quantities
const static double xRange = 12;
const static double h_0 = xRange / N_of_Divisions; 

const static double Psi_left = 1.0e-3; // left boundary condition 
const static double Psi_right = 0.0; // right boundary condition

double* CalculateEigenEnergies(double* E_pot) {
    double* Psi = new double[N_of_Divisions+1];
    double* EigenEnergies = new double[N_max+1];
    Psi[0] = Psi_left;
    Psi[1] = Psi_left + 1.0e-3; // Add arbitrary small value

    int N_quantum; //N_quantum is Energy Quantum Number
    int Nodes_plus; // Number of nodes (+1) in wavefunction
    double K_square; // Square of wave vector
    double E_lowerLimit = 0.0; // Eigen-energy must be positive 
    double E_upperLimit = 10.0;
    int End_sign = -1;
    bool Limits_are_defined = false;
    double Normalization_coefficient;
    double E_trial;

    // Calculations

    for (N_quantum=1; N_quantum <= N_max; ++N_quantum) {
        Limits_are_defined = false;

        while (Limits_are_defined == false) {
            Nodes_plus = 0;
            E_trial = E_upperLimit;

            for (int i=2; i <= N_of_Divisions; ++i) {
                K_square = 2.0*(E_trial - E_pot[i]);
                Psi[i] = 2.0*Psi[i-1]*(1.0 - (5.0*h_0*h_0*K_square / 12.0)) /(1.0 + (h_0*h_0*K_square/12.0))-Psi[i-2];
                
                if (Psi[i]*Psi[i-1] < 0) {
                    ++Nodes_plus;
                }
            }

            if (E_upperLimit < E_lowerLimit) {
                E_upperLimit = MAX(2*E_upperLimit, -2*E_upperLimit);
            }

            if (Nodes_plus > N_quantum) {
                E_upperLimit *= 0.7;
            } else if (Nodes_plus < N_quantum) {
                E_upperLimit *= 2.0;
            } else {
                Limits_are_defined = true; 
            }
        }

        End_sign = -End_sign;

        while ((E_upperLimit - E_lowerLimit) > Epsilon) {
            E_trial = (E_upperLimit + E_lowerLimit) / 2.0;
            
            for (int i=2; i <= N_of_Divisions; ++i) {
                K_square = 2.0*(E_trial - E_pot[i]);
                Psi[i] = 2.0*Psi[i-1] * (1.0 - (5.0*h_0*h_0*K_square / 12.0)) / (1.0 + (h_0*h_0*K_square/12.0)) - Psi[i-2]; 
            }

            if (End_sign*Psi[N_of_Divisions] > Psi_right) {
                E_lowerLimit = E_trial;
            } else {
                E_upperLimit = E_trial;
            }
        }

        E_trial = (E_upperLimit+E_lowerLimit)/2;
        EigenEnergies[N_quantum] = E_trial;
        E_upperLimit = E_trial;
        E_lowerLimit = E_trial;
    }
    
    return EigenEnergies;
}

double* CalculatePotentialEnergy() {
    double* E_pot = new double[N_of_Divisions+1]; 
    double dist;
    
    for (int i = 0; i <= N_of_Divisions; ++i) {
        dist = i*h_0 - 0.5*xRange;
        E_pot[i] = 0.5*dist*dist;
    }

    return E_pot;
}
