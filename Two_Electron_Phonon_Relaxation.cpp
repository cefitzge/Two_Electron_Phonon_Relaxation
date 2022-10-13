#include <iostream>
#include <cmath>
#include <string>
#include <tgmath.h>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_dawson.h>
//#include <boost/chrono.hpp>
//#include <Eigen/Dense>


using namespace std;
//using namespace Eigen;

double pi = 3.14159;
double hbar = 6.626 * pow(10, -34) / (2 * pi);
//double epsilon_r = 11.68;
double e = 1.602 * pow(10.0, -19);
//double k_const = 8.9876 * pow(10.0, 9);
double m = 9.11 * pow(10, -31);
double m_perp = .198 * m;
double m_z = .92 * m;
double omega0 = .004 * e / hbar; //8 meV
double E_VS = .58 * e / 1000; // valley splitting 
double g = 1.998;
double r = 4.8 * pow(10.0, -9); //dipole size, 4.8 nm
double rho = 2200; //silicon mass density 2200 kg/m^3 in SiO2 (Peihao's paper)
double v_t = 3750; //3750, 5420
double v_l = 5900; //5900, 9330

double Rashba = 45;
double Dressel = 0;

double xi_d = 5.0 * e; //dilation deformation 5 eV
double xi_u = 8.77 * e; //uniaxial sheer deformation 8.77 eV



double Zeeman_Energy(double B) {
    double ub = e * hbar / (2 * m);
    return g * ub * B;
}

double sin_gamma_sq(double B) {
    double e3 = abs(E_VS - Zeeman_Energy(B));
    double delta23 = 2.0*r * m_perp * (E_VS) * (Dressel + Rashba) / (sqrt(2) * hbar);
    double e3_tilde = sqrt(e3 * e3 + delta23 * delta23);
    return (sqrt(e3 * e3 + delta23 * delta23) + e3) / (2 * sqrt(e3 * e3 + delta23 * delta23));
}

double cos_gamma_sq(double B) {
    double e3 = abs(E_VS - Zeeman_Energy(B));
    double delta23 = 2.0*r * m_perp * (E_VS) * (Dressel + Rashba) / (sqrt(2) * hbar);
    double e3_tilde = sqrt(e3 * e3 + delta23 * delta23);
    return (sqrt(e3 * e3 + delta23 * delta23) - e3) / (2 * sqrt(e3 * e3 + delta23 * delta23));
}

double pure_transition(double B) {
    double Energy_diff = 1.0;
    if (E_VS >= Zeeman_Energy(B)) {
        Energy_diff = (E_VS - Zeeman_Energy(B));
    }
    else Energy_diff = (Zeeman_Energy(B) - E_VS);
    double pre = pow(Energy_diff, 5) * r * r / (4.0 * pi * rho * pow(hbar, 6)); 
    double omega_z = Zeeman_Energy(B) / hbar;

    double pre_l = pre / pow(v_l, 7);
    double integrals_l = ((4.0 / 3.0) * xi_d * xi_d + (8.0 / 15.0) * xi_d * xi_u + (4.0 / 35.0) * xi_u * xi_u);
    double rate_l = pre_l * (integrals_l);

    double pre_t = pre / pow(v_t, 7);
    double integrals_t = (xi_u * xi_u * (16.0 / 105.0));
    double rate_t = pre_t * (integrals_t);

    return (rate_l + rate_t);

}

int main() {
    double data_points = 120;
    double beginning = 2.0;
    double ending = 8.0;

    ofstream myfile;
    myfile.open("relaxation.txt");
    for (double i = beginning * data_points; i <= ending * data_points; i++) {
        double B = i / data_points;
        if (pure_transition(B) * cos_gamma_sq(B) * sin_gamma_sq(B) > 5.0) {
            myfile << B << " ";
            myfile << pure_transition(B) * cos_gamma_sq(B) * sin_gamma_sq(B) << endl;
        }
    }
    myfile.close();
    return 0;
}