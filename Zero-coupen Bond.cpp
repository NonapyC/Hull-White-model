// -------------------------------------------------------------------------------------
//  Comparison between Monte-Carlo simulation and analytical solution of Zero-Coupen bond price
//  under Hull-white model spot rate.
// -------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include "../header/random2.hpp"
#include "../header/stopwatch.hpp"
#include "../LibMK4/integ.hpp"
using namespace LibMK4;
using namespace std;

constexpr int NI = 100000; //# of event
constexpr int T_max = 500;
constexpr double R_0 = 0.05; //initial spot rate [/year]
constexpr double dt = 0.01;//[year]
double T_manki = T_max * dt;//maturity [year]

// -------------------------------------------------------------------------------------
// Parameter
// -------------------------------------------------------------------------------------
double Parameter_phi(double t){ //parametar
    return 0.0001;
}

double Parameter_a(double t){ //parametar
    return 0.1;
}

double Parameter_sigma(double t){ //parametar
    return 0.01;
}

double R(int T_max){ //Sample path of spot rate
    double R = R_0;
    for (int t = 0; t < T_max; t++) {
        double rand = rand_normal( 0.0 , 1.0 / sqrt(dt) );
        R += ( Parameter_phi( t * dt ) - Parameter_a( t * dt ) * R ) * dt + rand * Parameter_sigma( t * dt ) * dt;
    }
    return R;
}

double Monte_Carlo_Zero_Coupon_Bonds(){ //monte-carlo simulation
    double ave = 0.0;
    for (int i = 0; i < NI; i++) {
        double r_sum = 0.0;
        double R = R_0;
        for (int t = 0; t < T_max; t++) {
            double rand = rand_normal( 0.0 , 1.0 / sqrt(dt) );
            R += ( Parameter_phi( t * dt ) - Parameter_a( t * dt ) * R ) * dt + rand * Parameter_sigma( t * dt ) * dt;
            r_sum += R * dt;
        }
        ave += exp( - r_sum );
    } //全イベント終了
    return ave /= NI;
}

// -------------------------------------------------------------------------------------
// Analytic solution
// -------------------------------------------------------------------------------------
double H_2(double s){
    return ( 1.0 - exp( - Parameter_a(s) * ( T_manki - s ) ) ) / Parameter_a(s);
}

double H1_integ(double s){
    return - Parameter_phi(s) * H_2(s) + 0.5 * sqr( Parameter_sigma(s) * H_2(s) );
}

double H_1(double t){
    return exp( qgaus(H1_integ, t, T_manki) );
}

double Bond_Price_anal(double t){
    return H_1(t) * exp( - R(t) * H_2(t) );
}

// -------------------------------------------------------------------------------------
// main関数
// -------------------------------------------------------------------------------------
int main(){
    StopWatch stw;
    stw.Reset();
    MT::init_genrand((unsigned)time(NULL)); /*seed生成*/
    Hull_White_sample_Path(T_max);

    cout << "Zero-Coupen Bond Price (Numerical): B(t=0,T) =  " << Monte_Carlo_Zero_Coupon_Bonds() << endl;
    cout << "Zero-Coupen Bond Price (Analytic) : B(t=0,T) =  " << Bond_Price_anal(0) << endl;

    return 0;
}
