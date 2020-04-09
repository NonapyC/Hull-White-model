// -------------------------------------------------------------------------------------
//  Hull-white modelを用いた割引債価格の数値計算(モンテカルロ)と解析解との比較
//
//  Comparison between Monte-Carlo simulation and anaritic solution of Zero-Coupen bond price
//  under Hull-white modeling spot rate.
// -------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include "header/nonapy2.hpp"
#include "header/random2.hpp"
#include "header/stopwatch.hpp"
#include "LibMK4/integ.hpp"
using namespace LibMK4;
using namespace LibNS2;
using namespace std;

constexpr int NI = 300000; //イベント回数
constexpr int T_max = 400;
constexpr double R_0 = 0.05; //スポットレート初期値 [/year]
constexpr double dt = 0.1;//時間の刻み幅 [year]
double T_manki = T_max * dt;//満期 [year]

// -------------------------------------------------------------------------------------
// Parameter
// -------------------------------------------------------------------------------------
double Parameter_phi(double t){ //パラメーター
    return 0.005;
}

double Parameter_a(double t){ //パラメーター
    return 0.01;
}

double Parameter_sigma(double t){ //パラメーター
    return 0.01;
}

double R(int T_max){ //Sample Path of spot rate
    double R = R_0;
    for (int t = 0; t < T_max; t++) {
        double rand = rand_normal( 0.0 , 1.0 / sqrt(dt) );
        R += ( Parameter_phi( t * dt ) - Parameter_a( t * dt ) * R ) * dt + rand * Parameter_sigma( t * dt ) * dt;
    }
    return R;
}

double Monte_Carlo_Zero_Coupon_Bonds(){ //割引債価格をモンテカルロシミュレーションで求める
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

double Bond_Price_anal(double t){ //時刻0における債権価格の理論値
    return H_1(t) * exp( - R(t) * H_2(t) );
}

double Hull_White_sample_Path(int T_max){ //スポットレートサンプルパス (グラフ描画用)
    ofstream fs;
    fs.open("test.txt");
    double R = R_0;
    for (int t = 0; t < T_max; t++) {
        double rand = rand_normal( 0.0 , 1.0 / sqrt(dt) );
        R += ( Parameter_phi( t * dt ) - Parameter_a( t * dt ) * R ) * dt + rand * Parameter_sigma( t * dt ) * dt;
        fs <<  t * dt << " " << R << endl;
    }
    fs.close();
    return 0;
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

    FILE *pp;
    pp = popen("python3","w");
    fprintf(pp, "import numpy as np \n");
    fprintf(pp, "import matplotlib.pyplot as plt \n");
    fprintf(pp, "plt.style.use('ggplot') \n");
    fprintf(pp, "data = np.loadtxt('test.txt')\n");
    fprintf(pp, "x1 = data[:,0]\n");
    fprintf(pp, "y1 = data[:,1]\n");
    fprintf(pp, "plt.plot(x1,y1,label='Spot Rate : R(t)') \n");
    fprintf(pp, "plt.xlabel('t') \n");
    fprintf(pp, "plt.minorticks_on() \n");
    fprintf(pp, "plt.legend() \n");
    fprintf(pp, "plt.show() \n");
    pclose(pp);
    return 0;
}
