#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
using namespace std;

//path-dependent option with Monte-Carlo simulation

double callOption(double S0, double  sigma, double  r, double  T, double  X, int K, int N){
    //N is number of simulations for the Monte-Carlo
    //K is number of path-time-steps
    static mt19937 rng;
    normal_distribution<> ND(0., 1.);
    double sum=0;

    for (int n=0; n<N; n++){
        double dt = T/K;
        vector<double> path(K+1);
        path[0] = S0;
        //create path
        for(int i=1; i<K+1; i++){
            double phi = ND(rng);
            path[i] = path[i-1]*exp((r - 0.5*sigma*sigma)*dt + phi*sigma*sqrt(dt));
        }
        //calculate A
        double A=0.;
        for(int i=1; i<K+1; i++){
            A += fabs(path[i]-path[i-1]);
        }
        sum += max(A-X,0.);
    }

    return exp(-r*T)*sum/N;
}


int main(){
    double  X=150., T=1., S0=104.81, sigma =0.4, r=0.03;
    int K=20, N=100;
    cout<< "Call Option for N=100 : " << callOption(S0,sigma,r,T,X,K,N)<<endl;
    N=200;
    cout<< "Call Option for N=200 : " << callOption(S0,sigma,r,T,X,K,N)<<endl;
    N=500;
    cout<< "Call Option for N=500 : " << callOption(S0,sigma,r,T,X,K,N)<<endl;


}