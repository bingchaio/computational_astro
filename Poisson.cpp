#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fftw++.h>
#include<Array.h>

#define PI M_PI
#define N 9

using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;

double rho(double x, double y){
    return sin(2*PI*x)+sin(2*PI*y);
}

int main(){
	fftw::maxthreads = 1;
    double Lx = 1.0, Ly = 1.0;
    double dx = Lx/(N), dy = Ly/(N);
    
	array2<Complex> rho_x(N,N,sizeof(Complex));
	array2<Complex> phi_x(N,N,sizeof(Complex));
	array2<Complex> phi_k(N,N,sizeof(Complex));
	array2<Complex> rho_k(N,N,sizeof(Complex));
	for(int i = 0 ; i<N ; i++){
		for(int j = 0 ; j<N ; j++){
			// rho_x(i,j) = rho(dx*i, dy*j);
            if(i==0 && j==0) rho_x(i,j) = 3;
            //else if(i==5 && j==6) rho_x(i,j) = 2;
            //if(i==N/2 && j==N/2) rho_x(i,j) = 3;
            else rho_x(i,j) = 0;
        }
	}

	cout << "rho_x = \n" << rho_x << endl;
	
	fft2d Forward(N, N, -1, rho_x, rho_k);
	fft2d Backward(N, N, 1, phi_k, phi_x);
	
	Forward.fft0(rho_x, rho_k);
	cout << "rho_k = \n" << rho_k << endl;
    
    for(int i = 0 ; i<N ; i++){
		for(int j = 0 ; j<N ; j++){
            phi_k(i,j) = -dx*dx*rho_k(i,j) / ( 4.0 * ( pow(sin(PI*i/N),2.0) + pow(sin(PI*j/N),2.0) ) );
            if(i==0 && j==0) phi_k(i,j) = 0;
		}
	}
    
	Backward.fft0Normalized(phi_k, phi_x);
	
    cout << "phi_x = \n";
    for(int i = 0 ; i<N ; i++){
        for(int j = 0 ; j<N ; j++){
            cout << real(phi_x(i,j)) << "\t";
        }
        cout << endl;
    }
    cout << endl;
    
	/* Output results */
    array2<double> Lap(N,N,sizeof(Complex));
    Lap *= 0;
	for(int i = 0 ; i<N ; i++){
        for(int j = 0 ; j<N ; j++){
            Lap(i,j) = real(phi_x((i+1)%N,j)+phi_x((N+i-1)%N,j)+phi_x(i,(j+1)%N)+phi_x(i,(N+j-1)%N)-4.0*phi_x(i,j))/(dx*dx);
        }
    }
    cout << "g_x = \n" << Lap << endl;
    array2<double> Fx(N,N,sizeof(Complex));
    Lap *= 0;
    for(int i = 0 ; i<N ; i++){
        for(int j = 0 ; j<N ; j++){
            Fx(i,j) = real(phi_x((i+1)%N,j)-phi_x((N+i-1)%N,j))/dx;
        }
    }
    cout << "Lap phi_x = \n" << Fx << endl;
    
    return 0;
}
