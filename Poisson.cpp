#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fftw++.h>
#include<Array.h>

#define PI M_PI
#define N 6

using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;

double rho(double x, double y){
	return 1/(1+x+y*y);
}

int main(){
	fftw::maxthreads = 1;
    double Lx = 1.0, Ly = 1.0;
    double dx = Lx/(N), dy = Ly/(N); 
	double M = 0;

	array2<Complex> rho_x(N,N,sizeof(Complex));
	array2<Complex> phi_x(N,N,sizeof(Complex));
	array2<Complex> phi_k(N,N,sizeof(Complex));
	array2<Complex> rho_k(N,N,sizeof(Complex));
	fft2d Forward(N, N, -1, rho_x, rho_k);
	fft2d Backward(N, N, 1, phi_k, phi_x);
	
	for(int i = 1 ; i<N-1 ; i++){
		for(int j = 1 ; j<N-1 ; j++){
			if(i==N/2&&j==N/2) rho_x(i,j) = 1;//rho(dx*i, dy*j);
			M += real(rho_x(i,j));
		}
	}
	M /= (4*N-4);
	for(int i = 0 ; i<N ; i++) for(int j = 0 ; j<N ; j++) if(i==0||i==N-1||j==0||j==N-1) rho_x(i,j) = -M;

	cout << "rho_x = \n";
	cout << rho_x << endl;	
	
	Forward.fft0(rho_x, rho_k);
	cout << "rho_k = FFT(rho_x) = \n";
	cout << rho_k << endl;
   	
    for(int i = 0 ; i<N ; i++){
		for(int j = 0 ; j<N ; j++){
            phi_k(i,j) = -dx*dx*rho_k(i,j) / ( 4.0 * ( pow(sin(PI*min(i,N-i)/N),2.0) + pow(sin(PI*min(j,N-j)/N),2.0) ) ) ;
			if(i==0 && j==0) phi_k(i,j) = 0;
		}
	}
	cout << "phi_k = \n";
	cout << phi_k << endl;
    
	Backward.fft0Normalized(phi_k, phi_x);
	
    cout << "phi_x = IFFT(phi_k) = \n";
    cout << phi_x << endl;

	/* Output results */
    array2<Complex> Lap(N,N,sizeof(Complex));
    Lap *= 0;
	for(int i = 0 ; i<N ; i++){
        for(int j = 0 ; j<N ; j++){
            Lap(i,j) = (phi_x((i+1)%N,j)+phi_x((N+i-1)%N,j)+phi_x(i,(j+1)%N)+phi_x(i,(N+j-1)%N)-4.0*phi_x(i,j))/(dx*dx);
        }
    }
    cout << "Lap phi_x = \n";
	cout << Lap << endl;

	cout << "M = " << M << endl;
	
    return 0;
}
