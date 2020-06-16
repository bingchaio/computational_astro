#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fftw++.h>
#include<Array.h>

#define PI M_PI
#define N 4

using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;

double rho(double x, double y, double z){
	return 1/(1+x+y*y+z*z*z);
}

int main(){
	fftw::maxthreads = 1;
    double Lx = 1.0, Ly = 1.0, Lz = 1.0;
    double dx = Lx/N, dy = Ly/N, dz = Lz/N; 
	double M = 0;

	array3<Complex> rho_x(N,N,N,sizeof(Complex));
	array3<Complex> phi_x(N,N,N,sizeof(Complex));
	array3<Complex> phi_k(N,N,N,sizeof(Complex));
	array3<Complex> rho_k(N,N,N,sizeof(Complex));
	fft3d Forward(N, N, N, -1, rho_x, rho_k);
	fft3d Backward(N, N, N, 1, phi_k, phi_x);
	
	for(int i = 1 ; i<N-1 ; i++){
		for(int j = 1 ; j<N-1 ; j++){
			for(int k = 1 ; k<N-1 ; k++){
				if(i==N/2&&j==N/2&&k==N/2) rho_x(i,j,k) = 1; //rho(dx*i, dy*j, dz*k);
				else rho_x(i,j,k) = 0;
				M += real(rho_x(i,j,k));
			}
		}
	}
	M /= (pow(N,3)-pow(N-2,3));
	for(int i = 0 ; i<N ; i++) for(int j = 0 ; j<N ; j++) for(int k = 0 ; k<N ; k++) if(i==0||i==N-1||j==0||j==N-1||k==0||k==N-1) rho_x(i,j,k) = -M;

	// cout << "rho_x = \n";
	// cout << rho_x << endl;	
	
	Forward.fft0(rho_x, rho_k);
	// cout << "rho_k = FFT(rho_x) = \n";
	// cout << rho_k << endl;
   	
    for(int i = 0 ; i<N ; i++){
		for(int j = 0 ; j<N ; j++){
			for(int k = 0 ; k<N ; k++){
            	phi_k(i,j,k) = -dx*dx*rho_k(i,j,k) / ( 4.0 * ( pow(sin(PI*min(i,N-i)/N),2.0) + pow(sin(PI*min(j,N-j)/N),2.0) + pow(sin(PI*min(k,N-k)/N),2.0) ) ) ;
				if(i==0 && j==0 && k==0) phi_k(i,j,k) = 0;
			}
		}
	}
	// cout << "phi_k = \n";
    
	Backward.fft0Normalized(phi_k, phi_x);
	
    // cout << "phi_x = IFFT(phi_k) = \n";

	/* Output results */
    array3<Complex> Lap(N,N,N,sizeof(Complex));
	for(int i = 0 ; i<N ; i++){
        for(int j = 0 ; j<N ; j++){
			for(int k = 0 ; k<N ; k++){
            	Lap(i,j,k) = (phi_x((i+1)%N,j,k)+phi_x((N+i-1)%N,j,k)+phi_x(i,(j+1)%N,k)+phi_x(i,(N+j-1)%N,k)+phi_x(i,j,(k+1)%N)+phi_x(i,j,(N+k-1)%N)-6.0*phi_x(i,j,k))/(dx*dx);
			}
        }
    }
    // cout << "Lap phi_x = \n";
	// cout << Lap << endl;

	cout << "M = " << M << endl;
	
    return 0;
}
