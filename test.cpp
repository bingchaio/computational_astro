#include<iostream>
#include<cstdlib>
#include<cmath>
#include"Eigen/Dense"
#include"unsupported/Eigen/FFT"

int main(){
	int N = 100;
	double dx = 0.1;
	double *Phi_x = new double [N];
	double *Phi_p = new double [N];
	for(int i = 0 ; i<N ; i++){
		double x = i*dx;
		Phi_x[i] = sin(x);
		printf("%.3f \t %e\n", x, Phi_x[i]);
	}
	
	return 0;
}
