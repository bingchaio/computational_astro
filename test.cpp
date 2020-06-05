/*
 * ==============
 * PARTICLE MESH
 * ==============
 */

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<ctime>

int main(){
	/* Variables */
	double G = 1.0;
	double Lx = 1.0, Ly = 1.0, Lz = 1.0; // domain size of 3D box
	int Nx = 16, Ny = 16, Nz = 16; // # of grid points
	double dx = Lx/(Nx-1), dy = Ly/(Ny-1), dz = Lz/(Nz-1);
	double t = 0.;
	double dt = 0.1; // time step
	int n = 1000; // # of particles
	double m = 1.0; // particle mass
	double PDx = 0.2, PDy = 0.2, PDz = 0.2;
	double *x = new double [n];
	double *y = new double [n];
	double *z = new double [n];
	double *vx = new double [n];
	double *vy = new double [n];
	double *vz = new double [n];
	double ***U = new double ** [N]; // Dirichlet B.C.
	

	srand(time(NULL));
	/* Initialization */
	for(int i = 0 ; i<n ; i++){
		x[i] = Lx * .5*PDx*(rand()/(double)RAND_MAX-0.5);
		y[i] = Ly * .5*PDy*(rand()/(double)RAND_MAX-0.5);
		z[i] = Lz * .5*PDz*(rand()/(double)RAND_MAX-0.5);
		double r0 = pow(PDx,2) + pow(PDy,2) + pow(PDz,2);
		double v0 = sqrt(G*m/r0);
	}
	for(int i = 0 ; i<Nx ; i++){
		for(int j = 0 ; j<Ny ; j++){
			for(int k = 0 ; k<Nz ; k++){
				U[i][j][k] = 0.;
			}
		}
	}
	

	return 0;
	
}
