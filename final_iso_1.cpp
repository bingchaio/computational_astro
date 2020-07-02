// ==============
// PARTICLE MESH
// ==============

#include<fftw++.h>
#include "Array.h"
#include <ctime>
#include <omp.h>

#define PI M_PI

using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;

//--------------------------------------------------------mode selection----------------------------------------------
int mesh_mode = 2;  // 0: NGP ; 1: CIC ; 2: TSC
int force_mode = 2; // 0: NGP ; 1: CIC ; 2: TSC
int OI_mode = 0;    //Orbit integration mode. 0: DKD 1:KDK 2:fourth-order symplectic integrator 3:RK4  4:Hermite

//-----------------------------------------------------------constants-------------------------------------------------
double G = 1.0;                                  // gravitational constant
double Lx = 1.0, Ly = 1.0, Lz = 1.0;             // domain size of 3D box
int N = 64;                                      // # of grid points
int Nx = N, Ny = N, Nz = N;
double dx = Lx / (Nx-1), dy = Ly / (Ny-1), dz = Lz / (Nz-1); // spatial resolution
int n = 2;                                       // # of particles
double m = 1.0;                                  // particle mass
double t = 0.0;                                  // time
double PDx = 0.2, PDy = 0.2, PDz = 0.2;          // size of particle clumps
double dt = 0.1*sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2))/sqrt(n*G*m/sqrt(pow(PDx, 2) + pow(PDy, 2) + pow(PDz, 2))); //time steps
double t_end = dt*901.0;                         // ending time                             
double vmax = 1.0;                               // initial maximal velocity weight
double time_elapsed = 0.0;                       // elapsed time
struct timeval start, ending;                    // starting and ending time
const int NThread = 4;                           // number of threads
array3<double>  rho_x(2*Nx,2*Ny,2*Nz,sizeof(Complex)); // rho_x for fft
array3<double>  phi_x(2*Nx,2*Ny,2*Nz,sizeof(Complex)); // phi_x for fft
array3<Complex> phi_k(2*Nx,2*Ny,Nz+1,sizeof(Complex)); // phi_k for fft
array3<Complex> rho_k(2*Nx,2*Ny,Nz+1,sizeof(Complex)); // rho_k for fft
array3<double>  R_x(2*Nx,2*Ny,2*Nz,sizeof(Complex));   // Isolated BC symmetric discrete Green’s function
array3<Complex> R_k(2*Nx,2*Ny,Nz+1,sizeof(Complex));   // FT of Green function

//----------------------------------------------------------functions------------------------------------------------
//Particle Force Interpolation Function
void Get_Force_of_Particle(array3<double> *phi_x, double x, double y, double z, double & F_x, double & F_y, double & F_z, int mode) {
    int X_grid, Y_grid, Z_grid; //grid positions of particles
    if (mode == 0) {

	//grid positions of particles (left grid)
        X_grid = int( x / dx);  
        Y_grid = int( y / dy);
        Z_grid = int( z / dz);

	//choose the nearest grid
        if (abs(x - X_grid * dx) > abs(x - (X_grid + 1) * dx)) X_grid++; 
        if (abs(y - Y_grid * dy) > abs(y - (Y_grid + 1) * dy)) Y_grid++;
        if (abs(z - Z_grid * dz) > abs(z - (Z_grid + 1) * dz)) Z_grid++;

	//exclude the particles at the boundary
        if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx) && (Y_grid<Ny) && (Z_grid<Nz)){

            //calculate the force by using first-order difference of potential
	    
            F_x = -( - (*phi_x)((X_grid + Nx - 1)%Nx,Y_grid,Z_grid) * 0.5 + (*phi_x)((X_grid + 1)%Nx,Y_grid,Z_grid) * 0.5 )/dx; 
            F_y = -( - (*phi_x)(X_grid,(Y_grid + Ny - 1)%Ny,Z_grid) * 0.5 + (*phi_x)(X_grid,(Y_grid + 1)%Ny,Z_grid) * 0.5 )/dy;
            F_z = -( - (*phi_x)(X_grid,Y_grid,(Z_grid + Nz - 1)%Nz) * 0.5 + (*phi_x)(X_grid,Y_grid,(Z_grid + 1)%Nz) * 0.5 )/dz;
            
	    //calculate the force by using second-order difference of potential
	    /*	    
            F_x = -( (*phi_x)((X_grid + Nx - 2)%Nx,Y_grid,Z_grid) / 12. - (*phi_x)((X_grid + Nx - 1)%Nx,Y_grid,Z_grid) * (2. / 3.) + 
                     (*phi_x)((X_grid + 1)%Nx,Y_grid,Z_grid) * (2. / 3.) - (*phi_x)((X_grid + 2)%Nx,Y_grid,Z_grid) * (1. / 12.) )/dx;
            F_y = -( (*phi_x)(X_grid,(Y_grid - 2 + Ny)%Ny,Z_grid) / 12. - (*phi_x)(X_grid,(Y_grid + Ny - 1)%Ny,Z_grid) * (2. / 3.) +
                     (*phi_x)(X_grid,(Y_grid + 1)%Ny,Z_grid) * (2. / 3.) - (*phi_x)(X_grid,(Y_grid + 2)%Ny,Z_grid) * (1. / 12.) )/dy;
            F_z = -( (*phi_x)(X_grid,Y_grid,(Z_grid + Nz - 2)%Nz) / 12. - (*phi_x)(X_grid,Y_grid,(Z_grid + Nz - 1)%Nz) * (2. / 3.) +
                     (*phi_x)(X_grid,Y_grid,(Z_grid + 1)%Nz) * (2. / 3.) - (*phi_x)(X_grid,Y_grid,(Z_grid + 2)%Nz) * (1. / 12.) )/dz;
	    */
        }
    } else if (mode == 1) {

        double f; //the weigting factor

	//grid positions of particles (left grid)
        X_grid = int( x/ dx);
        Y_grid = int( y/ dy);
        Z_grid = int( z/ dz);

	//exclude the particles at the boundary
	if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
            for (int i = X_grid; i <= X_grid + 1; i++) {
                for (int j = Y_grid; j <= Y_grid + 1; j++) {
                    for (int k = Z_grid; k <= Z_grid + 1; k++) {

			//calculate the weigting factor by CIC
                        f = (1.0 - abs(x - i * dx) / dx) * (1.0 - abs(y - j * dy) / dy) * (1.0 - abs(z - k * dz) / dz);
			
			//calculate the force by using first-order difference of potential
			
			F_x -= f * (-(*phi_x)((i + Nx - 1)%Nx,j,k) * 0.5 + (*phi_x)((i + 1)%Nx,j,k) * 0.5)/dx;
                        F_y -= f * (-(*phi_x)(i,(j + Ny - 1)%Ny,k) * 0.5 + (*phi_x)(i,(j + 1)%Ny,k) * 0.5)/dy;
                        F_z -= f * (-(*phi_x)(i,j,(k + Nz - 1)%Nz) * 0.5 + (*phi_x)(i,j,(k + 1)%Nz) * 0.5)/dz;
			
			//calculate the force by using second-order difference of potential
			/*
                        F_x -= f * ((*phi_x)((i + Nx - 2)%Nx,j,k) / 12. - (*phi_x)((i + Nx - 1)%Nx,j,k) * (2. / 3.) +
                                    (*phi_x)((i + 1)%Nx,j,k) * (2. / 3.) - (*phi_x)((i + 2)%Nx,j,k) * (1. / 12.))/dx;
                        F_y -= f * ((*phi_x)(i,(j + Ny - 2)%Ny,k) / 12. - (*phi_x)(i,(j + Ny - 1)%Ny,k) * (2. / 3.) +
                                    (*phi_x)(i,(j + 1)%Ny,k) * (2. / 3.) - (*phi_x)(i,(j + 2)%Ny,k) * (1. / 12.))/dy;
                        F_z -= f * ((*phi_x)(i,j,(k + Nz - 2)%Nz) / 12. - (*phi_x)(i,j,(k + Nz - 1)%Nz) * (2. / 3.) +
                                    (*phi_x)(i,j,(k + 1)%Nz) * (2. / 3.) - (*phi_x)(i,j,(k + 2)%Nz) * (1. / 12.))/dz;
			*/
                    }
                }
	    }
        }
    } else if (mode == 2){

        double fx, fy, fz, f; //the weigting factor

	//grid positions of particles (left grid)
        X_grid = int( x / dx);
        Y_grid = int( y / dy);
        Z_grid = int( z / dz);

	//exclude the particles at the boundary
	if ((X_grid>0) && (Y_grid>0) && (Z_grid>0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
            for (int i = X_grid - 1; i <= X_grid + 1; i++) {
                if (i == X_grid) fx = 0.75 - pow(x - i * dx, 2) / pow(dx, 2); //calculate the weigting factor by TSC
                else fx = 0.5 * pow(1.5 - abs(x - i * dx) / dx, 2);
                for (int j = Y_grid - 1; j <= Y_grid + 1; j++) {
                    if (j == Y_grid) fy = 0.75 - pow(y - j * dy, 2) / pow(dy, 2); //calculate the weigting factor by TSC
                    else fy = 0.5 * pow(1.5 - abs(y - j * dy) / dy, 2);
                    for (int k = Z_grid - 1; k <= Z_grid + 1; k++) {
                        if (k == Z_grid) fz = 0.75 - pow(z - k * dz, 2) / pow(dz, 2); //calculate the weigting factor by TSC
                        else fz = 0.5 * pow(1.5 - abs(z - k * dz) / dz, 2);
                        f = fx * fy * fz;

			//calculate the force by using first-order difference of potential
                        F_x -= f * (-(*phi_x)((i + Nx - 1)%Nx,j,k) * 0.5 + (*phi_x)((i + 1)%Nx,j,k) * 0.5)/dx;
                        F_y -= f * (-(*phi_x)(i,(j + Ny - 1)%Ny,k) * 0.5 + (*phi_x)(i,(j + 1)%Ny,k) * 0.5)/dy;
                        F_z -= f * (-(*phi_x)(i,j,(k + Nz - 1)%Nz) * 0.5 + (*phi_x)(i,j,(k + 1)%Nz) * 0.5)/dz;
			
			//calculate the force by using second-order difference of potential
			/*
                        F_x -= f * ((*phi_x)((i + Nx - 2)%Nx,j,k) / 12. - (*phi_x)((i + Nx - 1)%Nx,j,k) * (2. / 3.) +
                                    (*phi_x)((i + 1)%Nx,j,k) * (2. / 3.) - (*phi_x)((i + 2)%Nx,j,k) * (1. / 12.))/dx;
                        F_y -= f * ((*phi_x)(i,(j + Ny - 2)%Ny,k) / 12. - (*phi_x)(i,(j + Ny - 1)%Ny,k) * (2. / 3.) +
                                    (*phi_x)(i,(j + 1)%Ny,k) * (2. / 3.) - (*phi_x)(i,(j + 2)%Ny,k) * (1. / 12.))/dy;
                        F_z -= f * ((*phi_x)(i,j,(k + Nz - 2)%Nz) / 12. - (*phi_x)(i,j,(k + Nz - 1)%Nz) * (2. / 3.) +
                                    (*phi_x)(i,j,(k + 1)%Nz) * (2. / 3.) - (*phi_x)(i,j,(k + 2)%Nz) * (1. / 12.))/dz;
			*/
			
                	}
		}
            }
        }
    }
}

//Poisson Solver (FFT)
void FFT(array3<double> *rho_x,array3<double> *phi_x,array3<Complex> *R_k){
    //fftw::maxthreads = get_max_threads();

    gettimeofday(&start, NULL);

    rcfft3d Forward(2*Nx, 2*Ny, 2*Nz, (*rho_x), rho_k);
    crfft3d Backward(2*Nx, 2*Ny, 2*Nz, phi_k, (*phi_x));

    	// fourier transform
    	
    	Forward.fft((*rho_x), rho_k);
	
    	// calculate the potential in k space
    
    	for(int i = 0 ; i<2*Nx ; i++){
		for(int j = 0 ; j<2*Ny ; j++){
	   		for(int k = 0 ; k<Nz+1 ; k++){
				phi_k(i,j,k) = (*R_k)(i,j,k)*rho_k(i,j,k);	
			}
      		}
    	}
        
        // inverse fourier transform
    	
    	Backward.fftNormalized(phi_k, (*phi_x));

	for(int i = 0 ; i<Nx ; i++){
		for(int j = 0 ; j<Ny ; j++){
	   		for(int k = 0 ; k<Nz ; k++){
				(*phi_x)(i,j,k) = (*phi_x)(i,j,k)/Nx/Ny/Nz;	
			}
      		}
    	}

    gettimeofday(&ending, NULL); 
    float delta = ((ending.tv_sec  - start.tv_sec) * 1000000u + ending.tv_usec - start.tv_usec) / 1.e6;
    time_elapsed += 1.0*(delta);
}

//Particle Mesh function
void mesh(array3<double> *rho_x, double *x, double *y, double *z, int mode) {

    int X_grid, Y_grid, Z_grid; //grid positions of particles

    if (mode == 0) {
        for (int p = 0; p < n; p++) {

	    //grid positions of particles (left grid)
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dy);
            Z_grid = int((z[p]) / dz); 

	    //choose the nearest grid
            if (abs(x[p] - X_grid * dx) > abs(x[p] - (X_grid + 1) * dx)) X_grid++;
            if (abs(y[p] - Y_grid * dy) > abs(y[p] - (Y_grid + 1) * dy)) Y_grid++;
            if (abs(z[p] - Z_grid * dz) > abs(z[p] - (Z_grid + 1) * dz)) Z_grid++;

	    // set the density by NGP
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx) && (Y_grid<Ny) && (Z_grid<Nz)) (*rho_x)(X_grid,Y_grid,Z_grid) += m / (dx * dy * dz);
        }
    } else if (mode == 1) {
        for (int p = 0; p < n; p++) {

	    //grid positions of particles (left grid)
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dy);
            Z_grid = int((z[p]) / dz);
	    
	    //exclude the particles at the boundary
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
                for (int i = X_grid; i <= X_grid + 1; i++) {
                    for (int j = Y_grid; j <= Y_grid + 1; j++) {
                        for (int k = Z_grid; k <= Z_grid + 1; k++) {
			     // set the density by CIC
                             (*rho_x)(i,j,k) += m * (1.0 - abs(x[p] - i * dx) / dx) * (1.0 - abs(y[p] - j * dy) / dy) * (1.0 - abs(z[p] - k * dz) / dz) / (dx * dy * dz);
                        }
                    }
                }
            }
        }
    } else if (mode == 2){

        double fx, fy, fz; //the weigting factor

        for (int p = 0; p < n; p++) {

	    //grid positions of particles (left grid)
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dy);
            Z_grid = int((z[p]) / dz);

	    //exclude the particles at the boundary
            if ((X_grid>0) && (Y_grid>0) && (Z_grid>0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
                for (int i = X_grid - 1; i <= X_grid + 1; i++) {
                    if (i == X_grid) fx = 0.75 - pow(x[p] - i * dx, 2) / pow(dx, 2); //calculate the weigting factor by TSC
                    else fx = 0.5 * pow(1.5 - abs(x[p] - i * dx) / dx, 2);
                    for (int j = Y_grid - 1; j <= Y_grid + 1; j++) {
                        if (j == Y_grid) fy = 0.75 - pow(y[p] - j * dy, 2) / pow(dy, 2); //calculate the weigting factor by TSC
                        else fy = 0.5 * pow(1.5 - abs(y[p] - j * dy) / dy, 2);
                        for (int k = Z_grid - 1; k <= Z_grid + 1; k++) {
                            if (k == Z_grid) fz = 0.75 - pow(z[p] - k * dz, 2) / pow(dz, 2); //calculate the weigting factor by TSC
                            else fz = 0.5 * pow(1.5 - abs(z[p] - k * dz) / dz, 2);

			    // set the density by TSC
                            (*rho_x)(i,j,k) += m * fx * fy * fz / (dx * dy * dz);
                        }
              	    }
            	}
            }
        }
    }
}

void i_density(array3<double> *rho_x, double *x, double *y, double *z, int mode) {
 
    int X_grid, Y_grid, Z_grid; //grid positions of particles
    if (mode == 0) {
        for (int p = 0; p < n; p++) {

	    //grid positions of particles (left grid)
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dy);
            Z_grid = int((z[p]) / dz); 

	    //choose the nearest grid
            if (abs(x[p] - X_grid * dx) > abs(x[p] - (X_grid + 1) * dx)) X_grid++;
            if (abs(y[p] - Y_grid * dy) > abs(y[p] - (Y_grid + 1) * dy)) Y_grid++;
            if (abs(z[p] - Z_grid * dz) > abs(z[p] - (Z_grid + 1) * dz)) Z_grid++;

	    // set the density by NGP
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx) && (Y_grid<Ny) && (Z_grid<Nz)) (*rho_x)(X_grid,Y_grid,Z_grid) = 0.0;
        }
    } else if (mode == 1) {
        for (int p = 0; p < n; p++) {

	    //grid positions of particles (left grid)
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dy);
            Z_grid = int((z[p]) / dz);
	    
	    //exclude the particles at the boundary
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
                for (int i = X_grid; i <= X_grid + 1; i++) {
                    for (int j = Y_grid; j <= Y_grid + 1; j++) {
                        for (int k = Z_grid; k <= Z_grid + 1; k++) {
                             (*rho_x)(i,j,k) = 0.0;
                        }
                    }
                }
            }
        }
    } else if (mode == 2){

        double fx, fy, fz; //the weigting factor

        for (int p = 0; p < n; p++) {

	    //grid positions of particles (left grid)
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dy);
            Z_grid = int((z[p]) / dz);

	    //exclude the particles at the boundary
            if ((X_grid>0) && (Y_grid>0) && (Z_grid>0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
                for (int i = X_grid - 1; i <= X_grid + 1; i++) {
                    for (int j = Y_grid - 1; j <= Y_grid + 1; j++) {
                        for (int k = Z_grid - 1; k <= Z_grid + 1; k++) {
                            (*rho_x)(i,j,k) = 0.0;
                        }
              	    }
            	}
            }
        }
    }
}

double Get_Kenergy(double *x, double *y, double *z, double *vx, double *vy, double *vz){
    double KE = 0.0;
    for(int p = 0 ; p<n ; p++){
        KE += 0.5 * m * (vx[p]*vx[p]+vy[p]*vy[p]+vz[p]*vz[p]); //kinetic energy
    }
    return KE;
}

double Get_Uenergy(double *x, double *y, double *z, double *vx, double *vy, double *vz){
    double UE = 0.0;
    for(int p = 0 ; p<n ; p++){
        for(int q = p+1 ; q<n ; q++){
            UE -= G*m*m/sqrt(pow(x[p]-x[q],2)+pow(y[p]-y[q],2)+pow(z[p]-z[q],2)); //potential energy
        }
    }
    return UE;
}

int main() {
    /* Variables */
    double * x = new double[n];          //positions of the particles
    double * y = new double[n];
    double * z = new double[n];
    double * vx = new double[n];         //velocities of the particles
    double * vy = new double[n];
    double * vz = new double[n];
    double * F_x = new double[n];        // forces of the particles
    double * F_y = new double[n];
    double * F_z = new double[n];
    
    srand(time(NULL));
    /* Initialization */
    //Random distribution
    double r0 = pow(pow(PDx, 2) + pow(PDy, 2) + pow(PDz, 2),0.5); //mean distance
    double v0 = vmax*sqrt(2* G * n * m / r0 );                    //Virial speed
    for (int i = 0; i < n; i++) {
        x[i] = PDx * (rand() / (double) RAND_MAX -0.5) + Lx/2;
        y[i] = PDy * (rand() / (double) RAND_MAX -0.5) + Ly/2;
        z[i] = PDz * (rand() / (double) RAND_MAX -0.5) + Lz/2;
        vx[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) *2.0;
        vy[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) *2.0;
        vz[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) *2.0;
    }
    
    x[0] = 0.6;
    y[0] = 0.5;
    z[0] = 0.5;
    x[1] = 0.4;
    y[1] = 0.5;
    z[1] = 0.5;
    vx[0] = 0.0;
    vy[0] = 0.5*sqrt(1.0/0.2);
    vz[0] = 0.0;
    vx[1] = 0.0;
    vy[1] = -0.5*sqrt(1.0/0.2);
    vz[1] = 0.0;
    
    printf("isolated N = %d mesh mode = %d orbit mode = %d NThread = %d dt = %.3e r0 = 0.2\n",N,mesh_mode,OI_mode,NThread,dt);

    //initialize R
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {	    
            for (int k = 0; k < Nz; k++) {
		rho_x(i,j,k) = 0.0;
		rho_x(2*Nx-i-1,j,k) = 0.0; 
		rho_x(i,2*Ny-j-1,k) = 0.0;
		rho_x(i,j,2*Nz-k-1) = 0.0;
		rho_x(2*Nx-i-1,2*Ny-j-1,k) = 0.0;
		rho_x(2*Nx-i-1,j,2*Nz-k-1) = 0.0;
		rho_x(i,2*Ny-j-1,2*Nz-k-1) = 0.0;
		rho_x(2*Nx-i-1,2*Ny-j-1,2*Nz-k-1) = 0.0;
		R_x(i,j,k) = -G/pow((pow(i*dx,2)+pow(j*dy,2)+pow(k*dz,2)),0.5);
		R_x(2*Nx-i-1,j,k) = -G/pow((pow((i+1)*dx,2)+pow(j*dy,2)+pow(k*dz,2)),0.5); 
		R_x(i,2*Ny-j-1,k) = -G/pow((pow(i*dx,2)+pow((j+1)*dy,2)+pow(k*dz,2)),0.5);
		R_x(i,j,2*Nz-k-1) = -G/pow((pow(i*dx,2)+pow(j*dy,2)+pow((k+1)*dz,2)),0.5);
		R_x(2*Nx-i-1,2*Ny-j-1,k) = -G/pow((pow((i+1)*dx,2)+pow((j+1)*dy,2)+pow(k*dz,2)),0.5);
		R_x(2*Nx-i-1,j,2*Nz-k-1) = -G/pow((pow((i+1)*dx,2)+pow(j*dy,2)+pow((k+1)*dz,2)),0.5);
		R_x(i,2*Ny-j-1,2*Nz-k-1) = -G/pow((pow(i*dx,2)+pow((j+1)*dy,2)+pow((k+1)*dz,2)),0.5);
		R_x(2*Nx-i-1,2*Ny-j-1,2*Nz-k-1) = -G/pow((pow((i+1)*dx,2)+pow((j+1)*dy,2)+pow((k+1)*dz,2)),0.5);
	    }
	}
    }	

    R_x(0,0,0) = 0.0;
    R_x(Nx,Ny,Nz) = 0.0;
    //FT of Green function
    rcfft3d Forward(2*Nx, 2*Ny, 2*Nz, R_x, R_k);
    Forward.fft(R_x, R_k);

    mesh(&rho_x, x, y, z, mesh_mode);
    FFT(&rho_x,&phi_x,&R_k);
    i_density(&rho_x, x, y, z, mesh_mode);
    for (int i = 0; i < n; i++) {
        F_x[i] = 0.0; 
	F_y[i] = 0.0; 
	F_z[i] = 0.0;
        Get_Force_of_Particle(&phi_x, x[i], y[i], z[i], F_x[i], F_y[i], F_z[i], force_mode);
    }

    while (t <= t_end) {

        // check conservation
        double Px = 0, Py = 0, Pz = 0;
        double X = 0, Y = 0, Z = 0;
        double M = 0;
        int n_in = 0;
	
        if((int)(t/dt)%(100)==0){
	    FILE *den_output;
            char fname[10];
            int t_out = (t/dt);
            sprintf(fname,"density_%04d", (t_out));
            den_output = fopen(fname,"w");
            cout << "================================\n";
            for(int p = 0 ; p<n ; p++){
                Px += m*vx[p];
                Py += m*vy[p];
                Pz += m*vz[p];
                if(x[p]>0&&x[p]<Lx&&y[p]>0&&y[p]<Ly&&z[p]>0&&z[p]<Lz) n_in++;
            }
	    mesh(&rho_x, x, y, z, mesh_mode);
    	    FFT(&rho_x,&phi_x,&R_k);
            for(int i = 0 ; i<Nx ; i++) for(int j = 0 ; j<Ny ; j++) for(int k = 0 ; k<Nz ; k++) M += rho_x(i,j,k)*dx*dy*dz;
            i_density(&rho_x, x, y, z, mesh_mode);
            printf("t = %.3f\n", t);
            printf("Px = %.12f \t Py = %.12f \t Pz = %.12f\tphi(0.5,0.5,0.5) = %.3f\n", Px, Py, Pz, phi_x(Nx/2,Ny/2,Nz/2));
            printf("n_in = %d\tM = %.3f\tKE = %.3e\tUE = %.3e\tE = %.3e\tt=%.6f\n", n_in, M,Get_Kenergy(x,y,z,vx,vy,vz),Get_Uenergy(x,y,z,vx,vy,vz),Get_Kenergy(x,y,z,vx,vy,vz)+Get_Uenergy(x,y,z,vx,vy,vz),time_elapsed);
	    for (int i = 0; i < n; i++) fprintf (den_output, "%g  %g  %g   \n",x[i], y[i], z[i] );
            fclose(den_output);
        }    
     
        //DKD
        if (OI_mode == 0) {

            //drift: update position by 0.5*dt
	    
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * 0.5 * dt;
                  y[i] += vy[i] * 0.5 * dt;
                  z[i] += vz[i] * 0.5 * dt;
            }

            //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
	    
            mesh(&rho_x, x, y, z, mesh_mode);
    	    FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);

	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x[i], y[i], z[i], F_x[i], F_y[i], F_z[i], force_mode);
                  vx[i] += F_x[i] / m * dt;
                  vy[i] += F_y[i] / m * dt;
                  vz[i] += F_z[i] / m * dt;
            }

            //drift: update position by 0.5*dt
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * 0.5 * dt;
                  y[i] += vy[i] * 0.5 * dt;
                  z[i] += vz[i] * 0.5 * dt;
            }    
        }
        
        //KDK
        else if (OI_mode == 1) { 
	    
	    //kick	    
	    for (int i = 0; i < n; i++) {
                  vx[i] += F_x[i] / m * 0.5 * dt;
                  vy[i] += F_y[i] / m * 0.5 * dt;
                  vz[i] += F_z[i] / m * 0.5 * dt;
            }

            //drift: update position by dt
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * dt;
                  y[i] += vy[i] * dt;
                  z[i] += vz[i] * dt;
            }

            //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
            mesh(&rho_x, x, y, z, mesh_mode);
    	    FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);

	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x[i], y[i], z[i], F_x[i], F_y[i], F_z[i], force_mode);
                  vx[i] += F_x[i] / m * 0.5 * dt;
                  vy[i] += F_y[i] / m * 0.5 * dt;
                  vz[i] += F_z[i] / m * 0.5 * dt;
            }
        }
        
        //fourth-order symplectic integration
        else if (OI_mode == 2) {
            //fourth-order symplectic coefficients
            double w1 = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));
            double w0 = 1.0 - 2.0 * w1;
            double c1 = w1 / 2.0;
            double c2 = (w1 + w0) / 2.0;
            double c3 = c2;
            double c4 = c1;
            double d1 = w1;
            double d2 = w0;
            double d3 = w1;
            
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * c1 * dt;
                  y[i] += vy[i] * c1 * dt;
                  z[i] += vz[i] * c1 * dt;
            }
            
            mesh(&rho_x, x, y, z, mesh_mode);
    	    FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);

	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x[i], y[i], z[i], F_x[i], F_y[i], F_z[i], force_mode);
                  vx[i] += F_x[i] / m * d1 * dt;
                  vy[i] += F_y[i] / m * d1 * dt;
                  vz[i] += F_z[i] / m * d1 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * c2 * dt;
                  y[i] += vy[i] * c2 * dt;
                  z[i] += vz[i] * c2 * dt;
            }
            
            mesh(&rho_x, x, y, z, mesh_mode);
    	    FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);

	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x[i], y[i], z[i], F_x[i], F_y[i], F_z[i], force_mode);
                  vx[i] += F_x[i] / m * d2 * dt;
                  vy[i] += F_y[i] / m * d2 * dt;
                  vz[i] += F_z[i] / m * d2 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * c3 * dt;
                  y[i] += vy[i] * c3 * dt;
                  z[i] += vz[i] * c3 * dt;
            }
            
            mesh(&rho_x, x, y, z, mesh_mode);
    	    FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);

	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x[i], y[i], z[i], F_x[i], F_y[i], F_z[i], force_mode);
                  vx[i] += F_x[i] / m * d3 * dt;
                  vy[i] += F_y[i] / m * d3 * dt;
                  vz[i] += F_z[i] / m * d3 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                  x[i] += vx[i] * c4 * dt;
                  y[i] += vy[i] * c4 * dt;
                  z[i] += vz[i] * c4 * dt;
            }
        }
        
        //RK4 mode
        else if (OI_mode == 3) {
            double ** kr1  = new double * [n];
            double ** kv1  = new double * [n];
            double ** kr2  = new double * [n];
            double ** kv2  = new double * [n];
            double ** kr3  = new double * [n];
            double ** kv3  = new double * [n];
            double ** kr4  = new double * [n];
            double ** kv4  = new double * [n];
            double * x_tmp = new double[n]; //positions of the particles
            double * y_tmp = new double[n];
            double * z_tmp = new double[n];
            for (int i = 0; i < n; i++) {
                  kr1[i] = new double[3];
                  kv1[i] = new double[3];
                  kr2[i] = new double[3];
                  kv2[i] = new double[3];
                  kr3[i] = new double[3];
                  kv3[i] = new double[3];
                  kr4[i] = new double[3];
                  kv4[i] = new double[3];
                  x_tmp[i] = x[i];
                  y_tmp[i] = y[i];
                  z_tmp[i] = z[i];
            }
	    mesh(&rho_x, x_tmp, y_tmp, z_tmp, mesh_mode);
            FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);
	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
               	  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x_tmp[i], y_tmp[i], z_tmp[i], F_x[i], F_y[i], F_z[i], force_mode);
                  kr1[i][0] = vx[i];
                  kr1[i][1] = vy[i];
                  kr1[i][2] = vz[i];
                  kv1[i][0] = F_x[i]/m;
                  kv1[i][1] = F_y[i]/m;
               	  kv1[i][2] = F_z[i]/m;
                  kr2[i][0] = kr1[i][0] + kv1[i][0] * 0.5 * dt;
                  kr2[i][1] = kr1[i][1] + kv1[i][1] * 0.5 * dt;
                  kr2[i][2] = kr1[i][2] + kv1[i][2] * 0.5 * dt;
                  x_tmp[i] += kr1[i][0]*0.5*dt;
                  y_tmp[i] += kr1[i][1]*0.5*dt;
                  z_tmp[i] += kr1[i][2]*0.5*dt;
            }
            mesh(&rho_x, x_tmp, y_tmp, z_tmp, mesh_mode);
            FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);
	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x_tmp[i], y_tmp[i], z_tmp[i], F_x[i], F_y[i], F_z[i], force_mode);
                  kv2[i][0] = F_x[i]/m;
                  kv2[i][1] = F_y[i]/m;
                  kv2[i][2] = F_z[i]/m;
                  kr3[i][0] = kr1[i][0] + kv2[i][0] * 0.5 * dt;
                  kr3[i][1] = kr1[i][1] + kv2[i][1] * 0.5 * dt;
                  kr3[i][2] = kr1[i][2] + kv2[i][2] * 0.5 * dt;
                  x_tmp[i]  = x[i] + kr2[i][0]*0.5*dt;
                  y_tmp[i]  = y[i] + kr2[i][1]*0.5*dt;
                  z_tmp[i]  = z[i] + kr2[i][2]*0.5*dt;
            }
	    mesh(&rho_x, x_tmp, y_tmp, z_tmp, mesh_mode);
            FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);
	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x_tmp[i], y_tmp[i], z_tmp[i], F_x[i], F_y[i], F_z[i], force_mode);
                  kv3[i][0] = F_x[i]/m;
                  kv3[i][1] = F_y[i]/m;
                  kv3[i][2] = F_z[i]/m;
                  kr4[i][0] = kr1[i][0] + kv3[i][0] * dt;
                  kr4[i][1] = kr1[i][1] + kv3[i][1] * dt;
                  kr4[i][2] = kr1[i][2] + kv3[i][2] * dt;
                  x_tmp[i]  = x[i] + kr3[i][0]*dt;
                  y_tmp[i]  = y[i] + kr3[i][1]*dt;
                  z_tmp[i]  = z[i] + kr3[i][2]*dt;
            }
	    mesh(&rho_x, x_tmp, y_tmp, z_tmp, mesh_mode);
            FFT(&rho_x,&phi_x,&R_k);
	    i_density(&rho_x, x, y, z, mesh_mode);
	    omp_set_num_threads( NThread );
    	    #  pragma omp parallel for
            for (int i = 0; i < n; i++) {
                  F_x[i] = 0.0; 
		  F_y[i] = 0.0; 
		  F_z[i] = 0.0;
                  Get_Force_of_Particle(&phi_x, x_tmp[i], y_tmp[i], z_tmp[i], F_x[i], F_y[i], F_z[i], force_mode);
                  kv4[i][0] = F_x[i]/m;
                  kv4[i][1] = F_y[i]/m;
                  kv4[i][2] = F_z[i]/m;
            }
            for (int i = 0; i < n; i++) {
                  x[i] += (kr1[i][0] + 2.0 * kr2[i][0] + 2.0 * kr3[i][0] + 1.0 * kr4[i][0]) / 6.0 * dt;
                  y[i] += (kr1[i][1] + 2.0 * kr2[i][1] + 2.0 * kr3[i][1] + 1.0 * kr4[i][1]) / 6.0 * dt;
                  z[i] += (kr1[i][2] + 2.0 * kr2[i][2] + 2.0 * kr3[i][2] + 1.0 * kr4[i][2]) / 6.0 * dt;
                  vx[i] += (kv1[i][0] + 2.0 * kv2[i][0] + 2.0 * kv3[i][0] + 1.0 * kv4[i][0]) / 6.0 * dt;
                  vy[i] += (kv1[i][1] + 2.0 * kv2[i][1] + 2.0 * kv3[i][1] + 1.0 * kv4[i][1]) / 6.0 * dt;
                  vz[i] += (kv1[i][2] + 2.0 * kv2[i][2] + 2.0 * kv3[i][2] + 1.0 * kv4[i][2]) / 6.0 * dt;
            }
        }
        
	t += dt;

    }
    
    return EXIT_SUCCESS;
}



