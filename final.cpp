// ==============
// PARTICLE MESH
// ==============

#include<fftw++.h>
#include "Array.h"

#define PI M_PI

using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;

//--------------------------------------------------------mode selection----------------------------------------------
int mesh_mode = 2; // 0: NGP ; 1: CIC ; 2: TSC
int force_mode = 2; // 0: NGP ; 1: CIC ; 2: TSC
int OI_mode = 2; //Orbit integration mode. 0: DKD 1:KDK 2:fourth-order symplectic integrator 3:RK4  4:Hermite

//-----------------------------------------------------------constants-------------------------------------------------
double G = 1.0; // gravitational constant
double Lx = 1.0, Ly = 1.0, Lz = 1.0; // domain size of 3D box
int N = 128; // # of grid points
int Nx = N, Ny = N, Nz = N;
double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz; // spatial resolution
int n = 1000; // # of particles
double m = 1.0; // particle mass

//----------------------------------------------------------functions------------------------------------------------
//Particle Force Interpolation Function
void Get_Force_of_Particle(double *** U, double x, double y, double z, double & F_x, double & F_y, double & F_z, int mode) {
	int X_grid, Y_grid, Z_grid;
    if (mode == 0) {
        X_grid = int( x / dx);
        Y_grid = int( y / dy);
        Z_grid = int( z / dz);
        if (abs(x - X_grid * dx) > abs(x - (X_grid + 1) * dx)) X_grid++;
        if (abs(y - Y_grid * dy) > abs(y - (Y_grid + 1) * dy)) Y_grid++;
        if (abs(z - Z_grid * dz) > abs(z - (Z_grid + 1) * dz)) Z_grid++;
        if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx) && (Y_grid<Ny) && (Z_grid<Nz)){
            
            F_x = -( - U[(X_grid + Nx - 1)%Nx][Y_grid][Z_grid] * 0.5 + U[(X_grid + 1)%Nx][Y_grid][Z_grid] * 0.5 );
            F_y = -( - U[X_grid][(Y_grid + Ny - 1)%Ny][Z_grid] * 0.5 + U[X_grid][(Y_grid + 1)%Ny][Z_grid] * 0.5 );
            F_z = -( - U[X_grid][Y_grid][(Z_grid + Nz - 1)%Nz] * 0.5 + U[X_grid][Y_grid][(Z_grid + 1)%Nz] * 0.5 );
            /*
            F_x = -( U[(X_grid + Nx - 2)%Nx][Y_grid][Z_grid] / 12. - U[(X_grid + Nx - 1)%Nx][Y_grid][Z_grid] * (2. / 3.) +
                     U[(X_grid + 1)%Nx][Y_grid][Z_grid] * (2. / 3.) - U[(X_grid + 2)%Nx][Y_grid][Z_grid] * (1. / 12.) );
            F_y = -( U[X_grid][(Y_grid - 2 + Ny)%Ny][Z_grid] / 12. - U[X_grid][(Y_grid + Ny - 1)%Ny][Z_grid] * (2. / 3.) +
                     U[X_grid][(Y_grid + 1)%Ny][Z_grid] * (2. / 3.) - U[X_grid][(Y_grid + 2)%Ny][Z_grid] * (1. / 12.) );
            F_z = -( U[X_grid][Y_grid][(Z_grid + Nz - 2)%Nz] / 12. - U[X_grid][Y_grid][(Z_grid + Nz - 1)%Nz] * (2. / 3.) +
                     U[X_grid][Y_grid][(Z_grid + 1)%Nz] * (2. / 3.) - U[X_grid][Y_grid][(Z_grid + 2)%Nz] * (1. / 12.) );
	    */
        }
    } else if (mode == 1) {
        double f;
        X_grid = int( x/ dx);
        Y_grid = int( y/ dx);
        Z_grid = int( z/ dx);
	if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx-1) && (Y_grid<Ny-1) && (Z_grid<Nz-1)){
            for (int i = X_grid; i <= X_grid + 1; i++) {
                for (int j = Y_grid; j <= Y_grid + 1; j++) {
                    for (int k = Z_grid; k <= Z_grid + 1; k++) {
                        f = (1.0 - abs(x - i * dx) / dx) * (1.0 - abs(y - j * dy) / dy) * (1.0 - abs(z - k * dz) / dz);

			F_x -= f * (-U[(i + Nx - 1)%Nx][j][k] * 0.5 + U[(i + 1)%Nx][j][k] * 0.5);
                        F_y -= f * (-U[i][(j + Ny - 1)%Ny][k] * 0.5 + U[i][(j + 1)%Ny][k] * 0.5);
                        F_z -= f * (-U[i][j][(k + Nz - 1)%Nz] * 0.5 + U[i][j][(k + 1)%Nz] * 0.5);
			/*
                        F_x -= f * (U[(i + Nx - 2)%Nx][j][k] / 12. - U[(i + Nx - 1)%Nx][j][k] * (2. / 3.) +
                                    U[(i + 1)%Nx][j][k] * (2. / 3.) - U[(i + 2)%Nx][j][k] * (1. / 12.));
                        F_y -= f * (U[i][(j + Ny - 2)%Ny][k] / 12. - U[i][(j + Ny - 1)%Ny][k] * (2. / 3.) +
                                    U[i][(j + 1)%Ny][k] * (2. / 3.) - U[i][(j + 2)%Ny][k] * (1. / 12.));
                        F_z -= f * (U[i][j][(k + Nz - 2)%Nz] / 12. - U[i][j][(k + Nz - 1)%Nz] * (2. / 3.) +
                                    U[i][j][(k + 1)%Nz] * (2. / 3.) - U[i][j][(k + 2)%Nz] * (1. / 12.));
			*/
                    }
                }
	    }
        }
    } else if (mode == 2){
        double fx, fy, fz, f;
        X_grid = int( x / dx);
        Y_grid = int( y / dx);
        Z_grid = int( z / dx);
	if ((X_grid>=1) && (Y_grid>=1) && (Z_grid>=1) && (X_grid<Nx-1) && (Y_grid<Ny-1) && (Z_grid<Nz-1)){
            for (int i = X_grid - 1; i <= X_grid + 1; i++) {
                if (i == X_grid) fx = 0.75 - pow(x - i * dx, 2) / pow(dx, 2);
                else fx = 0.5 * pow(1.5 - abs(x - i * dx) / dx, 2);
                for (int j = Y_grid - 1; j <= Y_grid + 1; j++) {
                    if (j == Y_grid) fy = 0.75 - pow(y - j * dy, 2) / pow(dy, 2);
                    else fy = 0.5 * pow(1.5 - abs(y - j * dy) / dy, 2);
                    for (int k = Z_grid - 1; k <= Z_grid + 1; k++) {
                        if (k == Z_grid) fz = 0.75 - pow(z - k * dz, 2) / pow(dz, 2);
                        else fz = 0.5 * pow(1.5 - abs(z - k * dz) / dz, 2);
                        f = fx * fy * fz;
                        F_x -= f * (-U[(i + Nx - 1)%Nx][j][k] * 0.5 + U[(i + 1)%Nx][j][k] * 0.5);
                        F_y -= f * (-U[i][(j + Ny - 1)%Ny][k] * 0.5 + U[i][(j + 1)%Ny][k] * 0.5);
                        F_z -= f * (-U[i][j][(k + Nz - 1)%Nz] * 0.5 + U[i][j][(k + 1)%Nz] * 0.5);
			/*
                        F_x -= f * (U[(i + Nx - 2)%Nx][j][k] / 12. - U[(i + Nx - 1)%Nx][j][k] * (2. / 3.) +
                                    U[(i + 1)%Nx][j][k] * (2. / 3.) - U[(i + 2)%Nx][j][k] * (1. / 12.));
                        F_y -= f * (U[i][(j + Ny - 2)%Ny][k] / 12. - U[i][(j + Ny - 1)%Ny][k] * (2. / 3.) +
                                    U[i][(j + 1)%Ny][k] * (2. / 3.) - U[i][(j + 2)%Ny][k] * (1. / 12.));
                        F_z -= f * (U[i][j][(k + Nz - 2)%Nz] / 12. - U[i][j][(k + Nz - 1)%Nz] * (2. / 3.) +
                                    U[i][j][(k + 1)%Nz] * (2. / 3.) - U[i][j][(k + 2)%Nz] * (1. / 12.));
			*/
                	}
		}
            }
        }
    }
}

//Poisson Solver (FFT)
void FFT(double ***rho,double ***U){
    fftw::maxthreads = 1;
    array3<Complex> rho_x(Nx,Ny,Nz,sizeof(Complex));
    array3<Complex> phi_x(Nx,Ny,Nz,sizeof(Complex));
    array3<Complex> phi_k(Nx,Ny,Nz,sizeof(Complex));
    array3<Complex> rho_k(Nx,Ny,Nz,sizeof(Complex));
    fft3d Forward(Nx, Ny, Nz, -1, rho_x, rho_k);
    fft3d Backward(Nx, Ny, Nz, 1, phi_k, phi_x);
    double M = 0; // total mass
    /*
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                rho_x(i,j,k) = rho[i][j][k];
                M += rho[i][j][k];
            }
        }
    }
    for(int i = 0 ; i<Nx ; i++) for(int j = 0 ; j<Ny ; j++) for(int k = 0 ; k<Nz ; k++) if(i==0||i==Nx-1||j==0||j==Ny-1||k==0||k==Nz-1) rho_x(i,j,k) = -M/(Nx*Ny*Nz-(Nx-2)*(Ny-2)*(Nz-2));
    */
    Forward.fft0(rho_x, rho_k);
	
    for(int i = 0 ; i<Nx ; i++){
		for(int j = 0 ; j<Ny ; j++){
			for(int k = 0 ; k<Nz ; k++){
                phi_k(i,j,k) = -4*PI*G*dx*dy*dz*rho_k(i,j,k) / ( 4.0 * ( pow(sin(PI*min(i,Nx-i)/Nx),2.0) + pow(sin(PI*min(j,Ny-j)/Ny),2.0) + pow(sin(PI*min(k,Nz-k)/Nz),2.0) ) ) ;
				if(i==0 && j==0 && k==0) phi_k(i,j,k) = 0;
			}
		}
    }
	
    Backward.fft0Normalized(phi_k, phi_x);

    for(int i = 0 ; i<Nx ; i++){
        for(int j = 0 ; j<Ny ; j++){
            for(int k = 0 ; k<Nz ; k++){
                U[i][j][k] = real(phi_x(i,j,k));
            }
        }
    }
}

//Particle Mesh function
void mesh(double ***rho, double *x, double *y, double *z, int mode) {
    int X_grid, Y_grid, Z_grid;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                rho[i][j][k] = 0.;
            }
        }
    }
    if (mode == 0) {
        for (int p = 0; p < n; p++) {
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dx);
            Z_grid = int((z[p]) / dx);
            if (abs(x[p] - X_grid * dx) > abs(x[p] - (X_grid + 1) * dx)) X_grid++;
            if (abs(y[p] - Y_grid * dy) > abs(y[p] - (Y_grid + 1) * dy)) Y_grid++;
            if (abs(z[p] - Z_grid * dz) > abs(z[p] - (Z_grid + 1) * dz)) Z_grid++;
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx) && (Y_grid<Ny) && (Z_grid<Nz)) rho[X_grid][Y_grid][Z_grid] += m / (dx * dy * dz);
        }
    } else if (mode == 1) {
        for (int p = 0; p < n; p++) {
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dx);
            Z_grid = int((z[p]) / dx);
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid+1<Nx) && (Y_grid+1<Ny) && (Z_grid+1<Nz)){
                for (int i = X_grid; i <= X_grid + 1; i++) {
                    for (int j = Y_grid; j <= Y_grid + 1; j++) {
                        for (int k = Z_grid; k <= Z_grid + 1; k++) {
                            rho[i][j][k] += m * (1.0 - abs(x[p] - i * dx) / dx) * (1.0 - abs(y[p] - j * dy) / dy) * (1.0 - abs(z[p] - k * dz) / dz) / (dx * dy * dz);
                        }
                    }
                }
            }
        }
    } else if (mode == 2){
        double fx, fy, fz;
        for (int p = 0; p < n; p++) {
            X_grid = int((x[p]) / dx);
            Y_grid = int((y[p]) / dx);
            Z_grid = int((z[p]) / dx);
            if ((X_grid>0) && (Y_grid>0) && (Z_grid>0) && (X_grid+2<Nx) && (Y_grid+2<Ny) && (Z_grid+2<Nz)){
                for (int i = X_grid - 1; i <= X_grid + 1; i++) {
                    if (i == X_grid) fx = 0.75 - pow(x[p] - i * dx, 2) / pow(dx, 2);
                    else fx = 0.5 * pow(1.5 - abs(x[p] - i * dx) / dx, 2);
                    for (int j = Y_grid - 1; j <= Y_grid + 1; j++) {
                        if (j == Y_grid) fy = 0.75 - pow(y[p] - j * dy, 2) / pow(dy, 2);
                        else fy = 0.5 * pow(1.5 - abs(y[p] - j * dy) / dy, 2);
                        for (int k = Z_grid - 1; k <= Z_grid + 1; k++) {
                            if (k == Z_grid) fz = 0.75 - pow(z[p] - k * dz, 2) / pow(dz, 2);
                            else fz = 0.5 * pow(1.5 - abs(z[p] - k * dz) / dz, 2);
                            rho[i][j][k] += m * fx * fy * fz / (dx * dy * dz);
                        }
              	    }
            	}
            }
        }
    }
}

double Get_Energy(double *x, double *y, double *z, double *vx, double *vy, double *vz){
    double E = 0;
    for(int p = 0 ; p<n ; p++){
        E += 0.5 * m * (vx[p]*vx[p]+vy[p]*vy[p]+vz[p]*vz[p]);
        for(int q = p+1 ; q<n ; q++){
            E += -G*m*m/sqrt(pow(x[p]-x[q],2)+pow(y[p]-y[q],2)+pow(z[p]-z[q],2));
        }
    }
    return E;
}

int main() {
    /* Variables */
    double t = 0.0; //time
    double t_end = 10.0; //ending time
    double dt = 0.01; // time step
    double PDx = 0.1, PDy = 0.1, PDz = 0.1; //size of particle clumps
    double * x = new double[n]; //positions of the particles
    double * y = new double[n];
    double * z = new double[n];
    double * vx = new double[n]; //velocities of the particles
    double * vy = new double[n];
    double * vz = new double[n];
    double *** rho = new double ** [Nx]; // mass density
    double *** U = new double ** [Nx]; // Periodic B.C.
    
    srand(time(NULL));
    /* Initialization */
    //Random distribution
    for (int i = 0; i < n; i++) {
        x[i] = Lx * PDx * (rand() / (double) RAND_MAX -0.5) + Lx/2;
        y[i] = Ly * PDy * (rand() / (double) RAND_MAX -0.5) + Ly/2;
        z[i] = Lz * PDz * (rand() / (double) RAND_MAX -0.5) + Lz/2;
        double r0 = pow(pow(PDx, 2) + pow(PDy, 2) + pow(PDz, 2),0.5);
        double v0 = sqrt(G * m / r0);
        vx[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) / 10.;
        vy[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) / 10.;
        vz[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) / 10.;
    }

    //set U = 0.0
    for (int i = 0; i < Nx; i++) {
        rho[i] = new double * [Ny];
        U[i] = new double * [Ny];
        for (int j = 0; j < Ny; j++) {
            rho[i][j] = new double[Nz];
            U[i][j] = new double[Nz];
            for (int k = 0; k < Nz; k++) {
                U[i][j][k] = 0.;
            }
        }
    }
    
    while (t <= t_end) {
	
        //DKD
        //Starting to calculate force on partilces
        if (OI_mode == 0) {
            //drift: update position by 0.5*dt
            for (int i = 0; i < n; i++) {
                x[i] += vx[i] * 0.5 * dt;
                y[i] += vy[i] * 0.5 * dt;
                z[i] += vz[i] * 0.5 * dt;
            }
            //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
            mesh(rho, x, y, z, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                vx[i] += F_x / m * dt;
                vy[i] += F_y / m * dt;
                vz[i] += F_z / m * dt;
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
            mesh(rho, x, y, z, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                vx[i] += F_x / m * 0.5 * dt;
                vy[i] += F_y / m * 0.5 * dt;
                vz[i] += F_z / m * 0.5 * dt;
            }
            //drift: update position by dt
            for (int i = 0; i < n; i++) {
                x[i] += vx[i] * dt;
                y[i] += vy[i] * dt;
                z[i] += vz[i] * dt;
            }
            //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
            mesh(rho, x, y, z, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                vx[i] += F_x / m * 0.5 * dt;
                vy[i] += F_y / m * 0.5 * dt;
                vz[i] += F_z / m * 0.5 * dt;
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
            
            mesh(rho, x, y, z, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                vx[i] += F_x / m * d1 * dt;
                vy[i] += F_y / m * d1 * dt;
                vz[i] += F_z / m * d1 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                x[i] += vx[i] * c2 * dt;
                y[i] += vy[i] * c2 * dt;
                z[i] += vz[i] * c2 * dt;
            }
            
            mesh(rho, x, y, z, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                vx[i] += F_x / m * d2 * dt;
                vy[i] += F_y / m * d2 * dt;
                vz[i] += F_z / m * d2 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                x[i] += vx[i] * c3 * dt;
                y[i] += vy[i] * c3 * dt;
                z[i] += vz[i] * c3 * dt;
            }
            
            mesh(rho, x, y, z, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                vx[i] += F_x / m * d3 * dt;
                vy[i] += F_y / m * d3 * dt;
                vz[i] += F_z / m * d3 * dt;
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
            mesh(rho, x_tmp, y_tmp, z_tmp, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x[i], y[i], z[i], F_x, F_y, F_z, force_mode);
                kr1[i][0] = vx[i];
                kr1[i][1] = vy[i];
                kr1[i][2] = vy[i];
                kv1[i][0] = F_x/m;
                kv1[i][1] = F_y/m;
                kv1[i][2] = F_z/m;
                kr2[i][0] = kr1[i][0] + kv1[i][0] * 0.5 * dt;
                kr2[i][1] = kr1[i][1] + kv1[i][1] * 0.5 * dt;
                kr2[i][2] = kr1[i][1] + kv1[i][1] * 0.5 * dt;
                x_tmp[i] += kr1[i][0]*0.5*dt;
                y_tmp[i] += kr1[i][1]*0.5*dt;
                z_tmp[i] += kr1[i][2]*0.5*dt;
            }
            mesh(rho, x_tmp, y_tmp, z_tmp, mesh_mode);
            FFT(rho,U);
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Force_of_Particle(U, x_tmp[i], y_tmp[i], z_tmp[i], F_x, F_y, F_z, force_mode);
                kv2[i][0] = F_x/m;
                kv2[i][1] = F_y/m;
                kv2[i][2] = F_z/m;
                kr3[i][0] = kr1[i][0] + kv2[i][0] * 0.5 * dt;
                kr3[i][1] = kr1[i][1] + kv2[i][1] * 0.5 * dt;
                kr3[i][2] = kr1[i][1] + kv2[i][1] * 0.5 * dt;
                x_tmp[i]  = x[i] + kr2[i][0]*0.5*dt;
                y_tmp[i]  = y[i] + kr2[i][1]*0.5*dt;
                z_tmp[i]  = z[i] + kr2[i][2]*0.5*dt;
            }
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                mesh(rho, x_tmp, y_tmp, z_tmp, mesh_mode);
                FFT(rho,U);
                Get_Force_of_Particle(U, x_tmp[i], y_tmp[i], z_tmp[i], F_x, F_y, F_z, force_mode);
                kv3[i][0] = F_x/m;
                kv3[i][1] = F_y/m;
                kv3[i][2] = F_z/m;
                kr4[i][0] = kr1[i][0] + kv3[i][0] * dt;
                kr4[i][1] = kr1[i][1] + kv3[i][1] * dt;
                kr4[i][2] = kr1[i][1] + kv3[i][1] * dt;
                x_tmp[i]  = x[i] + kr3[i][0]*0.5*dt;
                y_tmp[i]  = y[i] + kr3[i][1]*0.5*dt;
                z_tmp[i]  = z[i] + kr3[i][2]*0.5*dt;
            }
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                mesh(rho, x_tmp, y_tmp, z_tmp, mesh_mode);
                FFT(rho,U);
                Get_Force_of_Particle(U, x_tmp[i], y_tmp[i], z_tmp[i], F_x, F_y, F_z, force_mode);
                kv4[i][0] = F_x/m;
                kv4[i][1] = F_y/m;
                kv4[i][2] = F_z/m;
            }
            for (int i = 0; i < n; i++) {
                x[i] += (kr1[i][0] + 2.0 * kr2[i][0] + 2.0 * kr3[i][0] + 1.0 * kr4[i][0]) / 6.0 * dt;
                y[i] += (kr1[i][1] + 2.0 * kr2[i][1] + 2.0 * kr3[i][1] + 1.0 * kr4[i][1]) / 6.0 * dt;
                z[i] += (kr1[i][1] + 2.0 * kr2[i][2] + 2.0 * kr3[i][2] + 1.0 * kr4[i][2]) / 6.0 * dt;
                vx[i] += (kv1[i][0] + 2. * kv2[i][0] + 2.0 * kv3[i][0] + 1.0 * kv4[i][0]) / 6.0 * dt;
                vy[i] += (kv1[i][1] + 2. * kv2[i][1] + 2.0 * kv3[i][1] + 1.0 * kv4[i][1]) / 6.0 * dt;
                vz[i] += (kv1[i][1] + 2. * kv2[i][1] + 2.0 * kv3[i][1] + 1.0 * kv4[i][1]) / 6.0 * dt;
            }
        }
        
        // check conservation
        double Px = 0, Py = 0, Pz = 0;
        double X = 0, Y = 0, Z = 0;
        double M = 0;
        int n_in = 0;
        if((int)(t/dt)%100==0){
	    FILE *den_output;
            char fname[100];
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
            for(int i = 0 ; i<Nx ; i++) for(int j = 0 ; j<Ny ; j++) for(int k = 0 ; k<Nz ; k++) M += rho[i][j][k]*dx*dy*dz;
            printf("t = %.3f\n", t);
            printf("Px = %.3f \t Py = %.3f \t Pz = %.3f\tphi(0.5,0.5,0.5) = %.3f\n", Px, Py, Pz, U[Nx/2][Ny/2][Nz/2]);
            printf("n_in = %d\tM = %.3f\tE = %.3f\n", n_in, M, Get_Energy(x,y,z,vx,vy,vz));
	    for (int i = 0; i < n; i++) fprintf (den_output, "%g  %g  %g   \n",x[i], y[i], z[i] );
            fclose(den_output);
        }     
	t += dt;

    }
    return EXIT_SUCCESS;
}


