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
int OI_mode = 2;    //Orbit integration mode. 0: DKD 1:KDK 2:fourth-order symplectic integrator 3:RK4  4:Hermite

//-----------------------------------------------------------constants-------------------------------------------------
double G = 1.0;                                  // gravitational constant
double Lx = 1.0, Ly = 1.0, Lz = 1.0;             // domain size of 3D box
int N = 128;                                     // # of grid points
int Nx = N, Ny = N, Nz = N;
double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz; // spatial resolution
int n = 3;                                    // # of particles
double m = 1.0;                                  // particle mass
double t = 0.0;                                  // time
double PDx = 0.2, PDy = 0.2, PDz = 0.2;          // size of particle clumps
double dt = 0.1*sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2))/sqrt(n*G*m/sqrt(pow(PDx, 2) + pow(PDy, 2) + pow(PDz, 2))); //time steps
double t_end = dt*900.0;                         // ending time                             
double vmax = 1.0;                               // initial maximal velocity weight


//----------------------------------------------------------functions------------------------------------------------
//Particle Force Interpolation Function
//
void Get_Direct_Force_of_Particle(double * X, double *Y, double *Z, double x, double y, double z, double & F_x, double & F_y, double & F_z) {
    
    for (int i=0;i<n;i++){
	double d_15 = pow(pow(X[i]-x,2)+pow(Y[i]-y,2)+pow(Z[i]-z,2),1.5);
	if (d_15!=0.0){
	   F_x += G*m*m/d_15*(X[i]-x);
	   F_y += G*m*m/d_15*(Y[i]-y);
	   F_z += G*m*m/d_15*(Z[i]-z);
	}
    }

}

//Particle Mesh function
void mesh(double ***rho, double *x, double *y, double *z, int mode) {

    int X_grid, Y_grid, Z_grid; //grid positions of particles
    //initialize rho

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                rho[i][j][k] = 0.;
            }
        }
    }

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
            if ((X_grid>=0) && (Y_grid>=0) && (Z_grid>=0) && (X_grid<Nx) && (Y_grid<Ny) && (Z_grid<Nz)) rho[X_grid][Y_grid][Z_grid] += m / (dx * dy * dz);
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
                             rho[i][j][k] += m * (1.0 - abs(x[p] - i * dx) / dx) * (1.0 - abs(y[p] - j * dy) / dy) * (1.0 - abs(z[p] - k * dz) / dz) / (dx * dy * dz);
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
                            rho[i][j][k] += m * fx * fy * fz / (dx * dy * dz);
                        }
              	    }
            	}
            }
        }
    }
}


double Get_Energy(double *x, double *y, double *z, double *vx, double *vy, double *vz){
    double E = 0.0;
    for(int p = 0 ; p<n ; p++){
        E += 0.5 * m * (vx[p]*vx[p]+vy[p]*vy[p]+vz[p]*vz[p]); //kinetic energy
        for(int q = p+1 ; q<n ; q++){
            E += -G*m*m/sqrt(pow(x[p]-x[q],2)+pow(y[p]-y[q],2)+pow(z[p]-z[q],2)); //potential energy
        }
    }
    return E;
}

int main() {
    /* Variables */
    double t = 0.0; //time
    double * x = new double[n]; //positions of the particles
    double * y = new double[n];
    double * z = new double[n];
    double * vx = new double[n]; //velocities of the particles
    double * vy = new double[n];
    double * vz = new double[n];
    double *** rho = new double ** [Nx]; // mass density
    int frame = 0;
 
    //initialize rho, U and W
    for (int i = 0; i < Nx; i++) {
        rho[i] = new double * [Ny];
        for (int j = 0; j < Ny; j++) {
            rho[i][j] = new double[Nz];
        }
    }
    
    srand(time(NULL));
    /* Initialization */
    //Random distribution
    double r0 = pow(pow(PDx, 2) + pow(PDy, 2) + pow(PDz, 2),0.5); //mean distance
    double v0 = vmax*sqrt(G * m / r0 / 2);                        //Virial speed
    for (int i = 0; i < n; i++) {
        x[i] = PDx * (rand() / (double) RAND_MAX -0.5) + 0.5;
        y[i] = PDy * (rand() / (double) RAND_MAX -0.5) + 0.5;
        z[i] = PDz * (rand() / (double) RAND_MAX -0.5) + 0.5;
        vx[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) *2.0;
        vy[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) *2.0;
        vz[i] = v0 * ( rand() / (double) RAND_MAX - 0.5) *2.0;
    }
    if(n==2){
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
    }
    if(n==3){
        x[0] = 0.6;
        y[0] = 0.5;
        z[0] = 0.5;
        x[1] = 0.4;
        y[1] = 0.5;
        z[1] = 0.5;
        x[2] = 0.7;
        y[2] = 0.3;
        z[2] = 0.5;
        vx[0] = 0.0;
        vy[0] = 0.5*sqrt(1.0/0.2);
        vz[0] = 0.0;
        vx[1] = 0.0;
        vy[1] = -0.5*sqrt(1.0/0.2);
        vz[1] = 0.0;
        vx[2] = 0.; vy[2] = 0; vz[2] = 0;
    }
    
    printf("particle number = %d orbit mode = %d dt = %.3e\n particle size = %.2f vmax = %.3f\n",n ,OI_mode,dt,PDx,v0);
    
    while (t <= t_end) {
     
        // check conservation
        if((int)(t/dt)%10==0){
            mesh(rho, x, y, z, mesh_mode);
            double Px = 0, Py = 0, Pz = 0;
            double X = 0, Y = 0, Z = 0;
            double M = 0;
            double E = Get_Energy(x,y,z,vx,vy,vz);
            int n_in = 0;
            
            char fname[100], name[100];
            sprintf(fname,"./output/position_%04d", frame);
            sprintf(name,"./output/density_%04d", frame);
            //FILE *position_output = fopen(fname,"w");
            FILE *density_output = fopen(name,"w");
            
            cout << "================================\n";
            for(int p = 0 ; p<n ; p++){
                Px += m*vx[p];
                Py += m*vy[p];
                Pz += m*vz[p];
                if(x[p]>0&&x[p]<Lx&&y[p]>0&&y[p]<Ly&&z[p]>0&&z[p]<Lz) n_in++;
            }
            for(int i = 0 ; i<Nx ; i++) for(int j = 0 ; j<Ny ; j++) for(int k = 0 ; k<Nz ; k++) M += rho[i][j][k]*dx*dy*dz;
            printf("t = %.3f\n", t);
            printf("Px = %.3f \t Py = %.3f \t Pz = %.3f\n", Px, Py, Pz);
            printf("n_in = %d\tM = %.3f\tE = %.3f\n", n_in, M, E);
            
            fprintf(density_output, "%.4f\t%e\t%d\n", t, E, n_in);
            //for (int i = 0; i < n; i++) fprintf(position_output, "%g  %g  %g   \n",x[i], y[i], z[i] );
            for(int i = 0 ; i<Nx ; i++) for(int j = 0 ; j<Ny ; j++) for(int k = 0 ; k<Nz ; k++) if(rho[i][j][k]!=0) fprintf(density_output, "%d\t%d\t%d\t%e\n", i, j, k, rho[i][j][k]);
            //fclose(position_output);
            fclose(density_output);
            
            frame++;
        }

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
	    
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
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
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
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
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
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
                        
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
		vx[i] += F_x / m * d1 * dt;
                vy[i] += F_y / m * d1 * dt;
                vz[i] += F_z / m * d1 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                x[i] += vx[i] * c2 * dt;
                y[i] += vy[i] * c2 * dt;
                z[i] += vz[i] * c2 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
                vx[i] += F_x / m * d2 * dt;
                vy[i] += F_y / m * d2 * dt;
                vz[i] += F_z / m * d2 * dt;
            }
            
            for (int i = 0; i < n; i++) {
                x[i] += vx[i] * c3 * dt;
                y[i] += vy[i] * c3 * dt;
                z[i] += vz[i] * c3 * dt;
            }

            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
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

            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
                kr1[i][0] = vx[i];
                kr1[i][1] = vy[i];
                kr1[i][2] = vz[i];
                kv1[i][0] = F_x/m;
                kv1[i][1] = F_y/m;
                kv1[i][2] = F_z/m;
                kr2[i][0] = kr1[i][0] + kv1[i][0] * 0.5 * dt;
                kr2[i][1] = kr1[i][1] + kv1[i][1] * 0.5 * dt;
                kr2[i][2] = kr1[i][2] + kv1[i][2] * 0.5 * dt;
                x_tmp[i] += kr1[i][0]*0.5*dt;
                y_tmp[i] += kr1[i][1]*0.5*dt;
                z_tmp[i] += kr1[i][2]*0.5*dt;
            }
           
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
                kv2[i][0] = F_x/m;
                kv2[i][1] = F_y/m;
                kv2[i][2] = F_z/m;
                kr3[i][0] = kr1[i][0] + kv2[i][0] * 0.5 * dt;
                kr3[i][1] = kr1[i][1] + kv2[i][1] * 0.5 * dt;
                kr3[i][2] = kr1[i][2] + kv2[i][2] * 0.5 * dt;
                x_tmp[i]  = x[i] + kr2[i][0]*0.5*dt;
                y_tmp[i]  = y[i] + kr2[i][1]*0.5*dt;
                z_tmp[i]  = z[i] + kr2[i][2]*0.5*dt;
            }
	    
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
                kv3[i][0] = F_x/m;
                kv3[i][1] = F_y/m;
                kv3[i][2] = F_z/m;
                kr4[i][0] = kr1[i][0] + kv3[i][0] * dt;
                kr4[i][1] = kr1[i][1] + kv3[i][1] * dt;
                kr4[i][2] = kr1[i][2] + kv3[i][2] * dt;
                x_tmp[i]  = x[i] + kr3[i][0]*0.5*dt;
                y_tmp[i]  = y[i] + kr3[i][1]*0.5*dt;
                z_tmp[i]  = z[i] + kr3[i][2]*0.5*dt;
            }
            
            for (int i = 0; i < n; i++) {
                double F_x=0.0, F_y=0.0, F_z=0.0;
                Get_Direct_Force_of_Particle(x, y, z, x[i], y[i], z[i], F_x, F_y, F_z);
                kv4[i][0] = F_x/m;
                kv4[i][1] = F_y/m;
                kv4[i][2] = F_z/m;
            }
            for (int i = 0; i < n; i++) {
                x[i] += (kr1[i][0] + 2.0 * kr2[i][0] + 2.0 * kr3[i][0] + 1.0 * kr4[i][0]) / 6.0 * dt;
                y[i] += (kr1[i][1] + 2.0 * kr2[i][1] + 2.0 * kr3[i][1] + 1.0 * kr4[i][1]) / 6.0 * dt;
                z[i] += (kr1[i][2] + 2.0 * kr2[i][2] + 2.0 * kr3[i][2] + 1.0 * kr4[i][2]) / 6.0 * dt;
                vx[i] += (kv1[i][0] + 2. * kv2[i][0] + 2.0 * kv3[i][0] + 1.0 * kv4[i][0]) / 6.0 * dt;
                vy[i] += (kv1[i][1] + 2. * kv2[i][1] + 2.0 * kv3[i][1] + 1.0 * kv4[i][1]) / 6.0 * dt;
                vz[i] += (kv1[i][2] + 2. * kv2[i][2] + 2.0 * kv3[i][2] + 1.0 * kv4[i][2]) / 6.0 * dt;
            }
        }
        
	t += dt;

    }
    
    return EXIT_SUCCESS;
}



