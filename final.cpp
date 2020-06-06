
// ==============
// PARTICLE MESH
// ==============

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<ctime>

//--------------------------------------------------------mode selection----------------------------------------------
int mesh_mode = 0; // 0: NGP ; 1: CIC ; 2: TSC
int OI_mode   = 0; //Orbit integration mode. 0: DKD 1:KDK 2:fourth-order symplectic integrator 3:RK4 $:Hermite

//-----------------------------------------------------------constants-------------------------------------------------
double G = 1.0; // gravitational constant
double Lx = 1.0, Ly = 1.0, Lz = 1.0; // domain size of 3D box
int Nx = 3, Ny = 3, Nz = 3; // # of grid points
double dx = Lx/(Nx-1), dy = Ly/(Ny-1), dz = Lz/(Nz-1); // spatial resolution
int n = 1000; // # of particles
double m = 1.0; // particle mass

//----------------------------------------------------------functions------------------------------------------------
//Particle Force Interpolation Function
void Get_Force_of_Particle(double ***U, double x, double y, double z, double &F_x, double &F_y, double &F_z, int mode){

        int X_grid, Y_grid, Z_grid;
	if( mode == 0 ){
            	X_grid = int( (x + Lx/2)/dx );
            	Y_grid = int( (y + Ly/2)/dx );
            	Z_grid = int( (z + Lz/2)/dx );
            	if( abs(x - X_grid*dx) > abs(x - (X_grid + 1)*dx) ) X_grid++;
            	if( abs(y - Y_grid*dy) > abs(y - (Y_grid + 1)*dy) ) Y_grid++;
            	if( abs(z - Z_grid*dz) > abs(z - (Z_grid + 1)*dz) ) Z_grid++;
            	F_x = U[X_grid-2][Y_grid][Z_grid]/12. -  U[X_grid-1][Y_grid][Z_grid]*(2./3.) 
		    + U[X_grid+1][Y_grid][Z_grid]*(2./3.) - U[X_grid+2][Y_grid][Z_grid]*(1./12.);
        	F_y = U[X_grid][Y_grid-2][Z_grid]/12. -  U[X_grid][Y_grid-1][Z_grid]*(2./3.)
		    + U[X_grid][Y_grid+1][Z_grid]*(2./3.) - U[X_grid][Y_grid+2][Z_grid]*(1./12.);
        	F_z = U[X_grid][Y_grid][Z_grid-2]/12. -  U[X_grid][Y_grid][Z_grid-1]*(2./3.) 
	            + U[X_grid][Y_grid][Z_grid+1]*(2./3.) - U[X_grid][Y_grid][Z_grid+2]*(1./12.);
    	}else if( mode == 1 ){
	    	double f;
            	X_grid = int( (x + Lx/2)/dx );
           	Y_grid = int( (y + Ly/2)/dx );
           	Z_grid = int( (z + Lz/2)/dx );
            	for(int i = X_grid ; i<=X_grid+1 ; i++){
                	for(int j = Y_grid ; j<=Y_grid+1 ; j++){
                    		for(int k = Z_grid ; k<=Z_grid+1 ; k++){
					f = ( 1.0 - abs(x + 0.5 - i*dx)/dx)*( 1.0-abs(y + 0.5 - j*dy)/dy )*( 1.0-abs(z + 0.5 - k*dz)/dz);
                       			F_x += f*(U[X_grid-2][Y_grid][Z_grid]/12. -  U[X_grid-1][Y_grid][Z_grid]*(2./3.) 
		    	        		+ U[X_grid+1][Y_grid][Z_grid]*(2./3.) - U[X_grid+2][Y_grid][Z_grid]*(1./12.));
					F_y += f*(U[X_grid][Y_grid-2][Z_grid]/12. -  U[X_grid][Y_grid-1][Z_grid]*(2./3.)
		    	        		+ U[X_grid][Y_grid+1][Z_grid]*(2./3.) - U[X_grid][Y_grid+2][Z_grid]*(1./12.));
					F_z += f*(U[X_grid][Y_grid][Z_grid-2]/12. -  U[X_grid][Y_grid][Z_grid-1]*(2./3.) 
	            				+ U[X_grid][Y_grid][Z_grid+1]*(2./3.) - U[X_grid][Y_grid][Z_grid+2]*(1./12.));
                    		}
                	}
        	}
   	}else if( mode == 2 ){
            	double fx, fy, fz, f;
           	X_grid = int( (x + Lx/2)/dx );
            	Y_grid = int( (y + Ly/2)/dx );
            	Z_grid = int( (z + Lz/2)/dx );
            	for(int i = X_grid ; i<=X_grid+1 ; i++){
                	if(i == X_grid) fx = 0.75 - pow( x + 0.5 - i*dx , 2 ) / pow(dx,2);
                	else fx = 0.5 * pow( 1.5 - abs( x + 0.5 - i*dx )/dx , 2);
               	 	for(int j = Y_grid ; j<=Y_grid+1 ; j++){
                    		if(j == Y_grid) fy = 0.75 - pow( y + 0.5 - j*dy , 2 ) / pow(dy,2);
                    		else fy = 0.5 * pow( 1.5 - abs( y + 0.5 - j*dy )/dy , 2);
                   		for(int k = Z_grid ; k<=Z_grid+1 ; k++){
                        		if(k == Z_grid) fz = 0.75 - pow( z + 0.5 - k*dz , 2 ) / pow(dz,2);
                        		else fz = 0.5 * pow( 1.5 - abs( z + 0.5 - k*dz )/dz , 2);
					f = fx*fy*fz;
                        		F_x += f*(U[X_grid-2][Y_grid][Z_grid]/12. -  U[X_grid-1][Y_grid][Z_grid]*(2./3.) 
		    	        		+ U[X_grid+1][Y_grid][Z_grid]*(2./3.) - U[X_grid+2][Y_grid][Z_grid]*(1./12.));
					F_y += f*(U[X_grid][Y_grid-2][Z_grid]/12. -  U[X_grid][Y_grid-1][Z_grid]*(2./3.)
		    	        		+ U[X_grid][Y_grid+1][Z_grid]*(2./3.) - U[X_grid][Y_grid+2][Z_grid]*(1./12.));
					F_z += f*(U[X_grid][Y_grid][Z_grid-2]/12. -  U[X_grid][Y_grid][Z_grid-1]*(2./3.) 
	            				+ U[X_grid][Y_grid][Z_grid+1]*(2./3.) - U[X_grid][Y_grid][Z_grid+2]*(1./12.));
                    }
                }
            }
        }
}

//Poisson Solver (FFT)

//Particle Mesh function
void mesh(double ***rho, double *x, double *y, double *z, int mode){
    int X_grid, Y_grid, Z_grid;
    for(int i = 0 ; i<Nx ; i++){
        for(int j = 0 ; j<Ny ; j++){
            for(int k = 0 ; k<Nz ; k++){
                rho[i][j][k] = 0.;
            }
        }
    }
    if( mode == 0 ){
        for(int p = 0 ; p<n ; p++){
            X_grid = int( (x[p] + Lx/2)/dx );
            Y_grid = int( (y[p] + Ly/2)/dx );
            Z_grid = int( (z[p] + Lz/2)/dx );
            if( abs(x[p] - X_grid*dx) > abs(x[p] - (X_grid + 1)*dx) ) X_grid++;
            if( abs(x[p] - Y_grid*dy) > abs(x[p] - (Y_grid + 1)*dy) ) Y_grid++;
            if( abs(x[p] - Z_grid*dz) > abs(x[p] - (Z_grid + 1)*dz) ) Z_grid++;
            rho[X_grid][Y_grid][Z_grid] += m/(dx*dy*dz);
        }
    }else if( mode == 1 ){
        for(int p = 0 ; p<n ; p++){
            X_grid = int( (x[p] + Lx/2)/dx );
            Y_grid = int( (y[p] + Ly/2)/dx );
            Z_grid = int( (z[p] + Lz/2)/dx );
            for(int i = X_grid ; i<=X_grid+1 ; i++){
                for(int j = Y_grid ; j<=Y_grid+1 ; j++){
                    for(int k = Z_grid ; k<=Z_grid+1 ; k++){
                        rho[i][j][k] += m*( 1.0 - abs(x[p] + 0.5 - i*dx)/dx)*( 1.0-abs(y[p] + 0.5 - j*dy)/dy )*( 1.0-abs(z[p] + 0.5 - k*dz)/dz) /(dx*dy*dz);
                    }
                }
            }
        }
    }else if( mode == 2 ){
        double fx, fy, fz;
        for(int p = 0 ; p<n ; p++){
            X_grid = int( (x[p] + Lx/2)/dx );
            Y_grid = int( (y[p] + Ly/2)/dx );
            Z_grid = int( (z[p] + Lz/2)/dx );
            for(int i = X_grid-1 ; i<=X_grid+1 ; i++){
                if(i == X_grid) fx = 0.75 - pow( x[p] + 0.5 - i*dx , 2 ) / pow(dx,2);
                else fx = 0.5 * pow( 1.5 - abs(x[p] + 0.5 - i*dx)/dx , 2);
                for(int j = Y_grid-1 ; j<=Y_grid+1 ; j++){
                    if(j == Y_grid) fy = 0.75 - pow( y[p] + 0.5 - j*dy , 2 ) / pow(dy,2);
                    else fy = 0.5 * pow( 1.5 - abs(y[p] + 0.5 - j*dy)/dy , 2);
                    for(int k = Z_grid-1 ; k<=Z_grid+1 ; k++){
                        if(k == Z_grid) fz = 0.75 - pow( z[p] + 0.5 - k*dz , 2 ) / pow(dz,2);
                        else fz = 0.5 * pow( 1.5 - abs(z[p] + 0.5 - k*dz)/dz , 2);
                        rho[i][j][k] += m*fx*fy*fz/(dx*dy*dz);
                    }
                }
            }
        }
    }
}


int main(){
	/* Variables */
    double t = 0.0; //time
    double t_end = 1.0; //ending time
    double dt = 0.1; // time step
    double PDx = 0.2, PDy = 0.2, PDz = 0.2; //size of particle clumps
    double *x = new double [n]; //positions of the particles
    double *y = new double [n];
    double *z = new double [n];
    double *vx = new double [n]; //velocities of the particles
    double *vy = new double [n];
    double *vz = new double [n];
    double ***rho = new double ** [Nx]; // mass density
    double ***U = new double ** [Nx]; // Periodic B.C.
    
	srand(time(NULL));
	/* Initialization */
        //Random distribution
	for(int i = 0 ; i<n ; i++){
		x[i] = Lx * .5*PDx*(rand()/(double)RAND_MAX-0.5);
		y[i] = Ly * .5*PDy*(rand()/(double)RAND_MAX-0.5);
		z[i] = Lz * .5*PDz*(rand()/(double)RAND_MAX-0.5);
		double r0 = pow(PDx,2) + pow(PDy,2) + pow(PDz,2);
		double v0 = sqrt(G*m/r0);
		vx[i] = v0 * (2*rand()/(double)RAND_MAX-1.0);
		vy[i] = v0 * (2*rand()/(double)RAND_MAX-1.0);
		vz[i] = v0 * (2*rand()/(double)RAND_MAX-1.0);
	}
        //set U = 0.0
	for(int i = 0 ; i<Nx ; i++){
		rho[i] = new double *[Ny];
		U[i] = new double *[Ny];
		for(int j = 0 ; j<Ny ; j++){
			rho[i][j] = new double [Nz];
			U[i][j] = new double [Nz];
			for(int k = 0 ; k<Nz ; k++){
				U[i][j][k] = 0.;
			}
		}
	}
    
    double M = 0; // total mass
    mesh(rho,x,y,z,0);
    for(int i = 0 ; i<Nx ; i++){
        for(int j = 0 ; j<Ny ; j++){
            for(int k = 0 ; k<Nz ; k++){
                M += rho[i][j][k]*dx*dy*dz;
            }
        }
    }
    std::cout << M << std::endl;
    mesh(rho,x,y,z,1);
   	M = 0;
   	for(int i = 0 ; i<Nx ; i++){
        for(int j = 0 ; j<Ny ; j++){
            for(int k = 0 ; k<Nz ; k++){
                M += rho[i][j][k]*dx*dy*dz;
            }
        }
    }
    std::cout << M << std::endl;
	mesh(rho,x,y,z,2);
	M = 0;
    for(int i = 0 ; i<Nx ; i++){
        for(int j = 0 ; j<Ny ; j++){
            for(int k = 0 ; k<Nz ; k++){
                M += rho[i][j][k]*dx*dy*dz;
            }
        }
    }
    std::cout << M << std::endl;

    while (t<=t_end){
    //DKD
    //Starting to calculate force on partilces
    if (OI_mode == 0){
        //drift: update position by 0.5*dt
	for(int i = 0 ; i<n ; i++){          
           x[i] += vx[i]*0.5*dt;
           y[i] += vy[i]*0.5*dt;
           z[i] += vz[i]*0.5*dt;
	}
        //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
	mesh(rho,x,y,z,mesh_mode);
	for(int i = 0 ; i<n ; i++){
	    double F_x, F_y, F_z;
	    Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	    vx[i] += F_x/m*dt;
	    vy[i] += F_y/m*dt;
	    vz[i] += F_z/m*dt;
	}
        //drift: update position by 0.5*dt
	for(int i = 0 ; i<n ; i++){
           x[i] += vx[i]*0.5*dt;
           y[i] += vy[i]*0.5*dt;
           z[i] += vz[i]*0.5*dt;
        }
    }

    //KDK
    if (OI_mode == 1){
	//kick
	mesh(rho,x,y,z,mesh_mode);
        for(int i = 0 ; i<n ; i++){
	  double F_x, F_y, F_z;
	  Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	  vx[i] += F_x/m*0.5*dt;
	  vy[i] += F_y/m*0.5*dt;
	  vz[i] += F_z/m*0.5*dt;
	}
	//drift: update position by dt
	for(int i = 0 ; i<n ; i++){
	   x[i] += vx[i]*dt;
           y[i] += vy[i]*dt;
           z[i] += vz[i]*dt;
	}
        //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
	mesh(rho,x,y,z,mesh_mode);
	for(int i = 0 ; i<n ; i++){	   
	   double F_x, F_y, F_z;
	   Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	   vx[i] += F_x/m*0.5*dt;
	   vy[i] += F_y/m*0.5*dt;
	   vz[i] += F_z/m*0.5*dt;
        }
    }

    //fourth-order symplectic integration
    if (OI_mode == 2){
        //fourth-order symplectic coefficients
        double w1 = 1.0/( 2.0 - pow(2.0,1.0/3.0) );
        double w0 = 1.0 - 2.0*w1;
        double c1 = w1/2.0;
        double c2 = ( w1+w0 )/2.0;
        double c3 = c2;
        double c4 = c1; 
        double d1 = w1;
        double d2 = w0;
        double d3 = w1;

        for(int i = 0 ; i<n ; i++){
	   x[i] += vx[i]*c1*dt;
           y[i] += vy[i]*c1*dt;
           z[i] += vz[i]*c1*dt;
	}
	
	mesh(rho,x,y,z,mesh_mode);
        for(int i = 0 ; i<n ; i++){
	  double F_x, F_y, F_z;
	  Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	  vx[i] += F_x/m*d1*dt;
	  vy[i] += F_y/m*d1*dt;
	  vz[i] += F_z/m*d1*dt;
	}

	for(int i = 0 ; i<n ; i++){
	   x[i] += vx[i]*c2*dt;
           y[i] += vy[i]*c2*dt;
           z[i] += vz[i]*c2*dt;
	}

	mesh(rho,x,y,z,mesh_mode);
        for(int i = 0 ; i<n ; i++){
	  double F_x, F_y, F_z;
	  Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	  vx[i] += F_x/m*d2*dt;
	  vy[i] += F_y/m*d2*dt;
	  vz[i] += F_z/m*d2*dt;
	}

        for(int i = 0 ; i<n ; i++){
	   x[i] += vx[i]*c3*dt;
           y[i] += vy[i]*c3*dt;
           z[i] += vz[i]*c3*dt;
	}

	mesh(rho,x,y,z,mesh_mode);
        for(int i = 0 ; i<n ; i++){
	  double F_x, F_y, F_z;
	  Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	  vx[i] += F_x/m*d3*dt;
	  vy[i] += F_y/m*d3*dt;
	  vz[i] += F_z/m*d3*dt;
	}

	for(int i = 0 ; i<n ; i++){
	   x[i] += vx[i]*c4*dt;
           y[i] += vy[i]*c4*dt;
           z[i] += vz[i]*c4*dt;
	}
    }

    //RK4 mode
    if (OI_mode == 3){
       double **kr1 = new double *[n];
       double **kv1 = new double *[n];
       double **kr2 = new double *[n];
       double **kv2 = new double *[n];
       double **kr3 = new double *[n];
       double **kv3 = new double *[n];
       double **kr4 = new double *[n];
       double **kv4 = new double *[n];
       
       for(int i = 0 ; i<n ; i++){
          kr1[i] = new double [3];
	  kv1[i] = new double [3];
	  kr2[i] = new double [3];
          kv2[i] = new double [3];
	  kr3[i] = new double [3];
          kv3[i] = new double [3];
	  kr4[i] = new double [3];
          kv4[i] = new double [3];
	}	  
       for(int i = 0 ; i<n ; i++){
          
          double F_x, F_y, F_z;
	  mesh(rho,x,y,z,mesh_mode);
          Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z,mesh_mode);
	  kr1[i][0] = vx[i];
          kr1[i][1] = vy[i];
	  kr1[i][2] = vy[i];
          kv1[i][0] = F_x;
          kv1[i][1] = F_y;
	  kv1[i][2] = F_z;
	  kr2[i][0] = kr1[i][0] + kv1[i][0]*0.5*dt;
          kr2[i][1] = kr1[i][1] + kv1[i][1]*0.5*dt;
	  kr2[i][2] = kr1[i][1] + kv1[i][1]*0.5*dt;

       }
       for(int i = 0 ; i<n ; i++){
          
          double F_x, F_y, F_z;
	  mesh(rho,x,y,z,mesh_mode);
          Get_Force_of_Particle(U,x[i]+kr1[i][0]*0.5*dt,y[i]+kr1[i][1]*0.5*dt,z[i]+kr1[i][2]*0.5*dt,F_x,F_y,F_z,mesh_mode);
          kv2[i][0] = F_x;
          kv2[i][1] = F_y;
          kv2[i][2] = F_z;
	  kr3[i][0] = kr1[i][0] + kv2[i][0]*0.5*dt;
          kr3[i][1] = kr1[i][1] + kv2[i][1]*0.5*dt;
	  kr3[i][2] = kr1[i][1] + kv2[i][1]*0.5*dt;
       }
       for(int i = 0 ; i<n ; i++){

          double F_x, F_y, F_z;
	  mesh(rho,x,y,z,mesh_mode);
          Get_Force_of_Particle(U,x[i]+kr2[i][0]*0.5*dt,y[i]+kr2[i][1]*0.5*dt,z[i]+kr2[i][2]*0.5*dt,F_x,F_y,F_z,mesh_mode);
          kv3[i][0] = F_x;
          kv3[i][1] = F_y;
          kv3[i][2] = F_z;
          kr4[i][0] = kr1[i][0] + kv3[i][0]*dt;
          kr4[i][1] = kr1[i][1] + kv3[i][1]*dt;
          kr4[i][2] = kr1[i][1] + kv3[i][1]*dt;
       }
       for(int i = 0 ; i<n ; i++){

          double F_x, F_y, F_z;
	  mesh(rho,x,y,z,mesh_mode);
          Get_Force_of_Particle(U,x[i]+kr3[i][0]*dt,y[i]+kr3[i][1]*dt,z[i]+kr3[i][2]*dt,F_x,F_y,F_z,mesh_mode);
          kv4[i][0] = F_x;
          kv4[i][1] = F_y;
          kv4[i][2] = F_z;
       }
       for(int i = 0 ; i<n ; i++){
	  x[i] += (kr1[i][0]+2.0*kr2[i][0]+2.0*kr3[i][0]+1.0*kr4[i][0])/6.0*dt;
          y[i] += (kr1[i][1]+2.0*kr2[i][1]+2.0*kr3[i][1]+1.0*kr4[i][1])/6.0*dt;
          z[i] += (kr1[i][1]+2.0*kr2[i][2]+2.0*kr3[i][2]+1.0*kr4[i][2])/6.0*dt;

	  vx[i] += (kv1[i][0]+2.*kv2[i][0]+2.0*kv3[i][0]+1.0*kv4[i][0])/6.0*dt;
          vy[i] += (kv1[i][1]+2.*kv2[i][1]+2.0*kv3[i][1]+1.0*kv4[i][1])/6.0*dt;
          vz[i] += (kv1[i][1]+2.*kv2[i][1]+2.0*kv3[i][1]+1.0*kv4[i][1])/6.0*dt;
	}
   }
	
   t += dt;
   }   
   return EXIT_SUCCESS;
}


