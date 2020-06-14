// ==============
// PARTICLE MESH
// ==============

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<ctime>

<<<<<<< HEAD
int mesh_mode = 0; // 0: NGP ; 1: CIC ; 2: TSC
=======
//--------------------------------------------------------mode selection----------------------------------------------
int mesh_mode = 0; // 0: NGP ; 1: CIC ; 2: TSC
int OI_mode   = 0; //Orbit integration mode. 0: DKD 1:KDK 2:RK4

//-----------------------------------------------------------constants-------------------------------------------------
>>>>>>> Herby
double G = 1.0;
double Lx = 1.0, Ly = 1.0, Lz = 1.0; // domain size of 3D box
int Nx = 3, Ny = 3, Nz = 3; // # of grid points
double dx = Lx/(Nx-1), dy = Ly/(Ny-1), dz = Lz/(Nz-1);
int n = 1000; // # of particles
double m = 1.0; // particle mass

<<<<<<< HEAD
=======
//----------------------------------------------------------functions------------------------------------------------
void Get_Force_of_Particle(double ***F, double x, double y, double z, double &F_x, double &F_y, double &F_z);
>>>>>>> Herby
void mesh(double ***rho, double *x, double *y, double *z, int mode){
    int I, J, K;
    for(int i = 0 ; i<Nx ; i++){
        for(int j = 0 ; j<Ny ; j++){
            for(int k = 0 ; k<Nz ; k++){
                rho[i][j][k] = 0.;
            }
        }
    }
    if(mode==0){
        for(int p = 0 ; p<n ; p++){
            I = int((x[p]+0.5)/dx);
            J = int((y[p]+0.5)/dy);
            K = int((z[p]+0.5)/dz);
            if(abs(x[p]-I*dx)>abs(x[p]-(I+1)*dx)) I++;
            if(abs(x[p]-J*dy)>abs(x[p]-(J+1)*dy)) J++;
            if(abs(x[p]-K*dz)>abs(x[p]-(K+1)*dz)) K++;
            rho[I][J][K] += m/(dx*dy*dz);
        }
    }else if(mode==1){
        for(int p = 0 ; p<n ; p++){
            I = int((x[p]+0.5)/dx);
            J = int((y[p]+0.5)/dy);
            K = int((z[p]+0.5)/dz);
            for(int i = I ; i<=I+1 ; i++){
                for(int j = J ; j<=J+1 ; j++){
                    for(int k = K ; k<=K+1 ; k++){
                        rho[i][j][k] += m*(1.0-abs(x[p]+0.5-i*dx)/dx)*(1.0-abs(y[p]+0.5-j*dy)/dy)*(1.0-abs(z[p]+0.5-k*dz)/dz)/(dx*dy*dz);
                    }
                }
            }
        }
    }else if(mode==2){
        double fx, fy, fz;
        for(int p = 0 ; p<n ; p++){
            I = int((x[p]+0.5+.5*dx)/dx);
            J = int((y[p]+0.5+.5*dy)/dy);
            K = int((z[p]+0.5+.5*dz)/dz);
            for(int i = I-1 ; i<=I+1 ; i++){
                if(i==I) fx = 0.75 - pow( x[p]+0.5-i*dx , 2 ) / pow(dx,2);
                else fx = 0.5 * pow( 1.5-abs(x[p]+0.5-i*dx)/dx , 2);
                for(int j = J-1 ; j<=J+1 ; j++){
                    if(j==J) fy = 0.75 - pow( y[p]+0.5-j*dy , 2 ) / pow(dy,2);
                    else fy = 0.5 * pow( 1.5-abs(y[p]+0.5-j*dy)/dy , 2);
                    for(int k = K-1 ; k<=K+1 ; k++){
                        if(k==K) fz = 0.75 - pow( z[p]+0.5-k*dz , 2 ) / pow(dz,2);
                        else fz = 0.5 * pow( 1.5-abs(z[p]+0.5-k*dz)/dz , 2);
                        rho[i][j][k] += m*fx*fy*fz/(dx*dy*dz);
                    }
                }
            }
        }
    }
}


int main(){
	/* Variables */
    double t = 0.;
    double dt = 0.1; // time step
    double PDx = 0.2, PDy = 0.2, PDz = 0.2;
    double *x = new double [n];
    double *y = new double [n];
    double *z = new double [n];
    double *vx = new double [n];
    double *vy = new double [n];
    double *vz = new double [n];
    double ***rho = new double ** [Nx]; // mass density
    double ***U = new double ** [Nx]; // Dirichlet B.C.
    
	srand(time(NULL));
	/* Initialization */
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
<<<<<<< HEAD

	return 0;
}
=======
    
    //Starting to calculate force on partilces
    if (OI_mode == 0){ //DKD mode
	for(int i = 0 ; i<n ; i++){
           //drift: update position by 0.5*dt
           x[i] += vx[i]*0.5*dt;
           y[i] += vy[i]*0.5*dt;
           z[i] += vz[i]*0.5*dt;
	}
        //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
	for(int i = 0 ; i<n ; i++){ 
	    double F_x, F_y, F_z;
	    Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z);
	    vx[i] += F_x/m*dt;
	    vy[i] += F_y/m*dt;
	    vz[i] += F_z/m*dt;
	}
	for(int i = 0 ; i<n ; i++){
           //drift: update position by 0.5*dt
           x[i] += vx[i]*0.5*dt;
           y[i] += vy[i]*0.5*dt;
           z[i] += vz[i]*0.5*dt;
        }
    }

    if (OI_mode == 1){//KDK mode
       for(int i = 0 ; i<n ; i++){
	  //kick
	  double F_x, F_y, F_z;
	  Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z);
	  vx[i] += 0.5*F_x/m*dt;
          vy[i] += 0.5*F_y/m*dt;
          vz[i] += 0.5*F_z/m*dt;
	}
	//drift: update position by dt
	for(int i = 0 ; i<n ; i++){
	   x[i] += vx[i]*0.5*dt;
           y[i] += vy[i]*0.5*dt;
           z[i] += vz[i]*0.5*dt;
	}
	for(int i = 0 ; i<n ; i++){
	   //kick: calculate a(t+0.5*dt) and use that to update velocity by dt
	   double F_x, F_y, F_z;
	   Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z);
	   vx[i] += 0.5*F_x/m*dt;
           vy[i] += 0.5*F_y/m*dt;
           vz[i] += 0.5*F_z/m;
        }
    }
    if (OI_mode == 2){//RK4 mode
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
          Get_Force_of_Particle(U,x[i],y[i],z[i],F_x,F_y,F_z);
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
          Get_Force_of_Particle(U,x[i]+kr1[i][0]*0.5*dt,y[i]+kr1[i][1]*0.5*dt,z[i]+kr1[i][2]*0.5*dt,F_x,F_y,F_z);
          kv2[i][0] = F_x;
          kv2[i][1] = F_y;
          kv2[i][2] = F_z;
	  kr3[i][0] = kr1[i][0] + kv2[i][0]*0.5*dt;
          kr3[i][1] = kr1[i][1] + kv2[i][1]*0.5*dt;
	  kr3[i][2] = kr1[i][1] + kv2[i][1]*0.5*dt;
       }
       for(int i = 0 ; i<n ; i++){

          double F_x, F_y, F_z;
          Get_Force_of_Particle(U,x[i]+kr2[i][0]*0.5*dt,y[i]+kr2[i][1]*0.5*dt,z[i]+kr2[i][2]*0.5*dt,F_x,F_y,F_z);
          kv3[i][0] = F_x;
          kv3[i][1] = F_y;
          kv3[i][2] = F_z;
          kr4[i][0] = kr1[i][0] + kv3[i][0]*dt;
          kr4[i][1] = kr1[i][1] + kv3[i][1]*dt;
          kr4[i][2] = kr1[i][1] + kv3[i][1]*dt;
       }
       for(int i = 0 ; i<n ; i++){

          double F_x, F_y, F_z;
          Get_Force_of_Particle(U,x[i]+kr3[i][0]*dt,y[i]+kr3[i][1]*dt,z[i]+kr3[i][2]*dt,F_x,F_y,F_z);
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
	   
	return 0;
}

void Get_Force_of_Particle(double ***U, double x, double y, double z, double &F_x, double &F_y, double &F_z){

        int X_grid = (x/dx);
        int Y_grid = (y/dx);
        int Z_grid = (z/dx);
        printf("%d\n",X_grid);
        F_x = U[X_grid-2][Y_grid][Z_grid]/12. -  U[X_grid-1][Y_grid][Z_grid]*(2./3.)
                + U[X_grid+1][Y_grid][Z_grid]*(2./3.) - U[X_grid+2][Y_grid][Z_grid]*(1./12.);
        F_y = U[X_grid][Y_grid-2][Z_grid]/12. -  U[X_grid][Y_grid-1][Z_grid]*(2./3.)
                + U[X_grid][Y_grid+1][Z_grid]*(2./3.) - U[X_grid][Y_grid+2][Z_grid]*(1./12.);
        F_z = U[X_grid][Y_grid][Z_grid-2]/12. -  U[X_grid][Y_grid][Z_grid-1]*(2./3.)
                + U[X_grid][Y_grid][Z_grid+1]*(2./3.) - U[X_grid][Y_grid][Z_grid+2]*(1./12.);
	//Herby: This part is still not finished
}

>>>>>>> Herby
