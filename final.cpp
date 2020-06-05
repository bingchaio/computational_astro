#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <math.h>

//Particle Mesh

int main( int argc, char *argv[] )
{
// constant      
   // 3-D computational domain size L (3 real numbers). Assuming a 3D box.
   // number of equally spaced sampling points in 3D N (3 natural numbers)
   // time of update dt (real number)
   // number of particles n (natural number)
   const double pi = 3.1415926535897932384626433832795; //pi
   //Note: it will be easier to assume the box has the same length on each side.

// derived constant
   // spatial resolution d (3 real numbers, length/number of equally sampling points)

// initial condition
   //set t = 0.0
   //set the array of each particle, including its mass, 3D positions and velocities 
   //Note: it will be easier to assume each particle has the same mass 
   //randomly distributed the particles' positions and velocities as initial conditions
   //Note: Using the Virial theorem  to give an estimate on the velocities.
   //set the 3D array of potential u
   //initialize the potential

// define the 3D interpolating function

// define the 3D particle mesh function
   // Note: to get the potential on each grid by assigning a mass distribution function to each particle
   // 1. Nearest-Grid-Point (NGP)
   // delta function mass distribution
   // Let the position of the particle be xp, the position of the grid be xi.
   // If (xi-xp) and (yi-yp) and (zi-zp) are all smaller and equal than dx/2, dy/2, dz/2, then set the density of i grid to be mp/(dxdydz)
   // And sum over all the particles

   // 2. Cloud-In-Cell (CIC)
   // constant mass distribution
   // Let the position of the particle be xp, the position of the grid be xi.
   // Set the 3d width of the distribution as delta x,y,z.
   // 
    
3
   
//update
   while (cri>=1.0e-14 || t/dt>=1.0e6)
   {
	// update using SOR
	cri = 0.0;
	if (opt == 0)
	{
        #  pragma omp parallel
        {
  	    #  pragma omp for collapse (2)
	    for (int i=2; i<N-1 ;i+=2)
   	    {
                for (int j=2; j<N-1 ;j+=2)
                {
	            u_in[i][j]     = u[i][j];
		    u_in[i-1][j-1] = u[i-1][j-1];
		    u[i][j]        = u[i][j] + dt*( u[i+1][j] + u[i-1][j] 
                                   + u[i][j+1] + u[i][j-1] - 4.0*u[i][j]-den[i][j]*pow(dx,2.0))/pow(dx,2.0);
		    u[i-1][j-1]    = u[i-1][j-1] + dt*( u[i][j-1] + u[i-2][j-1] 
                                   + u[i-1][j] + u[i-1][j-2] - 4.0*u[i-1][j-1]-den[i-1][j-1]*pow(dx,2.0))/pow(dx,2.0);
		    du[i][j]       = abs(u[i][j]-u_in[i][j])/abs(u[i][j])/dt*pow(dx,2.0);
		    du[i-1][j-1]   = abs(u[i-1][j-1]-u_in[i-1][j-1])/abs(u[i-1][j-1])/dt*pow(dx,2.0);
                }
            }
	    #  pragma omp for collapse (2)
	    for (int i=1; i<N-2 ;i+=2)
   	    {
                for (int j=2; j<N-1 ;j+=2)
                {
	            u_in[i][j]     = u[i][j];
		    u_in[i+1][j-1] = u[i+1][j-1];
		    u[i][j]        = u[i][j] + dt*( u[i+1][j] + u[i-1][j] 
                                   + u[i][j+1] + u[i][j-1] - 4.0*u[i][j]-den[i][j]*pow(dx,2.0))/pow(dx,2.0);
		    u[i+1][j-1]    = u[i+1][j-1] + dt*( u[i+2][j-1] + u[i][j-1] 
                                   + u[i+1][j] + u[i+1][j-2] - 4.0*u[i+1][j-1]-den[i+1][j-1]*pow(dx,2.0))/pow(dx,2.0);
		    du[i][j]       = abs(u[i][j]-u_in[i][j])/abs(u[i][j])/dt*pow(dx,2.0);
		    du[i+1][j-1]   = abs(u[i+1][j-1]-u_in[i+1][j-1])/abs(u[i+1][j-1])/dt*pow(dx,2.0);
                }
            }
	}
        }
	if (opt == 1)
	{
	    #  pragma omp parallel for collapse (2)
	    for (int i=1; i<N-1 ;i++)
   	    {
                for (int j=1; j<N-1 ;j++)
                {
	            u_in[i][j]     = u[i][j];
		    u[i][j]        = u[i][j] + dt*( u[i+1][j] + u[i-1][j] 
                                   + u[i][j+1] + u[i][j-1] - 4.0*u[i][j]-den[i][j]*pow(dx,2.0))/pow(dx,2.0);
		    du[i][j]       = abs(u[i][j]-u_in[i][j])/abs(u[i][j])/dt*pow(dx,2.0);
                }
            }
	}
	if (opt == 2)
	{
	    for (int i=1; i<N-1 ;i++)
   	    {
                for (int j=1; j<N-1 ;j++)
                {
	            u_in[i][j]     = u[i][j];
		    u[i][j]        = u[i][j] + dt*( u[i+1][j] + u[i-1][j] 
                                   + u[i][j+1] + u[i][j-1] - 4.0*u[i][j]-den[i][j]*pow(dx,2.0))/pow(dx,2.0);
		    du[i][j]       = abs(u[i][j]-u_in[i][j])/abs(u[i][j])/dt*pow(dx,2.0);
                }
            }
	}
	for (int i=1; i<N-1 ;i++)
   	{
            for (int j=1; j<N-1 ;j++)
            {
	        cri += du[i][j];
            }
        }
	cri = cri/pow(double(N),2.0);
        t += dt;
	count++;
	if (count%1000 == 0)
	{
	     printf("%f %4e\n",t/dt, cri);
	}
   }
   double **du_err;
   du_err = new double *[N];
   for (int i=0;i<N;i++)
   {
       du_err[i] = new double [N];
   }
   for (int i=1; i<N-1 ;i++)
   {
       for (int j=1; j<N-1 ;j++)
       {
	   du_err[i][j] = abs(u[i][j]-u_ref[i][j]);
       }
   }
 
   double err = 0.0;
   for (int i=1; i<N-1 ;i++)
   {
       for (int j=1; j<N-1 ;j++)
       {
	   err += du_err[i][j];
       }
   }
   err = err/pow(double(N),2.0);	
   printf("%f %4e",t/dt, err);

   return EXIT_SUCCESS;
}
