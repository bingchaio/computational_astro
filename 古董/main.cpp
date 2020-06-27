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
   // If (xi-xp) and (yi-yp) and (zi-zp) are all smaller than dx/2, dy/2, dz/2, then set the density of i grid to be mp/(dxdydz)
   // And sum over all the particles

   // 2. Cloud-In-Cell (CIC)
   // constant mass distribution
   // Let the position of the particle be xp, the position of the grid be xi.
   // Set the 3d width of the distribution as delta x,y,z.
    

   //Herby: Writting for test
}

