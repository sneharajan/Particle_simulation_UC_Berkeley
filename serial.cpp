#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

#define IDX(i,j,NB) ((i)*(NB)+(j))
#define OUTBOUND(i,j,NB) ((i) < 0 || (i) >= (NB) || (j) < 0 || (j) >= (NB))
#define cutoff  0.01

inline int getBinID(particle_t & particle, int NB)
{
    return floor(particle.x / cutoff) * NB + floor(particle.y / cutoff); 
}

extern double size;

void clear_bins(int NB2, int * binMax)
{
    for (int k = 0; k < NB2; k ++)
    {
        binMax[k] = 0;
    }
}

// note: bins is passed by reference, not by value.
void assign_bins(int NB, int n, particle_t * particles, int ** bins, int * binMax)
{
    clear_bins(NB*NB, binMax);
    for (int k = 0; k < n; k ++)
    {
        int binid = getBinID(particles[k],NB);
        bins[binid][binMax[binid]++] = k;
    }
}

void apply_force_binned(int NB, int ** bins, int * binMax, particle_t * particles, particle_t &particle)
{
    int dstbinid = getBinID(particle, NB);
    int i = dstbinid / NB;
    int j = dstbinid % NB;
    int istart = max(i-1,0), iend = min(i+2,NB),
        jstart = max(j-1,0), jend = min(j+2,NB);
    for (int _i = istart; _i < iend; _i ++)
    { 
        for (int _j = jstart; _j < jend; _j ++)
        {
            int * neighborbin = bins[IDX(_i,_j,NB)];
            int pnum = binMax[IDX(_i,_j,NB)];
            for (int src = 0; src < pnum; src ++)
            { 
                apply_force(particle, particles[neighborbin[src]]);
            } 
        }
    }

}


int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n ); 
    init_particles( n, particles );

    // create bins: we will use this to store the index of objects
    // It is highly inefficient: we may well end up with a large vector
    int NB = floor( size / cutoff) + 1; // NB is the number of bins per edge

    int P_PER_NB = 100; // magical number: we assume that this is enough for each bin
    int * bins[NB*NB];
    int binMax[NB*NB];
    for (int i = 0; i < NB*NB; i ++) {
        bins[i] = new int[P_PER_NB];
        binMax[i] = 0;
    }

    assign_bins(NB, n, particles, bins, binMax);

    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    // this for loop can't be parallelized because future steps depend on
    // previous ones!
    for( int step = 0; step < NSTEPS; step++ )
    {
        // compute forces
        // we need to use blocks
        // first, initialize
        for (int k = 0; k < n; k ++)
        {
            particles[k].ax = particles[k].ay = 0;
            apply_force_binned(NB, bins, binMax, particles, particles[k]);
        }

        //
        //  move particles

        for( int i = 0; i < n; i++ )
        {
            move( particles[i] );
        }

        // update bins again
        assign_bins(NB, n, particles, bins, binMax);

        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    for (int i = 0; i < NB*NB; i ++) {
        delete bins[i];
    }

    if( fsave )
        fclose( fsave );

    return 0;
}
