/************************************************************************
 *
 * g2munge.c - Munges Gadget2 SnapFormat=1 IC files with 
 *             accelerations in desirable ways.  For use with 
 *             single files only, though extension to multiple files
 *             would be trivial to implement
 *
 * Copyright(c) 2014 Kevin Croker
 * Contains bits of code originally appearing in read_snapshot.c, Copyright(c) Vorkel Springel
 *    included with Gadget2 public distribution
 *
 * GPL v3
 * 
 * Functions:
 *  + merge : combine two snapshots
 *  + stats : give [min,max] of the physical quantites
 *  + xlate : Euclidian translation of a snapshot
 *  + dump  : Dump a list of positions in ASCII
 *  + ascii : Read positions in ASCII
 *  + grid  : Dump a grid of uniformly spaced squirticles
 *  + rand  : Randomly place random particles
 *  - scale : mechanically scale a dataset (is this useful?)
 *            Let l'/l = R
 *            Let k be the order of the potential (as in U(ar) = a^k U(r) )
 *            Then Landau doth spaketh (L1.10)
 *               v'/v = R^(k/2)
 *               E'/E = R^k
 *               M'/M (angular momentum) = R^(1+k/2)
 *            (NOTE: k = -1 for Newtonian gravity)
 *  ~ strip : remove everything EXCEPT a type of particle from the snapshot
 *  + flatten : set all particle types to a single particle type
 *  + binit : bin the particles within radius ranges of whatever
 *  + ident : makes an identical copy of the input (for debugging really)
 *  + burn  : imprint a density onto an initial condition assuming spherical symmetry
 *  + mass  : allows you to set particle type masses
 *
 * Consider eventual sick visualizations using OpenGL:https://sites.google.com/site/...
 *   ...dlampetest/python/vectorized-particle-system-and-geometry-shaders
 *
 *************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_rng.h>
#include "tkgadget2.h"
#include "burns.h"
 
// CAREFUL: Returns an allocated object, you have to free it explicitly!
// I think r is the origin...
gsl_histogram *binit(snapshot *snap, float radius, int nbins, float *r) {

  gsl_histogram *h;
  particle_data *P;

  h = gsl_histogram_alloc(nbins);
  gsl_histogram_set_ranges_uniform(h, 0.0, radius);

  for(P = snap->P + snap->NumPart; P > snap->P; --P)
    gsl_histogram_increment(h, mag3(r, P->Pos)); 

  return h;
}

gsl_histogram *binitk(snapshot *snap, float radius, int nbins, float *r) {

  gsl_histogram *h;
  particle_data *P;
  int res;
  double *ranges, *ksums;
  int k;

  h = gsl_histogram_alloc(nbins);
  
  ranges = (double *)malloc(sizeof(double *)*(nbins+1));
  ksums = (double *)malloc(sizeof(double *)*(nbins+1));
  
  //
  // From solving the recurrence relation
  //

  // Calculate the partial sums
  //  ksums[0] = 0;
  // for(k = 1; k <= nbins; ++k)
  // ksums[k] = ksums[k-1] + pow(k, 1.0/3.0);

  // Set up the ranges
  ranges[0] = 0.0;
  for(k = 1; k <= nbins; ++k)
    ranges[k] = pow(k, 1.0/3.0) * radius / pow(nbins, 1.0/3.0);

  // Clean up and assign ranges
  //free(ksums);
  gsl_histogram_set_ranges(h, ranges, nbins+1);
  free(ranges);

  //gsl_histogram_set_ranges_uniform(h, 0.0, radius);

  for(P = snap->P + snap->NumPart; P > snap->P; --P)
    gsl_histogram_increment(h, mag3(r, P->Pos)); 

  return h;
}

void dumpr(snapshot *snapA) {

  int total, k, i;
  particle_data *P;

  // Do Ids respect particle type??
  // Dump out the particles in different groups for cool colorization options
  total = 0;
  printf("# g2munge: %d\n", snapA->NumPart);
  
  for(i = 0; i < 6; ++i) {
    
    // Output empty groups so that when we plot in gnuplot, groups 
    // can be consistently indexed within the file
    //if(snapA->head.npart[i] == 0)
    //  continue;
    
    printf("\n# Group %d: %d\n", i, snapA->head.npart[i]);
    for(k = 0; k < snapA->head.npart[i]; ++k) {
	
      P = snapA->sortedP[total];
      //
      // yes, its slow to do this check every time.  If necessary cleave into two loops
      //
      printf("%u %f %f %f %f %f %f %f %f\n",
	     P->Id,
	     P->Pos[0], P->Pos[1], P->Pos[2], 
	     P->Vel[0], P->Vel[1], P->Vel[2], 
	     snapA->head.mass[i] ? snapA->head.mass[i] : P->Mass, 
	     !i ? P->U : -1);
      ++total;
    }
  }
}

void extrema(snapshot *snap, float *eP, float *eV, float *eA) {

  particle_data *P;
  int k;

  // XXX won't work, minima need to be initiall +NAN, maxima need to be -NAN
  for(k = 0; k < 3; ++k) {
    
    // Set maxima
    eP[2*k] = -INFINITY;
    eP[2*k + 1] = INFINITY;
    
    eV[2*k] = -INFINITY;
    eV[2*k + 1] = INFINITY;
    
    eA[2*k] = -INFINITY;
    eA[2*k + 1] = INFINITY;
  }
 
  for(P = snap->P + snap->NumPart; P > snap->P; --P) {
      
    // Do min-max
    for(k = 0; k < 3; ++k) {
      P->Pos[k] > eP[2*k] ? eP[2*k] = P->Pos[k] : (P->Pos[k] < eP[2*k + 1] ? eP[2*k+1] = P->Pos[k] : 0);
      P->Vel[k] > eV[2*k] ? eV[2*k] = P->Vel[k] : (P->Vel[k] < eV[2*k + 1] ? eV[2*k+1] = P->Vel[k] : 0);
      P->Acc[k] > eA[2*k] ? eA[2*k] = P->Acc[k] : (P->Acc[k] < eA[2*k + 1] ? eA[2*k+1] = P->Acc[k] : 0);
    }
  }
}

const char *usage = "Usage: %s <snapshot filename | - (stdin)> <command string>\n";

int main(int argc, char **argv)
{
  char input_fname[200];
  int type, snapshot_number, files;
  float r[3], cm[3], L[3], E, vsq, totalMass;
  float maxPos[6], maxVel[6], maxAcc[6];
  int i, sanity, k;
  snapshot snapA, snapB, snapC;
  particle_data *P, *pA, *pB, *pC;
  snapshot *writeout;
  profileFunc burnf;
  gsl_histogram *strips, *realbins;
  gsl_rng *entropy;
  float Rmax, e, density;
  size_t index;
  float params[10];
  particle_data **pIter, **pIter2;
  int numprev;

  // Load in the parameters
  if(argc < 3) {

    fprintf(stderr, usage, argv[0]);
    exit(1);
  }
 
  // Hardcoding the numfiles are we?
  files = 1;			/* number of files per snapshot */

  // This will always be the source unless ascii
  // in which case these parameters are ignored
  //
  // If we don't load a snapshot, zero out snapA
  //
  if(strcmp(argv[2], "ascii") && strcmp(argv[2], "grid") && strcmp(argv[2], "rand") && strcmp(argv[2], "nodetest")) {
    
    // Stop doing this name basename shit, makes it really hard to automate
    sprintf(input_fname, "%s", argv[1]);
    load_snapshot(&snapA, input_fname, files);
  }
  else
    memset(&snapA, 0, sizeof(snapshot));

  // Initially, don't write anything out
  writeout = NULL;

  // unit_conversion();		/* optional stuff */

  // Parse the command string
  if(!strcmp(argv[2], "burn")) {

    // 3      4          
    // burn <resolution> [args]
    // tophatburn: [args] = <cut> <0|1>
    //
    // Prunes particles until a region (0, radius) is overdense by some specified function
    // Function is specified by function pointer
    // Presently assumes that the input snapshot's distribution is uniform
    // but this isn't hard to generalize (use binit)
    //
    if(argc < 5) {
      
      fprintf(stderr, "Usage: burn <resolution> <mask> [args]\nPresent masks are:\n\ttophat\n\tthin\n");
      exit(2);
    }
    
    // Get the resolution
    e = atof(argv[3]);
    fprintf(stderr, "burn/NOTIFY: read resolution of %f\n", e);

    // Determine extent of the dataset
    // This will upperbound the radius
    extrema(&snapA, maxPos, maxVel, maxAcc);
    
    // Check negative extent as well!
    // Note, I wrote extrema to produce [max,min] which is retarded...
    for(i = 0; i < 3; ++i) {
      
      if(abs(maxPos[i+1]) > maxPos[i])
	maxPos[i] = abs(maxPos[i+1]);
    }

    Rmax = gsl_hypot3(maxPos[0], maxPos[2], maxPos[4]);
    fprintf(stderr, "burn/NOTIFY: computed Rmax of %f\n", Rmax);

    // Go one shell further out
    //Rmax += e;

    // Define your function
    // Note the function is 1 where things are to be kept
    // Note the function is 0 where things are to be completely cut
    if(!strcmp(argv[4], "thin")) {
      
      if(argc != 7) {
	
	fprintf(stderr, "Usage: burn/thin: <radius> <factor between 0.0 and 1.0>\n");
	exit(3);
      }

      params[0] = atof(argv[5]);
      params[1] = atof(argv[6]);
      
      fprintf(stderr, "burn/thin/NOTIFY: thinning beyond %f by a factor of %f\n", params[0], params[1]);
      burnf = thinBurn;
    }
    else if(!strcmp(argv[4], "tophat")) {
      
      if(argc != 7) {

	fprintf(stderr, "Usage: burn/tophat: <radius> <0/1: 1 to invert the tophat>\n");
	exit(3);
      }

      params[0] = atof(argv[5]);
      params[1] = atof(argv[6]);

      burnf = tophatBurn;
    }

    // Cut into shells
    strips = gsl_histogram_alloc(ceil(Rmax/e));
    gsl_histogram_set_ranges_uniform(strips, 0.0, Rmax+e);

    // Definitely need to make a first pass and compute the density in the regions!
    memset(r, 0, sizeof(float)*3);
    realbins = binit(&snapA, Rmax+e, ceil(Rmax/e), r);

    // Make sure to floor, or else we won't prune enough later on...
    for(i = 0; i < ceil(Rmax/e); ++i)
      strips->bin[i] = floor ( (1.0 - (*burnf)(i*e, params)) * realbins->bin[i]);

#undef DEBUG

#ifdef DEBUG
    // Do the bins sum to numpart?
    for(i = 0, sanity = 0; i < ceil(Rmax/e); ++i)
      sanity += (int)realbins->bin[i];
    
    fprintf(stderr, "burn: realbins contain %d particles (of %d)\n", sanity, snapA.NumPart);

    // Print off the deltas
    for(i = 0; i < ceil(Rmax/e); ++i)
      fprintf(stderr, "burn: delta[%d] = %d\n", i, (int)(realbins->bin[i] - strips->bin[i]));
#endif

    // We don't need realbins anymore, free it
    gsl_histogram_free(realbins);

    // Use the sorted array and shuffle it
    // XXX Not seeding...
    entropy = gsl_rng_alloc(gsl_rng_taus);
    gsl_ran_shuffle(entropy, snapA.sortedP, snapA.NumPart, sizeof(particle_data *));
    gsl_rng_free(entropy);

#ifdef DEBUG
    // Do an N(N-1)/2 check of the consistency of the pointer list
    for(pIter = snapA.sortedP;
    	pIter < snapA.sortedP + snapA.NumPart;
    	++pIter) {

      for(pIter2 = pIter+1;
    	  pIter2 < snapA.sortedP + snapA.NumPart;
    	  ++pIter2) {

    	// Check for duplicate
    	if(*pIter == *pIter2)
    	  fprintf(stderr, "burn: DUPLICATE IN SORTED LIST at %d\n", (int)(pIter2 - snapA.sortedP)/sizeof(particle_data *));
      }
    }
#endif

    // Keep track of actually pruned count
    k = 0;
    sanity = 0;
    
    // Prune
    for(pIter = snapA.sortedP, i = 0; 
	pIter < snapA.sortedP + snapA.NumPart; 
	++pIter, ++i) {

      // GO FETCH
      P = *pIter;

      e = mag3(P->Pos, r);
      gsl_histogram_find(strips, e, &index);

      //      if(gsl_histogram_get(strips, index) > 0.0) {
      //      fprintf(stderr, "burn: BEFORE strips[%d] = %d\n", index, (int)strips->bin[index]); 

      if(strips->bin[index] > 0.0) {
	// Writeout expects header statistics to honestly reflect eventual P array
	--snapA.head.npart[P->Type];

	// Was it massive?
	// (XXX possible) note that checking the mass directly won't work
	// if they are flagged massive in the header, but still given no mass here
	// Header valus is what matters.  Question is, were they particles with
	// a dynamic mass spec.  The answer is yes if the header says 0
	if(snapA.head.mass[P->Type] == 0.0)
	  --snapA.Nmassive;
	
	// Was it gassive?
	if(P->Type == 0)
	  --snapA.Ngas;

	// Flag for prune
	P->Type = -1;
	
	// Record that we pruned
	//gsl_histogram_accumulate(strips, e, -1.0);
	strips->bin[index] -= 1.0;

	//fprintf(stderr, "burn: AFTER strips[%d] = %d\n", index, (int)strips->bin[index]);
	// Update the count of prunes
	++k;
      }
#ifdef DEBUG
      else {
	
	// Particle was NOT removed, output its sex life
	fprintf(stderr, "burn: particle UNPRUNED @ (%f, %f, %f) with mag %f\n", P->Pos[0], P->Pos[1], P->Pos[2], e);
	fprintf(stderr, "burn:\t index %d, remaining here %d\n", index, (int)strips->bin[index]);
	++sanity;
      }
#endif
    }
    
#ifdef DEBUG
    fprintf(stderr, "burn: processed %d (of %d) particles\n\tskipped %d particles, pruned %d particles, sum %d\n", i, snapA.NumPart, sanity, k, sanity + k);

    // Did we empty all the strips?
    sanity = 0;
    for(i = 0; i < strips->n; ++i) {
      fprintf(stderr, "burn: strips->bin[%d] = %d\n", i, (int)strips->bin[i]);
      sanity += (int)strips->bin[i];
    }
    fprintf(stderr, "burn: strip residuals %d (should be zero)\n", sanity);
    fflush(stderr);
#endif

    gsl_histogram_free(strips);

    // Note that snapA.NumPart still contains the number of objects allocated within the P array
    if( !(pA = (particle_data *)malloc(sizeof(particle_data)*(snapA.NumPart - k)))) {

      fprintf(stderr, "burnr/DEATH: could not allocate prune buffer, shitting out the mouth...\n");
      exit(3);
    }

    // Would be better to just avoid these when writing out, instead of doing ridiculous memory copy
    // operations, but this is more straightforward right now and does not risk breaking working code
    // also, doing the writeout means comparing repeatedly over very many passes over the P struct
    // (one for each section)...
    // And id's have to be reset
    P = snapA.P;
    pC = pA;
    for(++P, i = 0; P <= snapA.P + snapA.NumPart; ++P) {

      if(P->Type < 0)
	continue;
      else {
	P->Id = ++i;
	memcpy(pC++, P, sizeof(particle_data));
      }
    }
  
    // Now, replace the P array
    free(snapA.P + 1);
    snapA.P = pA - 1;
    
#ifdef DEBUG
    fprintf(stderr, "burn/DEBUG: particles before %d\n", snapA.NumPart);
#endif

    // Now update the global number of particles
    snapA.NumPart -= k; 
   
#ifdef DEBUG
    fprintf(stderr, "burn/DEBUG: particles after %d\n", snapA.NumPart);
#endif
    
    // And update the npartTotal array which we always seem to forget about
    memcpy(snapA.head.npartTotal, snapA.head.npart, sizeof(int)*6);

    // Flag for writeout
    writeout = &snapA;

  }
  else if(!strcmp(argv[2], "retype")) {

    sanity = atoi(argv[3]);
    i = atoi(argv[4]);
    k = atoi(argv[5]);
    P = snapA.P + 1 + snapA.NumPart;
    
    fprintf(stderr, "retype: shifting %d particles of type %d into type %d...\n", sanity, i, k);
    while(--P > snapA.P && sanity > 0) {

      if(P->Type == i) {
	P->Type = k;
	sanity--;

	// Update headers
	snapA.head.npart[i]--;
	snapA.head.npartTotal[i]--;
	snapA.head.npart[k]++;
	snapA.head.npartTotal[k]++;
      }
    }

    fprintf(stderr, "retype: %d particles unshifted\n", sanity);

    // Update mass table
    snapA.head.mass[k] = snapA.head.mass[i];

    // Flag for writeout
    writeout = &snapA;
  }

  else if(!strcmp(argv[2], "rand")) {

    //
    // rand <#> <max l> <mass> [randoms?]
    // 
    // Produce a random arrangement of particles in a cube side length [0,l] 
    // with N particles
    //

    if(argc != 6 && argc != 8) {

      fprintf(stderr, "Usage: rand <#> <max l> <mass> [<random mod> <random base>]\n");
      exit(1);
    }
    
    snapA.NumPart = atoi(argv[3]);
    Rmax = atof(argv[4]);

    if( !(snapA.P = (particle_data *)malloc(sizeof(particle_data)*snapA.NumPart))) {
      
      exit(2);
    }
    
    fprintf(stderr, "rand: placing %d particles randomly in a box of size %f\n", snapA.NumPart, Rmax);

    // Zero it all out, and 1 index it 
    memset(snapA.P, 0, sizeof(particle_data)*snapA.NumPart);
    P = snapA.P;
    --snapA.P;
    sanity = 1;
    
    // KC 10/13/15
    // Don't forget to seed, dummy
    entropy = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(entropy, time(NULL));

    for(i = 0; i < snapA.NumPart; ++i) {

      P->Pos[0] = Rmax * gsl_rng_uniform(entropy);
      P->Pos[1] = Rmax * gsl_rng_uniform(entropy);
      P->Pos[2] = Rmax * gsl_rng_uniform(entropy);

      // Note that uniform types can be created with mod = 1, base = desired type 
      if(argc > 6) {
	P->Type = (rand() % atoi(argv[6])) + atoi(argv[7]);
      }
      else
	P->Type = 1;
      
      // Accumulate here in case of randoms
      snapA.head.npartTotal[P->Type]++;
      snapA.head.npart[P->Type]++;
      
      (P++)->Id = sanity++;
    }
    
    //gsl_rng_free(entropy);

    // Set other things in the header that need to be set
    for(i = 0; i < 6; ++i)
      snapA.head.mass[i] = atof(argv[5]);
    
    // Flag for writeout
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "nodetest")) {

    //
    // nodetest <#> <max l> <mass>
    //
    // Produce a nested layout of particles such that their locations fall in the center
    // of an oct tree of depth
    //
    snapA.NumPart = atoi(argv[3]);
    Rmax = atof(argv[4]);
    
    if(snapA.NumPart < 0)
      exit(1);

    if( !(snapA.P = (particle_data *)malloc(sizeof(particle_data)*snapA.NumPart))) {
      
      exit(2);
    }
    
    // Zero it all out, and 1 index it 
    memset(snapA.P, 0, sizeof(particle_data)*snapA.NumPart);

    i = snapA.NumPart;
    P = snapA.P;
    while(P != snapA.P + i) {
          
      P->Pos[0] = Rmax / 2.0;
      P->Pos[1] = Rmax / 2.0;
      P->Pos[2] = Rmax / 2.0;
      P->Type = 0;
      P->Id = (P - snapA.P) + 1;
      Rmax /= 2.0;
      
      ++P;
     }
    snapA.head.npartTotal[1] = i;
    snapA.head.npart[1] = i;
    snapA.head.mass[1] = 1.0;
   
    --snapA.P;
    
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "grid")) {

    //
    // grid <# density> <max l> <mass> [randoms?]
    // 
    // Produce a grid arrangement of particles in a cube side length [0,l] 
    // with N particles
    //
    
    // Initialize
    memset(&snapA, 0, sizeof(snapshot));

    if(argc != 6 && argc != 8) {

      fprintf(stderr, "Usage: grid <# density> <max l> <mass> [<random mod> <random base>]\n");
      exit(1);
    }

    // Initialize
    memset(&snapA, 0, sizeof(snapshot));
    
    // Determine the increment
    Rmax = atof(argv[4]);

    i = ceil(pow(atof(argv[3]) * gsl_pow_3(Rmax), 1.0/3.0));
    snapA.NumPart = gsl_pow_3(i);
    
    fprintf(stderr, "grid: populating with %d particles\n", snapA.NumPart);

    e = Rmax / i;
    fprintf(stderr, "grid: Rmax = %f, linear step is %e\n", Rmax, e);

    if( !(snapA.P = (particle_data *)malloc(sizeof(particle_data)*snapA.NumPart))) {
      
      exit(2);
    }
    
    // Zero it all out, and 1 index it 
    memset(snapA.P, 0, sizeof(particle_data)*snapA.NumPart);
    --snapA.P;

    
    // Create the particles
    P = snapA.P + 1;
    sanity = 1;
    
    //entropy = gsl_rng_alloc(gsl_rng_taus);
    // Check that the types are in the right range
    if(atoi(argv[6]) - 1 + atoi(argv[7]) > 5) {

      fprintf(stderr, "grid: ERROR %d + rand() %% %d exceeds permissable type of [0,5]\n", atoi(argv[7]), atoi(argv[6]));
      exit(3);
    }
    // XXX
    // This code makes all the types right next to each other!  I guess this was good for error testing!
    for(r[0] = 0.0; r[0] < Rmax; r[0] += e) {
      
      for(r[1] = 0.0; r[1] < Rmax; r[1] += e) {

	for(r[2] = 0.0; r[2] < Rmax; r[2] += e) {

	  // Bail if we've made enough
	  if(sanity > snapA.NumPart) {
	    
	    // Output particle deficit
	    fprintf(stderr, "grid: final dimensions filled (%f, %f, %f)\n", r[0], r[1], r[2]);
    	    r[1] = r[0] = Rmax + 1.0;
	    break;
	  }

	  memcpy(P->Pos, r, sizeof(float)*3);
	  if(argc > 6) {
	    P->Type = (rand() % atoi(argv[6])) + atoi(argv[7]);
	  }
	  else
	     P->Type = 1;
	  
	  // Accumulate here in case of randoms
	  snapA.head.npartTotal[P->Type]++;
	  snapA.head.npart[P->Type]++;

	  (P++)->Id = sanity++;
	}
      }
    }
    
    
    // Set other things in the header that need to be set
    for(i = 0; i < 6; ++i)
      snapA.head.mass[i] = atof(argv[5]);
      
    // Flag for writeout
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "binit")) {

    //
    // binit <radius> <bins> [x] [y] [z]
    //   
    // Note binit() returns the histogram because its useful, so we must explicitly free it
    if(argc == 5) {
     
      memset(r, 0, sizeof(float)*3);
      strips = binit(&snapA, atof(argv[3]), atoi(argv[4]), r);
    }
    else if(argc == 8) {
      
      r[0] = atof(argv[5]);
      r[1] = atof(argv[6]);
      r[2] = atof(argv[7]);
      strips = binit(&snapA, atof(argv[3]), atoi(argv[4]), r); 
    }
    else {
      
      fprintf(stderr, "Usage: binit <radius> <bins> [x] [y] [z]\n");
      exit(2);
    }

    // Go over the histogram and output the bins
    gsl_histogram_fprintf(stdout, strips, "%e", "%e");
    gsl_histogram_free(strips);
  }
  else if(!strcmp(argv[2], "binitk")) {

    //
    // binit <radius> <bins> [x] [y] [z]
    //   
    // Note binit() returns the histogram because its useful, so we must explicitly free it
    if(argc == 5) {
     
      memset(r, 0, sizeof(float)*3);
      strips = binitk(&snapA, atof(argv[3]), atoi(argv[4]), r);
    }
    else if(argc == 8) {
      
      r[0] = atof(argv[5]);
      r[1] = atof(argv[6]);
      r[2] = atof(argv[7]);
      strips = binitk(&snapA, atof(argv[3]), atoi(argv[4]), r); 
    }
    else {
      
      fprintf(stderr, "Usage: binitk <radius> <bins> [x] [y] [z]\n");
      exit(2);
    }

    // Go over the histogram and output the bins
    gsl_histogram_fprintf(stdout, strips, "%e", "%e");
    gsl_histogram_free(strips);
  }
  else if(!strcmp(argv[2], "ident")) {
    
    // Just write out what we read in
    writeout=&snapA;
  }
  else if(!strcmp(argv[2], "strip")) {

    // XXX Not yet implemented
    if(argc != 4) {

      fprintf(stderr, "Usage: strip <type to keep>\n");
      exit(2);
    }
  
    // Shallow copy snapshot A into B
    memcpy(&snapB, &snapA, sizeof(snapshot));
    
    // Allocate a new particle_data struct
    snapB.P = (particle_data *)malloc(sizeof(particle_data)*snapA.NumPart);
    if(!snapB.P) {

      fprintf(stderr, "strip: failed to allocate buffer, dying\n");
      exit(3);
    }

    // Zero out the particle count in the target
    snapB.NumPart = 0;
      
    // Signal to write out
    writeout = &snapB;
  }
  else if(!strcmp(argv[2], "flatten")) {

    if(argc != 4) {

      fprintf(stderr, "Usage: flatten <new type>\n");
      exit(2);
    }

    i = atoi(argv[3]);
    if(i < 0 || i > 5) {

      fprintf(stderr, "flatter/ERROR: type must be in (0,5), given %d\n", i);
      exit(4);
    }
    
    // Reset individual particle numbers
    memset(snapA.head.npart, 0, sizeof(int)*6);
    snapA.Nmassive = 0;
    snapA.Ngas = 0;
    memset(snapA.head.npartTotal, 0, sizeof(int)*6);

    // Tag all particles as new type i
    snapA.head.npart[i] = snapA.NumPart;
    snapA.head.npartTotal[i] = snapA.NumPart;

    // Did we become gas?
    if(i == 0)
      snapA.Ngas = snapA.NumPart;
      
    // Did we become massive?
    if(snapA.head.mass[i] == 0)
      snapA.Nmassive = snapA.NumPart;
    
    // Change each individual particle
    P = snapA.P;
    for(++P; P <= snapA.P + snapA.NumPart; ++P) {

      // Should really be a proper double comparison...
      if(snapA.head.mass[P->Type] > 0.0) {

	// The original type originally had uniform mass.  Does the destination have the same property?
	// If not, we stash the previous uniform mass into the new per-partice mass
	if(snapA.head.mass[i] > 0.0);
	else
	  P->Mass = snapA.head.mass[P->Type];
      }	 
      
      // XXX If we became a gas, we should probably acquire a temperature...
      if(i == 0)
	fprintf(stderr, "flatten/WARNING: becoming a gas, but U will not be set to anything sane!\n");
      
      // Change type
      P->Type = i;
    }

    // Flag to write it out
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "ascii")) {

    //
    // ascii
    // 
    // creates an IC and outputs to stdout for ascii values given in the dump format
    //
    // # g2munge: N\n
    // \n
    // # Group x: N_x\n
    // %u %f %f %f %f %f %f %f %f\n
    // ...
    // \n
    // # Group x+1: N_{x+1}\n
    // ...
    //
    // Where the floats appear as: id Xpos Ypos Zpos Xvel Yvel Zvel mass enerugi
    // i.e. delimited by space, with newlines denoting next particle
    //
    // (NOTE: this is also the format that is produced during a dump, so a dump can be read back 
    // into Gadget format if desired.)
    //
   
    // Signal to the abuser
    fprintf(stderr, "ascii/NOTIFY: ignoring listed input basefile and snapshot...\n");

    // Attempt to parse out the total particle counts
    if(!scanf("# g2munge: %d\n", &snapA.NumPart)) {
      
      fprintf(stderr, "ascii/ERROR: could not determine total particle count, header line bad?\n");
      exit(4);
    }
    else
      fprintf(stderr, "ascii/NOTIFY: searching for %d particles...\n", snapA.NumPart);
    
    // Allocate stuffs
    if( !(snapA.P = (particle_data *)malloc(snapA.NumPart * sizeof(particle_data)))) {

      fprintf(stderr, "ascii/ERROR: could not allocate my buffs!  death\n");
      exit(3);
    }
    
    // Adjust for weird indexing
    P = snapA.P;
    --snapA.P;
    
    // Number of previously read particles
    numprev = 0;
    
    while(!feof(stdin)) {

      // Read in a block "\n"
      if(scanf("\n# Group %d: %d\n", &k, &i) < 2) {
	
	fprintf(stderr, "ascii/ERROR: could not determine (easily) number or type\n");
	exit(4);
      }
      
      if(k < 0 || k > 5) {

	fprintf(stderr, "ascii/ERROR: unknown type %d\n", k);
	exit(4);
      }

      fprintf(stderr, "ascii/NOTIFY: Found type %d, reading %d of them...\n", k, i);
      
      // Set things
      snapA.head.npart[k] = i;
      snapA.head.npartTotal[k] = i;
      snapA.head.num_files = 1;

      // Gas check
      if(k == 0)
	snapA.Ngas = i;

      // Reuse variables:
      //  type = a flag for whether they were all the same mass
      type = 1;
      
      while(i-- > 0) {

	// Directly read em' in
	if( scanf("%u %f %f %f %f %f %f %f %f\n", 
		  &(P->Id), 
		  P->Pos, (P->Pos)+1, (P->Pos)+2,
		  P->Vel, (P->Vel)+1, (P->Vel)+2, 
		  &(P->Mass), &(P->U)) < 8 ) {
	  
	  fprintf(stderr, "ascii/ERROR: unable to read expected particle %d from the end\n", snapA.head.npart[k]-i);
	  exit(4);
	}
	
	// See if we should just tag these masses in the header
	// Short-circuit if we've already detected discrepancy
	if( (P > snapA.P + numprev) && type && (P-1)->Mass != P->Mass)
	  type = 0;
	
	// Got one, so go to the next particle
	++P;
      }

      // This is how many we read on this round
      numprev = P - snapA.P;

      // If we read anything for these guys...
      if(snapA.head.npart[k] > 0) {
	
	// If they were all the same type, set uniform mass in the header.  If
	// not, update the number of specific massive particles in the header
	if(type) {
	  snapA.head.mass[k] = (P-1)->Mass;
	  fprintf(stderr, "ascii/NOTIFY: masses for type %d set uniformly to %f\n", k, (P-1)->Mass);
	}
	else {
	  snapA.Nmassive += snapA.head.npart[k];
	  fprintf(stderr, "ascii/NOTIFY: particles of type %d have distinct masses\n", k);
	}
      }
    }
    
    // Flag for writeout
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "xlate")) {

    if(argc != 6) {

      fprintf(stderr, "Usage: xlate <x> <y> <z>\n");
      exit(2);
    }
          
    r[0] = atof(argv[3]);
    r[1] = atof(argv[4]);
    r[2] = atof(argv[5]);

    // Just do it here
    P = snapA.P;
    for(++P; P <= snapA.P + snapA.NumPart; ++P) {
      
      // For some reason, things are not stored as gsl_vectors.  Okay.
      P->Pos[0] += r[0];
      P->Pos[1] += r[1];
      P->Pos[2] += r[2];
    }
    
    // Signal to write it out
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "dump"))
    dumpr(&snapA);
  else if(!strcmp(argv[2], "merge")) {

    if(argc != 4) {
      
      fprintf(stderr, "Usage: merge <snapshot B filename>\n");
      exit(2);
    }
    else {
      
      // Hardcoding numfiles again, load up the B snapshot
      files = 1;			/* number of files per snapshot */
      //sprintf(input_fname, "%s_%03d", argv[3], atoi(argv[4]));
      load_snapshot(&snapB, argv[3], files);
  
      // Notify that we are preferentially using snapB's parameters
      fprintf(stderr, "merge/NOTIFY: merging into snapshot B, keeping B parameters as best can...\n");
      memcpy(&snapC.head, &snapB.head, sizeof(header));
      sanity = 0;

      // Adjust some header stuff
      snapC.NumPart = snapA.NumPart + snapB.NumPart;
      snapC.Ngas = snapA.Ngas + snapB.Ngas;
      snapC.Nmassive = snapA.Nmassive + snapB.Nmassive;

      // Increase BoxSize if necessary
      if(snapA.head.BoxSize > snapB.head.BoxSize) {
	fprintf(stderr, "merge/NOTIFY: enlargening box size from %f to %f\n", snapB.head.BoxSize, snapA.head.BoxSize);
	snapC.head.BoxSize = snapA.head.BoxSize;
      }
      
      // XXX: Allocate target memory!
      if(! (snapC.P = (particle_data *)malloc(sizeof(particle_data) * snapC.NumPart))) {
     	fprintf(stderr, "Could not allocate merge buffer, dying\n");
	exit(3);
      }

      // Be silly (but consistent)
      --snapC.P;

      // Note that particles are organized by type in each snapshot, so we can't just
      // block copy shit all at once, we must create a third snapshot.
      // and copy over the blocks one by one
      //
      // Zero index like normal people
      pC = snapC.P + 1;
      pA = snapA.P + 1;
      pB = snapB.P + 1;
      
      // Go through and increment all the ids of snapB by the number of particles in snapA
      for(i = 1; i <= snapB.NumPart; ++i)
	snapB.P[i].Id += snapA.NumPart;

      for(i = 0; i < 6; ++i) {

	// XXX HERE IS A BUG SOMEWHERE: MERGE ONLY GETS TYPES CORRECT IF DONE IN
	// SOME SPECIFIC ORDER.  SHOULD NOT MATTER
	//
	// This loop indexes over particle types
	// each header contains a number of particles of each given type
	memcpy(pC, pA, snapA.head.npart[i]*sizeof(particle_data));
	pC += snapA.head.npart[i];
	pA += snapA.head.npart[i];
	
	memcpy(pC, pB, snapB.head.npart[i]*sizeof(particle_data));
	pC += snapB.head.npart[i];
	pB += snapB.head.npart[i];

	// Set total particle count
	snapC.head.npart[i] = snapA.head.npart[i] + snapB.head.npart[i];

	// Check for the masses being the same
	if(snapB.head.mass[i] != snapA.head.mass[i])
	  fprintf(stderr, "merge/WARNING: mass %d NOT EQUAL (%f != %f), using %f\n",
		  i, snapA.head.mass[i], snapB.head.mass[i], snapB.head.mass[i]);
	
	// Declare the mass as that from snapB
	snapC.head.mass[i] = snapB.head.mass[i];
	     
	// Sanity check?
	sanity += snapC.head.npart[i];
      }
      
      // Sanity check
      if(sanity != snapC.NumPart) {

	fprintf(stderr, "merge/BUG: num parts %d not equal to computed num whole %d\n", snapC.NumPart, sanity);
	exit(4);
      }
      
      // Update nPartTotal!
      memcpy(&snapC.head.npartTotal, &snapC.head.npart, sizeof(int)*6);

      // XXX Should we care about particles that are now on top of each other?
      
      // Signal to write out
      writeout = &snapC;
    }
  }
  else if(!strcmp(argv[2], "mass")) {

    // 
    // mass
    //

    if(argc < 5) {
      
      fprintf(stderr, "Usage: mass <type> <value>\n");
      exit(3);
    }

    i = atoi(argv[3]);
    if(i < 0 || i > 5) {

      fprintf(stderr, "mass/ERROR: invalid type %d\n", i);
      exit(4);
    }

    // Update the global mass
    e = atof(argv[4]);
    snapA.head.mass[i] = e; 
    
    // Update Nmassive if we raped it
    if(e > 0.0)
      snapA.Nmassive -= snapA.head.npart[i];

    // Flag for writeout
    writeout = &snapA;
  }
  else if(!strcmp(argv[2], "stats")) {

    //
    // stats
    //
    memset(L, 0, sizeof(float)*3);
    memset(cm, 0, sizeof(float)*3);
    E = 0.0;
    totalMass = 0.0;
    
    // Go ahead and compute Newtonian energies and angular momenta (about the origin)
    P = snapA.P;
    printf("# Redshift: %e\n\n", snapA.head.redshift);
    printf("# Scalefactor: %e\n\n", snapA.head.time); 
    printf("# NumPart: %d\n\n", snapA.NumPart);
    for(i = 0; i < 6; ++i)
      printf("# npart[%d]: %d\n# mass[%d]: %e\n\n", i, snapA.head.npart[i], i, snapA.head.mass[i]);

    // Get the extrema
    extrema(&snapA, maxPos, maxVel, maxAcc);

    /* // Accumulate angular momentum */
    /* lm(&L, P); */
    /* } */

    /* // Finish up */
    /* for(k = 0; k < 3; ++k) */
    /*   cm[k] /= totalMass; */
    
    /* printf("# KE | Lx | Ly | Lz | CMx | CMy | CMz |\n"); */
    /* printf("%f %f %f %f %f %f %f\n", E, L[0], L[1], L[2], cm[0], cm[1], cm[2]); */
    /* printf("\n"); */
    printf("# Output in x_min x_max v_min v_max a_min a_max\n# For X, Y, Z\n"); 
    for(i =0; i < 3; ++i)
      printf("%e %e %e %e %e %e\n", maxPos[2*i + 1], maxPos[2*i], maxVel[2*i + 1], maxVel[2*i], maxAcc[2*i+1], maxAcc[2*i]);
  }
  else
    fprintf(stderr, "Unknown command: %s\n", argv[2]);
  
  if(writeout)
    dumpBinarySnapshot(writeout);
}

/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int unit_conversion(snapshot *snap)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;
  particle_data *P;

  /* physical constants in cgs units */
  GRAVITY = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm = 3.085678e21;	/*  code length unit in cm/h */
  UnitMass_in_g = 1.989e43;	/*  code mass unit in g/h */
  UnitVelocity_in_cm_per_s = 1.0e5;

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);


  Xh = 0.76;			/* mass fraction of hydrogen */
  HubbleParam = 0.65;

  P = snap->P;
  
  for(i = 1; i <= snap->NumPart; i++)
    {
      if(P[i].Type == 0)	/* gas particle */
	{
	  // XXX Don't trust me anymore!!
	  //MeanWeight = 4.0 / (3 * Xh + 1 + 4 * Xh * P[i].Ne) * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u = P[i].U * UnitEnergy_in_cgs / UnitMass_in_g;

	  gamma = 5.0 / 3;

	  /* get temperature in Kelvin */

	  P[i].Temp = MeanWeight / BOLTZMANN * (gamma - 1) * u;
	}
    }
}
