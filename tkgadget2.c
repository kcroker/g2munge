#include "tkgadget2.h"
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

// Convenience function that avoids overflow
double mag3(float *q1, float *q2) {

  return gsl_hypot3(q1[0] - q2[0], q1[1] - q2[1], q1[2] - q2[2]); 
}

size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      fprintf(stderr, "I/O error (fwrite) on has occured.\n");
      fflush(stdout);
      //      endrun(777);
    }
  return nwritten;
}

// Compute an angular momentum with respect to (0,0,0)
void lm(float *L, particle_data *P) {
  
  L[0] += P->Mass * (P->Pos[1] * P->Vel[2] - P->Pos[2] * P->Vel[1]);
  L[1] += P->Mass * (-P->Pos[0] * P->Vel[2] + P->Pos[2] * P->Vel[0]);
  L[2] += P->Mass * (P->Pos[0] * P->Vel[1] - P->Pos[1] * P->Vel[0]);
   
}

/*
 * Dumps a snapshot to stdout in binary format, Gadget snapshot style #1.
 * Note that when we write out, we need to make 6 passes now to get
 * the types correct... (because that is how they are read back in)
 *
 * Not the most efficient, but suffices for present purposes
 * 
*/
void dumpBinarySnapshot(snapshot *snap) {

  particle_data *P;
  size_t blklen;
  int NumPart, N_gas;
  int i, dicks, k;
  int masscount;
  int dummy[3];

  P = snap->P;
  NumPart = snap->NumPart;

#ifdef DEBUG 
  fprintf(stderr, "Header says %d gas particles\n", snap->head.npart[0]); 
#endif

  // Probably encodes the length of the next block
#define BLKLEN my_fwrite(&blklen, sizeof(int), 1, stdout);

  // Write out the header
  blklen=sizeof(header);
  BLKLEN;
  my_fwrite(&snap->head, sizeof(header), 1, stdout);
  BLKLEN;

  blklen=NumPart*3*sizeof(float);

  BLKLEN;
  for(k = 0; k < 6; ++k) {

    for(i = 1; i <= NumPart; ++i) {
      
      if(P[i].Type == k)
	my_fwrite(P[i].Pos, sizeof(float)*3, 1, stdout);
    }
  }

  BLKLEN;

  BLKLEN;
  for(k = 0; k < 6; ++k) {
    
    for(i=1;i<=NumPart;i++) {
      
      if(P[i].Type == k)
	my_fwrite(P[i].Vel,sizeof(float)*3,1, stdout);
    }
  }
  BLKLEN;
 
  // Ids are tracked from 1
  // Never need to worry about the previous values, as these are IC
  // and so won't be used to mangle something in the middle of a run
  //
  // Perhaps the ID code is used in Father computations???
  blklen=NumPart*sizeof(int);
  BLKLEN;
  for(k = 0; k < 6; ++k) {
    
    for(i=1;i<=NumPart;i++) {

      if(P[i].Type == k)
	my_fwrite(&P[i].Id,sizeof(int),1,stdout);
    }
  }
  BLKLEN;

  blklen=snap->Nmassive * sizeof(float);
  // Note N massive is again sorted by type
  // Again, we need to bracket this region if there are massive particles
  if(snap->Nmassive) {
    
    BLKLEN;
    dicks = 0;
    for(k = 0; k < 6; ++k) {

      // If some type does not have variable masses, just skip it
      if(snap->head.mass[k] != 0)
	continue;

      // Do a pass over all particles, looking for this type
      for(i=1;i<=NumPart;i++) {

	if(P[i].Type == k) {
	  ++dicks;
	  my_fwrite(&P[i].Mass,sizeof(float),1,stdout);
	}
      }
    }
    BLKLEN;
    
    #ifdef DEBUG
    fprintf(stderr, "Wrote out %d masses, expected to write %d masses\n", dicks, snap->Nmassive);
    #endif
  }
  
  N_gas = snap->head.npart[0];
  if(N_gas > 0) {
    
    blklen=N_gas*sizeof(float);
    BLKLEN;
    
    for(i=1;i<=N_gas;i++) {
      // This must be deprecated or something...
      //dummy[0]=P[i].U;
      // Yeah, from looking at U, also culling from the same data structure now
      my_fwrite(&P[i].U, sizeof(float), 1, stdout);
      #ifdef DEBUG 
      //      fprintf(stderr, "Wrote out energy[%d]: %f\n", i, P[i].U);
      #endif
    }
    BLKLEN;

    /* // Unnecessary */
    /* //blklen=N_gas*sizeof(float);  /\* added density  *\/ */
    /* BLKLEN; */
    /* for(i=1;i<=N_gas;i++) { */
    /*   dummy[0]=P[i].Rho; */
    /*   my_fwrite(dummy,sizeof(float),1,stdout); */
    /* } */
    /* BLKLEN; */

      
    /* // XXX How does this line do anything?  THere was no #endif, or code to exec on that line... */
    
    /* //#istdoutef COOLING */
    /* // Maybe this is happening even without cooling */
    /* //if(snap->head.flag_cooling) { */
	
    /* //  blklen=N_gas*sizeof(float);  /\* electron abundance *\/ */
    /* BLKLEN; */
    /* for(i=1;i<=N_gas;i++) { */
    /*   dummy[0]= P[i].Hsml; */
    /*   fprintf(stderr, "P[i = %d].Hsml = %f\n", i, P[i].Hsml); */
    /*   my_fwrite(dummy,sizeof(float),1,stdout); */
    /* } */
    /* BLKLEN; */
    /* //} */
  }

  // Not required for initial onditions
  // Now add my accelerations, bitch
  // blklen = NumPart*3*sizeof(float);
  //BLKLEN;
  //for(i=1;i<=NumPart;i++) {
  //  my_fwrite(P[i].Acc,sizeof(float),3,stdout);
  // }
  //BLKLEN;
  /// XXX not in yet!!
}

/* 
 * This routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */

int load_snapshot(snapshot *snap, char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, j, k, dummy, ntot_withmasses;
  int t, n, off, pc, pc_new, pc_sph;
  particle_data *P;
  unsigned int *buffer;

  // Zero out the snapshot
  memset(snap, 0, sizeof(snapshot));
  
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new) {
    if(files > 1)
      sprintf(buf, "%s.%d", fname, i);
    else {
      if(!strcmp(fname, "-")) {
	fprintf(stderr, "g2munge/NOTIFY: null filename given, reading from stdin...\n");
	fd = stdin;
      }
      else
	sprintf(buf, "%s", fname);
    }

    if(strcmp(fname, "-") && !(fd = fopen(buf, "r"))) {
      fprintf(stderr, "can't open file `%s`\n", buf);
      exit(0);
    }
    
    fflush(stdout);

    // This code should suffice to read in the header...
    SKIP;
    fprintf(stderr, "Expecting HEADER block of length: %d (of %lu expected)\n", dummy, sizeof(header));
    fread(&snap->head, sizeof(header), 1, fd);
    SKIP;

    // Why is any of this code necessary?
    /* if(files == 1) */
    /* 	{ */
    /* 	  // This initializes NumPart = 0... */
    /* 	  for(k = 0, snap->NumPart = 0, ntot_withmasses = 0; k < 6; k++) */
    /* 	    snap->NumPart += snap->head.npart[k]; */
    /* 	  snap->Ngas = snap->head.npart[0]; */
    /* 	} */
    /* else */
    /* 	{ */
    /* 	  for(k = 0, snap->NumPart = 0, ntot_withmasses = 0; k < 6; k++) */
    /* 	    snap->NumPart += snap->head.npartTotal[k]; */
    /* 	  snap->Ngas = snap->head.npartTotal[0]; */
    /* 	} */

    // Set convenience variables for easy indexing
    snap->Ngas = snap->head.npart[0];
    #ifdef DEBUG 
    fprintf(stderr, "Nmassive before: %d\n", snap->Nmassive);
    #endif
    for(k = 0; k < 6; ++k) {
	
      snap->NumPart += snap->head.npart[k];
      if(snap->head.mass[k] == 0.0)
	snap->Nmassive += snap->head.npart[k];
    }
    #ifdef DEBUG 
    fprintf(stderr, "Nmassive after: %d\n", snap->Nmassive);
    #endif

    // The number of particles has been determined for this simulation
    if(i == 0) {
      
      // Just allocate here explicitly, no need to call another function
      if(!(snap->P = malloc(sizeof(particle_data) * snap->NumPart))) {
	
	fprintf(stderr, "Could not allocate %ld bytes for %s snaphot, dying\n", 
		sizeof(particle_data)*snap->NumPart, 
		buf);
	exit(2);
      }
      
      // They index from 1 for ... the science?
      --snap->P;

      // Assign shortcut so I don't need to mangle the code below
      P = snap->P;

      // Allocate the sorted list, which we will then just drop data into
      // This is zero indexed!
      if( !(snap->sortedP = (particle_data **)malloc(sizeof(particle_data *) * snap->NumPart))) {
	
	fprintf(stderr, "Could not allocate data for the sorted order, dyning\n");
	exit(2);
      }
    }

    SKIP;
    fprintf(stderr, "Expecting POS block of length: %d\nThis gives %lu entries (of %u expected)\n", dummy, dummy/(3*sizeof(float)), snap->NumPart);
     
    for(k = 0, pc_new = pc; k < 6; k++)	{
      for(n = 0; n < snap->head.npart[k]; n++) {
	fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	pc_new++;
      }
    }
    SKIP;

    SKIP;
    fprintf(stderr, "Expecting VEL block of length: %d\nThis gives %lu entries (of %u expected)\n", dummy, dummy/(3*sizeof(float)), snap->NumPart);
    for(k = 0, pc_new = pc; k < 6; k++)	{
      for(n = 0; n < snap->head.npart[k]; n++) {
	fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	pc_new++;
      }
    }
    SKIP;

    // Now I can look at the block
    SKIP;
    fprintf(stderr, "Expecting ID block of length: %d\nThis gives %u ids (of %u expected)\n", dummy, dummy/4, snap->NumPart);
    for(k = 0, pc_new = pc; k < 6; k++)	{
      for(n = 0; n < snap->head.npart[k]; n++)	{
	fread(&P[pc_new].Id, sizeof(unsigned int), 1, fd);
	if(P[pc_new].Id > snap->NumPart)
	  fprintf(stderr, "load/WARNING: id %d for %dth particle in block exceeds particle count\n", 
		  P[pc_new].Id, pc_new);
	else {
	  
	  // This creates a thing that will fuck FORTRAN-style optimization 
	  #ifdef DEBUG 
	  //fprintf(stderr, "Read id %d for particle %d.\n", P[pc_new].Id, pc_new);
	  #endif
	  snap->sortedP[P[pc_new].Id - 1] = &P[pc_new];
	  // So convenient though
	}
	pc_new++;
      }
    }
    SKIP;

    // This portion is bracketed
    // If there are massive particles, we need to skip 
    // because we need to read in the next block
    // But if not, we still want to initialize by looping over things

    if(snap->Nmassive > 0)
      SKIP;
    
    fprintf(stderr, "Expecting MASS block of length: %d\nThis gives %lu masses (of %d expected)\n", dummy, dummy/sizeof(float), snap->Nmassive);
   
    for(k = 0, pc_new = pc; k < 6; k++)	{
      for(n = 0; n < snap->head.npart[k]; n++) {
	P[pc_new].Type = k;
	
	// If mass is defined to be zero in the header, that means read it from the file.
	// Gotcha.
	#ifdef DEBUG 
	//fprintf(stderr, "head.mass[%d]=%f\n", k, snap->head.mass[k]);
	#endif
	if(snap->head.mass[k] == 0)
	  fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	else
	  P[pc_new].Mass = snap->head.mass[k];
	
	#ifdef DEBUG 
	//fprintf(stderr, "Mass read for type %d[%d] is: %f\n", k, n, P[pc_new].Mass); 
	#endif
	pc_new++;
      }
    }
    
    if(snap->Nmassive > 0)
      SKIP;

    // At this point, pc_new is set to the next location to make a particle
    // To keep from fucking that up, we start using a new variable pc_sph
    
    // This is the gas code...
    // It looks like they slammed the two previous arrays
    // into one.  The first particles in each block
    // area always the gas particles.... 
    if(snap->head.npart[0] > 0) {
      SKIP;
      fprintf(stderr, "Expecting ENERGY block of length: %u\nThis gives %lu energies (of %u expected)\n", dummy, dummy/sizeof(float), snap->head.npart[0]);
   
      for(n = 0, pc_sph = pc; n < snap->head.npart[0]; n++) {
	fread(&P[pc_sph].U, sizeof(float), 1, fd);
	#ifdef DEBUG 
	//	fprintf(stderr, "Read energy[%d]: %f\n", n+1, P[pc_sph].U);
	#endif
	pc_sph++;
      }
      SKIP;

      /* SKIP; */
      /* for(n = 0, pc_sph = pc; n < snap->head.npart[0]; n++) { */
      /* 	fread(&P[pc_sph].Rho, sizeof(float), 1, fd); */
      /* 	pc_sph++; */
      /* } */
      /* SKIP; */

      /* // Are the Ne the smoothing lengths? */
      /* // THEY ARE NOW! */
      /* // LOL Ne is no longer present.  */
      /* // Now its hsnl */
      /* //if(snap->head.flag_cooling) */
      /* //  { */
      /* SKIP; */
      /* for(n = 0, pc_sph = pc; n < snap->head.npart[0]; n++) { */
      /* 	fread(&P[pc_sph].Hsml, sizeof(float), 1, fd); */
      /* 	//fprintf(stderr, "P[pc_sph = %d].Hsnl = %f\n", pc_sph, P[pc_sph].Hsnl); */
      /* 	pc_sph++; */
      /* } */
      /* SKIP; */
    }
    
    // XXX where is smoothing length in the above code?
    // Okay, so the above was SPH stuff... if potential is turned off in Makefile, then we should be okay
    // Load accelerations in the stype 1 snapshot format
    if(feof(fd))
      fprintf(stderr, "load_snapshot/NOTIFY: No acceleration data found\n");
    else {
      
      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++) {
	for(n = 0; n < snap->head.npart[k]; n++)	{
	  // fread(P[pc_new].Acc, sizeof(float), 3, fd); is equivalent.....
	  // extra referencing and defreferencing is silly.
	  fread(&P[pc_new].Acc[0], sizeof(float), 3, fd);
	  pc_new++;
	}
      }
      SKIP;
    }

    // head.npartTotal was never being set!
    memcpy(snap->head.npartTotal, snap->head.npart, sizeof(int)*6);

    // Again, pc_new should be set to the right place to begin adding particles
    // from the next file...
    fclose(fd);
  }
}
