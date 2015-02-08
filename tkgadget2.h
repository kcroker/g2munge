#include <stdio.h>

#define PI 3.1415926

// Define a profileFunc function pointer type
typedef float (*profileFunc)(float, float *);

typedef struct {
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} header;

// Extended to include the acceleration data
typedef struct
{
  float Pos[3];
  float Vel[3];
  float Acc[3];
  unsigned int Id;
  float Mass;
  int Type;

  // Looks like Temp is just used as a convenience
  // relating U to some unit change
  float Temp;
  float Rho, U, Hsml;
  
} particle_data;

// Group everything together in a reasonable way
typedef struct
{

  header head;
  particle_data *P;
  particle_data **sortedP;
  int NumPart;
  int Ngas;
  int Nmassive;

} snapshot;

// Function prototypes
int load_snapshot(snapshot *snap, char *fname, int files);
void dumbBinarySnapshot(snapshot *snap);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
double mag3(float *q1, float *q2);
