
// first change the working directory:
// cd /Users/sophiemccoy/Dropbox/Documents/Manuscripts/Cissell/Mat\ Metapop
// then look at file list in the directory to check:
// ls

// To compile locally:
// gcc -Wall -I/usr/local/include -L/usr/local/lib  Mat_model_ecc_edit.c -o TestECC -lgsl -lgslcblas -lm -O3 -g
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ecissell/gsl/lib/
export LD_LIBRARY_PATH
// The inputs in order are:
// SimulationFileName NSub DisturbFile GRate GridSize TimeSteps pDispersal randomseedgrid StartGrid.txt (test with .asc)
// Substrate IDs: (Mat(=1), Sediment (=0)
// Final grid size in mm - 5000000000 (250,000 x 20,000)
// GridSize input argument is one side of a square...
// EX: ./TestECC 2 testmatDP.txt 0.0 500 100 0.0 123 testGrid.txt

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <fcntl.h>

int PrintMatrix(FILE * Where, gsl_matrix * M){
  int NR = M->size1;
  int NC = M->size2;
  int i, j;
  for (i = 0; i < NR; i++){
    for (j = 0; j < NC; j++){
		fprintf(Where, "%2.2f ", gsl_matrix_get(M, i, j));
    }
    fprintf(Where, "\n");
  }
  fprintf(Where, "\n");
  return 0;
}

int PrintGrid(FILE * Where, gsl_matrix * M){
  int NR = M->size1;
  int NC = M->size2;
  int i, j;
  for (i = 0; i < NR; i++){
    for (j = 0; j < NC; j++){
		fprintf(Where, "%d ", (int) gsl_matrix_get(M, i, j));
    }
    fprintf(Where, "\n");
  }
  fprintf(Where, "\n");
  return 0;
}

int PrintAb(FILE * Where, gsl_vector * M){
  int NR = M->size;
  int i;
  for (i = 0; i < NR; i++){
	  fprintf(Where, "%d ", (int) gsl_vector_get(M, i));
  }
  fprintf(Where, "\n");
  return 0;
}

int main(int argc, char *argv[]){

    //fprintf(stderr, "Number of arguments %d \n\n\n\n", argc);
  int NSub = atoi(argv[1]); // number of substrate types
  char * FileDP = argv[2]; // file storing the disturbance probabilities
  double GRate = atof(argv[3]); // mat local growth rate probability
  int GridSize = atoi(argv[4]); // size of the grid
  int TimeSteps = atoi(argv[5]); // number of steps (GridSize^2 predation
								 // events for each TimeStep)
  double ProbDispersal = atof(argv[6]); // 0.0 = use only
												// metacommunity; 1.0,
												// use only local abundances
  int seed = atoi(argv[7]); // random seed

  // read the disturbance probability vector
  gsl_vector * DP = gsl_vector_calloc(3);

  FILE * F;
  F = fopen(FileDP, "rb");
  gsl_vector_fscanf(F, DP);
  fclose(F);

  // Random number generator
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set (r, seed);

  // initialize the grid
  // initialize with empty space
  gsl_matrix * Space = gsl_matrix_calloc(GridSize, GridSize);

  // Initialize abundances
  gsl_vector * Ab = gsl_vector_calloc(NSub);

  int x, y, i, j, k, x2, y2, neigh, spfocal, spneigh;
    if (argc == 9){
        // read the grid from a file
        F = fopen(argv[8], "rb"); // the 9th argument (count starts with 0) is the file of the grid to read in
        gsl_matrix_fscanf(F, Space);
        fclose(F);
        for (i = 0; i < GridSize; i++){
            for (j = 0; j < GridSize; j++){
                k = (int) gsl_matrix_get(Space, i, j);
                gsl_vector_set(Ab, k, gsl_vector_get(Ab, k) + 1);
            }
        }
    }
    // Initialize STDOutput file
    //int file = open("Output/OutAbundance.txt", O_WRONLY | O_CREAT, 0777);
  // main loop
  for (i = 0; i < TimeSteps; i++){
	  if (i % 1 == 0) PrintAb(stdout, Ab);
	  for (j = 0; j < (GridSize * GridSize); j++){
		  //PrintAb(stderr, Ab);
		  // pick a random cell
		  x = gsl_rng_uniform_int(r, GridSize);
		  y = gsl_rng_uniform_int(r, GridSize);
		  // pick a random neighbor
		  neigh = gsl_rng_uniform_int(r, 8);
		  switch(neigh){
		  case 0:
			  x2 = x - 1;
			  y2 = y - 1;
			  break;
		  case 1:
			  x2 = x;
			  y2 = y - 1;
			  break;
		  case 2:
			  x2 = x + 1;
			  y2 = y - 1;
			  break;
		  case 3:
			  x2 = x - 1;
			  y2 = y;
			  break;
		  case 4:
			  x2 = x + 1;
			  y2 = y;
			  break;
		  case 5:
			  x2 = x - 1;
			  y2 = y + 1;
			  break;
		  case 6:
			  x2 = x ;
			  y2 = y + 1;
			  break;
		  case 7:
			  x2 = x + 1;
			  y2 = y + 1;
			  break;
		  }
		  //fprintf(stderr, "%d %d -> %d %d\n", x, y, x2, y2);
      // remove boundaries
		  if (x2 == (GridSize)) x2 = 0;
		  if (y2 == (GridSize)) y2 = 0;
		  if (x2 == -1) x2 = GridSize - 1;
		  if (y2 == -1) y2 = GridSize - 1;
		  // store the ID of focal and neighbor
		  spfocal = (int) gsl_matrix_get(Space, x, y);
		  spneigh = (int) gsl_matrix_get(Space, x2, y2);

		  // if the focal cell was empty, we can potentially disperse into it

		  if (spfocal == 0){
			  double DispProb = gsl_rng_uniform(r);
					  if (DispProb < ProbDispersal){
						  gsl_matrix_set(Space, x, y, 1);
						  // change abundances
						  gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) - 1);
						  gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) + 1);
					  }
				  }
		  // if the focal species was nonempty, we can disturb
      // we will add all disturbances as separate discrete logical events --
        // with probability read-in from a vector set as arg [2] --
        // if one of the disturbance events happen, none else will happen --
        // i.e. only one discrete disturbance event can occur to any given mat focal cell per time step then break
        //
      // We can also think about adding in a bounded random number generator for dist probs
        // to model heterogeneity in predation pressure
        // Then import final asc raster into Q to get landscape stats and compare
        // heterogenous vs. constant predation pressure effects on landscape patterns
		  int HasSomethingHappened = 0;
      if (spfocal == 1 && spneigh == 0){ // Start with a growth
        if(gsl_rng_uniform(r) < GRate){
          // the mat grows
          gsl_matrix_set(Space, x2, y2, 1);
          // change abundances
          gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) - 1);
          gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) + 1);
          spneigh=1;
        }
      }
		  if (spfocal == 1){ // First we remove via goatfish disturbance
			  double DisturbProbG = gsl_rng_uniform(r);
			  if (DisturbProbG < gsl_vector_get(DP, 0)){
				  gsl_matrix_set(Space, x, y, 0);
				  gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
				  gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
				  HasSomethingHappened = 1;
			  }
		  }
		  if (spfocal == 1 && HasSomethingHappened == 0){ // Next we remove via fish
        double DisturbProbF = gsl_rng_uniform(r);
			  if (DisturbProbF < gsl_vector_get(DP, 1)){
				  gsl_matrix_set(Space, x, y, 0);
				  gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
				  gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
				  HasSomethingHappened = 1;
			  }
		  }
      if (spfocal == 1 && HasSomethingHappened == 0){ // Next we remove via senescence (wholesale microscale processes
                                        //including viral predation and localized senescence)
        double DisturbProbS = gsl_rng_uniform(r);
			  if (DisturbProbS < gsl_vector_get(DP, 2)){
				  gsl_matrix_set(Space, x, y, 0);
				  gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
				  gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
				  HasSomethingHappened = 1;
        }
      }
  }
}
// save output
char OutFileName[1000];
FILE * OutFile;
sprintf(OutFileName, "Output/Out-%d-%s-%f-%d-%d-%0.2f-%d-Abundances.txt", NSub, FileDP, GRate, GridSize, TimeSteps, ProbDispersal, seed);
OutFile = fopen(OutFileName, "w");
PrintAb(OutFile, Ab);
fclose(OutFile);

//sprintf(OutFileName, "Output/Out-%d-%s-%f-%d-%d-%0.2f-%d-Grid.txt", NSub, FileDP, GRate, GridSize, TimeSteps, ProbDispersal, seed);
//OutFile = fopen(OutFileName, "w");
//PrintGrid(OutFile, Space);
//fclose(OutFile);

// free memory
gsl_vector_free(DP);
gsl_matrix_free(Space);
gsl_vector_free(Ab);
gsl_rng_free (r);
return 0;
}
