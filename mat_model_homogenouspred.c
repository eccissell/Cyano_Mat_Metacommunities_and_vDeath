
// first change the working directory:
// cd /Users/sophiemccoy/Dropbox/Documents/Manuscripts/Cissell/Mat\ Metapop
// then look at file list in the directory to check:
// ls

// To compile locally:
// gcc -Wall -I/usr/local/include -L/usr/local/lib  Mat_model_ecc_edit.c -o TestECC -lgsl -lgslcblas -lm -O3 -g
// The inputs in order are:
// SimulationFileName NSub DisturbFile GRate GridSize TimeSteps pDispersal randomseedgrid StartGrid.txt (test with .asc)
// Substrate IDs: (Mat(=1), Sediment (=0)
// Final grid size in mm - 5000000000 (250,000 x 20,000)
// GridSize input argument is one side of a square...
// EX: ./TestECC 2 testmatDP.txt 0.0 500 100 0.0 123 testGrid.txt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>         // For memset
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <sys/types.h>

// Define medium levels for parameters
#define MEDIUM_DISTURBANCE_PROB_GOATFISH 0.21
#define MEDIUM_DISTURBANCE_PROB_FISH 0.28
#define MEDIUM_DISTURBANCE_PROB_SENESCENCE 0.21
#define MEDIUM_DISPERSAL 0.001
#define MEDIUM_GROWTH_RATE 0.5

// Function prototypes
void get_neighbor(int *x2, int *y2, int x, int y, int GridSize, int direction);
void output_grid(gsl_matrix *Space, int GridSize, int seed, int t);

// Main function
int main(int argc, char *argv[]) {
    if (argc < 8) {
        fprintf(stderr, "Usage: %s NSub GRate GridSize TimeSteps ProbDispersal seed DisturbanceProbFile [GridFile]\n", argv[0]);
        return 1;
    }

    // Parse input arguments
    int NSub = atoi(argv[1]);                  // Number of substrate types
    double GRate = atof(argv[2]);              // Growth rate probability
    int GridSize = atoi(argv[3]);              // Grid size (length of one side)
    int TimeSteps = atoi(argv[4]);             // Number of timesteps
    double ProbDispersal = atof(argv[5]);      // Dispersal probability
    int seed = atoi(argv[6]);                  // Random seed
    const char *DisturbanceProbFile = argv[7]; // Disturbance probabilities file
    const char *GridFile = (argc > 8) ? argv[8] : NULL; // Optional grid file

    // Initialize RNG
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, seed);

    // Initialize grid and abundance vector
    gsl_matrix *Space = gsl_matrix_calloc(GridSize, GridSize);
    gsl_vector *Ab = gsl_vector_calloc(NSub);

    // Load grid from file (if provided)
    if (GridFile) {
        FILE *F = fopen(GridFile, "rb");
        if (F == NULL) {
            fprintf(stderr, "Error: Could not open grid file %s\n", GridFile);
            return 1;
        }
        gsl_matrix_fscanf(F, Space);
        fclose(F);
    } else {
        // Initialize Space with zeros (sediment) if no grid file is provided
        gsl_matrix_set_all(Space, 0);
    }

    // Calculate initial abundances
    int initial_mat_count = 0;
    int initial_sediment_count = 0;
    for (int i = 0; i < GridSize; i++) {
        for (int j = 0; j < GridSize; j++) {
            int cell = gsl_matrix_get(Space, i, j);
            if (cell == 1) {
                initial_mat_count++;
            } else if (cell == 0) {
                initial_sediment_count++;
            }
        }
    }
    gsl_vector_set(Ab, 0, initial_sediment_count);
    gsl_vector_set(Ab, 1, initial_mat_count);

    // Read disturbance probabilities from file
    gsl_vector *DisturbanceProbs = gsl_vector_calloc(3); // Assuming three disturbance types
    FILE *F = fopen(DisturbanceProbFile, "rb");
    if (F == NULL) {
        fprintf(stderr, "Error: Could not open disturbance probabilities file %s\n", DisturbanceProbFile);
        return 1;
    }
    gsl_vector_fscanf(F, DisturbanceProbs);
    fclose(F);

    // Extract the probabilities
    double disturb_prob_goatfish = gsl_vector_get(DisturbanceProbs, 0);
    double disturb_prob_fish = gsl_vector_get(DisturbanceProbs, 1);
    double disturb_prob_senescence = gsl_vector_get(DisturbanceProbs, 2);

    // Create grids directory if it doesn't exist
    struct stat st = {0};
    if (stat("grids", &st) == -1) {
        mkdir("grids", 0700);
    }

    // Determine if we should output grids based on parameters and seed value
    int isMediumDisturbance = (disturb_prob_goatfish == MEDIUM_DISTURBANCE_PROB_GOATFISH) &&
                              (disturb_prob_fish == MEDIUM_DISTURBANCE_PROB_FISH) &&
                              (disturb_prob_senescence == MEDIUM_DISTURBANCE_PROB_SENESCENCE);
    int isMediumDispersal = (ProbDispersal == MEDIUM_DISPERSAL);
    int isMediumGrowthRate = (GRate == MEDIUM_GROWTH_RATE);
    int isSeedInSet = (seed % 10 == 0) && (seed >= 10) && (seed <= 100);
    int shouldOutputGrids = isMediumDisturbance && isMediumDispersal && isMediumGrowthRate && isSeedInSet;

    // Open output file for writing results
    FILE *outputFile = fopen("output_abundance.csv", "w");
    if (outputFile == NULL) {
        fprintf(stderr, "Error: Could not open output file.\n");
        return 1;
    }

    // Write header and initial abundances
    fprintf(outputFile, "TimeStep,MatAbundance,SedimentAbundance\n");
    fprintf(outputFile, "0,%d,%d\n", initial_mat_count, initial_sediment_count);

    printf("Time 0: MatAbundance = %d, SedimentAbundance = %d\n", initial_mat_count, initial_sediment_count);

    // Main simulation loop
    for (int t = 1; t <= TimeSteps; t++) {
        // Determine if mats and sediments are present
        int mat_exists = (gsl_vector_get(Ab, 1) > 0);
        int sediment_exists = (gsl_vector_get(Ab, 0) > 0);

        // Loop over each cell in the grid
        for (int x = 0; x < GridSize; x++) {
            for (int y = 0; y < GridSize; y++) {
                int spfocal = gsl_matrix_get(Space, x, y);

                // Dispersal logic
                if (spfocal == 0 && mat_exists) {
                    double DispProb = gsl_rng_uniform(rng);
                    if (DispProb < ProbDispersal) {
                        gsl_matrix_set(Space, x, y, 1);
                        // Update abundances immediately
                        gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) - 1);
                        gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) + 1);
                    }
                }

                // Growth logic
                else if (spfocal == 1) {
                    int x2, y2;
                    get_neighbor(&x2, &y2, x, y, GridSize, gsl_rng_uniform_int(rng, 8));
                    int spneigh = gsl_matrix_get(Space, x2, y2);

                    if (spneigh == 0) {
                        double growth_rand = gsl_rng_uniform(rng);
                        if (growth_rand < GRate) {
                            gsl_matrix_set(Space, x2, y2, 1);
                            // Update abundances immediately
                            gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) - 1);
                            gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) + 1);
                        }
                    }
                }

                // Disturbance logic
                int HasSomethingHappened = 0;
                if (spfocal == 1) {
                    double disturb_rand;
                    // Goatfish disturbance
                    if (!HasSomethingHappened) {
                        disturb_rand = gsl_rng_uniform(rng);
                        if (disturb_rand < disturb_prob_goatfish) {
                            gsl_matrix_set(Space, x, y, 0);
                            // Update abundances immediately
                            gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
                            gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
                            HasSomethingHappened = 1;
                        }
                    }
                    // Fish disturbance
                    if (!HasSomethingHappened) {
                        disturb_rand = gsl_rng_uniform(rng);
                        if (disturb_rand < disturb_prob_fish) {
                            gsl_matrix_set(Space, x, y, 0);
                            // Update abundances immediately
                            gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
                            gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
                            HasSomethingHappened = 1;
                        }
                    }
                    // Senescence disturbance
                    if (!HasSomethingHappened) {
                        disturb_rand = gsl_rng_uniform(rng);
                        if (disturb_rand < disturb_prob_senescence) {
                            gsl_matrix_set(Space, x, y, 0);
                            // Update abundances immediately
                            gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
                            gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
                            HasSomethingHappened = 1;
                        }
                    }
                }

                // Check for negative abundances
                if (gsl_vector_get(Ab, 0) < 0 || gsl_vector_get(Ab, 1) < 0) {
                    fprintf(stderr, "Error: Negative abundance detected at time %d, position (%d,%d)\n", t, x, y);
                    fprintf(stderr, "MatAbundance: %g, SedimentAbundance: %g\n", gsl_vector_get(Ab, 1), gsl_vector_get(Ab, 0));
                    exit(1);
                }
            }
        }

        // Output the grid if conditions are met
        if (shouldOutputGrids) {
            output_grid(Space, GridSize, seed, t);
        }

        // Write abundances for the current timestep
        fprintf(outputFile, "%d,%g,%g\n", t, gsl_vector_get(Ab, 1), gsl_vector_get(Ab, 0));
        printf("Time %d: MatAbundance = %g, SedimentAbundance = %g\n", t, gsl_vector_get(Ab, 1), gsl_vector_get(Ab, 0));
    }

    // Free memory and close file
    gsl_matrix_free(Space);
    gsl_vector_free(Ab);
    gsl_vector_free(DisturbanceProbs);
    gsl_rng_free(rng);
    fclose(outputFile);

    return 0;
}

// Function to calculate neighbor coordinates using modular arithmetic
void get_neighbor(int *x2, int *y2, int x, int y, int GridSize, int direction) {
    int dx[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int dy[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    *x2 = (x + dx[direction] + GridSize) % GridSize;
    *y2 = (y + dy[direction] + GridSize) % GridSize;
}

// Function to output the grid to a file
void output_grid(gsl_matrix *Space, int GridSize, int seed, int t) {
    // Construct the filename
    char filename[256];
    snprintf(filename, sizeof(filename),
             "grids/grid_seed_%d_timestep_%d.txt",
             seed, t);

    // Open the file for writing
    FILE *gridFile = fopen(filename, "w");
    if (gridFile == NULL) {
        fprintf(stderr, "Error: Could not open grid file %s for writing.\n", filename);
        exit(1);
    }

    // Write the grid to the file
    for (int i = 0; i < GridSize; i++) {
        for (int j = 0; j < GridSize; j++) {
            int cell = gsl_matrix_get(Space, i, j);
            fprintf(gridFile, "%d", cell);
            if (j < GridSize - 1) {
                fprintf(gridFile, ",");
            }
        }
        fprintf(gridFile, "\n");
    }

    fclose(gridFile);
}
