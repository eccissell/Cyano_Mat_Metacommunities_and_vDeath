#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <sys/types.h>

// Define constants
#define MEDIUM_DISTURBANCE_PROB_GOATFISH 0.21
#define MEDIUM_DISTURBANCE_PROB_FISH 0.28
#define MEDIUM_DISTURBANCE_PROB_SENESCENCE 0.21
#define HOTSPOT_RADIUS 80          // Radius of hotspots
#define COLDSPOT_RADIUS 80         // Radius of coldspots
#define HOTSPOT_INTENSITY 0.2      // Added disturbance probability in hotspots
#define COLDSPOT_INTENSITY 0.2     // Reduced disturbance probability in coldspots
#define NUM_HOTSPOTS 78            // Number of hotspots
#define NUM_COLDSPOTS 78           // Number of coldspots
#define MEDIUM_DISPERSAL 0.001
#define MEDIUM_GROWTH_RATE 0.5

// Function prototypes
void get_neighbor(int *x2, int *y2, int x, int y, int GridSize, int direction);
double calculate_distance(int x1, int y1, int x2, int y2, int GridSize);
void output_grid(gsl_matrix *Space, int GridSize, int seed, int t);

int main(int argc, char *argv[]) {
    if (argc < 7) {
        fprintf(stderr, "Usage: %s GRate GridSize TimeSteps ProbDispersal seed [GridFile]\n", argv[0]);
        return 1;
    }

    // Parse input arguments
    double GRate = atof(argv[1]);              // Growth rate probability
    int GridSize = atoi(argv[2]);              // Grid size (length of one side)
    int TimeSteps = atoi(argv[3]);             // Number of timesteps
    double ProbDispersal = atof(argv[4]);      // Dispersal probability
    int seed = atoi(argv[5]);                  // Random seed
    const char *GridFile = (argc > 6) ? argv[6] : NULL; // Optional grid file

    // Initialize RNG
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, seed);

    // Initialize grid and abundance vector
    gsl_matrix *Space = gsl_matrix_calloc(GridSize, GridSize);
    gsl_vector *Ab = gsl_vector_calloc(2); // 0 = Sediment, 1 = Mat

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
        gsl_matrix_set_all(Space, 0); // Initialize Space with zeros
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

    // Create grids directory if it doesn't exist
    struct stat st = {0};
    if (stat("grids", &st) == -1) {
        mkdir("grids", 0700);
    }

    // Determine if we should output grids based on parameters and seed value
    int isMediumDispersal = (ProbDispersal == MEDIUM_DISPERSAL);
    int isMediumGrowthRate = (GRate == MEDIUM_GROWTH_RATE);
    int isSeedInSet = (seed % 10 == 0) && (seed >= 10) && (seed <= 100);
    int shouldOutputGrids = isMediumDispersal && isMediumGrowthRate && isSeedInSet;

    // Open output file for writing results
    FILE *outputFile = fopen("output_abundance.csv", "w");
    if (outputFile == NULL) {
        fprintf(stderr, "Error: Could not open output file.\n");
        return 1;
    }

    fprintf(outputFile, "TimeStep,MatAbundance,SedimentAbundance\n");
    fprintf(outputFile, "0,%d,%d\n", initial_mat_count, initial_sediment_count);
    printf("Time 0: MatAbundance = %d, SedimentAbundance = %d\n", initial_mat_count, initial_sediment_count);

    // Initialize hotspot and coldspot centers
    int hotspot_centers[NUM_HOTSPOTS][2];
    int coldspot_centers[NUM_COLDSPOTS][2];

    // Main simulation loop
    for (int t = 1; t <= TimeSteps; t++) {
        // Update hotspot and coldspot centers
        for (int h = 0; h < NUM_HOTSPOTS; h++) {
            hotspot_centers[h][0] = gsl_rng_uniform_int(rng, GridSize);
            hotspot_centers[h][1] = gsl_rng_uniform_int(rng, GridSize);
        }
        for (int c = 0; c < NUM_COLDSPOTS; c++) {
            coldspot_centers[c][0] = gsl_rng_uniform_int(rng, GridSize);
            coldspot_centers[c][1] = gsl_rng_uniform_int(rng, GridSize);
        }

        // Process each cell
        for (int x = 0; x < GridSize; x++) {
            for (int y = 0; y < GridSize; y++) {
                int spfocal = gsl_matrix_get(Space, x, y);

                // Dispersal logic
                if (spfocal == 0 && gsl_vector_get(Ab, 1) > 0) {
                    double DispProb = gsl_rng_uniform(rng);
                    if (DispProb < ProbDispersal) {
                        gsl_matrix_set(Space, x, y, 1);
                        gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) - 1);
                        gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) + 1);
                    }
                }

                // Growth logic
                else if (spfocal == 1) {
                    int x2, y2;
                    get_neighbor(&x2, &y2, x, y, GridSize, gsl_rng_uniform_int(rng, 8));
                    if (gsl_matrix_get(Space, x2, y2) == 0 && gsl_rng_uniform(rng) < GRate) {
                        gsl_matrix_set(Space, x2, y2, 1);
                        gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) - 1);
                        gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) + 1);
                    }
                }

                // Disturbance logic with hotspots and coldspots
                double disturb_prob_goatfish = MEDIUM_DISTURBANCE_PROB_GOATFISH;
                double disturb_prob_fish = MEDIUM_DISTURBANCE_PROB_FISH;
                double disturb_prob_senescence = MEDIUM_DISTURBANCE_PROB_SENESCENCE;

                // Check if the cell is in a hotspot
                int is_in_hotspot = 0;
                for (int h = 0; h < NUM_HOTSPOTS; h++) {
                    if (calculate_distance(x, y, hotspot_centers[h][0], hotspot_centers[h][1], GridSize) <= HOTSPOT_RADIUS) {
                        disturb_prob_goatfish += HOTSPOT_INTENSITY;
                        disturb_prob_fish += HOTSPOT_INTENSITY;
                        disturb_prob_senescence += HOTSPOT_INTENSITY;
                        is_in_hotspot = 1;
                        break;
                    }
                }

                // Check if the cell is in a coldspot
                if (!is_in_hotspot) {
                    for (int c = 0; c < NUM_COLDSPOTS; c++) {
                        if (calculate_distance(x, y, coldspot_centers[c][0], coldspot_centers[c][1], GridSize) <= COLDSPOT_RADIUS) {
                            disturb_prob_goatfish -= COLDSPOT_INTENSITY;
                            disturb_prob_fish -= COLDSPOT_INTENSITY;
                            disturb_prob_senescence -= COLDSPOT_INTENSITY;
                            break;
                        }
                    }
                }

                // Clamp probabilities
                disturb_prob_goatfish = fmax(0, fmin(1, disturb_prob_goatfish));
                disturb_prob_fish = fmax(0, fmin(1, disturb_prob_fish));
                disturb_prob_senescence = fmax(0, fmin(1, disturb_prob_senescence));

                // Apply disturbances
                if (spfocal == 1 && gsl_rng_uniform(rng) < disturb_prob_goatfish) {
                    gsl_matrix_set(Space, x, y, 0);
                    gsl_vector_set(Ab, 1, gsl_vector_get(Ab, 1) - 1);
                    gsl_vector_set(Ab, 0, gsl_vector_get(Ab, 0) + 1);
                }
            }
        }

        // Output the grid if conditions are met
        if (shouldOutputGrids) {
            output_grid(Space, GridSize, seed, t);
        }

        // Write abundances to output file
        fprintf(outputFile, "%d,%g,%g\n", t, gsl_vector_get(Ab, 1), gsl_vector_get(Ab, 0));

        // Print abundances to the console for debugging
        printf("Time %d: MatAbundance = %g, SedimentAbundance = %g\n",
               t, gsl_vector_get(Ab, 1), gsl_vector_get(Ab, 0));
    }

    gsl_matrix_free(Space);
    gsl_vector_free(Ab);
    gsl_rng_free(rng);
    fclose(outputFile);

    return 0;
}

// Helper function to calculate neighbor coordinates
void get_neighbor(int *x2, int *y2, int x, int y, int GridSize, int direction) {
    int dx[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
    int dy[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    *x2 = (x + dx[direction] + GridSize) % GridSize;
    *y2 = (y + dy[direction] + GridSize) % GridSize;
}

// Helper function to calculate distance
double calculate_distance(int x1, int y1, int x2, int y2, int GridSize) {
    int dx = abs(x1 - x2);
    int dy = abs(y1 - y2);
    if (dx > GridSize / 2) dx = GridSize - dx;
    if (dy > GridSize / 2) dy = GridSize - dy;
    return sqrt(dx * dx + dy * dy);
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
