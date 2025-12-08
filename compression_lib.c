#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include "wavedec.h" // Ensure this header contains your wave_object and wt_object structs

#define SIGNAL_LENGTH 16

#define NUM_LEVELS 4
#define TEST_WAVELET "db4"
#define NUM_CHANNELS 1
#define NUM_SAMPLES 533130


// Compute the Euclidean norm of all values in sparse_rep
double compute_energy(double **sparse_rep, int num_channels, int num_coefs)
{
    double energy = 0.0;
    for (int i = 0; i < num_channels; i++)
    {
        double row_energy = 0.0; // Compute energy per row
        for (int j = 0; j < num_coefs; j++)
        {
            row_energy += sparse_rep[i][j] * sparse_rep[i][j];
        }
        energy += row_energy; // Accumulate row-wise
    }
    return sqrt(energy); // Euclidean norm
}
// in the future, i'll unroll loops and align things for optimization
// later, we'll use the MVP API for fast convolutions and more

void read_eeg_signal(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error opening file.\n");
        exit(1);
    }
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        for (int j = 0; j < NUM_SAMPLES; j++)
        {
            if (fscanf(file, "%lf", &data[i][j]) != 1)
            {
                printf("Error reading channel %d sample %d\n", i, j);
                fclose(file);
                exit(2);
            }
        }
    }
    fclose(file);
}

// Main function
// could take in inputs double cr, data, signal_length, num_levels, num_channels
void compress()
{
    double cr = 0.001; // target compression ratio
    // copy test signal into data for now
    for (int i = 0; i < NUM_CHANNELS; i++) {
        for (int j = 0; j < SIGNAL_LENGTH; j++) {
            data[i][j] = test_signal[j];
        }
    }

    wave_object wave = wave_init("db4");
    wt_object wave_transform = wt_init(wave, "dwt", SIGNAL_LENGTH, NUM_LEVELS);

    double tolerance = 0.1 * cr;

    double means[NUM_CHANNELS];

    // grab mean of each channel of the data
    for (int i = 0; i < NUM_CHANNELS; i++) {
        double sum = 0.0;
        for (int j = 0; j < SIGNAL_LENGTH; j++) {
            sum += data[i][j];
        }
        double mean = sum / SIGNAL_LENGTH;
        means[i] = mean;  
        // demean the data
        for (int j = 0; j < SIGNAL_LENGTH; j++) {
            data[i][j] -= mean;
        }    
    }

    // now do dwt on each row of data
    for (int i = 0; i < NUM_CHANNELS; i++) {
        dwt(wave_transform, data[i]);
    }
    // print the length of wave_transform->dwt_coeff
    printf("DWT output length: %d\n", wave_transform->outlength);
    // print the wave_transform->dwt_coeff
    printf("DWT coefficients:\n");
    for (int i = 0; i < wave_transform->outlength; i++) {
        printf("%f ", wave_transform->dwt_coeff[i]);
    }
    printf("\n");

    // store each wave_transform->dwt_coeff into sparse_rep
    double **sparse_rep = (double **)malloc(NUM_CHANNELS * sizeof(double *));
    for (int i = 0; i < NUM_CHANNELS; i++) {
        sparse_rep[i] = (double *)malloc(wave_transform->outlength * sizeof(double));
        memcpy(sparse_rep[i], wave_transform->dwt_coeff, wave_transform->outlength * sizeof(double));
    }



    // this section of code sets up the binary search
    // mn is 1e-20 initially
    double mn = 1e-20;
    double mx = 0.0;
    // find max positive value in sparse_rep
    for (int i = 0; i < NUM_CHANNELS; i++) {
        for (int j = 0; j < wave_transform->outlength; j++) {
            if (sparse_rep[i][j] > mx) {
                mx = sparse_rep[i][j];  
            }
        }
    }
    double quant = mx;

    double original_energy = compute_energy(sparse_rep, NUM_CHANNELS, wave_transform->outlength);

    bool searching = true;
    double **temp_sparse_rep = (double **)malloc(NUM_CHANNELS * sizeof(double *));
    for (int i = 0; i < NUM_CHANNELS; i++)
    {
        temp_sparse_rep[i] = (double *)malloc(wave_transform->outlength * sizeof(double));
    }
    double q_temp;
    int num_nnz;
    int num_nnz_bits;
    double q_max;

    double bpp; // bits per pixel
    int iterations = 0;

    // make a thing of 'codewords' that represent the non zero coefficients
    int32_t codewords[NUM_CHANNELS][wave_transform->outlength]; // should be large enough
    while (searching)
    {
        // quantize sparse_rep with current quant value
        q_max = 0;
        printf("Current quant value: %g\n", quant);
        for (int i = 0; i < NUM_CHANNELS; i++) {
            for (int j = 0; j < wave_transform->outlength; j++) {
                q_temp = (int32_t)(sparse_rep[i][j]/quant);
                codewords[i][j] = q_temp;
                // print the codeword
                printf("%d ", codewords[i][j]);
                temp_sparse_rep[i][j] = q_temp * quant; 
                // q_max can be found here
                if (q_max < q_temp) {
                    q_max = q_temp;
                }
            }   
        }
        printf("\n");

        num_nnz_bits = 0; // reset these
        num_nnz = 0;
        for (int i = 0; i < NUM_CHANNELS; i++) {
            for (int j = 0; j < wave_transform->outlength; j++) {
                if (temp_sparse_rep[i][j] != 0) {
                    num_nnz += 1;
                }
            }   
        }
        num_nnz_bits = num_nnz * ceil(log2(q_max) + 1); // bits needed to represent each non zero coefficient
        bpp = (double) num_nnz_bits / (NUM_CHANNELS * wave_transform->outlength);
        // print num nnz bits and bpp
        printf("Num NNZ bits: %d, BPP: %f\n", num_nnz_bits, bpp);

        if (bpp > cr) {
            mn = quant;
            quant = (quant + mx) / 2.0;
        }
        if (bpp < cr) {
            mx = quant;
            quant = (quant + mn) / 2.0;
        }

        if (((cr - tolerance) < bpp) && (bpp < cr + tolerance)) {
            searching = false;
        }

        iterations += 1;
        
        if (iterations > 50) {
            searching = false; // prevent infinite loop
        }

        
    }

    // by now, num_nnz, bpp, quant are all set for the final temp_sparse_rep

    // print out my codewords
    for (int i = 0; i < NUM_CHANNELS; i++) {
        for (int j = 0; j < wave_transform->outlength; j++) {
            printf("%d ", codewords[i][j]);
        }
        printf("\n");
    }

    double sparsity = (double)num_nnz / (NUM_CHANNELS * wave_transform->outlength);

    double compressed_energy = compute_energy(temp_sparse_rep, NUM_CHANNELS, wave_transform->outlength);
    double energy_ratio = compressed_energy / original_energy;
/*
    printf("Sparsity: %g\n", sparsity);

    // Print bits per pixel
    printf("bpp: %g\n", bpp);

    // Print number of nonzero coefficients
    printf("#nnz: %d\n", num_nnz);

    // Print q_max
    printf("q_max: %g\n", q_max);

    // Print number of bits
    printf("Bits: %d\n", num_nnz_bits);

    // Print energy ratio and compression energy
    printf("Energy Ratio: %g\n", energy_ratio);
    printf("Original Energy: %g\n", original_energy);
    printf("Compressed Energy: %g\n", compressed_energy);

    // print sparse_rep
    printf("Final sparse representation after quantization:\n");
    for (int i = 0; i < NUM_CHANNELS; i++) {
        for (int j = 0; j < wave_transform->outlength; j++) {
            printf("%f ", temp_sparse_rep[i][j]);
        }
        printf("\n");
    }
*/


    // free allocated memory
    for (int i = 0; i < NUM_CHANNELS; i++) {
        free(sparse_rep[i]);
        free(temp_sparse_rep[i]);
    }
    free(sparse_rep);
    free(temp_sparse_rep);
}