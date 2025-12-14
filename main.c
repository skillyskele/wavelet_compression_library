#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include "wavedec.h" // Ensure this header contains your wave_object and wt_object structs
#include "compression_types.h"
#include "dwt.h"




#define NUM_LEVELS 5 // should be floor(log2(SIGNAL_LENGTH))
#define TEST_WAVELET "db4"
#define NUM_CHANNELS 1
#define NUM_SAMPLES 533130


wave_object_t wave_pool[MAX_WAVE_OBJECTS];
int wave_pool_used[MAX_WAVE_OBJECTS] = {0}; // Track usage


wt_object_t wt_pool[MAX_WT_OBJECTS];
int wt_pool_used[MAX_WT_OBJECTS] = {0}; // Track usage



// Compute the Euclidean norm of all values in sparse_rep
COEFFICIENT_TYPE compute_energy(COEFFICIENT_TYPE *sparse_rep, int num_channels, int num_coefs)
{
    COEFFICIENT_TYPE energy = 0.0;
    for (int i = 0; i < num_channels; i++)
    {
        COEFFICIENT_TYPE row_energy = 0.0; // Compute energy per row
        for (int j = 0; j < num_coefs; j++)
        {
            row_energy += sparse_rep[i * num_coefs + j] * sparse_rep[i * num_coefs + j];
        }
        energy += row_energy; // Accumulate row-wise
    }
    return sqrt(energy); // Euclidean norm
}
// in the future, i'll unroll loops and align things for optimization
// later, we'll use the MVP API for fast convolutions and more


void read_eeg_signal(const char *filename, int signal_length, int num_channels, SAMPLE_TYPE* data)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Error opening file.\n");
        exit(1);
    }
    for (int i = 0; i < num_channels; i++)
    {
        for (int j = 0; j < signal_length; j++)
        {
            if (fscanf(file, "%f", &data[i*signal_length + j]) != 1) // the format specifier should match SAMPLE_TYPE
            {
                printf("Error reading channel %d sample %d\n", i, j);
                fclose(file);
                exit(2);
            }
        }
    }
    fclose(file);
}


typedef struct
{
    COEFFICIENT_TYPE bpp;
    COEFFICIENT_TYPE quant;
    int num_nnz;
    int signal_length;
    COEFFICIENT_TYPE *sparse_rep; // pointer to compressed array
    COEFFICIENT_TYPE *wc;         // pointer to original coefficients
    size_t rep_size;              // size of arrays
    int *lengths;                 // lengths of each level, similar to bookkeeping in matlab
    SAMPLE_TYPE *means;           // means used for demeaning each channel
} CompressResult;


// could take in inputs double cr, data, signal_length, num_levels, num_channels
// signal length can be a macro
// but we still pass it in!
// data is passed in as a pointer,though it's a 2d array...use data[i * signal_length + j] to access elements
// num_levels can also be a macro (it's equal to floor(log2(signal_length)) (confirm with matlab)
// num_channels can also be a macro, but it's still passed in when we call it.
// the data comes in from a 12 bit resolution ADC, so it's uint16_t type
CompressResult compress(wave_object wave, wt_object wave_transform, COEFFICIENT_TYPE cr, SAMPLE_TYPE* data, int signal_length, int num_levels, int num_channels)
{


    COEFFICIENT_TYPE tolerance = 0.1 * cr;




    // grab mean of each channel of the data
    // save the means used to demean the data later
    SAMPLE_TYPE means[num_channels];
    for (int i = 0; i < num_channels; i++) {
        SAMPLE_TYPE sum = 0;
        for (int j = 0; j < signal_length; j++) {
            sum += data[i*signal_length + j];
        }
        SAMPLE_TYPE mean = sum / signal_length;
        means[i] = mean;
        // demean the data
        for (int j = 0; j < signal_length; j++) {
            data[i*signal_length + j] -= mean;
        }    
    }


    size_t output_length;
    dwt(wave_transform, &data[0]);
    output_length = wave_transform->outlength;
    COEFFICIENT_TYPE sparse_rep[num_channels * output_length];


    // Store channel 0's coefficients
    for (int j = 0; j < output_length; j++)
    {
        sparse_rep[0 * output_length + j] = wave_transform->dwt_coeff[j];
    }


    // Process and store remaining channels
    // for (int i = 1; i < num_channels; i++)
    // {
    //     dwt(wave_transform, &data[i * signal_length]);
    //     for (int j = 0; j < output_length; j++)
    //     {
    //         sparse_rep[i * output_length + j] = wave_transform->dwt_coeff[j];
    //     }
    // }   // COMMENTED OUT SINCE WE ONLY HAVE ONE CHANNEL




    // this section of code sets up the binary search
    // mn is 1e-20 initially
    COEFFICIENT_TYPE mn = 1e-20;
    COEFFICIENT_TYPE mx = 0.0;
    // find max positive value in sparse_rep
    for (int i = 0; i < num_channels; i++) {
        for (int j = 0; j < output_length; j++) {
            if (sparse_rep[i * output_length + j] > mx) {
                mx = sparse_rep[i * output_length + j];  
            }
        }
    }


    COEFFICIENT_TYPE quant = mx;


    //COEFFICIENT_TYPE original_energy = compute_energy(sparse_rep, num_channels, output_length);


    bool searching = true;
   
    // flatten the temp_sparse_rep
    COEFFICIENT_TYPE temp_sparse_rep[num_channels * output_length]; // flattened version
   
    COEFFICIENT_TYPE q_temp;
    int32_t q_max;
    int num_nnz;
    int num_nnz_bits;


    COEFFICIENT_TYPE bpp; // bits per pixel
    int iterations = 0;


    // make a thing of 'codewords' that represent the non zero coefficients
    int32_t codewords[num_channels * output_length]; // should be large enough
    while (searching)
    {
        // quantize sparse_rep with current quant value
        q_max = 0;
        for (int i = 0; i < num_channels; i++) {
            for (int j = 0; j < output_length; j++) {
                q_temp = (sparse_rep[i * output_length + j]/quant); // sparse_rep./quant
                codewords[i * output_length + j] = (int32_t) q_temp; // fix(sparse_rep./quant), or temp_sparse_rep / quant
                temp_sparse_rep[i * output_length + j] = codewords[i * output_length + j] * quant; // multiply int32_t by COEFFICIENT_TYPE (float) to get float
                // q_max can be found here
                if (q_max < (int32_t) q_temp) {
                    q_max = round(q_temp); // include math.h
                }
            }  
        }


        num_nnz_bits = 0; // reset these
        num_nnz = 0;
        for (int i = 0; i < num_channels; i++) {
            for (int j = 0; j < output_length; j++) {
                if (temp_sparse_rep[i * output_length + j] != 0) {
                    num_nnz += 1;
                }
            }  
        }
        num_nnz_bits = num_nnz * ceil(log2(q_max) + 1); // bits needed to represent each non zero coefficient
        bpp = (COEFFICIENT_TYPE) num_nnz_bits / (num_channels * output_length);
       
        if (bpp > cr) {
            mn = quant;
            quant = (quant + mx) / 2.0;
        }
        if (bpp < cr) {
            mx = quant;
            quant = (quant + mn) / 2.0;
        }


        if (((cr - tolerance) < bpp) && (bpp < (cr + tolerance))) {
            searching = false;
        }


        iterations += 1;
       
        if (iterations > 50) {
            searching = false; // prevent infinite loop
        }


       
    }


    //COEFFICIENT_TYPE sparsity = (COEFFICIENT_TYPE)num_nnz / (num_channels * output_length);


    //COEFFICIENT_TYPE compressed_energy = compute_energy(temp_sparse_rep, num_channels, output_length);
    //COEFFICIENT_TYPE energy_ratio = compressed_energy / original_energy;




    COEFFICIENT_TYPE *wc = malloc(num_channels * output_length * sizeof(COEFFICIENT_TYPE));
    COEFFICIENT_TYPE *sparse = malloc(num_channels * output_length * sizeof(COEFFICIENT_TYPE));
    int *lengths = malloc((num_levels + 2) * sizeof(int)); // +2 for final approximation and length[0]
    SAMPLE_TYPE *means_out = malloc(num_channels * sizeof(SAMPLE_TYPE));
    memcpy(wc, sparse_rep, num_channels * output_length * sizeof(COEFFICIENT_TYPE)); // original coefficients
    memcpy(sparse, temp_sparse_rep, num_channels * output_length * sizeof(COEFFICIENT_TYPE)); // compressed coefficients
    memcpy(lengths, wave_transform->length, (num_levels + 2) * sizeof(int));
    memcpy(means_out, means, num_channels * sizeof(SAMPLE_TYPE));
    CompressResult result = {bpp, quant, num_nnz, signal_length, sparse, wc, num_channels * output_length, lengths, means_out};
    return result;
}



int main()
{
    // signal length varies from 16 to 20000 by step size 32
    // num levels is floor(log2(signal_length))
    // num channels is 1 for now
    int num_levels;
    int num_channels = 1;
    COEFFICIENT_TYPE cr = 0.05;


   
    FILE *out = fopen("results/compression_results.csv", "w");
    fprintf(out, "Signal Length,Bits Per Pixel,Quantization Step Size,Sparse Representation Size,Wavelet Coefficients Size,Number of Non-Zero Coefficients\n");
    for (int signal_length = 16; signal_length <= 20000; signal_length += 32) {
        //printf("Processing signal length: %d\n", signal_length);
        num_levels = floor(log2(signal_length));
        SAMPLE_TYPE *data = malloc(num_channels * signal_length * sizeof(SAMPLE_TYPE));
        read_eeg_signal("eeg.txt", signal_length, num_channels, data);
        wave_object wave = wave_init("db4");
        // check if wave is NULL
        if (wave == NULL)
        {
            printf("Error initializing wave object.\n");
            exit(1);
        }


        wt_object wave_transform = wt_init(wave, "dwt", signal_length, num_levels);
        // check if wave_transform is NULL
        if (wave_transform == NULL)
        {
            printf("Error initializing wave transform object.\n");
            exit(1);
        }


        CompressResult result = compress(wave, wave_transform, cr, data, signal_length, num_levels, num_channels); // RESULT MUST CONTAIN QUANT VALUE, REMEMBER, WE ACTUALLY AIM TO SEND THE QUANT VALUE!!!!

        // print all result fields for debugging
        // printf("Bits Per Pixel: %f\n", result.bpp);
        // printf("Quantization Step Size: %f\n", result.quant);
        // printf("Number of Non-Zero Coefficients: %d\n", result.num_nnz);
        // printf("Signal Length: %d\n", result.signal_length);
        // printf("Sparse Representation Size: %zu\n", result.rep_size);
        // for (int i = 0; i < 10 && i < result.rep_size; i++) { // print first 10 values
        //     printf("Sparse Rep[%d]: %f\n", i, result.sparse_rep[i]);
        // }
        // for (int i = 0; i < 10 && i < result.rep_size; i++) { // print first 10 values
        //     printf("Wavelet Coef[%d]: %f\n", i, result.wc[i]);
        // }
        // for (int i = 0; i < num_levels + 2; i++) {
        //     printf("Length[%d]: %d\n", i, result.lengths[i]);
        // }
        // for (int i = 0; i < num_channels; i++) {
        //     printf("Mean[%d]: %f\n", i, result.means[i]);
        // }


        // Save detailed arrays to separate files for each signal_length
        char sparse_fname[64], wc_fname[64], quant_fname[64], lengths_fname[64], means_fname[64];
        sprintf(sparse_fname, "results/sparse_rep/sparse_rep_%d.txt", signal_length);
        sprintf(wc_fname, "results/wc/wc_%d.txt", signal_length);
        sprintf(quant_fname, "results/quant/quant_%d.txt", signal_length);
        sprintf(lengths_fname, "results/lengths/lengths_%d.txt", signal_length);
        sprintf(means_fname, "results/means/means_%d.txt", signal_length);


        FILE *f_sparse = fopen(sparse_fname, "w");
        FILE *f_wc = fopen(wc_fname, "w");
        FILE *f_quant = fopen(quant_fname, "w");
        FILE *f_lengths = fopen(lengths_fname, "w");
        FILE *f_means = fopen(means_fname, "w");
        for (size_t i = 0; i < result.rep_size; i++)
        {
            fprintf(f_sparse, "%f\n", result.sparse_rep[i]);
            fprintf(f_wc, "%f\n", result.wc[i]);
        }

        fprintf(f_quant, "%f\n", result.quant);

        for (int i = 0; i < num_levels + 2; i++) {
            fprintf(f_lengths, "%d\n", result.lengths[i]);
        }

        for (int i = 0; i < num_channels; i++) {
            fprintf(f_means, "%f\n", result.means[i]);
        }

        fclose(f_sparse);
        fclose(f_wc);
        fclose(f_quant);
        fclose(f_lengths);
        fclose(f_means);


        // now put other results in a csv
       
        fprintf(out, "%d,%f,%f,\"1x%zu float\",\"1x%zu float\",%d\n",
        result.signal_length, result.bpp, result.quant,
        result.rep_size, result.rep_size, result.num_nnz);


        free(data);
        free(result.sparse_rep);
        free(result.wc);
        free(result.lengths);
        free(result.means);
        wave_free(wave);
        wt_free(wave_transform);
    }
   
   
    fclose(out);
    return 0;
   
   
}
