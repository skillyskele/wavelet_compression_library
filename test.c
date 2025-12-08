#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wavedec.h" // Ensure this header contains your wave_object and wt_object structs

#define SIGNAL_LENGTH 16
#define NUM_LEVELS 3
#define TEST_WAVELET "db4"

// Custom test signal
static const double test_signal[SIGNAL_LENGTH] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

// Main function
int main()
{
    // Allocate testing wave_object (ensure wave_pool is well-defined elsewhere as static array)
    wave_object wave_test = wave_init(TEST_WAVELET); // Initialize with the db4 wavelet coefficients
    if (wave_test == NULL)
    {
        printf("Error initializing wavelet object.\n");
        return -1; // Exit if initialization fails
    }

    // Print wavelet filter coefficients - to verify initialization
    printf("Low-pass decomposition coefficients: \n");
    for (int i = 0; i < wave_test->lpd_len; i++)
    {
        printf("%f ", wave_test->lpd[i]);
    }
    printf("\n");

    printf("High-pass decomposition coefficients: \n");
    for (int i = 0; i < wave_test->hpd_len; i++)
    {
        printf("%f ", wave_test->hpd[i]);
    }
    printf("\n");

    printf("Low-pass reconstruction coefficients: \n");
    for (int i = 0; i < wave_test->lpr_len; i++)
    {
        printf("%f ", wave_test->lpr[i]);
    }
    printf("\n");

    printf("High-pass reconstruction coefficients: \n");
    for (int i = 0; i < wave_test->hpr_len; i++)
    {
        printf("%f ", wave_test->hpr[i]);
    }
    printf("\n");

    printf("Wavelet Name: %s\n", wave_test->wname);
    printf("Filter Length: %d\n\n", wave_test->filtlength);

    // Initialize DWT object
    wt_object wt_test = wt_init(wave_test, "dwt", SIGNAL_LENGTH, NUM_LEVELS);
    if (wt_test == NULL)
    {
        printf("Error initializing wavelet transform object.\n");
        return -1; // Exit if initialization fails
    }

    printf("DWT object initialized successfully.\n");
    printf("Input signal length: %d\n", wt_test->siglength);
    printf("Output signal length: %d\n", wt_test->outlength);
    printf("Number of decomposition levels: %d\n", wt_test->J);
    printf("Max iteration level: %d\n", wt_test->MaxIter);
    printf("Signal is even: %s\n\n", wt_test->even ? "Yes" : "No");

    // Print initialized array for testing
    printf("Output coefficient buffer initialized to zero:\n");
    for (int i = 0; i < wt_test->outlength; i++)
    {
        printf("%f ", wt_test->dwt_coeff[i]);
    }
    printf("\n");

    // Test filtcoef
    double lp1[wave_test->filtlength], hp1[wave_test->filtlength];
    double lp2[wave_test->filtlength], hp2[wave_test->filtlength];

    // Initialize filter coefficients from "db4"
    filtcoef(wave_test->filtlength, lp1, hp1, lp2, hp2);

    printf("Filters initialized:\n");
    printf("Low-pass decomposition filter (lp1): ");
    for (int i = 0; i < wave_test->filtlength; i++)
    {
        printf("%f ", lp1[i]);
    }
    printf("\n");

    printf("High-pass decomposition filter (hp1): ");
    for (int i = 0; i < wave_test->filtlength; i++)
    {
        printf("%f ", hp1[i]);
    }
    printf("\n");

    printf("Low-pass reconstruction filter (lp2): ");
    for (int i = 0; i < wave_test->filtlength; i++)
    {
        printf("%f ", lp2[i]);
    }
    printf("\n");

    printf("High-pass reconstruction filter (hp2): ");
    for (int i = 0; i < wave_test->filtlength; i++)
    {
        printf("%f ", hp2[i]);
    }
    printf("\n");

    // Test filter coefficients and transformations
    double reversed_signal[SIGNAL_LENGTH];
    copy_reverse(test_signal, SIGNAL_LENGTH, reversed_signal);

    printf("Testing reverse copy of input signal:\n");
    for (int i = 0; i < SIGNAL_LENGTH; i++)
    {
        printf("%f ", reversed_signal[i]);
    }
    printf("\n");

    double qmf_out[wave_test->filtlength];
    qmf_even(wave_test->lpd, wave_test->filtlength, qmf_out);

    printf("QMF (even transform) of low-pass decomposition filter:\n");
    for (int i = 0; i < wave_test->filtlength; i++)
    {
        printf("%f ", qmf_out[i]);
    }
    printf("\n");

    // alright it's time to test dwt baby

    dwt(wt_test, test_signal);
    printf("\nDWT completed successfully!\n");
    printf("DWT lengths (length[0..J], length[J+1]=orig): ");
    for (int k = 0; k <= NUM_LEVELS + 1; ++k)
    {
        printf("%d ", wt_test->length[k]);
    }
    printf("\n\n");

    // Print final approximation coefficients (stored at the beginning)
    printf("Final approximation coefficients (cA%d):\n", NUM_LEVELS);
    for (int i = 0; i < wt_test->length[0]; i++)
    {
        printf("%f ", wt_test->dwt_coeff[i]);
    }
    printf("\n\n");

    // Print detail coefficients per level
    // They are stored after the approximation coefficients
    int offset = wt_test->length[0]; // Start after approximation coefficients
    for (int k = 1; k <= wt_test->J; ++k)
    {
        printf("Detail coefficients level %d (cD%d):\n", wt_test->J - k + 1, wt_test->J - k + 1);
        for (int i = 0; i < wt_test->length[k]; i++)
        {
            printf("%f ", wt_test->dwt_coeff[offset + i]);
        }
        printf("\n\n");
        offset += wt_test->length[k]; // Move to next level
    }

    // Verify total length
    printf("Total output length: %d (expected: %d)\n", offset, wt_test->outlength);

    return 0; // Exit with success
}