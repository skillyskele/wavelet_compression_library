import pywt
import numpy as np

# Example signal
signal = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])

# Perform multi-level discrete wavelet decomposition
wavelet = 'db4'  # Specify wavelet type
level = 3        # Specify number of levels for decomposition
coeffs = pywt.wavedec(signal, wavelet, level=level, mode='sym')  # Using symmetric padding

# Print coefficients
print("Wavelet Coefficients (Approximation + Details):")
for i, coef in enumerate(coeffs):
    level_type = "Approximation" if i == 0 else f"Detail Level {i}"
    print(f"{level_type}: {coef}")