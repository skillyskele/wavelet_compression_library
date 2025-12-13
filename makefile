# Makefile for Wavelet Transform Project

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -g

# Linker flags (libraries)
LDFLAGS = -lm

# Source files
SRCS = main.c wavedec.c dwt.c

# Header files
HEADERS = wavedec.h dwt.h compression_types.h

# Output executable
OUTPUT = wave_test

# Rule to build the output
$(OUTPUT): $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) -o $(OUTPUT) $(SRCS) $(LDFLAGS)

# Clean rule to remove the executable and any object files
clean:
	rm -f $(OUTPUT) *.o