PROGRAM_NAME = jacobi
SOURCE_FILES = jacobi_main.cpp jacobi_serial.cpp jacobi_mpi.cpp matrix_io.cpp

COMPILER_FLAGS = -std=c++17 -Wall -Wextra -Wall -Werror \
								 -Wno-sign-compare \
								 -Wno-unused-result \
								 -Wno-cast-function-type \
								 -Wno-unused-variable \
								 -Wno-unused-parameter 

default: debug

release: $(SOURCE_FILES)
	mpicxx $(SOURCE_FILES) -O3 -o $(PROGRAM_NAME) ${COMPILER_FLAGS}

debug: $(SOURCE_FILES)
	mpicxx $(SOURCE_FILES) -g -o $(PROGRAM_NAME) ${COMPILER_FLAGS}

run: release
	mpirun -np 4 ./$(PROGRAM_NAME)

plot: plot.py
	python3 plot.py