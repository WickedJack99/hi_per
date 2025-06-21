PROGRAM_NAME = jacobi
SOURCE_FILES = main.cpp jacobi.cpp matrix_io.cpp

COMPILER_FLAGS = -std=c++17 -Wall -Wextra -Wall -Werror -Wno-sign-compare -Wno-unused-result -Wno-cast-function-type

default: debug

release: $(SOURCE_FILES)
	mpicxx $(SOURCE_FILES) -O3 -o $(PROGRAM_NAME) ${COMPILER_FLAGS}

debug: $(SOURCE_FILES)
	mpicxx $(SOURCE_FILES) -g -o $(PROGRAM_NAME) ${COMPILER_FLAGS}

run: release
	mpirun -np 4 ./$(PROGRAM_NAME) 3072 3e-6 100000

plot: plot.py
	python3 plot.py