PROGRAM_NAME = game_of_life
SOURCE_FILES = game_of_life.cpp \
							 utils.cpp \
							 patterns.cpp \
							 matrix_io.cpp \
							 super_grid.cpp

COMPILER_FLAGS = -std=c++17 -Wall -Wextra -Wall -Werror \
								 -Wno-sign-compare \
								 -Wno-unused-result \
								 -Wno-cast-function-type \
								 -Wno-unused-variable \
								 -Wno-unused-parameter 

default: debug

release: $(SOURCE_FILES)
	mpicxx $(SOURCE_FILES) main.cpp -O3 -o $(PROGRAM_NAME) ${COMPILER_FLAGS}

debug: $(SOURCE_FILES)
	mpicxx $(SOURCE_FILES) main.cpp -g -o $(PROGRAM_NAME) ${COMPILER_FLAGS}

test:
	mpicxx $(SOURCE_FILES) test.cpp -g -o game_of_life_test ${COMPILER_FLAGS}
	mpirun -np 16 --oversubscribe ./game_of_life_test

tests:
	mpicxx $(SOURCE_FILES) test.cpp -g -o game_of_life_test ${COMPILER_FLAGS}
	mpirun -np 1 --oversubscribe ./game_of_life_test

runserial:
	mpirun -np 1 ./$(PROGRAM_NAME) 1 1

run: release
	mpirun -np 16 --oversubscribe ./$(PROGRAM_NAME) 4 4

ani: animate.py
	python3 animate.py
