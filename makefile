#-----------------------------------#
#	branchingtime makefile			#
#	Makefile for GNU C++ Compiler	#
#-----------------------------------#
CC = g++
EXE = branch.exe
CFLAGS = -Wall -lgsl -lgslcblas -lm
FILES = branch_vec.cpp fastmethod.cpp generalmethods.cpp
PNG = `freetype-config --cflags` -I/usr/local/include  -L/usr/local/lib -lpng -lpngwriter -lz -lfreetype

build: $(FILES)
	$(CC) $(CFLAGS) -o $(EXE) $(FILES) 


parallel: $(FILES)
	$(CC) $(CFLAGS) -o $(EXE) $(FILES) -fopenmp 

clean:
	rm -f *.o

# run make profile, then execute program, then run gprof -a branch_vec gmon.out
profile:  $(FILES)
	$(CC) $(CFLAGS) -pg -g -o $(EXE) $(FILES)

