all:
	mpicxx -o main main.cpp -lblacs-openmpi -lscalapack-openmpi -lgfortran

run:
	mpirun -n 4 ./main

