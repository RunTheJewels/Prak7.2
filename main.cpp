#include "scalapack.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <mpi.h>


/* 
 * Входные данные n - кол-во атомов
 * а - вектор длины n-1 амплитуд взаимодействия между соседними атомами
 * w - вектор длины n частот атомов
 * k - размерность стока
 * Emin Emax - минимальные и максимальные значения количества возбуждений (в сумме с энергией стока)
 * в выбираемых состояниях
 * 
 * */
 
 
using namespace std;

struct Matrix {

	int         desc[9];
	complex_d * data;

	int localRows, localCols;

	Matrix();

	void alloc(int ctxt, int N);
	void print();

	~Matrix();
};

Matrix::~Matrix() {
	delete[] data;
}

Matrix::Matrix() {
	for (int i = 0; i < 9; i++)
		desc[i] = -1;
	data = NULL;
	localRows = localCols = 0;
}

void Matrix::alloc(int ctxt, int N) {

	if (data != NULL)
		delete[] data;

	int gridRows, gridCols, procRow, procCol, rootRow = 0, rootCol = 0;
	int NB = 1;
	int info;


	Cblacs_gridinfo(ctxt, &gridRows, &gridCols, &procRow, &procCol);

	localRows = numroc_(&N, &NB, &procRow, &rootRow, &gridRows);
	localCols = numroc_(&N, &NB, &procCol, &rootCol, &gridCols);

	data = new complex_d [localRows * localCols];
	descinit_(desc, &N, &N, &NB, &NB, &rootRow, &rootCol, &ctxt, &localRows, &info);
}

void Matrix::print() {

	complex_d buf;
	MPI_Status status;

	int ctxt = desc[1];

	int N = desc[2];

	int gridRows, gridCols, procRow, procCol;
	Cblacs_gridinfo(ctxt, &gridRows, &gridCols, &procRow, &procCol);

	if ((procRow == 0) && (procCol == 0))
			cout << "Размер гамильтониана = " << N << endl;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (((j-procCol) >= 0) && ((i-procRow) >= 0) && (((j-procCol) % gridCols) == 0) && 
			(((i-procRow) % gridRows) == 0)) {

				int locJ = (j-procCol)/gridCols;
				int locI = (i-procRow)/gridRows;

				MPI_Send(&data[locJ * localRows + locI], 1, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);

			}
			if ((procRow == 0) && (procCol == 0)) {
				MPI_Recv(&buf, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
				cout << buf.real() << "\t";
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}
		if ((procRow == 0) && (procCol == 0))
			cout << endl;

		MPI_Barrier(MPI_COMM_WORLD);
	}
}

int combinations(unsigned int n, unsigned int m)
{
    if (m > n)
        return 0;

    int res = n--;

    if (m == 0)
    	return 1;

    for (unsigned int i = 2; i < m + 1; ++i, --n)
        res = res * n / i;

    return res;
}

void hamiltonian(int n, double *a, double *w, int k, int Emin, int Emax, int ctxt, Matrix &H) {

	int N = 0;

	vector<int> HPhotons; 
	vector<int> HSizes;
	vector<int> HOffsets; 

	HOffsets.push_back(0);

	for (int i = Emax; i >= Emin; i--) {
		for (int j = i; j >= max(i - k, 0); j--) {

			if (n >= j) {

				HPhotons.push_back(j);

				int c = combinations(n, j);
				HSizes.push_back(c);
				N += c;
				HOffsets.push_back(N);
			}
		}
	}
	
	int gridRows, gridCols, procRow, procCol;
	Cblacs_gridinfo(ctxt, &gridRows, &gridCols, &procRow, &procCol);

	H.alloc(ctxt, N);

	for (int j = 0; j < H.localCols; j++) {
		for (int i = 0; i < H.localRows; i++) {
		
			int globalJ = procCol + gridCols * j;
			int globalI = procRow + gridRows * i;

			int blockI;
			int blockJ;

			for (int k = 0; k < HPhotons.size(); k++) {
		
				if ((HOffsets[k] <= globalI) && (HOffsets[k + 1] > globalI)) {
					blockI = k;
					break;
				}
			}
			for (int k = 0; k < HPhotons.size(); k++) {
		
				if ((HOffsets[k] <= globalJ) && (HOffsets[k + 1] > globalJ)) {
					blockJ = k;
					break;
				}
			}

			if (blockI == blockJ) {
				int offsetI = globalI - HOffsets[blockI];
				int offsetJ = globalJ - HOffsets[blockJ];

				vector<int> left(n);
				
				for (int k = 0; k < n; k++) {
					if (k >= (n - HPhotons[blockI]))
						left[k] = 1;
					else
						left[k] = 0;
				}

				vector<int> top = left;

				for (int k = 0; k < offsetI; k++) {
					next_permutation(left.begin(), left.end());
				}
				for (int k = 0; k < offsetJ; k++) {
					next_permutation(top.begin(), top.end());
				}

				if (offsetI == offsetJ) {
					H.data[j * H.localRows + i] = 0;
					for (int k = 0; k < left.size(); k++) {
						if (left[k] == 1)
							H.data[j * H.localRows + i] += w[k];
					}
				}
				else {
					H.data[j * H.localRows + i] = 0;
					int a_index = -1;
					for (int k = 0; k < n; k++) {

						if (left[k] != top[k]) {
							if (a_index == -1) {
								a_index = -2;
								continue;
							}

							if (a_index == -2) {
								if (left[k - 1] != left[k])
									a_index = k - 1;
								else
									break;
								continue;
							}

							if (a_index >= 0) {
								a_index = -1;
								break;
							}
						}
						else {
							if (a_index == -2)
								break;
						}
					}

					if (a_index >= 0)
						H.data[j * H.localRows + i] = a[a_index];
				}
			}
			else
				H.data[j * H.localRows + i] = 0;

		}
	}
}

int main(int argc, char ** argv) {

	MPI_Init(&argc, &argv);

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if ((sqrt(size)) * (sqrt(size)) != size) return -1;
	
	int context;
	Cblacs_get(-1, 0, &context);
	Cblacs_gridinit(&context, (char *) "Row", sqrt(size), sqrt(size));
	
	Matrix H;

	int n = 30;
	double a[] = {
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3,
		0.1, 0.2, 0.3
		};
	double w[] = {
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000,
		1, 10, 100, 1000
		};
	int k = 0;
	int Emin = 1;
	int Emax = 1;

	hamiltonian(n, a, w, k, Emin, Emax, context, H);

	H.print();

	Cblacs_gridexit(context);
	Cblacs_exit(0);

	return 0;
}
