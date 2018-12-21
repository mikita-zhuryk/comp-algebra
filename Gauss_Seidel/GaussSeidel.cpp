#include "GaussSeidel.h"

using namespace CMA;

GaussSeidel::GaussSeidel(size_t dim, int acc): Method(dim), accuracy(pow(10, acc)) { }

void GaussSeidel::solve(ostream& out) {
	x = b;
	Vector<double> newX(n, 0);
	int k = 1;
	do {
		if (k != 1) {
			x = newX;
		}
		for (int i = 0; i < n; ++i) {
			newX[i] = b[i];
			for (int j = 0; j < i; ++j) {
				newX[i] -= A[i][j] * newX[j];
			}
			for (int j = i + 1; j < n; ++j) {
				newX[i] -= A[i][j] * x[j];
			}
			newX[i] /= A[i][i];
		}
		++k;
	} while ((newX - x).norm() > accuracy);
	out << "Number of iterations: " << k << endl << endl;
}

double GaussSeidel::determinant() {
	return 0;
}