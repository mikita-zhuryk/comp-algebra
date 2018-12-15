#include "GaussSeidel.h"

using namespace CMA;

GaussSeidel::GaussSeidel(size_t dim): Method(dim) { }

void GaussSeidel::solve(ostream& out) {
	x = b;
	Vector<double> newX(x);
	double eps = pow(10, -5);
	int i = 0;
	this->printSolution(out);
	while (calcDeficiency().norm() > eps) {
		for (int i = 0; i < n; ++i) {
			newX[i] = b[i];
			for (int j = 0; j < i; ++j) {
				newX[i] -= A[i][j] * newX[j];
			}
			for (int j = i + 1; j < n; ++j) {
				newX[i] -= A[i][j] * x[j];
			}
			newX[i] /= A[i][i];
			x[i] = newX[i];
		}
		out << i + 1 << " iteration" << endl;
		this->printSolution(out);
		++i;
	}
}

double GaussSeidel::determinant() {
	return 0;
}