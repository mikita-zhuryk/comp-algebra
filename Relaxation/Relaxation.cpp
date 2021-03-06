#include "Relaxation.h"

using namespace CMA;

Relaxation::Relaxation(size_t dim, double w, int acc) : Method(dim), omega(w), accuracy(pow(10, acc)) { }

void Relaxation::solve(ostream& out) {
	b = A.transpose() * b;
	A = A.transpose() * A;
	x = b;
	int k = 1;
	Vector<double> temp(n, 0);
	do {
		if (k != 1) {
			x = temp;
		}
		for (size_t i = 0; i < n; ++i) {
			temp[i] = (1 - omega) * x[i] + omega * b[i] / A[i][i];
			for (size_t j = 0; j < i; ++j) {
				temp[i] -= omega * temp[j] * A[i][j] / A[i][i];
			}
			for (size_t j = i + 1; j < n; ++j) {
				temp[i] -= omega * x[j] * A[i][j] / A[i][i];
			}
		}
		++k;
	} while ((temp - x).norm() > accuracy);
	out << "Number of iterations: " << k << endl << endl;
}

double Relaxation::determinant() {
	return 0;
}