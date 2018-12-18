#include "IterJacobi.h"

using namespace CMA;

IterJacobi::IterJacobi(size_t dim, int acc) : EigenMethod(dim) {
	accuracy = pow(10, acc);
	T = Matrix<double>::identity(n, 1);
}

void IterJacobi::find(ostream& out) {
	buildDiag();
	Matrix<double> temp = T.transpose();
	for (size_t i = 0; i < n; ++i) {
		eigenValues[i] = A[i][i];
		eigenVectors.push_back(temp[i]);
	}
}

void IterJacobi::buildDiag() {
	pair<size_t, size_t> max = findMax();
	double cos = 0;
	double sin = 0;
	double mu = 0;
	Matrix<double> Ti;
	size_t i, j;
	do {
		i = max.first;
		j = max.second;
		Ti = Matrix<double>::identity(n, 1);
		mu = 2 * A[i][j] / (A[i][i] - A[j][j]);
		cos = sqrt((1 + 1 / sqrt(1 + pow(mu, 2))) / 2);
		sin = (mu / abs(mu)) * sqrt((1 - 1 / sqrt(1 + pow(mu, 2))) / 2);
		Ti[j][j] = Ti[i][i] = cos;
		Ti[i][j] = -sin;
		Ti[j][i] = sin;
		T *= Ti;
		A = Ti.transpose() * A * Ti;
		max = findMax();
	} while (abs(A[max.first][max.second]) > accuracy);
}

pair<size_t, size_t> IterJacobi::findMax() {
	double max = 0;
	pair<size_t, size_t> index(100, 100);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			if (abs(A[i][j]) > max) {
				max = abs(A[i][j]);
				index = pair<size_t, size_t>(i, j);
			}
		}
	}
	return index;
}