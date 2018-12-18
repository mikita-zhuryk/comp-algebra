#include "Reflections.h"

using namespace CMA;

Reflections::Reflections(size_t dim) : Method(dim) {}

Reflections::Reflections(Matrix<double> M, Vector<double> b) : Method(M, b) {}

void Reflections::solve(ostream& out) {
	Matrix<double> V;
	for (size_t i = 0; i < n - 1; ++i) {
		pickY(i);
		calcOmega();
		V = Matrix<double>::identity(n, 1) - omega.matrixMult(omega) * 2;
		A = V * A;
		b = V * b;
	}
	out << *this;
	getSolution();
}

void Reflections::pickY(size_t step) {
	Matrix<double> AT = A.transpose();
	y = AT[step];
	for (size_t i = 0; i < step; ++i) {
		y[i] = 0;
	}
	e = Vector<double>(n, 0.0);
	e[step] = 1;
}

void Reflections::calcOmega() {
	double alpha = sqrt(y.dot(y));
	Vector<double> changedY = y - e * alpha;
	double p = sqrt(2 * y.dot(changedY));
	omega = changedY / p;
}

void Reflections::getSolution() {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (size_t j = n - 1; j > i; --j) {
			x[i] -= x[j] * A[i][j];
		}
		x[i] /= A[i][i];
	}
}

double Reflections::determinant() {
	for (size_t i = 0; i < n; ++i) {
		detA *= A[i][i];
	}
	return detA;
}

