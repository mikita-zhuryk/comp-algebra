#include "Danilevsky.h"
#include <fstream>

using namespace CMA;

Danilevsky::Danilevsky(size_t dim): EigenMethod(dim) {
	S = Matrix<double>::identity(n, 1);
	invS = Matrix<double>::identity(n, 1);
}

void Danilevsky::solve(ostream& out) {
	Frobenius();
	collectPoly();
	findEigenVectors();
}

void Danilevsky::Frobenius() {
	Matrix<double> M;
	Matrix<double> invM;
	for (size_t k = n - 1; k > 0; --k) {
		M = Matrix<double>::identity(n, 1);
		invM = Matrix<double>::identity(n, 1);
		M[k - 1] = A[k] * (-1.0 / A[k][k - 1]);
		M[k - 1][k - 1] = 1.0 / A[k][k - 1];
		invM[k - 1] = A[k];
		invS = invM * invS;
		S = S * M;
		A = invM * A * M;
	}
}

void Danilevsky::collectPoly() {
	eigenPolynomial[0] = 1;
	for (size_t j = 0; j < n; ++j) {
		eigenPolynomial[j + 1] = A[0][j] * pow(-1, n % 2);
	}
}

void Danilevsky::findEigenVectors() {
	Vector<double> y;
	Vector<double> x;
	for (size_t i = 0; i < n; ++i) {
		y = buildY(eigenValues[i]);
		x = S * y;
		eigenVectors.push_back(x / x.norm());
	}
}

Vector<double> Danilevsky::buildY(double eigenValue) {
	Vector<double> y(n);
	for (size_t i = 0; i < n; ++i) {
		y[i] = pow(eigenValue, n - 1 - i);
	}
	return y;
}