#include "Sqrt_Method.h"

using namespace std;

Sqrt_Method::Sqrt_Method(size_t dim): Method(dim) {
	S = Matrix<double>(n, 0);
	y = Vector<double>(n, 0);
}

void Sqrt_Method::solve(ostream& out) {
	Matrix<double> AT = A.transpose();
	Vector<double> temp(n);
	A = AT * A;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			temp[i] += AT[i][j] * b[j];
		}
	}
	b = temp;
	calculateS();
	out << "S:" << endl << S << endl;
	solveForY();
	solveForX();
	findInverse();
}

void Sqrt_Method::calculateS() {
	S[1][1] = sqrt(A[1][1]);
	for (size_t j = 0; j < n; ++j) {
		S[1][j] = A[1][j] / S[1][1];
	}
	double squareSum = 0;
	double strangeSum = 0;
	for (size_t i = 0; i < n; ++i) {
		squareSum = 0;
		for (size_t k = 0; k < i; ++k) {
			squareSum += pow(S[k][i], 2);
		}
		S[i][i] = sqrt(A[i][i] - squareSum);
		for (size_t j = i + 1; j < n; ++j) {
			strangeSum = 0;
			for (size_t k = 0; k < i; ++k) {
				strangeSum += S[k][i] * S[k][j];
			}
			S[i][j] = (A[i][j] - strangeSum) / S[i][i];
		}
	}
}

void Sqrt_Method::solveForY() {
	y[0] = b[0] / S[0][0];
	double partialSum = 0;
	for (size_t i = 1; i < n; ++i) {
		partialSum = 0;
		for (size_t k = 0; k < i; ++k) {
			partialSum += S[k][i] * y[k];
		}
		y[i] = (b[i] - partialSum) / S[i][i];
	}
}

void Sqrt_Method::solveForX() {
	x[n - 1] = y[n - 1] / S[n - 1][n - 1];
	double partialSum = 0;
	for (int i = n - 2; i >= 0; --i) {
		partialSum = 0;
		for (size_t k = i + 1; k < n; ++k) {
			partialSum += S[i][k] * x[k];
		}
		x[i] = (y[i] - partialSum) / S[i][i];
	}
}

double Sqrt_Method::determinant() {
	for (size_t i = 0; i < n; ++i) {
		detA *= abs(S[i][i]);
	}
	return detA;
}

void Sqrt_Method::findInverse() {
	inverse[n - 1][n - 1] = 1 / pow(S[n - 1][n - 1], 2);
	double sum = 0;
	for (int i = n - 2; i >= 0; --i) {
		for (size_t j = i + 1; j < n; ++j) {
			sum = 0;
			for (size_t k = i + 1; k < n; ++k) {
				sum += S[i][k] * inverse[k][j];
			}
			inverse[i][j] = -sum / S[i][i];
			inverse[j][i] = inverse[i][j];
		}
		sum = 0;
		for (size_t k = i + 1; k < n; ++k) {
			sum += S[i][k] * inverse[i][k];
		}
		inverse[i][i] = (1 / S[i][i] - sum) / S[i][i];
	}
}

void Sqrt_Method::printInverse(ostream& out) {
	if (!out) {
		throw invalid_argument("Bad output stream in printInverse(ostream&).");
	}
	out << "(A^TA)^(-1):" << endl;
	out << inverse << endl;
}

void Sqrt_Method::printInverseDeficiency(ostream& out) {
	if (!out) {
		throw invalid_argument("Bad output stream in Sqrt_Method::printInverseDeficiency(ostream&).");
	}
	Matrix<double> AT = initial_A.transpose();
	out << "R:" << endl;
	out << Matrix<double>(n, 0) + 1 - AT * initial_A * inverse << endl;
}