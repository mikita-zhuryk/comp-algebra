#include "GaussMatrix.h"
#include <algorithm>
#include <iomanip>
#include <fstream>

void GaussMatrix::setDimensions(size_t dim) {
	if (dim <= 0) {
		throw invalid_argument("Matrix dimensions are less than or equal to 0.");
	}
	n = dim;
}

void GaussMatrix::swapRows(size_t k, size_t l) {
	if (k != l) {
		swap(A[k], A[l]);
		swap(b[k], b[l]);
		swapCount++;
	}
}

size_t GaussMatrix::findMaxInRows(UINT step) {
	size_t max_index = step;
	for (size_t i = step + 1; i < n; ++i) {
		if (abs(A[i][step]) > abs(A[max_index][step])) {
			max_index = i;
		}
	}
	return max_index;
}

GaussMatrix::GaussMatrix(size_t dim) {
	setDimensions(dim);
	A = Matrix<double>(n, 0);
	initial_A = Matrix<double>(n, 0);
	inverse = Matrix<double>(n, 0);
	inverse += 1;
	b = Vector<double>(n, 0);
	initial_b = Vector<double>(n, 0);
	x = Vector<double>(n, 0);
	swapCount = 0;
	detA = 0;
}

GaussMatrix::~GaussMatrix() {}

void GaussMatrix::makeUpperTriangular() {
	detA = A[0][0];
	for (size_t k = 0; k < n; ++k) {
		swapRows(k, findMaxInRows(k));
		b[k] /= A[k][k];
		for (size_t j = k + 1; j < n; ++j) {
			A[k][j] /= A[k][k];
			b[j] -= b[k] * A[j][k];
			inverse[j] -= inverse[k] * A[j][k];
			for (size_t i = k + 1; i < n; ++i) {
				A[i][j] -= A[k][j] * A[i][k];
			}
		}
		detA *= A[k][k];
		inverse[k] *= 1 / A[k][k];
		A[k][k] = 1;
	}
}

void GaussMatrix::printUpperTriangular(ostream& out) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			out << setw(OUTPUT_WIDTH) << 0 << " ";
		}
		for (size_t j = i; j < n; ++j) {
			out << setw(OUTPUT_WIDTH) << A[i][j] << " ";
		}
		out << setw(OUTPUT_WIDTH) << b[i] << endl;
	}
	out << endl;
}

void GaussMatrix::getSolution() {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (size_t j = n - 1; j > i; --j) {
			x[i] -= x[j] * A[i][j];
		}
		for (int j = i - 1; j >= 0; --j) {
			inverse[j] -= inverse[i] * A[j][i];
			cout << inverse << endl;
		}
	}
}

void GaussMatrix::printSolution(ostream& out) {
	out << "x = (" << x << ")" << endl;
	out << endl;
}

void GaussMatrix::deficiency(ostream& out) {
	Vector<double> deficiency(n);
	for (size_t i = 0; i < n; ++i) {
		deficiency[i] = -initial_b[i];
		for (size_t j = 0; j < n; ++j) {
			deficiency[i] = deficiency[i] + initial_A[i][j] * x[j];
		}
	}
	out << "r = (" << deficiency << ")" << endl;
	out << endl;
}

double GaussMatrix::determinant() const {
	return detA;
}

void GaussMatrix::printInverse(ostream& out) {
	out << inverse << endl;
}

void GaussMatrix::inverseDeficiency(ostream& out) {
	out << (Matrix<double>(n, 0) + 1) - inverse * initial_A << endl;
}

istream& operator>>(istream& in, GaussMatrix& obj) {
	for (size_t i = 0; i < obj.n; ++i) {
		for (size_t j = 0; j < obj.n; ++j) {
			in >> obj.A[i][j];
			obj.initial_A[i][j] = obj.A[i][j];
		}
		in >> obj.b[i];
		obj.initial_b[i] = obj.b[i];
	}
	return in;
}

ostream& operator<<(ostream& out, GaussMatrix& obj) {
	for (size_t i = 0; i < obj.n; ++i) {
		out << obj.A[i] << " ";
		out << setw(OUTPUT_WIDTH) << obj.b[i] << endl;
	}
	out << endl;
	return out;
}
