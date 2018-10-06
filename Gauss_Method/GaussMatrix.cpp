#include "GaussMatrix.h"
#include <iomanip>

GaussMatrix::GaussMatrix(size_t dim): Method(dim) {
	swapCount = 0;
}

GaussMatrix::~GaussMatrix() {}

void GaussMatrix::solve(ostream& out) {
	makeUpperTriangular();
	printUpperTriangular(out);
	getSolution();
	printSolution(out);
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

void GaussMatrix::makeUpperTriangular() {
	for (size_t k = 0; k < n; ++k) {
		swapRows(k, findMaxInRows(k));
		b[k] /= A[k][k];
		for (size_t j = k + 1; j < n; ++j) {
			A[k][j] /= A[k][k];
			b[j] -= b[k] * A[j][k];
			for (size_t i = k + 1; i < n; ++i) {
				A[i][j] -= A[k][j] * A[i][k];
			}
		}
		inverse[k] *= 1 / A[k][k];
		for (size_t i = k + 1; i < n; ++i) {
			inverse[i] -= inverse[k] * A[i][k];
		}
		detA *= A[k][k];
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
		}
	}
}

double GaussMatrix::determinant() const {
	return ((swapCount % 2) ? -detA : detA);
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
