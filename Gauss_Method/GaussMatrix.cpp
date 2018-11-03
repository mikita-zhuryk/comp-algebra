#include "GaussMatrix.h"
#include <iomanip>

GaussMatrix::GaussMatrix(size_t dim): Method(dim) {
	swapCount = 0;
}

GaussMatrix::~GaussMatrix() {}

void GaussMatrix::solve(ostream& out) {
	makeUpperTriangular(A, b, true);
	printUpperTriangular(out);
	getSolution(A, b);
}

void GaussMatrix::swapRows(size_t k, size_t l) {
	if (k != l) {
		swap(A[k], A[l]);
		swap(b[k], b[l]);
		swap(inverse[k], inverse[l]);
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

void GaussMatrix::makeUpperTriangular(Matrix<double>& A, Vector<double>& b, bool savePivot = false) {
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
		if (savePivot) {
			pivot.push_back(A[k][k]);
		}
		A[k][k] = 1;
	}
}

void GaussMatrix::printUpperTriangular(ostream& out) {
	if (!out) {
		throw invalid_argument("Bad output stream in printUpperTriangular(ostream&).");
	}
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

void GaussMatrix::getSolution(Matrix<double>& A, Vector<double>& b) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (size_t j = n - 1; j > i; --j) {
			x[i] -= x[j] * A[i][j];
		}
	}
}

double GaussMatrix::determinant() {
	for (size_t i = 0; i < n; ++i) {
		detA *= pivot[i];
	}
	return ((swapCount % 2) ? -detA : detA);
}

void GaussMatrix::findInverse() {
	Vector<double> temp(n, 0);
	for (size_t k = 0; k < n; ++k) {
		temp[k] = 1;
		A = initial_A;
		this->makeUpperTriangular(A, temp);
		this->getSolution(A, temp);
		for (size_t j = 0; j < n; ++j) {
			inverse[j][k] = x[j];
			temp[j] = 0;
		}
	}
}