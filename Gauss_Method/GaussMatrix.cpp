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
	b = Vector<double>(n, 0);
	x = Vector<double>(n, 0);
	swapCount = 0;
	det = 0;
}

GaussMatrix::~GaussMatrix() {}

void GaussMatrix::makeUpperTriangular() {
	det = A[0][0];
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
		det *= A[k][k];
		A[k][k] = 1;
	}
}

void GaussMatrix::printUpperTriangular(ostream& out) const {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << 0 << " ";
		}
		for (size_t j = i; j < n; ++j) {
			out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << A[i][j] << " ";
		}
		out << b[i] << endl;
	}
	out << endl;
}

void GaussMatrix::getSolution() {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = n - 1; j > i; --j) {
			x[i] -= x[j] * A[i][j];
		}
	}
}

void GaussMatrix::printSolution(ostream& out) const {
	out << "x = (";
	for (size_t i = 0; i < n - 1; ++i) {
		out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << x[i] << ", ";
	}
	out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << x[n - 1] << ")" << endl;
	out << endl;
}

void GaussMatrix::deficiency(ostream& out) const {
	GaussMatrix temp(5);
	ifstream fIn("input.txt");
	fIn >> temp;
	double* deficiency = new double[n];
	for (size_t i = 0; i < n; ++i) {
		deficiency[i] = -temp.b[i];
		for (size_t j = 0; j < n; ++j) {
			deficiency[i] = deficiency[i] + temp.A[i][j] * x[j];
		}
	}
	out << "r = (";
	for (size_t i = 0; i < n - 1; ++i) {
		out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << deficiency[i] << ", ";
	}
	out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << deficiency[n - 1] << ")" << endl;
	out << endl;
	delete[] deficiency;
}

double GaussMatrix::determinant() const {
	return det;
}

istream& operator>>(istream& in, GaussMatrix& obj) {
	for (size_t i = 0; i < obj.n; ++i) {
		for (size_t j = 0; j < obj.n; ++j) {
			in >> obj.A[i][j];
		}
		in >> obj.b[i];
	}
	return in;
}

ostream& operator<<(ostream& out, const GaussMatrix& obj) {
	for (size_t i = 0; i < obj.n; ++i) {
		for (size_t j = 0; j < obj.n; ++j) {
			out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << obj.A[i][j] << " ";
		}
		out << obj.b[i] << endl;
	}
	out << endl;
	return out;
}
