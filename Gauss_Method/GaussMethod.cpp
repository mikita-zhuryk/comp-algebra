#include "GaussMethod.h"
#include <iomanip>

using namespace CMA;

GaussMethod::GaussMethod(size_t dim): Method(dim) {
	swapCount = 0;
}

GaussMethod::~GaussMethod() {}

void GaussMethod::solve(ostream& out) {
	makeUpperTriangular();
	printUpperTriangular(out);
	getSolution(b);
}

void GaussMethod::swapRows(size_t k, size_t l) {
	if (k != l) {
		swap(A[k], A[l]);
		swap(b[k], b[l]);
		swap(inverse[k], inverse[l]);
		swapCount++;
	}
}

size_t GaussMethod::findMaxInRows(UINT step) {
	size_t max_index = step;
	for (size_t i = step + 1; i < n; ++i) {
		if (abs(A[i][step]) > abs(A[max_index][step])) {
			max_index = i;
		}
	}
	return max_index;
}

void GaussMethod::makeUpperTriangular() {
	for (size_t k = 0; k < n; ++k) {
		swapRows(k, findMaxInRows(k));
		b[k] /= A[k][k];
		inverse[k] /= A[k][k];
		pivot.push_back(A[k][k]);
		for (size_t j = k + 1; j < n; ++j) {
			A[k][j] /= A[k][k];
		}
		A[k][k] = 1;
		for (size_t j = k + 1; j < n; ++j) {
			inverse[j] -= inverse[k] * A[j][k];
			b[j] -= b[k] * A[j][k];
			A[j] -= A[k] * A[j][k];
		}
	}
	double eps = pow(10.0, -15);
	if ((abs(A[n - 1][n - 1]) <= eps) && (abs(b[n - 1]) >= eps)) {
		throw "No solutions";
	}
}

void GaussMethod::printUpperTriangular(ostream& out) {
	if (!out) {
		throw invalid_argument("Bad output stream in printUpperTriangular(ostream&).");
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			out << setw(OUTPUT_WIDTH) << A[i][j] << " ";
		}
		out << setw(OUTPUT_WIDTH) << b[i] << endl;
	}
	out << endl;
}

void GaussMethod::getSolution(Vector<double>& b) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (size_t j = n - 1; j > i; --j) {
			x[i] -= x[j] * A[i][j];
		}
	}
}

double GaussMethod::determinant() {
	for (size_t i = 0; i < n; ++i) {
		detA *= pivot[i];
	}
	return ((swapCount % 2) ? -detA : detA);
}

void GaussMethod::findInverse() {
	Vector<double> temp(n, 0);
	for (size_t k = 0; k < n; ++k) {
		for (size_t j = 0; j < n; ++j) {
			temp[j] = inverse[j][k];
		}
		getSolution(temp);
		for (size_t j = 0; j < n; ++j) {
			inverse[j][k] = x[j];
		}
	}
}