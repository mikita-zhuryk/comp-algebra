#include "GaussMatrix.h"
#include <algorithm>
#include <iomanip>

void GaussMatrix::setDimensions(size_t dim) {
	if (dim <= 0) {
		throw invalid_argument("Matrix dimensions less than or equal to 0.");
	}
	n = dim;
}

void GaussMatrix::swapRows(size_t k, size_t l) {
	if (k != l) {
		swap(values[k], values[l]);
		swap(b[k], b[l]);
		swap(x[k].first, x[l].first);
		swapCount++;
	}
}

size_t GaussMatrix::findMaxInRows(UINT step) {
	size_t max_index = step;
	for (size_t i = step + 1; i < n; ++i) {
		if (abs(values[i][step]) > abs(values[max_index][step])) {
			max_index = i;
		}
	}
	return max_index;
}

GaussMatrix::GaussMatrix(size_t dim) {
	setDimensions(dim);
	values = new double*[n];
	for (size_t i = 0; i < n; ++i) {
		values[i] = new double[n];
	}
	b.resize(n);
	x.resize(n);
	for (size_t i = 0; i < n; ++i) {
		b[i].first = i;
		b[i].second = 0;
		x[i].first = i;
		x[i].second = 0;
		for (size_t j = 0; j < n; ++j) {
			values[i][j] = 0.0;
		}
	}
	swapCount = 0;
}

void GaussMatrix::makeUpperTriangular() {
	for (size_t i = 0; i < n - 1; ++i) {
		swapRows(i, findMaxInRows(i));
		b[i].second /= values[i][i];
		b[i + 1].second -= b[i].second * values[i + 1][i];
		for (size_t j = i + 1; j < n; ++j) {
			values[i][j] /= values[i][i];
			values[i + 1][j] -= values[i][j] * values[i + 1][i];
		}
		values[i][i] = 1;
	}
	b[n - 1].second /= values[n - 1][n - 1];
	values[n - 1][n - 1] = 1;
}

void GaussMatrix::printUpperTriangular(ostream& out) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << 0 << " ";
		}
		for (size_t j = i; j < n; ++j) {
			out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << values[i][j] << " ";
		}
		out << b[i].second << endl;
	}
}

void GaussMatrix::getSolution() {
	double tempSum = 0;
	for (int i = n - 1; i >= 0; --i) {
		tempSum = 0;
		for (int j = n - 1; j > i; --j) {
			tempSum += values[i][j] * x[j].second;
		}
		x[i].second = b[i].second - tempSum;
	}
	//sort(x.begin(), x.end(), [](pair<size_t, double> x_i, pair<size_t, double> x_j) { return x_i.first > x_j.second; });
	sort(x.begin(), x.end(), IndexComparator());
}

void GaussMatrix::printSolution(ostream& out) {
	out << "x = (";
	for (size_t i = 0; i < n - 1; ++i) {
		out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << x[i].second << ", ";
	}
	out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << x[n - 1].second << ")";
}

istream& operator>>(istream& in, GaussMatrix& obj) {
	for (size_t i = 0; i < obj.n; ++i) {
		for (size_t j = 0; j < obj.n; ++j) {
			in >> obj.values[i][j];
		}
		in >> obj.b[i].second;
	}
	return in;
}
