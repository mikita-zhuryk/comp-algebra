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

GaussMatrix::~GaussMatrix() {
	for (size_t i = 0; i < n; ++i) {
		delete[] values[i];
	}
	delete[] values;
	b.~vector();
	x.~vector();
}

void GaussMatrix::makeUpperTriangular() {
	for (size_t k = 0; k < n; ++k) {
		swapRows(k, findMaxInRows(k));
		b[k].second /= values[k][k];
		for (size_t j = k + 1; j < n; ++j) {
			values[k][j] /= values[k][k];
			b[j].second -= b[k].second * values[j][k];
			for (size_t i = k + 1; i < n; ++i) {
				values[i][j] -= values[k][j] * values[i][k];
			}
		}
		values[k][k] = 1;
	}
}

void GaussMatrix::printUpperTriangular(ostream& out) const {
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
	for (int i = n - 1; i >= 0; --i) {
		x[i].second = b[i].second;
		for (int j = n - 1; j > i; --j) {
			x[i].second -= x[j].second * values[i][j];
		}
	}
	sort(x.begin(), x.end(), IndexComparator());
}

void GaussMatrix::printSolution(ostream& out) const {
	out << "x = (";
	for (size_t i = 0; i < n - 1; ++i) {
		out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << x[i].second << ", ";
	}
	out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << x[n - 1].second << ")";
}

void GaussMatrix::calculate() const {
	GaussMatrix temp(5);
	ifstream fIn("input.txt");
	fIn >> temp;
	double* deficiency = new double[n];
	for (size_t i = 0; i < n; ++i) {
		deficiency[i] = -temp.b[i].second;
		for (size_t j = 0; j < n; ++j) {
			deficiency[i] = deficiency[i] + temp.values[i][j] * x[j].second;
		}
	}
	delete[] deficiency;
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

ostream& operator<<(ostream& out, const GaussMatrix& obj) {
	for (size_t i = 0; i < obj.n; ++i) {
		for (size_t j = 0; j < obj.n; ++j) {
			out << setw(OUTPUT_WIDTH) << setprecision(DECIMAL_PRECISION) << obj.values[i][j] << " ";
		}
		out << obj.b[i].second << endl;
	}
	out << endl;
	return out;
}
