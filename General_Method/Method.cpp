#include "Method.h"

using namespace std;

Method::Method() {}

void Method::setDimensions(size_t dim) {
	if (dim <= 0) {
		throw invalid_argument("Matrix dimensions are less than or equal to 0.");
	}
	n = dim;
}

Method::Method(size_t dim) {
	setDimensions(dim);
	A = Matrix<double>(n, 0);
	initial_A = Matrix<double>(n, 0);
	inverse = Matrix<double>(n, 0);
	inverse += 1;
	b = Vector<double>(n, 0);
	initial_b = Vector<double>(n, 0);
	x = Vector<double>(n, 0);
	detA = 1;
}

Method::~Method() {}

void Method::printSolution(ostream& out) {
	out << "x = (" << x << ")" << endl;
	out << endl;
}

void Method::deficiency(ostream& out) {
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

double Method::determinant() const {
	return detA;
}

void Method::printInverse(ostream& out) {
	out << inverse << endl;
}

void Method::printInverseDeficiency(ostream& out) {
	out << (Matrix<double>(n, 0) + 1) - inverse * initial_A << endl;
}