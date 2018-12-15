#include "Pass_Method.h"

using namespace CMA;

PassMethod::PassMethod(size_t dim, bool makeTriDiag): Method(dim) {
	mid = dim / 2;
	alpha = vector<double>(dim);
	beta = vector<double>(dim);
	ksi = vector<double>(dim);
	eta = vector<double>(dim);
	delta = vector<double>(dim);
	diag = vector<double>(dim);
	upper = vector<double>(dim);
	lower = vector<double>(dim);
	if (makeTriDiag) {
		transformed = true;
	}
}

void PassMethod::solve(ostream& out) {
	if (transformed) {
		buildTriDiag();
		out << *this;
	}
	fillVectors();
	calcRight();
	calcLeft();
	printCoeffs(out);
	calcX();
}

void PassMethod::buildTriDiag() {
	for (size_t k = 1; k < n; ++k) {
		for (size_t j = k + 1; j < n; ++j) {
			b[j] -= b[k] * A[j][k - 1] / A[k][k - 1];
			A[j] -= A[k] * (A[j][k - 1] / A[k][k - 1]);
		}
	}
	for (int i = n - 2; i >= 0; --i) {
		for (int j = i - 1; j >= 0; --j) {
			b[j] -= b[i] * A[j][i + 1] / A[i][i + 1];
			A[j] -= A[i] * (A[j][i + 1] / A[i][i + 1]);
		}
	}
}

void PassMethod::fillVectors() {
	diag[0] = A[0][0];
	for (size_t i = 1; i < n; ++i) {
		diag[i] = A[i][i];
		lower[i] = -A[i][i - 1];
		upper[i - 1] = -A[i - 1][i];
	}
}

void PassMethod::calcRight() {
	alpha[1] = upper[0] / diag[0];
	beta[1] = b[0] / diag[0];
	double denum = diag[0];
	for (size_t i = 1; i < n - 1; ++i) {
		denum = diag[i] - alpha[i] * lower[i];
		alpha[i + 1] = upper[i] / denum;
		beta[i + 1] = (b[i] + beta[i] * lower[i]) / denum;
	}
}

void PassMethod::calcLeft() {
	ksi[n - 1] = lower[n - 1] / diag[n - 1];
	eta[n - 1] = b[n - 1] / diag[n - 1];
	double denum = diag[n - 1];
	for (size_t i = n - 2; i > mid; --i) {
		denum = diag[i] - ksi[i + 1] * upper[i];
		ksi[i] = lower[i] / denum;
		eta[i] = (b[i] + upper[i] * eta[i + 1]) / denum;
	}
}

void PassMethod::printCoeffs(ostream& out) {
	out << "Alpha: ";
	for (size_t i = 0; i < n; ++i) {
		if (alpha[i]) {
			out << alpha[i] << " ";
		}
	}
	out << endl;
	out << "Ksi: ";
	for (size_t i = 0; i < n; ++i) {
		if (ksi[i]) {
			out << ksi[i] << " ";
		}
	}
	out << endl << endl;
}

void PassMethod::calcX() {
	x[mid] = (beta[mid + 1] + alpha[mid + 1] * eta[mid + 1]) / (1 - alpha[mid + 1] * ksi[mid + 1]);
	for (int i = mid - 1; i >= 0; --i) {
		x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
	}
	for (size_t i = mid; i < n - 1; ++i) {
		x[i + 1] = ksi[i + 1] * x[i] + eta[i + 1];
	}
}

void PassMethod::printDelta(ostream& out) {
	out << Vector<double>(delta) << endl;
}

double PassMethod::determinant() {
	detA *= diag[0];
	for (size_t i = 1; i < n; ++i) {
		delta[i] = diag[i] - alpha[i] * lower[i];
		detA *= delta[i];
	}
	return detA;
}

