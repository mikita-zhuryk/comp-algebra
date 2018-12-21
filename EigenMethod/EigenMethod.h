#pragma once
#include "../Matrix/Matrix.h"
#include "../Vector/Vector.h"
#include <list>
#define DECIMAL_PRECISION 5

using namespace std;

namespace CMA {

	class EigenMethod {

		EigenMethod() {}

		void setDimensions(size_t dim) {
			if (dim <= 0) {
				throw invalid_argument("Matrix dimensions are less than or equal to 0.");
			}
			n = dim;
		}

	protected:

		Matrix<double> A;
		Matrix<double> initial_A;
		Vector<double> eigenPolynomial;
		Vector<double> eigenValues;
		list<Vector<double>> eigenVectors;
		size_t n;

	public:

		EigenMethod(size_t dim) {
			setDimensions(dim);
			A = Matrix<double>(n, 0);
			initial_A = Matrix<double>(n, 0);
			eigenPolynomial = Vector<double>(n + 1, 0);
			eigenValues = Vector<double>(n);
			eigenVectors = list<Vector<double>>();
		}

		~EigenMethod() {}

		void run(ostream& out, bool findInv = false) {
			out.setf(ios_base::fixed, ios_base::floatfield);
			out.precision(DECIMAL_PRECISION);
			find(out);
			printSolution(out);
			out.setf(ios::scientific, ios_base::floatfield);
			printDeficiency(out);
			printVecDeficiency(out);
			out.setf(ios::fixed, ios_base::floatfield);
		}

		friend istream& operator>>(istream& in, EigenMethod& obj) {
			if (!in) {
				throw invalid_argument("Bad input stream in operator>>(istream&, EigenMethod&).");
			}
			in.peek();
			if (in.eof()) {
				throw invalid_argument("Empty input file in operator>>(istream&, EigenMethod&).");
			}
			for (size_t i = 0; i < obj.n; ++i) {
				for (size_t j = 0; j < obj.n; ++j) {
					in >> obj.A[i][j];
					obj.initial_A[i][j] = obj.A[i][j];
				}
			}
			double temp;
			for (size_t i = 0; i < obj.n; ++i) {
				in >> obj.eigenValues[i];
			}
			return in;
		}

		friend ostream& operator<<(ostream& out, EigenMethod& obj) {
			if (!out) {
				throw invalid_argument("Bad output stream in operator<<(ostream&, EigenMethod&).");
			}
			obj.printSolution(out);
			return out;
		}

	protected:

		virtual void find(ostream&) = 0;

		virtual void printSolution(ostream& out) {
			if (!out) {
				throw invalid_argument("Bad output stream in printSolution(ostream&).");
			}
			out << A << endl;
			out << "Polynomial:\n" << eigenPolynomial << endl;
			out << "Eigenvalues:\n" << eigenValues << endl;
			out << "Eigenvectors: \n";
			for (auto vec : eigenVectors) {
				out << vec << endl;
			}
			out << endl;
		}

		virtual void printDeficiency(ostream& out) {
			Vector<double> def(n);
			for (size_t i = 0; i < n; ++i) {
				def[i] = singleValDeficiency(eigenValues[i]);
			}
			out << "Eigenvalues deficiency:\n" << def << endl;
		}

		double singleValDeficiency(double v) {
			double result = 0;
			for (size_t i = 0; i < n + 1; ++i) {
				result += eigenPolynomial[i] * pow(v, n - i);
			}
			return result;
		}

		virtual void printVecDeficiency(ostream& out) {
			list<Vector<double>> def;
			size_t i = 0;
			for (auto it = eigenVectors.begin(); it != eigenVectors.end(); ++it) {
				def.push_back(initial_A * *it - *it * eigenValues[i]);
				++i;
			}
			out << "Eigenvectors deficiency norms:\n";
			i = 1;
			for (auto it = def.begin(); it != def.end(); ++it) {
				out << i << ": " << (*it).norm() << endl;
				++i;
			}
		}

	};

}