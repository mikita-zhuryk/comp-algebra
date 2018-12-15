#include "../Matrix/Matrix.h"
#include "../Vector/Vector.h"
#pragma once
#define DECIMAL_PRECISION 5

using namespace std;

namespace CMA {

	class Method {

		Method() {}

		void setDimensions(size_t dim) {
			if (dim <= 0) {
				throw invalid_argument("Matrix dimensions are less than or equal to 0.");
			}
			n = dim;
		}

	protected:

		Matrix<double> A;
		Matrix<double> initial_A;
		Matrix<double> inverse;
		Vector<double> b;
		Vector<double> initial_b;
		Vector<double> x;
		double detA;
		size_t n;

	public:

		Method(size_t dim) {
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

		~Method() {}

		void run(ostream& out, bool findInv = false) {
			out.setf(ios_base::fixed, ios_base::floatfield);
			out.precision(DECIMAL_PRECISION);
			solve(out);
			printSolution(out);
			out.setf(ios_base::scientific, ios_base::floatfield);
			printDeficiency(out);
			out.setf(ios_base::fixed, ios_base::floatfield);
			out << "det = " << determinant() << endl << endl;
			if (findInv) {
				findInverse();
				printInverse(out);
				out.setf(ios_base::scientific, ios_base::floatfield);
				printInverseDeficiency(out);
			}
		}

		virtual void solve(ostream&) = 0;

		virtual void findInverse() {

		}

		void printSolution(ostream& out) {
			if (!out) {
				throw invalid_argument("Bad output stream in printSolution(ostream&).");
			}
			out << "x = (" << x << ")" << endl;
			out << endl;
		}

		Vector<double> calcDeficiency() {
			Vector<double> deficiency(n);
			for (size_t i = 0; i < n; ++i) {
				deficiency[i] = -initial_b[i];
				for (size_t j = 0; j < n; ++j) {
					deficiency[i] = deficiency[i] + initial_A[i][j] * x[j];
				}
			}
			return deficiency;
		}

		void printDeficiency(ostream& out) {
			if (!out) {
				throw invalid_argument("Bad output stream in printDeficiency(ostream&).");
			}
			auto deficiency = calcDeficiency();
			out << "r = (" << deficiency << ")" << endl;
			out << "||r|| = " << deficiency.norm() << endl;
			out << endl;
		}

		virtual double determinant() = 0;

		virtual void printInverse(ostream& out) {
			if (!out) {
				throw invalid_argument("Bad output stream in printInverse(ostream&).");
			}
			out << "A^(-1):" << endl;
			out << inverse << endl;
		}

		virtual void printInverseDeficiency(ostream& out) {
			if (!out) {
				throw invalid_argument("Bad output stream in printInverseDeficiency(ostream&).");
			}
			out << "R:" << endl;
			Matrix<double> R = (Matrix<double>(n, 0) + 1) - inverse * initial_A;
			out << R << endl;
			out << "||R|| = " << R.norm() << endl;
		}

		friend istream& operator>>(istream& in, Method& obj) {
			if (!in) {
				throw invalid_argument("Bad input stream in operator>>(istream&, Method&).");
			}
			in.peek();
			if (in.eof()) {
				throw invalid_argument("Empty input file in operator>>(istream&, Method&).");
			}
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

		friend ostream& operator<<(ostream& out, Method& obj) {
			if (!out) {
				throw invalid_argument("Bad output stream in operator<<(ostream&, Method&).");
			}
			for (size_t i = 0; i < obj.n; ++i) {
				out << obj.A[i] << " ";
				out << setw(OUTPUT_WIDTH) << obj.b[i] << endl;
			}
			out << endl;
			return out;
		}

	};

}