#include "GaussMatrix.h"
#include <fstream>
#define PATH_TO_MATRIX "../Matrix.txt"
#define DECIMAL_PRECISION 5
#define MATRIX_DIM 5

using namespace std;

int main() {
	ifstream fIn(PATH_TO_MATRIX);
	ofstream fOut("Output.txt");
	try {
		GaussMatrix gauss(MATRIX_DIM);
		fIn >> gauss;
		fOut << gauss;
		fOut.setf(ios_base::fixed, ios_base::floatfield);
		fOut.precision(DECIMAL_PRECISION);
		gauss.makeUpperTriangular();
		gauss.printUpperTriangular(fOut);
		gauss.getSolution();
		gauss.printSolution(fOut);
		fOut.setf(ios_base::scientific, ios_base::floatfield);
		gauss.deficiency(fOut);
		fOut.setf(ios_base::fixed, ios_base::floatfield);
		fOut << "detA = " << gauss.determinant() << endl << endl;
		gauss.printInverse(fOut);
		fOut.setf(ios_base::scientific, ios_base::floatfield);
		gauss.printInverseDeficiency(fOut);
	}
	catch (invalid_argument ia) {
		cerr << ia.what() << endl;
	}
	catch (out_of_range oor) {
		cerr << oor.what() << endl;
	}
	return 0;
}