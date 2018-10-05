#include "GaussMatrix.h"
#include <fstream>
#define DECIMAL_PRECISION 5

using namespace std;

int main() {
	ifstream fIn("Input.txt");
	ofstream fOut("Output.txt");
	try {
		GaussMatrix gauss(5);
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
		gauss.inverseDeficiency(fOut);
		system("pause");
	}
	catch (invalid_argument ia) {
		cerr << ia.what() << endl;
	}
	return 0;
}