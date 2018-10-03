#include "GaussMatrix.h"
#include <fstream>

using namespace std;

int main() {
	ifstream fIn("Input.txt");
	ofstream fOut("Output.txt");
	try {
		GaussMatrix gauss(5);
		fIn >> gauss;
		fOut << gauss;
		gauss.makeUpperTriangular();
		gauss.printUpperTriangular(fOut);
		gauss.getSolution();
		gauss.printSolution(fOut);
		gauss.deficiency(fOut);
		fOut << "detA = " << gauss.determinant() << endl << endl;
		system("pause");
	}
	catch (invalid_argument ia) {
		cerr << ia.what() << endl;
	}
	return 0;
}