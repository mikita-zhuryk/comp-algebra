#include "GaussMatrix.h"
#include <fstream>
#define PATH_TO_MATRIX "../Matrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("GaussOutput.txt");
		GaussMatrix gauss(MATRIX_DIM);
		fIn >> gauss;
		fOut << gauss;
		gauss.run(fOut);
	}
	catch (invalid_argument ia) {
		cerr << ia.what() << endl;
	}
	catch (out_of_range oor) {
		cerr << oor.what() << endl;
	}
	catch (...) {
		cerr << "Unhandled exception." << endl;
	}
	return 0;
}