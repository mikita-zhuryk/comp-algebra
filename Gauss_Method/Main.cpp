#include "GaussMethod.h"
#include <fstream>
#include <string>
#define PATH_TO_MATRIX "../Matrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("GaussOutput.txt");
		CMA::GaussMethod gauss(MATRIX_DIM);
		fIn >> gauss;
		fOut << gauss;
		gauss.run(fOut, true, true);
	}
	catch (invalid_argument ia) {
		cerr << ia.what() << endl;
	}
	catch (out_of_range oor) {
		cerr << oor.what() << endl;
	}
	catch (string s) {
		cerr << s << endl;
	}
	catch (...) {
		cerr << "Unhandled exception." << endl;
	}
	return 0;
}