#include "Pass_Method.h"
#include <fstream>
#define PATH_TO_MATRIX "../MatrixTriDiag.txt"
//#define PATH_TO_MATRIX "../Matrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("PassOutput.txt");
		CMA::PassMethod pass(MATRIX_DIM, false);
		//CMA::PassMethod pass(MATRIX_DIM);
		fIn >> pass;
		fOut << pass;
		pass.run(fOut);
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