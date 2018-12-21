#include "GradientDescent.h"
#include <fstream>
#define PATH_TO_MATRIX "../Matrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("GradientOutput.txt");
		CMA::GradientDescent grad(MATRIX_DIM, -5);
		fIn >> grad;
		fOut << grad;
		grad.run(fOut, false);
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