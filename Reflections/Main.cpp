#include "Reflections.h"
#include <fstream>
#define PATH_TO_MATRIX "../Matrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("ReflOutput.txt");
		CMA::Reflections refl(MATRIX_DIM);
		fIn >> refl;
		fOut << refl;
		refl.run(fOut, true, false);
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