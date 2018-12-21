#include "IterJacobi.h"
#include <fstream>
#define PATH_TO_MATRIX "../EigenMatrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("IterJacobiOutput.txt");
		CMA::IterJacobi ja(MATRIX_DIM, -5);
		fIn >> ja;
		ja.run(fOut);
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
	system("pause");
	return 0;
}