#include "Danilevsky.h"
#include <fstream>
#define PATH_TO_MATRIX "../EigenMatrix.txt"
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		ifstream fIn(PATH_TO_MATRIX);
		ofstream fOut("DanilevskyOutput.txt");
		CMA::Danilevsky dan(MATRIX_DIM);
		fIn >> dan;
		dan.run(fOut);
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