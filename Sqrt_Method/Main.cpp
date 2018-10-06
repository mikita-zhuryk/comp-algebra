#include "Sqrt_Method.h"
#include <fstream>
#define PATH_TO_MATRIX "../Matrix.txt"
#define DECIMAL_PRECISION 5
#define MATRIX_DIM 5

using namespace std;

int main() {
	try {
		Sqrt_Method sqrt(MATRIX_DIM);
		sqrt.printSolution(cout);
	}
	catch (...) {}
	return 0;
}