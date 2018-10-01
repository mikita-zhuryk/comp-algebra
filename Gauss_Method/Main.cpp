#include "GaussMatrix.h"
#include <fstream>

using namespace std;

int main() {
	ifstream fIn("Input.txt");
	ofstream fOut("Output.txt", 'a');
	try {
		GaussMatrix gauss(5);
		fIn >> gauss;
		gauss.makeUpperTriangular();
		gauss.printUpperTriangular(fOut);
		gauss.getSolution();
		gauss.printSolution(fOut);
		gauss.calculate();
		system("pause");
	}
	catch (invalid_argument ia) {
		cerr << ia.what() << endl;
	}
	return 0;
}