#include "GaussMatrix.h"
#include <fstream>

using namespace std;

int main() {
	ifstream fIn("Input.txt");
	ofstream fOut("Output.txt", 'a');
	GaussMatrix gauss(5);
	fIn >> gauss;
	gauss.makeUpperTriangular();
	gauss.printUpperTriangular(fOut);
	gauss.getSolution();
	gauss.printSolution(fOut);
	return 0;
}