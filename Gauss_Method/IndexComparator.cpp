#include "IndexComparator.h"

using namespace std;

bool IndexComparator::operator()(pair<size_t, double> x_i, pair<size_t, double> x_j) {
	return x_i.first > x_j.first;
}