#ifndef readPlot3D_H
#define readPlot3D_H

#include <fstream>
#include <vector>

bool 
readHeader(std::vector<int> &dim, std::ifstream &file) {
	if (file.good()) {
		file >> dim[0];
		file >> dim[1];
		file >> dim[2];

		return true;
	}
	else {
		return false;
	}
}


bool
readCoords(std::vector<double> &coords, const int& size, std::ifstream &file) {
	double d = 0;
	if (file.good()) {
		for (int i = 0; i < size; i++) {
			file >> d;
			coords[i] = d;
		}
		return true;
	}
	else {
		return false;
	}
}


bool
readPlot3D(std::vector<int> &dim,
           std::vector<double> &xs, 
		   std::vector<double> &ys, 
	       std::vector<double> &zs, std::ifstream &file) {
	/*
		Read in coordinates into x, y, z arrays respectively. 
	*/
	if (!file) {
		return false;
	}

    int size = 0;
	if (readHeader(dim, file)) {
        size = dim[0] * dim[1] * dim[2];
		xs.resize(size, 0);
		ys.resize(size, 0);
		zs.resize(size, 0);
	}
	else {
		return false;
	}
	if (
		readCoords(xs, size, file) &&
		readCoords(ys, size,file) &&
		readCoords(zs, size, file)) {
		
		return true;
	}
	else {
		return false;
	}
}


inline int 
IJKtoXYZ(const int I, const int J, const int K, 
    const int NI, const int NJ, const int NK) {
    int ndx = (NI * NJ * K) + (NI * J) + I;

	return ndx;
}

#endif