#include "readPlot3D.h"

#define MESH "box.x"

// Globals Variables
// IJK Dimensions
std::vector<int> dim(3,0);
// X, Y, Z Coordinates of IJK Points
std::vector<double> xs;
std::vector<double> ys;
std::vector<double> zs;
// Area Vectors (A0X, A0Y, A0Z, A0mag, A1X, ...) 
std::vector<double> AI;
std::vector<double> AJ;
std::vector<double> AK;
// Volumes
std::vector<double> V;
std::vector<bool> closedV;


// Read Plot3D Points
bool 
readMesh() {
	std::ifstream f;
	f.open(MESH, std::ifstream::in);
	if (readPlot3D(dim, xs, ys, zs, f)){
        f.close();
		return true;
	}
	else {
        f.close();
		return false;
	}
}


bool
computeArea(std::vector<double> &A, int ndx, int a, int b, int c, int d) {
    /*
        Computes the area of each face by taking one half the 
        cross product of the diagonals, then divides
        each component by the magnitude to normalize it.
        The area vector stores AX,AY,AZ,Amag serially.
    */
    double dx1 = xs[a] - xs[c];
    double dx2 = xs[b] - xs[d];
    
    double dy1 = ys[a] - ys[c];
    double dy2 = ys[b] - ys[d];
    
    double dz1 = zs[a] - zs[c];
    double dz2 = zs[b] - zs[d];

    double AX = 0.5 * ((dy1 * dz2) - (dy2 * dz1));
    double AY = 0.5 * ((dz1 * dx2) - (dx1 * dz2));
    double AZ = 0.5 * ((dx1 * dy2) - (dx2 * dy1));
    double Amag = sqrt((AX * AX) + (AY * AY) + (AZ * AZ));
    if (Amag < 1E-15) {
        printf("zero area\n");
        return false;
    }
    AX = AX / Amag;
    AY = AY / Amag;
    AZ = AZ / Amag;
    A[ndx] = AX;
    A[ndx + 1] = AY;
    A[ndx + 2] = AZ;
    A[ndx + 3] = Amag;

    return true;
}


bool 
computeFaceAreas() { 
    /*
      Iterate throught he I,J,K faces and calculate the 
      area of each one using the helper funciton computeArea.
    */
    int NI = dim[0];
    int NJ = dim[1];
    int NK = dim[2];
    int size = NI * NJ * NK;  
    int ndx = 0;
    // Indices for four face points
    int a, b, c, d;

    // Appropriately resize Area Vectors
    int I_size = 4 * (NI * (NJ - 1) * (NK - 1));
    int J_size = 4 * ((NI - 1) * NJ * (NK - 1));
    int K_size = 4 * ((NI - 1) * (NJ - 1) * NK);
    AI.resize(I_size, 0);
    AJ.resize(J_size, 0);
    AK.resize(K_size, 0);

    // Compute I-const normals...
    for (int K = 0; K < NK-1; ++K) {
        for (int J = 0; J < NJ-1; ++J) {
            for (int I = 0; I < NI; ++I) {
                a = IJKtoXYZ(I, J, K, NI, NJ, NK);
                b = IJKtoXYZ(I, J, K+1, NI, NJ, NK);
                c = IJKtoXYZ(I, J+1, K+1, NI, NJ, NK);
                d = IJKtoXYZ(I, J+1, K, NI, NJ, NK);
                
                computeArea(AI, ndx, a, b, c, d);
                ndx += 4;
            }
        }
    }

    // Compute J-const normals...
    ndx = 0;
    for (int K = 0; K < NK - 1; ++K) {
        for (int I = 0; I < NI - 1; ++I) {
            for (int J = 0; J < NJ; ++J) {
                a = IJKtoXYZ(I, J, K, NI, NJ, NK);
                b = IJKtoXYZ(I+1, J, K, NI, NJ, NK);
                c = IJKtoXYZ(I+1, J, K+1, NI, NJ, NK);
                d = IJKtoXYZ(I, J, K+1, NI, NJ, NK);

                computeArea(AJ, ndx, a, b, c, d);
                ndx += 4;
            }
        }
    }

    // Compute K-const normals...
    ndx = 0;
    for (int J = 0; J < NJ - 1; ++J) {
        for (int I = 0; I < NI - 1; ++I) {
            for (int K = 0; K < NK; ++K) {
                a = IJKtoXYZ(I, J, K, NI, NJ, NK);
                b = IJKtoXYZ(I + 1, J, K, NI, NJ, NK);
                c = IJKtoXYZ(I + 1, J + 1, K, NI, NJ, NK);
                d = IJKtoXYZ(I, J + 1, K, NI, NJ, NK);

                computeArea(AK, ndx, a, b, c, d);
                ndx += 4;
            }
        }
    }

    return true; 
}


bool
computeCentroids(std::vector<double> &cent, int I, int J, int K) {
    /*
        Computes the part centroid each face of the I, J, Kth cell.
        Note that not all components of the centroid are needed
        to compute the volume of the cell.
        The required components are laid out serially in the cent 
        vector (for efficiency, not necessarily clarity).
        Also checks if the volume is closed by summing the 
        area vectors of each face.
    */
    int NI = dim[0];
    int NJ = dim[1];
    int NK = dim[2];
    int a, b, c, d;
    double closed = 0;

    // I- Face
    a = IJKtoXYZ(I, J, K, NI, NJ, NK);
    b = IJKtoXYZ(I, J, K + 1, NI, NJ, NK);
    c = IJKtoXYZ(I, J + 1, K + 1, NI, NJ, NK);
    d = IJKtoXYZ(I, J + 1, K, NI, NJ, NK);
    cent[0] = (xs[a] + xs[b] + xs[c] + xs[d]) / 4.0;

    // I+ Face
    a = IJKtoXYZ(I + 1, J, K, NI, NJ, NK);
    b = IJKtoXYZ(I + 1, J, K + 1, NI, NJ, NK);
    c = IJKtoXYZ(I + 1, J + 1, K + 1, NI, NJ, NK);
    d = IJKtoXYZ(I + 1, J + 1, K, NI, NJ, NK);
    cent[1] = (xs[a] + xs[b] + xs[c] + xs[d]) / 4.0;

    // J- Face
    a = IJKtoXYZ(I, J, K, NI, NJ, NK);
    b = IJKtoXYZ(I + 1, J, K, NI, NJ, NK);
    c = IJKtoXYZ(I + 1, J, K + 1, NI, NJ, NK);
    d = IJKtoXYZ(I, J, K + 1, NI, NJ, NK);
    cent[2] = (ys[a] + ys[b] + ys[c] + ys[d]) / 4.0;
    cent[3] = (zs[a] + zs[b] + zs[c] + zs[d]) / 4.0;

    // J+ Face
    a = IJKtoXYZ(I, J + 1, K, NI, NJ, NK);
    b = IJKtoXYZ(I + 1, J + 1, K, NI, NJ, NK);
    c = IJKtoXYZ(I + 1, J + 1, K + 1, NI, NJ, NK);
    d = IJKtoXYZ(I, J + 1, K + 1, NI, NJ, NK);
    cent[4] = (ys[a] + ys[b] + ys[c] + ys[d]) / 4.0;
    cent[5] = (zs[a] + zs[b] + zs[c] + zs[d]) / 4.0;

    // K- Face
    a = IJKtoXYZ(I, J, K, NI, NJ, NK);
    b = IJKtoXYZ(I + 1, J, K, NI, NJ, NK);
    c = IJKtoXYZ(I + 1, J + 1, K, NI, NJ, NK);
    d = IJKtoXYZ(I, J + 1, K, NI, NJ, NK);
    cent[6] = (ys[a] + ys[b] + ys[c] + ys[d]) / 4.0;
    cent[7] = (zs[a] + zs[b] + zs[c] + zs[d]) / 4.0;

    // K+ Face
    a = IJKtoXYZ(I, J, K + 1, NI, NJ, NK);
    b = IJKtoXYZ(I + 1, J, K + 1, NI, NJ, NK);
    c = IJKtoXYZ(I + 1, J + 1, K + 1, NI, NJ, NK);
    d = IJKtoXYZ(I, J + 1, K + 1, NI, NJ, NK);
    cent[8] = (ys[a] + ys[b] + ys[c] + ys[d]) / 4.0;
    cent[9] = (zs[a] + zs[b] + zs[c] + zs[d]) / 4.0;
    
    // Check if Volume is closed by summing the area vectors of each face.
    int LoIndx = 4*I;
    int HiIndx = LoIndx + 4;
    int LoJndx = 4*((K * NJ * (NI - 1)) + (NJ * I) + J);
    int HiJndx = LoJndx + 4;
    int LoKndx = 4*((J*NK*(NI-1)) + (NK*I) + K);
    int HiKndx = LoKndx + 4;

    closed += AI[LoIndx + 3] * (AI[LoIndx] + AI[LoIndx + 1] + AI[LoIndx + 2]);
    closed -= AI[HiIndx + 3] * (AI[HiIndx] + AI[HiIndx + 1] + AI[HiIndx + 2]);
    closed += AJ[LoJndx + 3] * (AJ[LoJndx] + AJ[LoJndx + 1] + AJ[LoJndx + 2]);
    closed -= AJ[HiJndx + 3] * (AJ[HiJndx] + AJ[HiJndx + 1] + AJ[HiJndx + 2]);
    closed += AK[LoKndx + 3] * (AK[LoKndx] + AK[LoKndx + 1] + AK[LoKndx + 2]);
    closed -= AK[HiKndx + 3] * (AK[HiKndx] + AK[HiKndx + 1] + AK[HiKndx + 2]);

    if (abs(closed) < 1E-14) {
        return true;
    }
    else {
        return false;
    }
}


bool 
computeVolumes() { 
    /*
        Loops through cells computing the volume of each one,
        as well as checking if it is closed. 
        V = [(dy/dJ * dz/dK) - (dy/dk * dz/dJ)] / dI/dx 
    */
    int NI = dim[0];
    int NJ = dim[1];
    int NK = dim[2];
    int numCells = (NI - 1) * (NJ - 1) * (NK - 1);
    int numClosed = 0;
    double dy_J, dz_K, dy_K, dz_J, dI_x;
    V.resize(numCells, 0);
    closedV.resize(numCells + 1, false);

    // Info Pertaining to Centroid of Each Face
    std::vector<double> cents(10,0);

    int ndx = 0;
    for (int K = 0; K < NK - 1; K++) {
        for (int J = 0; J < NJ - 1; J++) {
            for (int I = 0; I < NI - 1; I++) {
                // Compute centroid of each face
                if (computeCentroids(cents, I, J, K)) {
                    // Volume is closed                    
                    closedV[ndx] = true;
                    numClosed++;
                }

                // Compute Volume
                dy_J = cents[4] - cents[2];
                dz_K = cents[9] - cents[7];
                dy_K = cents[8] - cents[6];
                dz_J = cents[5] - cents[3];

             
                dI_x = 1.0 / (cents[1] - cents[0]);
                V[ndx] = ((dy_J * dz_K) - (dy_K * dz_J)) / dI_x;

                ndx++;
            }
        }
    }  
    if (numCells == numClosed) {
        // All cells are closed
        closedV[ndx] = true;
    }
        
    return true; 
}


void 
printResults() {
    printf("Printing Face Areas...\n\n");
    printf("Printing I const Faces...\n");
    auto it = AI.begin();
    for (it; it != AI.end(); it += 4) {
        printf("Normal (%f, %f, %f), Area: %f\n", *it, *(it + 1), *(it + 2), *(it + 3));
    }
    printf("\n");

    printf("Printing J const Faces...\n");
    it = AJ.begin();
    for (it; it != AJ.end(); it += 4) {
        printf("Normal (%f, %f, %f), Area: %f\n", *it, *(it + 1), *(it + 2), *(it + 3));
    }
    printf("\n");

    printf("Printing K const Faces...\n");
    it = AK.begin();
    for (it; it != AK.end(); it += 4) {
        printf("Normal (%f, %f, %f), Area: %f\n", *it, *(it + 1), *(it + 2), *(it + 3));
    }
    printf("\n\n");

    printf("Printing Volumes...\n\n");
    it = V.begin();
    int ndx;
    for (it; it != V.end()-1; ++it) {
        ndx = it - V.begin();
        printf("Volume: %f  Closed: %d\n", *it, int(closedV[ndx]));
    }
    if (true == *(closedV.end() - 1)) {
        printf("\nAll Volumes are Closed!\n");
    }
    else {
        printf("\nCheck for Unclosed Volumes\n");
    }

}


/*
	Read Plot3D mesh.
	Calculate Face Areas and Normals
	Calculate Cell Volumes
*/
int main() {
	if (readMesh() && computeFaceAreas() && 
        computeVolumes()) {
        
        printResults();
		return 0;
	}
	else {
		return -1;
	}
}