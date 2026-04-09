#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
using namespace std;

#define F 5
#define cs_sq (1./3.)
#define wC (1./3.)
#define wNC (1./6.)

#define AR(x,y,f) (f+F*(x)+F*LX*(y))
#define ARU(x,y) (x+LX*y)

enum {_C, _N, _S, _W, _E};

const int LX = 128;
const int LY = 128;

double cx[F] = { 0,0,0,-1,1 };
double cy[F] = { 0,1,-1,0,0 };
double w[F] = { wC,wNC,wNC,wNC,wNC };

double *cells;
double *temp_cells;

double tau, uLB, vLB;

double fEq(int i, double phi, double u, double v) {
	return w[i] * phi * (1 + (cx[i] * u + cy[i] * v) / cs_sq);
}

void collideAndStream(double* c, double* tc) { //wskaxniki do twoch tablic
	int x, y, i, xp, yp, xm, ym;
	double fstar[F];
	double phi;

	for (y = 0;y < LY;y++) {
		for (x = 0;x < LX;x++) {
			//compute macroscopic variables
			phi = 0;
			for (i = 0;i < F;i++) {
				phi += c[AR(x, y, i)];
			}
			//collision step
			for (i = 0;i < F;i++) {
				fstar[i] = (1-1./tau)*c[AR(x, y, i)] + fEq(i, phi, uLB, vLB) / tau;
			}
			//streaming step & BC
			//warunki brzegowe periodyczne
			xp = (x + 1) % LX;
			xm = (x - 1 + LX) % LX;
			yp = (y + 1) % LY;
			ym = (y - 1 + LY) % LY;

			tc[AR(xp, y, _E)] = fstar[_E];
			tc[AR(xm, y, _W)] = fstar[_W];
			tc[AR(x, yp, _N)] = fstar[_N];
			tc[AR(x, ym, _S)] = fstar[_S];
			tc[AR(x, y, _C)] = fstar[_C];

		}
	}
}

void setInitialConditions(double phii, double ui, double vi, double* c) {
	int x, y, i;
	for (y = 0; y < LY; y++) {
		for (x = 0;x < LX;x++) {
			//w poźniejszej częsci set velocity to taylor-green vortex
			for (i = 0;i < F; i++) {
				c[AR(x, y, i)] = fEq(i, phii, ui, vi);
			}
		}
	}
	for (y = LY / 2 - 8;y < LY / 2 + 8; y++) {
		for (x = LX / 2 - 8;x < LX / 2 + 8;x++) {
			for (i = 0;i < F;i++) {
				c[AR(x, y, i)] = fEq(i, 1.0, ui, vi);
			}
		}
			
	}
		
}

void dumpStateVTK(const char* fname, double* cells) {
	FILE* fp;
	int x, y;

	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "2D-ADE data file \n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET RECTILINEAR_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", LX, LY, 1);
	fprintf(fp, "X_COORDINATES %d int\n", LX);
	for (x = 0; x < LX; x++)
		fprintf(fp, "%d ", x);
	fprintf(fp, "\n");
	fprintf(fp, "Y_COORDINATES %d int\n", LY);
	for (x = 0; x < LY; x++)
		fprintf(fp, "%d ", x);
	fprintf(fp, "\n");
	fprintf(fp, "Z_COORDINATES 1 int\n");
	fprintf(fp, "0\n");
	fprintf(fp, "POINT_DATA %d \n", LX * LY);
	fprintf(fp, "SCALARS temperature double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (y = 0; y < LY; y++)
		for (x = 0; x < LX; x++) {
			double rho = 0;
			int f;
			for (f = 0; f < F; f++) rho += cells[AR(x, y, f)];
			fprintf(fp, "%e \n", rho);
		}
	fprintf(fp, "\n");
	fclose(fp);
}

int main(int argc, char* argv[]) {
	int iter = 0;
	int maxIter = 1000;

	double Ma, Pe, N;	//dodajemy adwekcje

	Ma = 0.1;
	Pe = 100;
	N = LX;

	double phii = 0;
	double ui = Ma*sqrt(cs_sq);
	double vi = 0;
	double D = ui * N / Pe;
	uLB = ui;

	cells = (double*)malloc(LX * LY * F * sizeof(double));
	temp_cells = (double*)malloc(LX * LY * F * sizeof(double));

	tau = D / cs_sq + 0.5;	//D [m^2/s]?

	setInitialConditions(phii, ui, vi, cells);
	setInitialConditions(phii, ui, vi, temp_cells);
	dumpStateVTK("output_000.vtk", cells);
	dumpStateVTK("output_001.vtk", temp_cells);

	do {
		if (iter % 2)
			collideAndStream(temp_cells, cells);
		else
			collideAndStream(cells, temp_cells);
		iter++;
	} while (iter < maxIter);
	dumpStateVTK("output_final.vtk", cells);

	free(cells);
	free(temp_cells);

	return 0;
}