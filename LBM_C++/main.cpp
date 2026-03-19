// D2Q5 Lattice Boltzmann Method
#include <iostream>

#define F 5
#define cs_sq (1./3.)
#define wC (1./3.)
#define wNC (1./6.)

#define AR(x, y, f) (f + F*(x) + F*LX*(y))

enum {_C, _N, _S, _W, _E};

constexpr int LX = 128;
constexpr int LY = 128;

double cx[F] = {0, 0, 0, -1, 1};
double cy[F] = {0, 1, -1, 0, 0};
double w[F] = {wC, wNC, wNC, wNC, wNC};

double *cells;
double *temp_cells;

double tau, uLB, vLB;

double fEq(int i, double phi, double u, double v) {
    return w[i] * phi * (1 + (cx[i]*u + cy[i]*v)/cs_sq);
}

void collideAndStream(double *c, double *tc) {
    int x, y, i, xp, xm, yp, ym;
    double fstar[F];
    double phi;

    for (y = 0; y <LY; y++) {
        for (x = 0; x <LX; x++) {
            // Compute macroscopic variables
            phi = 0;
            for (i = 0; i < F; i++) {
                phi += c[AR(x, y, i)];
            }
            // Collision
            for (i = 0; i < F; i++) {
                fstar[i] = c[AR(x, y, i)] + 1./tau * (fEq(i, phi, uLB, vLB) - c[AR(x, y, i)]);
            }
            // Streaming & BC
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

void setInitialCondition(double phii, double ui, double vi, double *c) {
    int x, y, i;
    for (y = 0; y < LY; y++) {
        for (x = 0; x < LX; x++) {
            for (i = 0; i < F; i++) {
                c[AR(x, y, i)] = fEq(i, phii, ui, vi);
            }
        }
    }
    for (i = 0; i < F; i++) {
        c[AR(LX/2, LY/2, i)] = fEq(i, 1.0, ui, vi);
    }
}

int main(int argc, char *argv[]) {
    int iter = 0;
    int maxiter = 1000;

    double phii = 0;
    double ui = 0;
    double vi = 0;

    double D = 0.01;

    cells = new double[LX*LY*F];
    temp_cells = new double[LX*LY*F];

    tau = D/cs_sq + 0.5;

    setInitialCondition(phii, ui, vi, temp_cells);
    // Visualization

    do {
        if (iter % 2)
            collideAndStream(temp_cells, cells);
        else
            collideAndStream(cells, temp_cells);
        iter++;
    }while (iter < maxiter);

    //Visualization

    free(cells);
    free(temp_cells);

    return 0;
}