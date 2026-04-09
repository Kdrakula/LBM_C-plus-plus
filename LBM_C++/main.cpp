// D2Q5 Lattice Boltzmann Method
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <vector>

#define F 5
#define cs_sq (1./3.)
#define wC (1./3.)
#define wNC (1./6.)

#define AR(x,y,f) (f+F*(x)+F*LX*(y))

enum {_C, _N, _S, _W, _E};

constexpr int LX = 128;
constexpr int LY = 128;

double cx[F] = {0, 0, 0, -1, 1};
double cy[F] = {0, 1, -1, 0, 0};
double w[F] = {wC, wNC, wNC, wNC, wNC};

double tau, uLB, vLB;

double fEq(int i, double phi, double u, double v) {
	return w[i] * phi * (1 + (cx[i] * u + cy[i] * v) / cs_sq);
}

void collideAndStream(double *c, double *tc) {
    int xp, xm, yp, ym;
    double fstar[F];
    double phi;

    for (int y = 0; y < LY; y++) {
        for (int x = 0; x < LX; x++) {
            // Compute macroscopic variables
            phi = 0;
            for (int i = 0; i < F; i++) {
                phi += c[AR(x, y, i)];
            }

            // Collision
            for (int i = 0; i < F; i++) {
                fstar[i] = (1-1./tau)*c[AR(x, y, i)] + fEq(i, phi, uLB, vLB)/tau;
                // fstar[i] = c[AR(x, y, i)] + 1./tau * (fEq(i, phi, uLB, vLB) - c[AR(x, y, i)]);
            }

            // Streaming & Periodic Boundary Conditions
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
    for (int y = 0; y < LY; y++) {
        for (int x = 0; x < LX; x++) {
            for (int i = 0; i < F; i++) {
                c[AR(x, y, i)] = fEq(i, phii, ui, vi);
            }
        }
    }
    for (int y = LY / 2 - 8;y < LY / 2 + 8; y++) {
        for (int x = LX / 2 - 8;x < LX / 2 + 8;x++) {
            for (int i = 0;i < F;i++) {
                c[AR(x, y, i)] = fEq(i, 1.0, ui, vi);
            }
        }

    }
}

void savePPM(double *c, int iter, const std::string& output_dir) {
    std::ostringstream oss;
    oss << output_dir << "/frame_" << std::setw(4) << std::setfill('0') << iter << ".ppm";
    std::string filename = oss.str();
    std::ofstream out(filename);

    if (!out) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write PPM P3 Header
    out << "P3\n" << LX << " " << LY << "\n255\n";

    // Find max phi to dynamically scale the brightness
    double max_phi = 0.0;
    for (int y = 0; y < LY; y++) {
        for (int x = 0; x < LX; x++) {
            double phi = 0;
            for (int i = 0; i < F; i++) {
                phi += c[AR(x, y, i)];
            }
            if (phi > max_phi) max_phi = phi;
        }
    }
    if (max_phi == 0.0) max_phi = 1.0; // Prevent division by zero

    // Write pixel data (Grayscale mapping)
    for (int y = 0; y < LY; y++) {
        for (int x = 0; x < LX; x++) {
            double phi = 0;
            for (int i = 0; i < F; i++) {
                phi += c[AR(x, y, i)];
            }

            // Map phi to a 0-255 grayscale color
            int color = static_cast<int>((phi / max_phi) * 255.0);
            color = std::max(0, std::min(255, color)); // Clamp values

            out << color << " " << color << " " << color << " ";
        }
        out << "\n";
    }
    out.close();
    std::cout << "Saved " << filename << std::endl;
}

int main(int argc, char *argv[]) {
    int iter = 0;
    int maxiter = 1000;
    int save_interval = 10;

    std::string output_dir = "simulation_output";
    std::filesystem::create_directories(output_dir);

    double Ma = 0.1;
    double Pe = 100;
    double N = LX;

    double phii = 0;
    double ui = Ma*sqrt(cs_sq);
    double vi = 0;
    double D = ui*N/Pe;
    uLB = ui;

    std::vector<double> cells(LX * LY * F);
    std::vector<double> temp_cells(LX * LY * F);

    tau = D/cs_sq + 0.5;

    setInitialCondition(phii, ui, vi, temp_cells.data());

    savePPM(temp_cells.data(), iter, output_dir);

    do {
        if (iter % 2 == 0)
            collideAndStream(temp_cells.data(), cells.data());
        else
            collideAndStream(cells.data(), temp_cells.data());

        iter++;

        if (iter % save_interval == 0) {
            if (iter % 2 == 0)
                savePPM(temp_cells.data(), iter, output_dir);
            else
                savePPM(cells.data(), iter, output_dir);
        }

    } while (iter < maxiter);

    return 0;
}