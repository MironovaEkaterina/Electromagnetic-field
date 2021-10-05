#include <iostream>
#include <vector>
#include <cmath>

const double dx = 1E+9;
const double dy = 1E+10;
const double c = 3E+10;
const double ax = 0, bx = 5E+10, ay = 0, by = 5E+10;
#define M_PI 3.1415926535897932384626433832795

struct Field_characteristics {
public:
	std::vector<std::vector<std::vector<double>>> E;
	std::vector<std::vector<std::vector<double>>> B;
};

class Field {
public:
	Field_characteristics v;
	std::vector<std::vector<std::vector<double>>> J;
	int nx = (bx - ax) / dx + 1;
	int ny = (by - ay) / dy + 1;

	Field() {
		v.E.assign(3, std::vector<std::vector<double>>(nx, std::vector<double>(ny)));
		v.B.assign(3, std::vector<std::vector<double>>(nx, std::vector<double>(ny)));
		J.assign(3, std::vector<std::vector<double>>(nx, std::vector<double>(ny)));
		//начальные значения в узлах сетки и плотность тока
		for (int k = 0; k < 3; k++) 
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++){
						J[k][i][j] = 0;
						v.E[k][i][j] = 0;
						v.B[k][i][j] = 0;
				}

	};

	Field_characteristics Field::FDTD(double dt) { //обновление сеточных значений

	//обновление Е
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++) {
				v.E[0][i][j] = v.E[0][i][j] - 4.0 * M_PI * dt * J[0][i][j] + c * dt * (v.B[2][i][j + 1] - v.B[2][i][j - 1]) / (2.0 * dy);
				v.E[1][i][j] = v.E[1][i][j] - 4.0 * M_PI * dt * J[1][i][j] - c * dt * (v.B[2][i + 1][j] - v.B[2][i - 1][j]) / (2.0 * dx);
				v.E[2][i][j] = v.E[2][i][j] - 4.0 * M_PI * dt * J[2][i][j] + c * dt * ((v.B[1][i + 1][j] - v.B[1][i - 1][j]) / (2.0 * dx) - (v.B[0][i][j + 1] - v.B[0][i][j - 1]) / (2.0 * dy));
			}

		for (int i = 0; i < nx; i++) { // j==0
			v.E[0][i][0] = v.E[0][i][0] - 4.0 * M_PI * dt * J[0][i][0] + c * dt * (v.B[2][i][1] - v.B[2][i][ny - 1]) / (2.0 * dy);
			if ((i != 0)&&(i!=nx-1)) {
				v.E[1][i][0] = v.E[1][i][0] - 4.0 * M_PI * dt * J[1][i][0] - c * dt * (v.B[2][i + 1][0] - v.B[2][i - 1][0]) / (2.0 * dx);
				v.E[2][i][0] = v.E[2][i][0] - 4.0 * M_PI * dt * J[2][i][0] + c * dt * ((v.B[1][i + 1][0] - v.B[1][i - 1][0]) / (2.0 * dx) - (v.B[0][i][1] - v.B[0][i][ny - 1]) / (2.0 * dy));
			}
			if (i==0) {
				v.E[1][i][0] = v.E[1][i][0] - 4.0 * M_PI * dt * J[1][i][0] - c * dt * (v.B[2][i + 1][0] - v.B[2][nx - 1][0]) / (2.0 * dx);
				v.E[2][i][0] = v.E[2][i][0] - 4.0 * M_PI * dt * J[2][i][0] + c * dt * ((v.B[1][i + 1][0] - v.B[1][nx - 1][0]) / (2.0 * dx) - (v.B[0][i][1] - v.B[0][i][ny - 1]) / (2.0 * dy));
			} 
			if (i == nx - 1) {
				v.E[1][i][0] = v.E[1][i][0] - 4.0 * M_PI * dt * J[1][i][0] - c * dt * (v.B[2][0][0] - v.B[2][i - 1][0]) / (2.0 * dx);
				v.E[2][i][0] = v.E[2][i][0] - 4.0 * M_PI * dt * J[2][i][0] + c * dt * ((v.B[1][0][0] - v.B[1][i - 1][0]) / (2.0 * dx) - (v.B[0][i][1] - v.B[0][i][ny - 1]) / (2.0 * dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j==ny-1
			v.E[0][i][ny - 1] = v.E[0][i][ny - 1] - 4.0 * M_PI * dt * J[0][i][ny - 1] + c * dt * (v.B[2][i][0] - v.B[2][i][ny - 2]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.E[1][i][ny - 1] = v.E[1][i][ny - 1] - 4.0 * M_PI * dt * J[1][i][ny - 1] - c * dt * (v.B[2][i + 1][ny - 1] - v.B[2][i - 1][ny - 1]) / (2.0 * dx);
				v.E[2][i][ny - 1] = v.E[2][i][ny - 1] - 4.0 * M_PI * dt * J[2][i][ny - 1] + c * dt * ((v.B[1][i + 1][ny - 1] - v.B[1][i - 1][ny - 1]) / (2.0 * dx) - (v.B[0][i][0] - v.B[0][i][ny - 2]) / (2.0 * dy));
			}
			if (i == 0) {
				v.E[1][i][ny - 1] = v.E[1][i][ny - 1] - 4.0 * M_PI * dt * J[1][i][ny - 1] - c * dt * (v.B[2][i + 1][ny - 1] - v.B[2][nx - 1][ny - 1]) / (2.0 * dx);
				v.E[2][i][ny - 1] = v.E[2][i][ny - 1] - 4.0 * M_PI * dt * J[2][i][ny - 1] + c * dt * ((v.B[1][i + 1][ny - 1] - v.B[1][nx - 1][ny - 1]) / (2.0 * dx) - (v.B[0][i][0] - v.B[0][i][ny - 2]) / (2.0 * dy));
			}
			if (i == nx - 1) {
				v.E[1][i][ny - 1] = v.E[1][i][ny - 1] - 4.0 * M_PI * dt * J[1][i][ny - 1] - c * dt * (v.B[2][0][ny - 1] - v.B[2][i - 1][ny - 1]) / (2.0 * dx);
				v.E[2][i][ny - 1] = v.E[2][i][ny - 1] - 4.0 * M_PI * dt * J[2][i][ny - 1] + c * dt * ((v.B[1][0][ny - 1] - v.B[1][i - 1][ny - 1]) / (2.0 * dx) - (v.B[0][i][0] - v.B[0][i][ny - 2]) / (2.0 * dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) { // i==0
			v.E[0][0][j] = v.E[0][0][j] - 4.0 * M_PI * dt * J[0][0][j] + c * dt * (v.B[2][0][j + 1] - v.B[2][0][j - 1]) / (2.0 * dy);
			v.E[1][0][j] = v.E[1][0][j] - 4.0 * M_PI * dt * J[1][0][j] - c * dt * (v.B[2][1][j] - v.B[2][nx - 1][j]) / (2.0 * dx);
			v.E[2][0][j] = v.E[2][0][j] - 4.0 * M_PI * dt * J[2][0][j] + c * dt * ((v.B[1][1][j] - v.B[1][nx - 1][j]) / (2.0 * dx) - (v.B[0][0][j + 1] - v.B[0][0][j - 1]) / (2.0 * dy));  
		}

		for (int j = 1; j < ny - 1; j++) { // i==nx-1
			v.E[0][nx - 1][j] = v.E[0][nx - 1][j] - 4.0 * M_PI * dt * J[0][nx - 1][j] + c * dt * (v.B[2][nx - 1][j + 1] - v.B[2][nx - 1][j - 1]) / (2.0 * dy);
			v.E[1][nx - 1][j] = v.E[1][nx - 1][j] - 4.0 * M_PI * dt * J[1][nx - 1][j] - c * dt * (v.B[2][0][j] - v.B[2][nx - 2][j]) / (2.0 * dx);
			v.E[2][nx - 1][j] = v.E[2][nx - 1][j] - 4.0 * M_PI * dt * J[2][nx - 1][j] + c * dt * ((v.B[1][0][j] - v.B[1][nx - 2][j]) / (2.0 * dx) - (v.B[0][nx - 1][j + 1] - v.B[0][nx - 1][j - 1]) / (2.0 * dy));  
		}

	//обновление В
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++){
				v.B[0][i][j] = v.B[0][i][j] - c * dt * (v.E[2][i][j + 1] - v.E[2][i][j - 1]) / (2.0 * dy);
				v.B[1][i][j] = v.B[1][i][j] + c * dt * (v.E[2][i + 1][j] - v.E[2][i - 1][j]) / (2.0 * dx);
				v.B[2][i][j] = v.B[2][i][j] - c * dt * ((v.E[1][i + 1][j] - v.E[1][i - 1][j]) / (2.0 * dx) - (v.E[0][i][j + 1] - v.E[0][i][j - 1]) / (2.0 * dy));
			}

		for (int i = 0; i < nx; i++) { // j==0
			v.B[0][i][0] = v.B[0][i][0] - c * dt * (v.E[2][i][1] - v.E[2][i][ny - 1]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.B[1][i][0] = v.B[1][i][0] + c * dt * (v.E[2][i + 1][0] - v.E[2][i - 1][0]) / (2.0 * dx);
				v.B[2][i][0] = v.B[2][i][0] - c * dt * ((v.E[1][i + 1][0] - v.E[1][i - 1][0]) / (2.0 * dx) - (v.E[0][i][1] - v.E[0][i][ny - 1]) / (2.0 * dy));
			}
			if (i == 0) {
				v.B[1][i][0] = v.B[1][i][0] + c * dt * (v.E[2][i + 1][0] - v.E[2][nx - 1][0]) / (2.0 * dx);
				v.B[2][i][0] = v.B[2][i][0] - c * dt * ((v.E[1][i + 1][0] - v.E[1][nx - 1][0]) / (2.0 * dx) - (v.E[0][i][1] - v.E[0][i][ny - 1]) / (2.0 * dy));
			}
			if (i == nx-1) {
				v.B[1][i][0] = v.B[1][i][0] + c * dt * (v.E[2][0][0] - v.E[2][i - 1][0]) / (2.0 * dx);
				v.B[2][i][0] = v.B[2][i][0] - c * dt * ((v.E[1][0][0] - v.E[1][i - 1][0]) / (2.0 * dx) - (v.E[0][i][1] - v.E[0][i][ny - 1]) / (2.0 * dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j==ny-1
			v.B[0][i][ny - 1] = v.B[0][i][ny - 1] - c * dt * (v.E[2][i][0] - v.E[2][i][ny - 2]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.B[1][i][ny - 1] = v.B[1][i][ny - 1] + c * dt * (v.E[2][i + 1][ny - 1] - v.E[2][i - 1][ny - 1]) / (2.0 * dx);
				v.B[2][i][ny - 1] = v.B[2][i][ny - 1] - c * dt * ((v.E[1][i + 1][ny - 1] - v.E[1][i - 1][ny - 1]) / (2.0 * dx) - (v.E[0][i][0] - v.E[0][i][ny - 2]) / (2.0 * dy));
			}
			if (i == 0) {
				v.B[1][i][ny - 1] = v.B[1][i][ny - 1] + c * dt * (v.E[2][i + 1][ny - 1] - v.E[2][nx - 1][ny - 1]) / (2.0 * dx);
				v.B[2][i][ny - 1] = v.B[2][i][ny - 1] - c * dt * ((v.E[1][i + 1][ny - 1] - v.E[1][nx - 1][ny - 1]) / (2.0 * dx) - (v.E[0][i][0] - v.E[0][i][ny - 2]) / (2.0 * dy));
			}
			if (i == nx - 1) {
				v.B[1][i][ny - 1] = v.B[1][i][ny - 1] + c * dt * (v.E[2][0][ny - 1] - v.E[2][i - 1][ny - 1]) / (2.0 * dx);
				v.B[2][i][ny - 1] = v.B[2][i][ny - 1] - c * dt * ((v.E[1][0][ny - 1] - v.E[1][i - 1][ny - 1]) / (2.0 * dx) - (v.E[0][i][0] - v.E[0][i][ny - 2]) / (2.0 * dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) {  // i==0
			v.B[0][0][j] = v.B[0][0][j] - c * dt * (v.E[2][0][j + 1] - v.E[2][0][j - 1]) / (2.0 * dy);
			v.B[1][0][j] = v.B[1][0][j] + c * dt * (v.E[2][1][j] - v.E[2][nx - 1][j]) / (2.0 * dx);
			v.B[2][0][j] = v.B[2][0][j] - c * dt * ((v.E[1][1][j] - v.E[1][nx - 1][j]) / (2.0 * dx) - (v.E[0][0][j + 1] - v.E[0][0][j - 1]) / (2.0 * dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i==nx-1
			v.B[0][nx - 1][j] = v.B[0][nx - 1][j] - c * dt * (v.E[2][nx - 1][j + 1] - v.E[2][nx - 1][j - 1]) / (2.0 * dy);
			v.B[1][nx - 1][j] = v.B[1][nx - 1][j] + c * dt * (v.E[2][0][j] - v.E[2][nx - 2][j]) / (2.0 * dx);
			v.B[2][nx - 1][j] = v.B[2][nx - 1][j] - c * dt * ((v.E[1][0][j] - v.E[1][nx - 2][j]) / (2.0 * dx) - (v.E[0][nx - 1][j + 1] - v.E[0][nx - 1][j - 1]) / (2.0 * dy));
		}

		return (*this).v;
	}
};
