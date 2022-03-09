#include <iostream>
#include <vector>
#include <cmath>


const double c = 3E+10;
#define M_PI 3.1415926535897932384626433832795

struct Field_characteristics {
public:
	std::vector<std::vector<std::vector<double>>> Ex;
	std::vector<std::vector<std::vector<double>>> Ey;
	std::vector<std::vector<std::vector<double>>> Ez;
	std::vector<std::vector<std::vector<double>>> Bx;
	std::vector<std::vector<std::vector<double>>> By;
	std::vector<std::vector<std::vector<double>>> Bz;
};

class Field {
public:
	Field_characteristics v;
	std::vector<std::vector<std::vector<double>>> Jx;
	std::vector<std::vector<std::vector<double>>> Jy;
	std::vector<std::vector<std::vector<double>>> Jz;
	double dx, dy, dz;
	int nx, ny,nz;
	double ax, bx, ay, by,az,bz;

	Field(double dx,double dy,double dz,double ax,double ay, double az, double bx, double by, double bz) {
		nx = (bx - ax) / dx + 1;
		ny = (by - ay) / dy + 1;
		nz = (bz - az) / dz + 1;
		v.Ex.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz,0)));
		v.Bx.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz,0)));
		Jx.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz,0)));
		v.Ey.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0)));
		v.By.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0)));
		Jy.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0)));
		v.Ez.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0)));
		v.Bz.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0)));
		Jz.assign(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0)));
	};

	void Field::border_update(double dt) { //обновление граничных значений
		//обновление граничных значений E
		for (int i = 0; i < nx; i++) { // j==0
			v.Ex[i][0][0] = v.Ex[i][0][0] - 4.0 * M_PI * dt * Jx[i][0][0] + c * dt * (v.Bz[i][1][0] - v.Bz[i][ny - 1][0]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.Ey[i][0][0] = v.Ey[i][0][0] - 4.0 * M_PI * dt * Jy[i][0][0] - c * dt * (v.Bz[i + 1][0][0] - v.Bz[i - 1][0][0]) / (2.0 * dx);
				v.Ez[i][0][0] = v.Ez[i][0][0] - 4.0 * M_PI * dt * Jz[i][0][0] + c * dt * ((v.By[i + 1][0][0] - v.By[i - 1][0][0]) / (2.0 * dx) - (v.Bx[i][1][0] - v.Bx[i][ny - 1][0]) / (2.0 * dy));
			}
			if (i == 0) {
				v.Ey[i][0][0] = v.Ey[i][0][0] - 4.0 * M_PI * dt * Jy[i][0][0] - c * dt * (v.Bz[i + 1][0][0] - v.Bz[nx - 1][0][0]) / (2.0 * dx);
				v.Ez[i][0][0] = v.Ez[i][0][0] - 4.0 * M_PI * dt * Jz[i][0][0] + c * dt * ((v.By[i + 1][0][0] - v.By[nx - 1][0][0]) / (2.0 * dx) - (v.Bx[i][1][0] - v.Bx[i][ny - 1][0]) / (2.0 * dy));
			}
			if (i == nx - 1) {
				v.Ey[i][0][0] = v.Ey[i][0][0] - 4.0 * M_PI * dt * Jy[i][0][0] - c * dt * (v.Bz[0][0][0] - v.Bz[i - 1][0][0]) / (2.0 * dx);
				v.Ez[i][0][0] = v.Ez[i][0][0] - 4.0 * M_PI * dt * Jz[i][0][0] + c * dt * ((v.By[0][0][0] - v.By[i - 1][0][0]) / (2.0 * dx) - (v.Bx[i][1][0] - v.Bx[i][ny - 1][0]) / (2.0 * dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j==ny-1
			v.Ex[i][ny - 1][0] = v.Ex[i][ny - 1][0] - 4.0 * M_PI * dt * Jx[i][ny - 1][0] + c * dt * (v.Bz[i][0][0] - v.Bz[i][ny - 2][0]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.Ey[i][ny - 1][0] = v.Ey[i][ny - 1][0] - 4.0 * M_PI * dt * Jy[i][ny - 1][0] - c * dt * (v.Bz[i + 1][ny - 1][0] - v.Bz[i - 1][ny - 1][0]) / (2.0 * dx);
				v.Ez[i][ny - 1][0] = v.Ez[i][ny - 1][0] - 4.0 * M_PI * dt * Jz[i][ny - 1][0] + c * dt * ((v.By[i + 1][ny - 1][0] - v.By[i - 1][ny - 1][0]) / (2.0 * dx) - (v.Bx[i][0][0] - v.Bx[i][ny - 2][0]) / (2.0 * dy));
			}
			if (i == 0) {
				v.Ey[i][ny - 1][0] = v.Ey[i][ny - 1][0] - 4.0 * M_PI * dt * Jy[i][ny - 1][0] - c * dt * (v.Bz[i + 1][ny - 1][0] - v.Bz[nx - 1][ny - 1][0]) / (2.0 * dx);
				v.Ez[i][ny - 1][0] = v.Ez[i][ny - 1][0] - 4.0 * M_PI * dt * Jz[i][ny - 1][0] + c * dt * ((v.By[i + 1][ny - 1][0] - v.By[nx - 1][ny - 1][0]) / (2.0 * dx) - (v.Bx[i][0][0] - v.Bx[i][ny - 2][0]) / (2.0 * dy));
			}
			if (i == nx - 1) {
				v.Ey[i][ny - 1][0] = v.Ey[i][ny - 1][0] - 4.0 * M_PI * dt * Jy[i][ny - 1][0] - c * dt * (v.Bz[0][ny - 1][0] - v.Bz[i - 1][ny - 1][0]) / (2.0 * dx);
				v.Ez[i][ny - 1][0] = v.Ez[i][ny - 1][0] - 4.0 * M_PI * dt * Jz[i][ny - 1][0] + c * dt * ((v.By[0][ny - 1][0] - v.By[i - 1][ny - 1][0]) / (2.0 * dx) - (v.Bx[i][0][0] - v.Bx[i][ny - 2][0]) / (2.0 * dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) { // i==0
			v.Ex[0][j][0] = v.Ex[0][j][0] - 4.0 * M_PI * dt * Jx[0][j][0] + c * dt * (v.Bz[0][j + 1][0] - v.Bz[0][j - 1][0]) / (2.0 * dy);
			v.Ey[0][j][0] = v.Ey[0][j][0] - 4.0 * M_PI * dt * Jy[0][j][0] - c * dt * (v.Bz[1][j][0] - v.Bz[nx - 1][j][0]) / (2.0 * dx);
			v.Ez[0][j][0] = v.Ez[0][j][0] - 4.0 * M_PI * dt * Jz[0][j][0] + c * dt * ((v.By[1][j][0] - v.By[nx - 1][j][0]) / (2.0 * dx) - (v.Bx[0][j + 1][0] - v.Bx[0][j - 1][0]) / (2.0 * dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i==nx-1
			v.Ex[nx - 1][j][0] = v.Ex[nx - 1][j][0] - 4.0 * M_PI * dt * Jx[nx - 1][j][0] + c * dt * (v.Bz[nx - 1][j + 1][0] - v.Bz[nx - 1][j - 1][0]) / (2.0 * dy);
			v.Ey[nx - 1][j][0] = v.Ey[nx - 1][j][0] - 4.0 * M_PI * dt * Jy[nx - 1][j][0] - c * dt * (v.Bz[0][j][0] - v.Bz[nx - 2][j][0]) / (2.0 * dx);
			v.Ez[nx - 1][j][0] = v.Ez[nx - 1][j][0] - 4.0 * M_PI * dt * Jz[nx - 1][j][0] + c * dt * ((v.By[0][j][0] - v.By[nx - 2][j][0]) / (2.0 * dx) - (v.Bx[nx - 1][j + 1][0] - v.Bx[nx - 1][j - 1][0]) / (2.0 * dy));
		}
		//обновление граничных значений B
		for (int i = 0; i < nx; i++) { // j==0
			v.Bx[i][0][0] = v.Bx[i][0][0] - c * dt * (v.Ez[i][1][0] - v.Ez[i][ny - 1][0]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.By[i][0][0] = v.By[i][0][0] + c * dt * (v.Ez[i + 1][0][0] - v.Ez[i - 1][0][0]) / (2.0 * dx);
				v.Bz[i][0][0] = v.Bz[i][0][0] - c * dt * ((v.Ey[i + 1][0][0] - v.Ey[i - 1][0][0]) / (2.0 * dx) - (v.Ex[i][1][0] - v.Ex[i][ny - 1][0]) / (2.0 * dy));
			}
			if (i == 0) {
				v.By[i][0][0] = v.By[i][0][0] + c * dt * (v.Ez[i + 1][0][0] - v.Ez[nx - 1][0][0]) / (2.0 * dx);
				v.Bz[i][0][0] = v.Bz[i][0][0] - c * dt * ((v.Ey[i + 1][0][0] - v.Ey[nx - 1][0][0]) / (2.0 * dx) - (v.Ex[i][1][0] - v.Ex[i][ny - 1][0]) / (2.0 * dy));
			}
			if (i == nx - 1) {
				v.By[i][0][0] = v.By[i][0][0] + c * dt * (v.Ez[0][0][0] - v.Ez[i - 1][0][0]) / (2.0 * dx);
				v.Bz[i][0][0] = v.Bz[i][0][0] - c * dt * ((v.Ey[0][0][0] - v.Ey[i - 1][0][0]) / (2.0 * dx) - (v.Ex[i][1][0] - v.Ex[i][ny - 1][0]) / (2.0 * dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j==ny-1
			v.Bx[i][ny - 1][0] = v.Bx[i][ny - 1][0] - c * dt * (v.Ez[i][0][0] - v.Ez[i][ny - 2][0]) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				v.By[i][ny - 1][0] = v.By[i][ny - 1][0] + c * dt * (v.Ez[i + 1][ny - 1][0] - v.Ez[i - 1][ny - 1][0]) / (2.0 * dx);
				v.Bz[i][ny - 1][0] = v.Bz[i][ny - 1][0] - c * dt * ((v.Ey[i + 1][ny - 1][0] - v.Ey[i - 1][ny - 1][0]) / (2.0 * dx) - (v.Ex[i][0][0] - v.Ex[i][ny - 2][0]) / (2.0 * dy));
			}
			if (i == 0) {
				v.By[i][ny - 1][0] = v.By[i][ny - 1][0] + c * dt * (v.Ez[i + 1][ny - 1][0] - v.Ez[nx - 1][ny - 1][0]) / (2.0 * dx);
				v.Bz[i][ny - 1][0] = v.Bz[i][ny - 1][0] - c * dt * ((v.Ey[i + 1][ny - 1][0] - v.Ey[nx - 1][ny - 1][0]) / (2.0 * dx) - (v.Ex[i][0][0] - v.Ex[i][ny - 2][0]) / (2.0 * dy));
			}
			if (i == nx - 1) {
				v.By[i][ny - 1][0] = v.By[i][ny - 1][0] + c * dt * (v.Ez[0][ny - 1][0] - v.Ez[i - 1][ny - 1][0]) / (2.0 * dx);
				v.Bz[i][ny - 1][0] = v.Bz[i][ny - 1][0] - c * dt * ((v.Ey[0][ny - 1][0] - v.Ey[i - 1][ny - 1][0]) / (2.0 * dx) - (v.Ex[i][0][0] - v.Ex[i][ny - 2][0]) / (2.0 * dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) {  // i==0
			v.Bx[0][j][0] = v.Bx[0][j][0] - c * dt * (v.Ez[0][j + 1][0] - v.Ez[0][j - 1][0]) / (2.0 * dy);
			v.By[0][j][0] = v.By[0][j][0] + c * dt * (v.Ez[1][j][0] - v.Ez[nx - 1][j][0]) / (2.0 * dx);
			v.Bz[0][j][0] = v.Bz[0][j][0] - c * dt * ((v.Ey[1][j][0] - v.Ey[nx - 1][j][0]) / (2.0 * dx) - (v.Ex[0][j + 1][0] - v.Ex[0][j - 1][0]) / (2.0 * dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i==nx-1
			v.Bx[nx - 1][j][0] = v.Bx[nx - 1][j][0] - c * dt * (v.Ez[nx - 1][j + 1][0] - v.Ez[nx - 1][j - 1][0]) / (2.0 * dy);
			v.By[nx - 1][j][0] = v.By[nx - 1][j][0] + c * dt * (v.Ez[0][j][0] - v.Ez[nx - 2][j][0]) / (2.0 * dx);
			v.Bz[nx - 1][j][0] = v.Bz[nx - 1][j][0] - c * dt * ((v.Ey[0][j][0] - v.Ey[nx - 2][j][0]) / (2.0 * dx) - (v.Ex[nx - 1][j + 1][0] - v.Ex[nx - 1][j - 1][0]) / (2.0 * dy));
		}

	};

	Field_characteristics Field::FDTD(double dt) { //обновление сеточных значений
	//обновление Е
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++) {
				v.Ex[i][j][0] = v.Ex[i][j][0] - 4.0 * M_PI * dt * Jx[i][j][0] + c * dt * (v.Bz[i][j + 1][0] - v.Bz[i][j - 1][0]) / (2.0 * dy);
				v.Ey[i][j][0] = v.Ey[i][j][0] - 4.0 * M_PI * dt * Jy[i][j][0] - c * dt * (v.Bz[i + 1][j][0] - v.Bz[i - 1][j][0]) / (2.0 * dx);
				v.Ez[i][j][0] = v.Ez[i][j][0] - 4.0 * M_PI * dt * Jz[i][j][0] + c * dt * ((v.By[i + 1][j][0] - v.By[i - 1][j][0]) / (2.0 * dx) - (v.Bx[i][j + 1][0] - v.Bx[i][j - 1][0]) / (2.0 * dy));
			}
		(*this).border_update(dt);
	//обновление В
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++){
				v.Bx[i][j][0] = v.Bx[i][j][0] - c * dt * (v.Ez[i][j + 1][0] - v.Ez[i][j - 1][0]) / (2.0 * dy);
				v.By[i][j][0] = v.By[i][j][0] + c * dt * (v.Ez[i + 1][j][0] - v.Ez[i - 1][j][0]) / (2.0 * dx);
				v.Bz[i][j][0] = v.Bz[i][j][0] - c * dt * ((v.Ey[i + 1][j][0] - v.Ey[i - 1][j][0]) / (2.0 * dx) - (v.Ex[i][j + 1][0] - v.Ex[i][j - 1][0]) / (2.0 * dy));
			}
		return (*this).v;
	}

	Field_characteristics Field::FDTD_with_shift(double dtB, double dtE) {

		for (int i = 0; i < nx -1; i++)
			for (int j = 0; j < ny - 1; j++)
				for (int k = 0; k < nz - 1; k++) {
					v.Bx[i][j][k] = v.Bx[i][j][k] + c * dtB * ((v.Ey[i][j][k + 1] - v.Ey[i][j][k]) / dz - (v.Ez[i][j + 1][k] - v.Ez[i][j][k]) / dy);
					v.By[i][j][k] = v.By[i][j][k] + c * dtB * ((v.Ez[i + 1][j][k] - v.Ez[i][j][k]) / dx - (v.Ex[i][j][k + 1] - v.Ex[i][j][k]) / dz);
					v.Bz[i][j][k] = v.Bz[i][j][k] + c * dtB * ((v.Ex[i][j + 1][k] - v.Ex[i][j][k]) / dy - (v.Ey[i + 1][j][k] - v.Ey[i][j][k]) / dx);
				}
		(*this).border_update_with_shift(dtB,dtE);
		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny; j++)
				for (int k = 1; k < nz; k++) {
					v.Ex[i][j][k] = v.Ex[i][j][k] - 4.0 * M_PI * dtE * Jx[i][j][k] + c * dtE * ((v.Bz[i][j][k] - v.Bz[i][j - 1][k]) / dy - (v.By[i][j][k] - v.By[i][j][k - 1]) / dz);
					v.Ey[i][j][k] = v.Ey[i][j][k] - 4.0 * M_PI * dtE * Jy[i][j][k] + c * dtE * ((v.Bx[i][j][k] - v.Bx[i][j][k - 1]) / dz - (v.Bz[i][j][k] - v.Bz[i - 1][j][k]) / dx);
					v.Ez[i][j][k] = v.Ez[i][j][k] - 4.0 * M_PI * dtE * Jz[i][j][k] + c * dtE * ((v.By[i][j][k] - v.By[i - 1][j][k]) / dx - (v.Bx[i][j][k] - v.Bx[i][j - 1][k]) / dy);
				}

		return (*this).v;
	}

	void Field::border_update_with_shift(double dtB, double dtE) {
		
		//границы В 
		for (int i = 0; i < nx; i++)  // j==ny-1
			for (int k = 0; k < nz; k++) {
				if ((k == nz - 1) && (i != nx - 1)) {
					v.Bx[i][ny - 1][k] = v.Bx[i][ny - 1][k] + c * dtB * ((v.Ey[i][ny - 1][0] - v.Ey[i][ny - 1][k]) / dz - (v.Ez[i][0][k] - v.Ez[i][ny - 1][k]) / dy);
					v.By[i][ny - 1][k] = v.By[i][ny - 1][k] + c * dtB * ((v.Ez[i + 1][ny - 1][k] - v.Ez[i][ny - 1][k]) / dx - (v.Ex[i][ny - 1][0] - v.Ex[i][ny - 1][k]) / dz);
					v.Bz[i][ny - 1][k] = v.Bz[i][ny - 1][k] + c * dtB * ((v.Ex[i][0][k] - v.Ex[i][ny - 1][k]) / dy - (v.Ey[i + 1][ny - 1][k] - v.Ey[i][ny - 1][k]) / dx);
				}
				if ((k != nz - 1) && (i == nx - 1)) {
					v.Bx[i][ny - 1][k] = v.Bx[i][ny - 1][k] + c * dtB * ((v.Ey[i][ny - 1][k + 1] - v.Ey[i][ny - 1][k]) / dz - (v.Ez[i][0][k] - v.Ez[i][ny - 1][k]) / dy);
					v.By[i][ny - 1][k] = v.By[i][ny - 1][k] + c * dtB * ((v.Ez[0][ny - 1][k] - v.Ez[i][ny - 1][k]) / dx - (v.Ex[i][ny - 1][k + 1] - v.Ex[i][ny - 1][k]) / dz);
					v.Bz[i][ny - 1][k] = v.Bz[i][ny - 1][k] + c * dtB * ((v.Ex[i][0][k] - v.Ex[i][ny - 1][k]) / dy - (v.Ey[0][ny - 1][k] - v.Ey[i][ny - 1][k]) / dx);
				}
				if ((k == nz - 1) && (i == nx - 1)) {
					v.Bx[i][ny - 1][k] = v.Bx[i][ny - 1][k] + c * dtB * ((v.Ey[i][ny - 1][0] - v.Ey[i][ny - 1][k]) / dz - (v.Ez[i][0][k] - v.Ez[i][ny - 1][k]) / dy);
					v.By[i][ny - 1][k] = v.By[i][ny - 1][k] + c * dtB * ((v.Ez[0][ny - 1][k] - v.Ez[i][ny - 1][k]) / dx - (v.Ex[i][ny - 1][0] - v.Ex[i][ny - 1][k]) / dz);
					v.Bz[i][ny - 1][k] = v.Bz[i][ny - 1][k] + c * dtB * ((v.Ex[i][0][k] - v.Ex[i][ny - 1][k]) / dy - (v.Ey[0][ny - 1][k] - v.Ey[i][ny - 1][k]) / dx);
				}
			}
		
		for (int j = 0; j < ny-1; j++) // i==nx-1
			for (int k = 0; k < nz; k++) {
				if (k == nz-1) {
					v.Bx[nx - 1][j][k] = v.Bx[nx - 1][j][k] + c * dtB * ((v.Ey[nx - 1][j][0] - v.Ey[nx - 1][j][k]) / dz - (v.Ez[nx - 1][j + 1][k] - v.Ez[nx - 1][j][k]) / dy);
					v.By[nx - 1][j][k] = v.By[nx - 1][j][k] + c * dtB * ((v.Ez[0][j][k] - v.Ez[nx - 1][j][k]) / dx - (v.Ex[nx - 1][j][0] - v.Ex[nx - 1][j][k]) / dz);
					v.Bz[nx - 1][j][k] = v.Bz[nx - 1][j][k] + c * dtB * ((v.Ex[nx - 1][j + 1][k] - v.Ex[nx - 1][j][k]) / dy - (v.Ey[0][j][k] - v.Ey[nx - 1][j][k]) / dx);
				}
				else {
					v.Bx[nx - 1][j][k] = v.Bx[nx - 1][j][k] + c * dtB * ((v.Ey[nx - 1][j][k + 1] - v.Ey[nx - 1][j][k]) / dz - (v.Ez[nx - 1][j + 1][k] - v.Ez[nx - 1][j][k]) / dy);
					v.By[nx - 1][j][k] = v.By[nx - 1][j][k] + c * dtB * ((v.Ez[0][j][k] - v.Ez[nx - 1][j][k]) / dx - (v.Ex[nx - 1][j][k + 1] - v.Ex[nx - 1][j][k]) / dz);
					v.Bz[nx - 1][j][k] = v.Bz[nx - 1][j][k] + c * dtB * ((v.Ex[nx - 1][j + 1][k] - v.Ex[nx - 1][j][k]) / dy - (v.Ey[0][j][k] - v.Ey[nx - 1][j][k]) / dx);
				}
			}

		for (int i = 0; i < nx-1; i++) // k==nz-1
			for (int j = 0; j < ny-1; j++) {
				v.Bx[i][j][nz - 1] = v.Bx[i][j][nz - 1] + c * dtB * ((v.Ey[i][j][0] - v.Ey[i][j][nz - 1]) / dz - (v.Ez[i][j + 1][nz - 1] - v.Ez[i][j][nz - 1]) / dy);
				v.By[i][j][nz - 1] = v.By[i][j][nz - 1] + c * dtB * ((v.Ez[i + 1][j][nz - 1] - v.Ez[i][j][nz - 1]) / dx - (v.Ex[i][j][0] - v.Ex[i][j][nz - 1]) / dz);
				v.Bz[i][j][nz - 1] = v.Bz[i][j][nz - 1] + c * dtB * ((v.Ex[i][j + 1][nz - 1] - v.Ex[i][j][nz - 1]) / dy - (v.Ey[i + 1][j][nz - 1] - v.Ey[i][j][nz - 1]) / dx);
			}
		
		
		//границы Е
		for (int i = 0; i < nx; i++) // j==0
			for (int k = 0; k < nz; k++) {
				if ((k == 0) && (i != 0)) {
					v.Ex[i][0][k] = v.Ex[i][0][k] - 4.0 * M_PI * dtE * Jx[i][0][k] + c * dtE * ((v.Bz[i][0][k] - v.Bz[i][ny - 1][k]) / dy - (v.By[i][0][k] - v.By[i][0][nz - 1]) / dz);
					v.Ey[i][0][k] = v.Ey[i][0][k] - 4.0 * M_PI * dtE * Jy[i][0][k] + c * dtE * ((v.Bx[i][0][k] - v.Bx[i][0][nz - 1]) / dz - (v.Bz[i][0][k] - v.Bz[i - 1][0][k]) / dx);
					v.Ez[i][0][k] = v.Ez[i][0][k] - 4.0 * M_PI * dtE * Jz[i][0][k] + c * dtE * ((v.By[i][0][k] - v.By[i - 1][0][k]) / dx - (v.Bx[i][0][k] - v.Bx[i][ny - 1][k]) / dy);
				}
				if ((k == 0) && (i == 0)) {
					v.Ex[i][0][k] = v.Ex[i][0][k] - 4.0 * M_PI * dtE * Jx[i][0][k] + c * dtE * ((v.Bz[i][0][k] - v.Bz[i][ny - 1][k]) / dy - (v.By[i][0][k] - v.By[i][0][nz - 1]) / dz);
					v.Ey[i][0][k] = v.Ey[i][0][k] - 4.0 * M_PI * dtE * Jy[i][0][k] + c * dtE * ((v.Bx[i][0][k] - v.Bx[i][0][nz - 1]) / dz - (v.Bz[i][0][k] - v.Bz[nx - 1][0][k]) / dx);
					v.Ez[i][0][k] = v.Ez[i][0][k] - 4.0 * M_PI * dtE * Jz[i][0][k] + c * dtE * ((v.By[i][0][k] - v.By[nx - 1][0][k]) / dx - (v.Bx[i][0][k] - v.Bx[i][ny - 1][k]) / dy);
				}
				if ((k != 0) && (i == 0)) {
					v.Ex[i][0][k] = v.Ex[i][0][k] - 4.0 * M_PI * dtE * Jx[i][0][k] + c * dtE * ((v.Bz[i][0][k] - v.Bz[i][ny - 1][k]) / dy - (v.By[i][0][k] - v.By[i][0][k - 1]) / dz);
					v.Ey[i][0][k] = v.Ey[i][0][k] - 4.0 * M_PI * dtE * Jy[i][0][k] + c * dtE * ((v.Bx[i][0][k] - v.Bx[i][0][k - 1]) / dz - (v.Bz[i][0][k] - v.Bz[nx - 1][0][k]) / dx);
					v.Ez[i][0][k] = v.Ez[i][0][k] - 4.0 * M_PI * dtE * Jz[i][0][k] + c * dtE * ((v.By[i][0][k] - v.By[nx - 1][0][k]) / dx - (v.Bx[i][0][k] - v.Bx[i][ny - 1][k]) / dy);
				}
			}

		for (int j = 1; j < ny; j++) // i==0
			for (int k = 0; k < nz; k++) {
				if (k == 0) {
					v.Ex[0][j][k] = v.Ex[0][j][k] - 4.0 * M_PI * dtE * Jx[0][j][k] + c * dtE * ((v.Bz[0][j][k] - v.Bz[0][j - 1][k]) / dy - (v.By[0][j][k] - v.By[0][j][nz - 1]) / dz);
					v.Ey[0][j][k] = v.Ey[0][j][k] - 4.0 * M_PI * dtE * Jy[0][j][k] + c * dtE * ((v.Bx[0][j][k] - v.Bx[0][j][nz - 1]) / dz - (v.Bz[0][j][k] - v.Bz[nx - 1][j][k]) / dx);
					v.Ez[0][j][k] = v.Ez[0][j][k] - 4.0 * M_PI * dtE * Jz[0][j][k] + c * dtE * ((v.By[0][j][k] - v.By[nx - 1][j][k]) / dx - (v.Bx[0][j][k] - v.Bx[0][j - 1][k]) / dy);
				}
				else {
					v.Ex[0][j][k] = v.Ex[0][j][k] - 4.0 * M_PI * dtE * Jx[0][j][k] + c * dtE * ((v.Bz[0][j][k] - v.Bz[0][j - 1][k]) / dy - (v.By[0][j][k] - v.By[0][j][k - 1]) / dz);
					v.Ey[0][j][k] = v.Ey[0][j][k] - 4.0 * M_PI * dtE * Jy[0][j][k] + c * dtE * ((v.Bx[0][j][k] - v.Bx[0][j][k - 1]) / dz - (v.Bz[0][j][k] - v.Bz[nx - 1][j][k]) / dx);
					v.Ez[0][j][k] = v.Ez[0][j][k] - 4.0 * M_PI * dtE * Jz[0][j][k] + c * dtE * ((v.By[0][j][k] - v.By[nx - 1][j][k]) / dx - (v.Bx[0][j][k] - v.Bx[0][j - 1][k]) / dy);
				}

			}

		for (int i = 1; i < nx; i++) // k==0
			for (int j = 1; j < ny; j++) {
				v.Ex[i][j][0] = v.Ex[i][j][0] - 4.0 * M_PI * dtE * Jx[i][j][0] + c * dtE * ((v.Bz[i][j][0] - v.Bz[i][j - 1][0]) / dy - (v.By[i][j][0] - v.By[i][j][nz - 1]) / dz);
				v.Ey[i][j][0] = v.Ey[i][j][0] - 4.0 * M_PI * dtE * Jy[i][j][0] + c * dtE * ((v.Bx[i][j][0] - v.Bx[i][j][nz - 1]) / dz - (v.Bz[i][j][0] - v.Bz[i - 1][j][0]) / dx);
				v.Ez[i][j][0] = v.Ez[i][j][0] - 4.0 * M_PI * dtE * Jz[i][j][0] + c * dtE * ((v.By[i][j][0] - v.By[i - 1][j][0]) / dx - (v.Bx[i][j][0] - v.Bx[i][j - 1][0]) / dy);
			}
	}
};
