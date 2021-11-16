#include "el_field.h"
#include <fstream>

int main() {
	//сферическая волна
	int n = 120, T = 16, m = 16;
	double d = c, dx = c, dy = c, dz = c, dt = 0.1, t1 = 50, t;
	double ax = -n * c / 2, ay = -n * c / 2, bx = n * c / 2, by = n * c / 2, az= -n * c / 2, bz= -n * c / 2;
	double x, y;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	int i = 0;
	int k = int(T / dt);

	res = f.FDTD_with_shift(dt / 2, dt);

	for (double t = dt; t < t1; t += dt) {
		if (i <= k) {
			t = dt * i;
			for (int l = 0; l < f.nx; l++)
				for (int j = 0; j < f.ny; j++) {
					x = ax + l * dx;
					y = ay + j * dy;
					if ((x >= -m * c / 4.0) && (x <= m * c / 4.0) && (y >= -m * c / 4.0) && (y <= m * c / 4.0))
						f.Jz[l][j][0] = sin(2 * M_PI * t / T) * pow(cos(2 * M_PI * x / (m * c)), 2) * pow(cos(2 * M_PI * y / (m * c)), 2);
				}
		}
		else {
			f.Jz[f.nx / 2][f.ny / 2][0] = 0.0;
		}
		res = f.FDTD_with_shift(dt,dt);
		i++;
	}

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	for (int l = 0; l < f.nx; l++)
		for (int k = 0; k < f.ny; k++) {
			out << l << ';' << k << ';' << res.Ez[l][k][0] << ';' << std::endl;

		}

	/*double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 0;
	double real_value;
	double dt = 0.01, t = 1, dx = 1E+9, dy = 1E+10, dz = 1;
	Field f(dx,dy,dz,ax,ay,az,bx,by,bz);
	Field_characteristics res;

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	//std::ofstream out2("D:\\final_el_field\\Electromagnetic-field\\result_real.txt", std::ios_base::out | std::ios_base::trunc);
	std::ofstream out3("D:\\final_el_field\\Electromagnetic-field\\begin.txt", std::ios_base::out | std::ios_base::trunc);

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
		}

	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);

	for (int j = 0; j < f.nx; j++) {
		out << res.Ey[j][1][0] << ';' << j * dx + ax << std::endl;
		//out2 << sin(2.0 * M_PI * (j * dx - c * t) / (bx - ax)) << ';' << j * dx + ax << std::endl;
		out3 << sin(2.0 * M_PI * (j * dx) / (bx - ax)) << ';' << j * dx + ax << std::endl;
	}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++){
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ex[" << i << "][" << j << "][0]=" << res.Ex[i][j][0] << "  " << "E_real=" << 0 <<"     ";
			std::cout << "Bx[" << i << "][" << j << "][0]=" << res.Bx[i][j][0] << "  " << "B_real=" << 0 << std::endl;
		}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ey[" << i << "][" << j << "][0]=" << res.Ey[i][j][0] << "  " << "E_real=" << real_value << "     ";
			std::cout << "By[" << i << "][" << j << "][0]=" << res.By[i][j][0] << "  " << "B_real=" << 0 << std::endl;
		}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ez[" << i << "][" << j << "][0]=" << res.Ez[i][j][0] << "  " << "E_real=" << 0 << "     ";
			std::cout << "Bz[" << i << "][" << j << "][0]=" << res.Bz[i][j][0] << "  " << "B_real=" << real_value << std::endl;
		}*/

	/*double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 0;
	double real_value;
	double dt = 0.01, t = 5, dx = 1E+9, dy = 1E+10, dz = 1E+10;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2.0 * M_PI * (dx * i+dx/2) / (bx - ax));
		}

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	std::ofstream out2("D:\\final_el_field\\Electromagnetic-field\\result_real.txt", std::ios_base::out | std::ios_base::trunc);
	//std::ofstream out3("D:\\final_el_field\\Electromagnetic-field\\begin.txt", std::ios_base::out | std::ios_base::trunc);
	res = f.FDTD_with_shift(dt/2, dt);

	for (double i = dt; i < t; i += dt)
		res = f.FDTD_with_shift(dt,dt);

	for (int j = 0; j < f.nx; j++) {
		out << res.Ey[j][1][0] << ';' << j * dx + ax << std::endl;
		out2 << sin(2.0 * M_PI * (j * dx - c * t) / (bx - ax)) << ';' << j * dx + ax << std::endl;
		//out3 << sin(2.0 * M_PI * (j * dx) / (bx - ax)) << ';' << j * dx + ax << std::endl;
	}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			std::cout << "Ex[" << i << "][" << j << "][0]=" << res.Ex[i][j][0] << "  " << "E_real=" << 0 << "     ";
			std::cout << "Bx[" << i << "][" << j << "][0]=" << res.Bx[i][j][0] << "  " << "B_real=" << 0 << std::endl;
		}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ey[" << i << "][" << j << "][0]=" << res.Ey[i][j][0] << "  " << "E_real=" << real_value << "     ";
			std::cout << "By[" << i << "][" << j << "][0]=" << res.By[i][j][0] << "  " << "B_real=" << 0 << std::endl;
		}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * (t-dt/2)) / (bx - ax));
			std::cout << "Ez[" << i << "][" << j << "][0]=" << res.Ez[i][j][0] << "  " << "E_real=" << 0 << "     ";
			std::cout << "Bz[" << i << "][" << j << "][0]=" << res.Bz[i][j][0] << "  " << "B_real=" << real_value << std::endl;
		}*/

}