#include "el_field.h"
#include <fstream>

int main() {
	//сферическая волна
	int n = 120, T = 16, m = 16;
	double d = c, dx = c, dy = c, dz = c, dt = 0.1, t1 = 40, t;
	double ax = -n * c / 2, ay = -n * c / 2, bx = n * c / 2, by = n * c / 2;
	double x, y;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	int i = 0;
	int k = int(T / dt);

	for (double t = 0; t < t1; t += dt) {
		if (i <= k) {
			t = dt * i;
			for (int l = 0; l < f.nx; l++)
				for (int j = 0; j < f.ny; j++) {
					x = ax + l * dx;
					y = ay + j * dy;
					if ((x >= -m * c / 4.0) && (x <= m * c / 4.0) && (y >= -m * c / 4.0) && (y <= m * c / 4.0))
						f.J[2][l][j] = sin(2 * M_PI * t / T) * pow(cos(2 * M_PI * x / (m * c)), 2) * pow(cos(2 * M_PI * y / (m * c)), 2);
				}
		}
		else {
			f.J[2][f.nx / 2][f.ny / 2] = 0.0;
		}
		res = f.FDTD(dt);
		i++;
	}

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	for (int l = 0; l < f.nx; l++)
		for (int k = 0; k < f.ny; k++) {
			out << l << ';' << k << ';' << res.E[2][l][k] << ';' << std::endl;

		}

	/*double ax = 0, bx = 5E+10, ay = 0, by = 5E+10;
	double real_value;
	double dt = 0.01, t = 1, dx =1E+9, dy = 1E+10;
	Field f(dx, dy,ax,ay,bx,by);
	Field_characteristics res;

	//std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	//std::ofstream out2("D:\\final_el_field\\Electromagnetic-field\\result_real.txt", std::ios_base::out | std::ios_base::trunc);
	//std::ofstream out3("D:\\final_el_field\\Electromagnetic-field\\begin.txt", std::ios_base::out | std::ios_base::trunc);

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
		}

	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);

	for (int j = 0; j < f.nx; j++) {
		out << res.E[1][j][3] << ';' << j * dx + ax << std::endl;
		out2 << sin(2.0 * M_PI * (j * dx - c * t) / (bx - ax)) << ';' << j * dx + ax << std::endl;
		//out3 << sin(2.0 * M_PI * (j * dx) / (bx - ax)) << ';' << j * dx + ax << std::endl;
	}

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < f.nx; j++)
			for (int k = 0; k < 1; k++) {
				if (i != 1)
					real_value = 0;
				else
					real_value = sin(2.0 * M_PI * (j * dx - c * t) / (bx - ax));
				std::cout << "E[" << i << "][" << j << "][" << k << "]=" << res.E[i][j][k] << "  " << "E_real=" << real_value << "    ";
				if (i != 2)
					real_value = 0;
				else
					real_value = sin(2.0 * M_PI * (j * dx - c * t) / (bx - ax));
				std::cout << "B[" << i << "][" << j << "][" << k << "]=" << res.B[i][j][k] << "  " << "B_real=" << real_value << std::endl;
			}*/

}