#include "el_field.h"
#include <fstream>
#define N 16

int main() {
	/*int n = 120, T = 16, m = 16;
	double d = c, dx = c, dy = c, dz = c, dt = 0.1, t1 = 70, t;
	double ax = -n * c / 2, ay = -n * c / 2, bx = n * c / 2, by = n * c / 2, az = -n * c / 2, bz = n * c / 2;
	double x, z;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz, (bx - ax) / dx + 1, 1, (bz - az) / dz + 1);
	int i = 1;
	int k = int(T / dt);
	int y = 0;

	for (double t = dt; t < t1; t += dt) {
		if (i <= k) {
			t = dt * i - dt / 2;
			for (int l = 0; l < f.nx; l++)
				for (int j = 0; j < f.nz; j++) {
					x = ax + l * dx;
					z = az + j * dz;
					if ((x >= -m * c / 4.0) && (x <= m * c / 4.0) && (z >= -m * c / 4.0) && (z <= m * c / 4.0))
						f.Jy(l, 0, j) = sin(2 * M_PI * t / T) * pow(cos(2 * M_PI * x / (m * c)), 2) * pow(cos(2 * M_PI * z / (m * c)), 2);
				}
		}
		else {
			f.Jy(f.nx / 2, 0, f.nz / 2) = 0.0;
		}
		f.PSTD(dt, dt);
		i++;
	}

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	for (int l = 0; l < f.nx; l++)
		for (int k = 0; k < f.nz; k++) {
			out << l << ';' << k << ';' << f.Ey(l, 0, k) << ';' << std::endl;

		}
	return 0;*/

	double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 1E+9;
	double real_value;
	double dt = 0.001, t = 0.37, dx = 1E+9, dy = 1E+9, dz = 1E+9;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz, (bx - ax) / dx + 1, (by - ay) / dy + 1, 1);

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	std::ofstream out2("D:\\final_el_field\\Electromagnetic-field\\result_real.txt", std::ios_base::out | std::ios_base::trunc);

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i, j, 0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i, j, 0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}

	f.PSATD(dt, t);

	for (int j = 0; j < f.ny; j++) {
		out << f.Ex(1, j, 0) << ';' << j * dy + ay << std::endl;
		out2 << sin(2.0 * M_PI * (j * dy - c * t) / (by - ay)) << ';' << j * dy + ay << std::endl;
	}
	
}
