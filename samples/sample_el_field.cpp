#include "el_field.h"
#include <fstream>

int main() {
	double real_value;
	Field f;
	double dt = 0.001, t=1;
	Field_characteristics res;

	/*std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	std::ofstream out2("D:\\final_el_field\\Electromagnetic-field\\result_real.txt", std::ios_base::out | std::ios_base::trunc);

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
	}
	


	for (int i = 0; i < 3; i++)
		for (int j = 0; j < f.nx; j++)
			for (int k = 0; k < f.ny; k++) {
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


	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2.0 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = sin(2.0 * M_PI * (dy * j) / (by - ay));
		}

	for (double i = 0; i < t; i+=dt)
		res = f.FDTD(dt);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < f.nx; j++)
			for (int k = 0; k < f.ny; k++) {
				if (i != 0)
					real_value = 0;
				else
					real_value = sin(2.0 * M_PI * (k * dy - c * t) / (by - ay));
				std::cout << "E[" << i << "][" << j << "][" << k << "]=" << res.E[i][j][k] << "  " << "E_real=" << real_value << "    ";
				if (i != 2)
					real_value = 0;
				else
					real_value = sin(2.0 * M_PI * (k * dy - c * t) / (by - ay));
				std::cout << "B[" << i << "][" << j << "][" << k << "]=" << res.B[i][j][k] << "  " << "B_real=" << real_value << std::endl;
			}
}