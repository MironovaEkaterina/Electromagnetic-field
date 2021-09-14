#include "el_field.h"

int main() {
	double real_value;
	Field f;
	double dt = 1E-10;
	Field_characteristics res;

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}

	//for (int n = 0; n < 100; n++)
	res = f.FDTD(dt);


	for (int i = 0; i < 3; i++)
		for (int j = 0; j < f.nx; j++)
			for (int k = 0; k < f.ny; k++) {
				if (i != 1)
					real_value = 0;
				else
					real_value = sin(2 * M_PI * (j * dx - c * dt) / (bx - ax));
				std::cout << "E[" << i << "][" << j << "][" << k << "]=" << res.E[i][j][k] << "  " << "E_real=" << real_value << std::endl;
				if (i != 2)
					real_value = 0;
				else
					real_value = sin(2 * M_PI * (j * dx - c * dt) / (bx - ax));
				std::cout << "B[" << i << "][" << j << "][" << k << "]=" << res.B[i][j][k] << "  " << "B_real=" << real_value << std::endl;
			}
}