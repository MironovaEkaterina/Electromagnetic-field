#include "el_field.h"
#include <fstream>
#include <complex.h>
#include <complex>
#include <fftw3.h>
#define N 16

int main() {

	double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 1E+9;
	double real_value;
	double dt = 0.001, t = 90, dx = 1E+9, dy = 1E+9, dz = 1E+9;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	std::ofstream out2("D:\\final_el_field\\Electromagnetic-field\\result_real.txt", std::ios_base::out | std::ios_base::trunc);
	//std::ofstream out3("D:\\final_el_field\\Electromagnetic-field\\begin.txt", std::ios_base::out | std::ios_base::trunc);

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) 
			for (int k = 0; k < f.nz; k++) {
			f.v.Ey[i][j][k] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][k] = sin(2.0 * M_PI * (dx * i) / (bx - ax));
		}

	double* Ex_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	double* Ey_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	double* Ez_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	double* Bx_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	double* By_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	double* Bz_for_dft = (double*)malloc(sizeof(double) * f.nx * f.ny * f.nz);
	
	for (int i = 0; i < f.nx * f.ny * f.nz; i++) {
		Ex_for_dft[i] = 0.0;
		Ez_for_dft[i] = 0.0;
		Bx_for_dft[i] = 0.0;
		By_for_dft[i] = 0.0;
	}
	
	int q = 0;
	for (int i=0; i<f.nx; i++)
		for (int j = 0; j < f.ny; j++)
			for (int k = 0; k < f.nz; k++) {
			Ey_for_dft[q] = f.v.Ey[i][j][k];
			Bz_for_dft[q] = f.v.Bz[i][j][k];
			q++;
		}

	fftw_complex* Ex_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	fftw_complex* Ey_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	fftw_complex* Ez_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	fftw_complex* Bx_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	fftw_complex* By_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);
	fftw_complex* Bz_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * f.nx * f.ny * f.nz);

	fftw_plan Ex, Ey, Ez, Bx, By, Bz;
	fftw_plan Ex_, Ey_, Ez_, Bx_, By_, Bz_;

	Ex = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Ex_for_dft, Ex_out, FFTW_ESTIMATE);
	Ey = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Ey_for_dft, Ey_out, FFTW_ESTIMATE);
	Ez = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Ez_for_dft, Ez_out, FFTW_ESTIMATE);
	Bx = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Bx_for_dft, Bx_out, FFTW_ESTIMATE);
	By = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, By_for_dft, By_out, FFTW_ESTIMATE);
	Bz = fftw_plan_dft_r2c_3d(f.nx, f.ny, f.nz, Bz_for_dft, Bz_out, FFTW_ESTIMATE);

	fftw_execute(Ex);
	fftw_execute(Ey);
	fftw_execute(Ez);
	fftw_execute(Bx);
	fftw_execute(By);
	fftw_execute(Bz);

	for (double t1 = 0; t1 < t; t1 += dt) {
		q = 0;
		double wx, wy, wz;
		for (int i_ = 0; i_ < f.nx; i_++)
			for (int j_ = 0; j_ < f.ny; j_++) 
			    for (int k_ = 0; k_ < f.nz; k_++){
					if (i_ <= f.nx / 2)
						wx = 2 * M_PI * i_ / (bx - ax);
					else
						wx = 2 * M_PI * (i_ - f.nx) / (bx - ax);

					if (j_ <= f.ny / 2)
						wy = 2 * M_PI * j_ / (by - ay);
					else
						wy = 2 * M_PI * (j_ - f.ny) / (by - ay);

					if (k_ <= f.nz / 2)
						wz = 2 * M_PI * k_ / (bz - az);
					else
						wz = 2 * M_PI * (k_ - f.nz) / (bz - az);

				Ex_out[q][0] += c * dt * (-wy * Bz_out[q][1] + wz * By_out[q][1]);
				Ex_out[q][1] += c * dt * (wy * Bz_out[q][0] - wz * By_out[q][0]);

				Ey_out[q][0] += c * dt * (wx * Bz_out[q][1] - wz * Bx_out[q][1]);
				Ey_out[q][1] += c * dt * (-wx * Bz_out[q][0] + wz * Bx_out[q][0]);

				Ez_out[q][0] += c * dt * (-wx * By_out[q][1] + wy * Bx_out[q][1]);
				Ez_out[q][1] += c * dt * (wx * By_out[q][0] - wy * Bx_out[q][0]);


				Bx_out[q][0] += c * dt * (wy * Ez_out[q][1] - wz * Ey_out[q][1]);
				Bx_out[q][1] += c * dt * (-wy * Ez_out[q][0] + wz * Ey_out[q][0]);

				By_out[q][0] += c * dt * (-wx * Ez_out[q][1] + wz * Ex_out[q][1]);
				By_out[q][1] += c * dt * (wx * Ez_out[q][0] - wz * Ex_out[q][0]);

				Bz_out[q][0] += c * dt * (wx * Ey_out[q][1] - wy * Ex_out[q][1]);
				Bz_out[q][1] += c * dt * (-wx * Ey_out[q][0] + wy * Ex_out[q][0]);
				q++;
			}

	}

	Ex_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Ex_out, Ex_for_dft, FFTW_ESTIMATE);
	Ey_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Ey_out, Ey_for_dft, FFTW_ESTIMATE);
	Ez_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Ez_out, Ez_for_dft, FFTW_ESTIMATE);
	Bx_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Bx_out, Bx_for_dft, FFTW_ESTIMATE);
	By_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, By_out, By_for_dft, FFTW_ESTIMATE);
	Bz_ = fftw_plan_dft_c2r_3d(f.nx, f.ny, f.nz, Bz_out, Bz_for_dft, FFTW_ESTIMATE);

	fftw_execute(Ex_);
	fftw_execute(Ey_);
	fftw_execute(Ez_);
	fftw_execute(Bx_);
	fftw_execute(By_);
	fftw_execute(Bz_);

	for (int i = 0; i < f.nx * f.ny * f.nz; i++) {
		Ex_for_dft[i] /= f.nx * f.ny * f.nz;
		Ey_for_dft[i] /= f.nx * f.ny * f.nz;
		Ez_for_dft[i] /= f.nx * f.ny * f.nz;
		Bx_for_dft[i] /= f.nx * f.ny * f.nz;
		By_for_dft[i] /= f.nx * f.ny * f.nz;
		Bz_for_dft[i] /= f.nx * f.ny * f.nz;
	}

	q = 0;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++)
			for (int k = 0; k < f.nz; k++) {
			f.v.Ex[i][j][k] = Ex_for_dft[q];
			f.v.Ey[i][j][k] = Ey_for_dft[q];
			f.v.Ez[i][j][k] = Ez_for_dft[q];
			f.v.Bx[i][j][k] = Bx_for_dft[q];
			f.v.By[i][j][k] = By_for_dft[q];
			f.v.Bz[i][j][k] = Bz_for_dft[q];
			q++;
		}

	for (int j = 0; j < f.nx; j++) {
		out << f.v.Ey[j][1][0] << ';' << j * dx + ax << std::endl;
		out2 << sin(2.0 * M_PI * (j * dx - c * t) / (bx - ax)) << ';' << j * dx + ax << std::endl;
		//out3 << sin(2.0 * M_PI * (j * dx) / (bx - ax)) << ';' << j * dx + ax << std::endl;
	}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ex[" << i << "][" << j << "][0]=" << f.v.Ex[i][j][0] << "  " << "E_real=" << 0 << "     ";
			std::cout << "Bx[" << i << "][" << j << "][0]=" << f.v.Bx[i][j][0] << "  " << "B_real=" << 0 << std::endl;
		}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ey[" << i << "][" << j << "][0]=" << f.v.Ey[i][j][0] << "  " << "E_real=" << real_value << "     ";
			std::cout << "By[" << i << "][" << j << "][0]=" << f.v.By[i][j][0] << "  " << "B_real=" << 0 << std::endl;
		}

	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			real_value = sin(2.0 * M_PI * (i * dx - c * t) / (bx - ax));
			std::cout << "Ez[" << i << "][" << j << "][0]=" << f.v.Ez[i][j][0] << "  " << "E_real=" << 0 << "     ";
			std::cout << "Bz[" << i << "][" << j << "][0]=" << f.v.Bz[i][j][0] << "  " << "B_real=" << real_value << std::endl;
		}

	fftw_destroy_plan(Ex);
	fftw_destroy_plan(Ey);
	fftw_destroy_plan(Ez);
	fftw_destroy_plan(Bx);
	fftw_destroy_plan(By);
	fftw_destroy_plan(Bz);

	fftw_destroy_plan(Ex_);
	fftw_destroy_plan(Ey_);
	fftw_destroy_plan(Ez_);
	fftw_destroy_plan(Bx_);
	fftw_destroy_plan(By_);
	fftw_destroy_plan(Bz_);

	fftw_free(Ex_out);
	fftw_free(Ey_out);
	fftw_free(Ez_out);
	fftw_free(Bx_out);
	fftw_free(By_out);
	fftw_free(Bz_out);

	free(Ex_for_dft);
	free(Ey_for_dft);
	free(Ez_for_dft);
	free(Bx_for_dft);
	free(By_for_dft);
	free(Bz_for_dft);
}




	/*//сферическая волна
	int n = 120, T = 16, m = 16;
	double d = c, dx = c, dy = c, dz = c, dt = 0.1, t1 = 60, t;
	double ax = -n * c / 2, ay = -n * c / 2, bx = n * c / 2, by = n * c / 2, az= -n * c / 2, bz= -n * c / 2;
	double x, y;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	int i = 0;
	int k = int(T / dt);

	//res = f.FDTD_with_shift(dt / 2, dt);

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
		//res = f.FDTD_with_shift(dt,dt);
		res = f.FDTD(dt);
		i++;
	}

	std::ofstream out("D:\\final_el_field\\Electromagnetic-field\\result.txt", std::ios_base::out | std::ios_base::trunc);
	for (int l = 0; l < f.nx; l++)
		for (int k = 0; k < f.ny; k++) {
			out << l << ';' << k << ';' << res.Ez[l][k][0] << ';' << std::endl;

		}
	return 0;*/

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
