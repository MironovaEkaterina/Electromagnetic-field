#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>

const double c = 3E+10;
#define M_PI 3.1415926535897932384626433832795

template <typename T>
class Array {
public:
	std::vector<T> array;
	int nx, ny, nz, size;

	Array() {
		nx = 0;
		ny = 0;
		nz = 0;
		size = 0;
	}

	Array(int _nx, int _ny, int _nz) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		size = nx * ny * nz;
		array.assign(size, 0.0);
	}

	void resize(int _nx, int _ny, int _nz) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		size = nx * ny * nz;
		array.assign(size, 0.0);
	}

	~Array(){
		array.clear();
	}

	T& operator() (int i, int j, int k){
		return array[i*ny*nz+j*nz+k];
	}

	T operator() (int i, int j, int k) const {
		return array[i * ny * nz + j * nz + k];
	}

	T& operator[](int i){   
		return array[i];
	}

	T* data() {
		return array.data();
	}
};

class Field {
public:
	Array<double> Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz;
	double dx, dy, dz;
	int nx, ny, nz, size;
	double ax, bx, ay, by, az, bz;

	Field() {}

	Field(double dx_, double dy_, double dz_, double ax_, double ay_, double az_, double bx_, double by_, double bz_, int nx_, int ny_, int nz_) {
		dx = dx_;
		dy = dy_;
		dz = dz_;
		ax = ax_;
		bx = bx_;
		ay = ay_;
		by = by_;
		az = az_;
		bz = bz_;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		size = nx * ny * nz;
		Ex.resize(nx, ny, nz);
		Bx.resize(nx, ny, nz);
		Jx.resize(nx, ny, nz);
		Ey.resize(nx, ny, nz);
		By.resize(nx, ny, nz);
		Jy.resize(nx, ny, nz);
		Ez.resize(nx, ny, nz);
		Bz.resize(nx, ny, nz);
		Jz.resize(nx, ny, nz);
	};

	void Field::Initialize(double dx_, double dy_, double dz_, double ax_, double ay_, double az_, double bx_, double by_, double bz_, int nx_, int ny_, int nz_) {
		dx = dx_;
		dy = dy_;
		dz = dz_;
		ax = ax_;
		bx = bx_;
		ay = ay_;
		by = by_;
		az = az_;
		bz = bz_;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		size = nx * ny * nz;
		Ex.resize(nx, ny, nz);
		Bx.resize(nx, ny, nz);
		Jx.resize(nx, ny, nz);
		Ey.resize(nx, ny, nz);
		By.resize(nx, ny, nz);
		Jy.resize(nx, ny, nz);
		Ez.resize(nx, ny, nz);
		Bz.resize(nx, ny, nz);
		Jz.resize(nx, ny, nz);
	}

	void Field::PSATD(double dt, double t) {
		int q = 0;
		std::complex<double> i(0.0, 1.0);
		double wx, wy, wz, w, c_, s;

		Array<std::complex<double>> Ex_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ey_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ez_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Bx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> By_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Bz_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jy_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jz_out(nx, ny, nz / 2 + 1);

		fftw_plan plan;

		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Ex.data(), (fftw_complex*)(Ex_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Ey.data(), (fftw_complex*)(Ey_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Ez.data(), (fftw_complex*)(Ez_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Bx.data(), (fftw_complex*)(Bx_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, By.data(), (fftw_complex*)(By_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Bz.data(), (fftw_complex*)(Bz_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Jx.data(), (fftw_complex*)(Jx_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Jy.data(), (fftw_complex*)(Jy_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Jz.data(), (fftw_complex*)(Jz_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);


		for (double t1 = 0; t1 < t; t1 += dt) {
			q = 0;
			for (int i_ = 0; i_ < nx; i_++)
				for (int j_ = 0; j_ < ny; j_++)
					for (int k_ = 0; k_ < nz / 2 + 1; k_++) {
						if (i_ <= nx / 2)
							wx = 2 * M_PI * i_ / (bx - ax);
						else
							wx = 2 * M_PI * (i_ - nx) / (bx - ax);

						if (j_ <= ny / 2)
							wy = 2 * M_PI * j_ / (by - ay);
						else
							wy = 2 * M_PI * (j_ - ny) / (by - ay);

						if (k_ <= nz / 2)
							wz = 2 * M_PI * k_ / (bz - az);
						else
							wz = 2 * M_PI * (k_ - nz) / (bz - az);

						w = sqrt(wx * wx + wy * wy + wz * wz);
						c_ = cos(w * c * dt / 2);
						s = sin(w * c * dt / 2);
						if (w == 0)
							w = 1; // чтобы не было деления на 0
						Bx_out[q] = c_ * Bx_out[q] - i * s * ((wy / w) * Ez_out[q] - (wz / w) * Ey_out[q]);
						By_out[q] = c_ * By_out[q] + i * s * ((wx / w) * Ez_out[q] - (wz / w) * Ex_out[q]);
						Bz_out[q] = c_ * Bz_out[q] - i * s * ((wx / w) * Ey_out[q] - (wy / w) * Ex_out[q]);

						Ex_out[q] = c_ * Ex_out[q] + i * s * ((wy / w) * Bz_out[q] - (wz / w) * By_out[q]) + (1 - c_) * wx/w * (Ex_out[q] * wx / w + Ey_out[q] * wy / w + Ez_out[q] * wz / w);
						Ey_out[q] = c_ * Ey_out[q] - i * s * ((wx / w) * Bz_out[q] - (wz / w) * Bx_out[q]) + (1 - c_) * wy/w * (Ex_out[q] * wx / w + Ey_out[q] * wy / w + Ez_out[q] * wz / w);
						Ez_out[q] = c_ * Ez_out[q] + i * s * ((wx / w) * By_out[q] - (wy / w) * Bx_out[q]) + (1 - c_) * wz/w * (Ex_out[q] * wx / w + Ey_out[q] * wy / w + Ez_out[q] * wz / w);
						q++;
					}
		}

		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ex_out.data()), Ex.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ey_out.data()), Ey.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ez_out.data()), Ez.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bx_out.data()), Bx.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(By_out.data()), By.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bz_out.data()), Bz.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jx_out.data()), Jx.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jy_out.data()), Jy.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jz_out.data()), Jz.data(), FFTW_ESTIMATE);
		fftw_execute(plan);

		for (int i = 0; i < nx * ny * nz; i++) {
			Ex[i] /= size;
			Ey[i] /= size;
			Ez[i] /= size;
			Bx[i] /= size;
			By[i] /= size;
			Bz[i] /= size;
			Jx[i] /= size;
			Jy[i] /= size;
			Jz[i] /= size;
		}

		fftw_destroy_plan(plan);
	};

	void Field::PSTD(double dt, double t) {
		int q = 0;
		std::complex<double> i(0.0, 1.0);
		double wx, wy, wz, dtE, dtB;

		Array<std::complex<double>> Ex_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ey_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Ez_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Bx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> By_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Bz_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jx_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jy_out(nx, ny, nz / 2 + 1);
		Array<std::complex<double>> Jz_out(nx, ny, nz / 2 + 1);

		fftw_plan plan;

		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Ex.data(), (fftw_complex*)(Ex_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Ey.data(), (fftw_complex*)(Ey_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Ez.data(), (fftw_complex*)(Ez_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Bx.data(), (fftw_complex*)(Bx_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, By.data(), (fftw_complex*)(By_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Bz.data(), (fftw_complex*)(Bz_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Jx.data(), (fftw_complex*)(Jx_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Jy.data(), (fftw_complex*)(Jy_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_r2c_3d(nx, ny, nz, Jz.data(), (fftw_complex*)(Jz_out.data()), FFTW_ESTIMATE);
		fftw_execute(plan);

		dtE = dt;
		dtB = dt / 2; //для сдвига по времени на -dt/2 сетки B на первой итерации метода

		for (double t1 = 0; t1 < t; t1 += dt) {
			q = 0;
			for (int i_ = 0; i_ < nx; i_++)
				for (int j_ = 0; j_ < ny; j_++)
					for (int k_ = 0; k_ < nz/2+1; k_++) {
						if (i_ <= nx / 2)
							wx = 2 * M_PI * i_ / (bx - ax);
						else
							wx = 2 * M_PI * (i_ - nx) / (bx - ax);

						if (j_ <= ny / 2)
							wy = 2 * M_PI * j_ / (by - ay);
						else
							wy = 2 * M_PI * (j_ - ny) / (by - ay);

						if (k_ <= nz / 2)
							wz = 2 * M_PI * k_ / (bz - az);
						else
							wz = 2 * M_PI * (k_ - nz) / (bz - az);

						Bx_out[q] -= i * c * dtB * (wy * Ez_out[q] - wz * Ey_out[q]);
						By_out[q] += i * c * dtB * (wx * Ez_out[q] - wz * Ex_out[q]);
						Bz_out[q] -= i * c * dtB * (wx * Ey_out[q] - wy * Ex_out[q]);

						Ex_out[q] += i * c * dtE * (wy * Bz_out[q] - wz * By_out[q]) - 4 * M_PI * dtE * Jx_out[q];
						Ey_out[q] += i * c * dtE * (-wx * Bz_out[q] + wz * Bx_out[q]) - 4 * M_PI * dtE * Jy_out[q];
						Ez_out[q] += i * c * dtE * (wx * By_out[q] - wy * Bx_out[q]) - 4 * M_PI * dtE * Jz_out[q];					
						q++;
					}
			dtE = dt;
			dtB = dt;
		}

		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ex_out.data()), Ex.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ey_out.data()), Ey.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Ez_out.data()), Ez.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bx_out.data()), Bx.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(By_out.data()), By.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Bz_out.data()), Bz.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jx_out.data()), Jx.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jy_out.data()), Jy.data(), FFTW_ESTIMATE);
		fftw_execute(plan);
		plan = fftw_plan_dft_c2r_3d(nx, ny, nz, (fftw_complex*)(Jz_out.data()), Jz.data(), FFTW_ESTIMATE);
		fftw_execute(plan);

		for (int i = 0; i <nx * ny * nz; i++) {
			Ex[i] /= size;
			Ey[i] /= size;
			Ez[i] /= size;
			Bx[i] /= size;
			By[i] /= size;
			Bz[i] /= size;
			Jx[i] /= size;
			Jy[i] /= size;
			Jz[i] /= size;
		}

		fftw_destroy_plan(plan);
	};

	void Field::border_update(double dt) { //обновление граничных значений
		//обновление граничных значений E
		for (int i = 0; i < nx; i++) { // j==0
			Ex(i, 0, 0) = Ex(i, 0, 0) - 4.0 * M_PI * dt * Jx(i, 0, 0) + c * dt * (Bz(i, 1, 0) - Bz(i, ny - 1, 0)) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				Ey(i, 0, 0) = Ey(i, 0, 0) - 4.0 * M_PI * dt * Jy(i, 0, 0) - c * dt * (Bz(i + 1, 0, 0) - Bz(i - 1, 0, 0)) / (2.0 * dx);
				Ez(i, 0, 0) = Ez(i, 0, 0) - 4.0 * M_PI * dt * Jz(i, 0, 0) + c * dt * ((By(i + 1, 0, 0) - By(i - 1, 0, 0)) / (2.0 * dx) - (Bx(i, 1, 0) - Bx(i, ny - 1, 0)) / (2.0 * dy));
			}
			if (i == 0) {
				Ey(i, 0, 0) = Ey(i, 0, 0) - 4.0 * M_PI * dt * Jy(i, 0, 0) - c * dt * (Bz(i + 1, 0, 0) - Bz(nx - 1, 0, 0)) / (2.0 * dx);
				Ez(i, 0, 0) = Ez(i, 0, 0) - 4.0 * M_PI * dt * Jz(i, 0, 0) + c * dt * ((By(i + 1, 0, 0) - By(nx - 1, 0, 0)) / (2.0 * dx) - (Bx(i, 1, 0) - Bx(i, ny - 1, 0)) / (2.0 * dy));
			}
			if (i == nx - 1) {
				Ey(i, 0, 0) = Ey(i, 0, 0) - 4.0 * M_PI * dt * Jy(i, 0, 0) - c * dt * (Bz(0, 0, 0) - Bz(i - 1, 0, 0)) / (2.0 * dx);
				Ez(i, 0, 0) = Ez(i, 0, 0) - 4.0 * M_PI * dt * Jz(i, 0, 0) + c * dt * ((By(0, 0, 0) - By(i - 1, 0, 0)) / (2.0 * dx) - (Bx(i, 1, 0) - Bx(i, ny - 1, 0)) / (2.0 * dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j==ny-1
			Ex(i, ny - 1, 0) = Ex(i, ny - 1, 0) - 4.0 * M_PI * dt * Jx(i, ny - 1, 0) + c * dt * (Bz(i, 0, 0) - Bz(i, ny - 2, 0)) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				Ey(i, ny - 1, 0) = Ey(i, ny - 1, 0) - 4.0 * M_PI * dt * Jy(i, ny - 1, 0) - c * dt * (Bz(i + 1, ny - 1, 0) - Bz(i - 1, ny - 1, 0)) / (2.0 * dx);
				Ez(i, ny - 1, 0) = Ez(i, ny - 1, 0) - 4.0 * M_PI * dt * Jz(i, ny - 1, 0) + c * dt * ((By(i + 1, ny - 1, 0) - By(i - 1, ny - 1, 0)) / (2.0 * dx) - (Bx(i, 0, 0) - Bx(i, ny - 2, 0)) / (2.0 * dy));
			}
			if (i == 0) {
				Ey(i, ny - 1, 0) = Ey(i, ny - 1, 0) - 4.0 * M_PI * dt * Jy(i, ny - 1, 0) - c * dt * (Bz(i + 1, ny - 1, 0) - Bz(nx - 1, ny - 1, 0)) / (2.0 * dx);
				Ez(i, ny - 1, 0) = Ez(i, ny - 1, 0) - 4.0 * M_PI * dt * Jz(i, ny - 1, 0) + c * dt * ((By(i + 1, ny - 1, 0) - By(nx - 1, ny - 1, 0)) / (2.0 * dx) - (Bx(i, 0, 0) - Bx(i, ny - 2, 0)) / (2.0 * dy));
			}
			if (i == nx - 1) {
				Ey(i, ny - 1, 0) = Ey(i, ny - 1, 0) - 4.0 * M_PI * dt * Jy(i, ny - 1, 0) - c * dt * (Bz(0, ny - 1, 0) - Bz(i - 1, ny - 1, 0)) / (2.0 * dx);
				Ez(i, ny - 1, 0) = Ez(i, ny - 1, 0) - 4.0 * M_PI * dt * Jz(i, ny - 1, 0) + c * dt * ((By(0, ny - 1, 0) - By(i - 1, ny - 1, 0)) / (2.0 * dx) - (Bx(i, 0, 0) - Bx(i, ny - 2, 0)) / (2.0 * dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) { // i==0
			Ex(0, j, 0) = Ex(0, j, 0) - 4.0 * M_PI * dt * Jx(0, j, 0) + c * dt * (Bz(0, j + 1, 0) - Bz(0, j - 1, 0)) / (2.0 * dy);
			Ey(0, j, 0) = Ey(0, j, 0) - 4.0 * M_PI * dt * Jy(0, j, 0) - c * dt * (Bz(1, j, 0) - Bz(nx - 1, j, 0)) / (2.0 * dx);
			Ez(0, j, 0) = Ez(0, j, 0) - 4.0 * M_PI * dt * Jz(0, j, 0) + c * dt * ((By(1, j, 0) - By(nx - 1, j, 0)) / (2.0 * dx) - (Bx(0, j + 1, 0) - Bx(0, j - 1, 0)) / (2.0 * dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i==nx-1
			Ex(nx - 1, j, 0) = Ex(nx - 1, j, 0) - 4.0 * M_PI * dt * Jx(nx - 1, j, 0) + c * dt * (Bz(nx - 1, j + 1, 0) - Bz(nx - 1, j - 1, 0)) / (2.0 * dy);
			Ey(nx - 1, j, 0) = Ey(nx - 1, j, 0) - 4.0 * M_PI * dt * Jy(nx - 1, j, 0) - c * dt * (Bz(0, j, 0) - Bz(nx - 2, j, 0)) / (2.0 * dx);
			Ez(nx - 1, j, 0) = Ez(nx - 1, j, 0) - 4.0 * M_PI * dt * Jz(nx - 1, j, 0) + c * dt * ((By(0, j, 0) - By(nx - 2, j, 0)) / (2.0 * dx) - (Bx(nx - 1, j + 1, 0) - Bx(nx - 1, j - 1, 0)) / (2.0 * dy));
		}
		//обновление граничных значений B
		for (int i = 0; i < nx; i++) { // j==0
			Bx(i, 0, 0) = Bx(i, 0, 0) - c * dt * (Ez(i, 1, 0) - Ez(i, ny - 1, 0)) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				By(i, 0, 0) = By(i, 0, 0) + c * dt * (Ez(i + 1, 0, 0) - Ez(i - 1, 0, 0)) / (2.0 * dx);
				Bz(i, 0, 0) = Bz(i, 0, 0) - c * dt * ((Ey(i + 1, 0, 0) - Ey(i - 1, 0, 0)) / (2.0 * dx) - (Ex(i, 1, 0) - Ex(i, ny - 1, 0)) / (2.0 * dy));
			}
			if (i == 0) {
				By(i, 0, 0) = By(i, 0, 0) + c * dt * (Ez(i + 1, 0, 0) - Ez(nx - 1, 0, 0)) / (2.0 * dx);
				Bz(i, 0, 0) = Bz(i, 0, 0) - c * dt * ((Ey(i + 1, 0, 0) - Ey(nx - 1, 0, 0)) / (2.0 * dx) - (Ex(i, 1, 0) - Ex(i, ny - 1, 0)) / (2.0 * dy));
			}
			if (i == nx - 1) {
				By(i, 0, 0) = By(i, 0, 0) + c * dt * (Ez(0, 0, 0) - Ez(i - 1, 0, 0)) / (2.0 * dx);
				Bz(i, 0, 0) = Bz(i, 0, 0) - c * dt * ((Ey(0, 0, 0) - Ey(i - 1, 0, 0)) / (2.0 * dx) - (Ex(i, 1, 0) - Ex(i, ny - 1, 0)) / (2.0 * dy));
			}
		}

		for (int i = 0; i < nx; i++) { // j==ny-1
			Bx(i, ny - 1, 0) = Bx(i, ny - 1, 0) - c * dt * (Ez(i, 0, 0) - Ez(i, ny - 2, 0)) / (2.0 * dy);
			if ((i != 0) && (i != nx - 1)) {
				By(i, ny - 1, 0) = By(i, ny - 1, 0) + c * dt * (Ez(i + 1, ny - 1, 0) - Ez(i - 1, ny - 1, 0)) / (2.0 * dx);
				Bz(i, ny - 1, 0) = Bz(i, ny - 1, 0) - c * dt * ((Ey(i + 1, ny - 1, 0) - Ey(i - 1, ny - 1, 0)) / (2.0 * dx) - (Ex(i, 0, 0) - Ex(i, ny - 2, 0)) / (2.0 * dy));
			}
			if (i == 0) {
				By(i, ny - 1, 0) = By(i, ny - 1, 0) + c * dt * (Ez(i + 1, ny - 1, 0) - Ez(nx - 1, ny - 1, 0)) / (2.0 * dx);
				Bz(i, ny - 1, 0) = Bz(i, ny - 1, 0) - c * dt * ((Ey(i + 1, ny - 1, 0) - Ey(nx - 1, ny - 1, 0)) / (2.0 * dx) - (Ex(i, 0, 0) - Ex(i, ny - 2, 0)) / (2.0 * dy));
			}
			if (i == nx - 1) {
				By(i, ny - 1, 0) = By(i, ny - 1, 0) + c * dt * (Ez(0, ny - 1, 0) - Ez(i - 1, ny - 1, 0)) / (2.0 * dx);
				Bz(i, ny - 1, 0) = Bz(i, ny - 1, 0) - c * dt * ((Ey(0, ny - 1, 0) - Ey(i - 1, ny - 1, 0)) / (2.0 * dx) - (Ex(i, 0, 0) - Ex(i, ny - 2, 0)) / (2.0 * dy));
			}
		}

		for (int j = 1; j < ny - 1; j++) {  // i==0
			Bx(0, j, 0) = Bx(0, j, 0) - c * dt * (Ez(0, j + 1, 0) - Ez(0, j - 1, 0)) / (2.0 * dy);
			By(0, j, 0) = By(0, j, 0) + c * dt * (Ez(1, j, 0) - Ez(nx - 1, j, 0)) / (2.0 * dx);
			Bz(0, j, 0) = Bz(0, j, 0) - c * dt * ((Ey(1, j, 0) - Ey(nx - 1, j, 0)) / (2.0 * dx) - (Ex(0, j + 1, 0) - Ex(0, j - 1, 0)) / (2.0 * dy));
		}

		for (int j = 1; j < ny - 1; j++) { // i==nx-1
			Bx(nx - 1, j, 0) = Bx(nx - 1, j, 0) - c * dt * (Ez(nx - 1, j + 1, 0) - Ez(nx - 1, j - 1, 0)) / (2.0 * dy);
			By(nx - 1, j, 0) = By(nx - 1, j, 0) + c * dt * (Ez(0, j, 0) - Ez(nx - 2, j, 0)) / (2.0 * dx);
			Bz(nx - 1, j, 0) = Bz(nx - 1, j, 0) - c * dt * ((Ey(0, j, 0) - Ey(nx - 2, j, 0)) / (2.0 * dx) - (Ex(nx - 1, j + 1, 0) - Ex(nx - 1, j - 1, 0)) / (2.0 * dy));
		}

	};

	void Field::FDTD(double dt) { //обновление сеточных значений
	//обновление Е
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++) {
				Ex(i, j, 0) = Ex(i, j, 0) - 4.0 * M_PI * dt * Jx(i, j, 0) + c * dt * (Bz(i, j + 1, 0) - Bz(i, j - 1, 0)) / (2.0 * dy);
				Ey(i, j, 0) = Ey(i, j, 0) - 4.0 * M_PI * dt * Jy(i, j, 0) - c * dt * (Bz(i + 1, j, 0) - Bz(i - 1, j, 0)) / (2.0 * dx);
				Ez(i, j, 0) = Ez(i, j, 0) - 4.0 * M_PI * dt * Jz(i, j, 0) + c * dt * ((By(i + 1, j, 0) - By(i - 1, j, 0)) / (2.0 * dx) - (Bx(i, j + 1, 0) - Bx(i, j - 1, 0)) / (2.0 * dy));
			}
		(*this).border_update(dt);
		//обновление В
		for (int i = 1; i < nx - 1; i++)
			for (int j = 1; j < ny - 1; j++) {
				Bx(i, j, 0) = Bx(i, j, 0) - c * dt * (Ez(i, j + 1, 0) - Ez(i, j - 1, 0)) / (2.0 * dy);
				By(i, j, 0) = By(i, j, 0) + c * dt * (Ez(i + 1, j, 0) - Ez(i - 1, j, 0)) / (2.0 * dx);
				Bz(i, j, 0) = Bz(i, j, 0) - c * dt * ((Ey(i + 1, j, 0) - Ey(i - 1, j, 0)) / (2.0 * dx) - (Ex(i, j + 1, 0) - Ex(i, j - 1, 0)) / (2.0 * dy));
			}
	}

	void Field::FDTD_with_shift(double dtB, double dtE) {

		for (int i = 0; i < nx - 1; i++)
			for (int j = 0; j < ny - 1; j++)
				for (int k = 0; k < nz - 1; k++) {
					Bx(i,j,k) = Bx(i,j,k) + c * dtB * ((Ey(i,j,k + 1) - Ey(i,j,k)) / dz - (Ez(i,j + 1,k) - Ez(i,j,k)) / dy);
					By(i,j,k) = By(i,j,k) + c * dtB * ((Ez(i + 1,j,k) - Ez(i,j,k)) / dx - (Ex(i,j,k + 1) - Ex(i,j,k)) / dz);
					Bz(i,j,k) = Bz(i,j,k) + c * dtB * ((Ex(i,j + 1,k) - Ex(i,j,k)) / dy - (Ey(i + 1,j,k) - Ey(i,j,k)) / dx);
				}
		(*this).border_update_with_shift(dtB, dtE);
		for (int i = 1; i < nx; i++)
			for (int j = 1; j < ny; j++)
				for (int k = 1; k < nz; k++) {
					Ex(i,j,k) = Ex(i,j,k) - 4.0 * M_PI * dtE * Jx(i,j,k) + c * dtE * ((Bz(i,j,k) - Bz(i,j - 1,k)) / dy - (By(i,j,k) - By(i,j,k - 1)) / dz);
					Ey(i,j,k) = Ey(i,j,k) - 4.0 * M_PI * dtE * Jy(i,j,k) + c * dtE * ((Bx(i,j,k) - Bx(i,j,k - 1)) / dz - (Bz(i,j,k) - Bz(i - 1,j,k)) / dx);
					Ez(i,j,k) = Ez(i,j,k) - 4.0 * M_PI * dtE * Jz(i,j,k) + c * dtE * ((By(i,j,k) - By(i - 1,j,k)) / dx - (Bx(i,j,k) - Bx(i,j - 1,k)) / dy);
				}
	}

	void Field::border_update_with_shift(double dtB, double dtE) {

		//границы В 
		for (int i = 0; i < nx; i++)  // j==ny-1
			for (int k = 0; k < nz; k++) {
				if ((k == nz - 1) && (i != nx - 1)) {
					Bx(i,ny - 1,k) = Bx(i,ny - 1,k) + c * dtB * ((Ey(i,ny - 1,0) - Ey(i,ny - 1,k)) / dz - (Ez(i,0,k) - Ez(i,ny - 1,k)) / dy);
					By(i,ny - 1,k) = By(i,ny - 1,k) + c * dtB * ((Ez(i + 1,ny - 1,k) - Ez(i,ny - 1,k)) / dx - (Ex(i,ny - 1,0) - Ex(i,ny - 1,k)) / dz);
					Bz(i,ny - 1,k) = Bz(i,ny - 1,k) + c * dtB * ((Ex(i,0,k) - Ex(i,ny - 1,k)) / dy - (Ey(i + 1,ny - 1,k) - Ey(i,ny - 1,k)) / dx);
				}
				if ((k != nz - 1) && (i == nx - 1)) {
					Bx(i,ny - 1,k) = Bx(i,ny - 1,k) + c * dtB * ((Ey(i,ny - 1,k + 1) - Ey(i,ny - 1,k)) / dz - (Ez(i,0,k) - Ez(i,ny - 1,k)) / dy);
					By(i,ny - 1,k) = By(i,ny - 1,k) + c * dtB * ((Ez(0,ny - 1,k) - Ez(i,ny - 1,k)) / dx - (Ex(i,ny - 1,k + 1) - Ex(i,ny - 1,k)) / dz);
					Bz(i,ny - 1,k) = Bz(i,ny - 1,k) + c * dtB * ((Ex(i,0,k) - Ex(i,ny - 1,k)) / dy - (Ey(0,ny - 1,k) - Ey(i,ny - 1,k)) / dx);
				}
				if ((k == nz - 1) && (i == nx - 1)) {
					Bx(i,ny - 1,k) = Bx(i,ny - 1,k) + c * dtB * ((Ey(i,ny - 1,0) - Ey(i,ny - 1,k)) / dz - (Ez(i,0,k) - Ez(i,ny - 1,k)) / dy);
					By(i,ny - 1,k) = By(i,ny - 1,k) + c * dtB * ((Ez(0,ny - 1,k) - Ez(i,ny - 1,k)) / dx - (Ex(i,ny - 1,0) - Ex(i,ny - 1,k)) / dz);
					Bz(i,ny - 1,k) = Bz(i,ny - 1,k) + c * dtB * ((Ex(i,0,k) - Ex(i,ny - 1,k)) / dy - (Ey(0,ny - 1,k) - Ey(i,ny - 1,k)) / dx);
				}
			}

		for (int j = 0; j < ny - 1; j++) // i==nx-1
			for (int k = 0; k < nz; k++) {
				if (k == nz - 1) {
					Bx(nx - 1,j,k) = Bx(nx - 1,j,k) + c * dtB * ((Ey(nx - 1,j,0) - Ey(nx - 1,j,k)) / dz - (Ez(nx - 1,j + 1,k) - Ez(nx - 1,j,k)) / dy);
					By(nx - 1,j,k) = By(nx - 1,j,k) + c * dtB * ((Ez(0,j,k) - Ez(nx - 1,j,k)) / dx - (Ex(nx - 1,j,0) - Ex(nx - 1,j,k)) / dz);
					Bz(nx - 1,j,k) = Bz(nx - 1,j,k) + c * dtB * ((Ex(nx - 1,j + 1,k) - Ex(nx - 1,j,k)) / dy - (Ey(0,j,k) - Ey(nx - 1,j,k)) / dx);
				}
				else {
					Bx(nx - 1,j,k) = Bx(nx - 1,j,k) + c * dtB * ((Ey(nx - 1,j,k + 1) - Ey(nx - 1,j,k)) / dz - (Ez(nx - 1,j + 1,k) - Ez(nx - 1,j,k)) / dy);
					By(nx - 1,j,k) = By(nx - 1,j,k) + c * dtB * ((Ez(0,j,k) - Ez(nx - 1,j,k)) / dx - (Ex(nx - 1,j,k + 1) - Ex(nx - 1,j,k)) / dz);
					Bz(nx - 1,j,k) = Bz(nx - 1,j,k) + c * dtB * ((Ex(nx - 1,j + 1,k) - Ex(nx - 1,j,k)) / dy - (Ey(0,j,k) - Ey(nx - 1,j,k)) / dx);
				}
			}

		for (int i = 0; i < nx - 1; i++) // k==nz-1
			for (int j = 0; j < ny - 1; j++) {
				Bx(i,j,nz - 1) = Bx(i,j,nz - 1) + c * dtB * ((Ey(i,j,0) - Ey(i,j,nz - 1)) / dz - (Ez(i,j + 1,nz - 1) - Ez(i,j,nz - 1)) / dy);
				By(i,j,nz - 1) = By(i,j,nz - 1) + c * dtB * ((Ez(i + 1,j,nz - 1) - Ez(i,j,nz - 1)) / dx - (Ex(i,j,0) - Ex(i,j,nz - 1)) / dz);
				Bz(i,j,nz - 1) = Bz(i,j,nz - 1) + c * dtB * ((Ex(i,j + 1,nz - 1) - Ex(i,j,nz - 1)) / dy - (Ey(i + 1,j,nz - 1) - Ey(i,j,nz - 1)) / dx);
			}


		//границы Е
		for (int i = 0; i < nx; i++) // j==0
			for (int k = 0; k < nz; k++) {
				if ((k == 0) && (i != 0)) {
					Ex(i,0,k) = Ex(i,0,k) - 4.0 * M_PI * dtE * Jx(i,0,k) + c * dtE * ((Bz(i,0,k) - Bz(i,ny - 1,k)) / dy - (By(i,0,k) - By(i,0,nz - 1)) / dz);
					Ey(i,0,k) = Ey(i,0,k) - 4.0 * M_PI * dtE * Jy(i,0,k) + c * dtE * ((Bx(i,0,k) - Bx(i,0,nz - 1)) / dz - (Bz(i,0,k) - Bz(i - 1,0,k)) / dx);
					Ez(i,0,k) = Ez(i,0,k) - 4.0 * M_PI * dtE * Jz(i,0,k) + c * dtE * ((By(i,0,k) - By(i - 1,0,k)) / dx - (Bx(i,0,k) - Bx(i,ny - 1,k)) / dy);
				}
				if ((k == 0) && (i == 0)) {
					Ex(i,0,k) = Ex(i,0,k) - 4.0 * M_PI * dtE * Jx(i,0,k) + c * dtE * ((Bz(i,0,k) - Bz(i,ny - 1,k)) / dy - (By(i,0,k) - By(i,0,nz - 1)) / dz);
					Ey(i,0,k) = Ey(i,0,k) - 4.0 * M_PI * dtE * Jy(i,0,k) + c * dtE * ((Bx(i,0,k) - Bx(i,0,nz - 1)) / dz - (Bz(i,0,k) - Bz(nx - 1,0,k)) / dx);
					Ez(i,0,k) = Ez(i,0,k) - 4.0 * M_PI * dtE * Jz(i,0,k) + c * dtE * ((By(i,0,k) - By(nx - 1,0,k)) / dx - (Bx(i,0,k) - Bx(i,ny - 1,k)) / dy);
				}
				if ((k != 0) && (i == 0)) {
					Ex(i,0,k) = Ex(i,0,k) - 4.0 * M_PI * dtE * Jx(i,0,k) + c * dtE * ((Bz(i,0,k) - Bz(i,ny - 1,k)) / dy - (By(i,0,k) - By(i,0,k - 1)) / dz);
					Ey(i,0,k) = Ey(i,0,k) - 4.0 * M_PI * dtE * Jy(i,0,k) + c * dtE * ((Bx(i,0,k) - Bx(i,0,k - 1)) / dz - (Bz(i,0,k) - Bz(nx - 1,0,k)) / dx);
					Ez(i,0,k) = Ez(i,0,k) - 4.0 * M_PI * dtE * Jz(i,0,k) + c * dtE * ((By(i,0,k) - By(nx - 1,0,k)) / dx - (Bx(i,0,k) - Bx(i,ny - 1,k)) / dy);
				}
			}

		for (int j = 1; j < ny; j++) // i==0
			for (int k = 0; k < nz; k++) {
				if (k == 0) {
					Ex(0,j,k) = Ex(0,j,k) - 4.0 * M_PI * dtE * Jx(0,j,k) + c * dtE * ((Bz(0,j,k) - Bz(0,j - 1,k)) / dy - (By(0,j,k) - By(0,j,nz - 1)) / dz);
					Ey(0,j,k) = Ey(0,j,k) - 4.0 * M_PI * dtE * Jy(0,j,k) + c * dtE * ((Bx(0,j,k) - Bx(0,j,nz - 1)) / dz - (Bz(0,j,k) - Bz(nx - 1,j,k)) / dx);
					Ez(0,j,k) = Ez(0,j,k) - 4.0 * M_PI * dtE * Jz(0,j,k) + c * dtE * ((By(0,j,k) - By(nx - 1,j,k)) / dx - (Bx(0,j,k) - Bx(0,j - 1,k)) / dy);
				}
				else {
					Ex(0,j,k) = Ex(0,j,k) - 4.0 * M_PI * dtE * Jx(0,j,k) + c * dtE * ((Bz(0,j,k) - Bz(0,j - 1,k)) / dy - (By(0,j,k) - By(0,j,k - 1)) / dz);
					Ey(0,j,k) = Ey(0,j,k) - 4.0 * M_PI * dtE * Jy(0,j,k) + c * dtE * ((Bx(0,j,k) - Bx(0,j,k - 1)) / dz - (Bz(0,j,k) - Bz(nx - 1,j,k)) / dx);
					Ez(0,j,k) = Ez(0,j,k) - 4.0 * M_PI * dtE * Jz(0,j,k) + c * dtE * ((By(0,j,k) - By(nx - 1,j,k)) / dx - (Bx(0,j,k) - Bx(0,j - 1,k)) / dy);
				}

			}

		for (int i = 1; i < nx; i++) // k==0
			for (int j = 1; j < ny; j++) {
				Ex(i,j,0) = Ex(i,j,0) - 4.0 * M_PI * dtE * Jx(i,j,0) + c * dtE * ((Bz(i,j,0) - Bz(i,j - 1,0)) / dy - (By(i,j,0) - By(i,j,nz - 1)) / dz);
				Ey(i,j,0) = Ey(i,j,0) - 4.0 * M_PI * dtE * Jy(i,j,0) + c * dtE * ((Bx(i,j,0) - Bx(i,j,nz - 1)) / dz - (Bz(i,j,0) - Bz(i - 1,j,0)) / dx);
				Ez(i,j,0) = Ez(i,j,0) - 4.0 * M_PI * dtE * Jz(i,j,0) + c * dtE * ((By(i,j,0) - By(i - 1,j,0)) / dx - (Bx(i,j,0) - Bx(i,j - 1,0)) / dy);
			}
	}
};
