#include "el_field.h"
#include <gtest.h>

const double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az=0,bz=0;

//Для тестов беру узел (3,4)

TEST(Field, FDTD_right_Ex_axis_x)
{
	double dt = 0.001, t = 0.1, real_value,e,dx = 1E+9, dy = 1E+10,dz= 1E+9;
	real_value = 0; 
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Ex(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e)<0.1);
}

TEST(Field, FDTD_right_Ey_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10,dz = 1E+9;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Ey(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Ez_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Ez(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Bx_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Bx(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_By_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.By(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Bz_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i,j,0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Bz(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}


TEST(Field, FDTD_right_Ex_axis_y) //тоже беру узел (3,4)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i,j,0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i,j,0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Ex(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Ey_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i,j,0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i,j,0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Ey(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Ez_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i,j,0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i,j,0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Ez(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Bx_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i,j,0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i,j,0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Bx(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_By_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i,j,0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i,j,0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.By(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_right_Bz_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1E+9;
	real_value = -sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ex(i,j,0) = sin(2 * M_PI * (dy * j) / (by - ay));
			f.Bz(i,j,0) = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	e = f.Bz(3,4,0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_with_shift_right_Ey_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	f.FDTD_with_shift(dt / 2, dt);
	for (double i = dt; i < t; i += dt)
		f.FDTD_with_shift(dt,dt);
	e = f.Ey(3, 4, 0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, FDTD_with_shift_right_Bz_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	real_value = sin(2 * M_PI * (3 * dx - c * (t-dt/2)) / (bx - ax));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.Ey(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.Bz(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	f.FDTD_with_shift(dt / 2, dt);
	for (double i = dt; i < t; i += dt)
		f.FDTD_with_shift(dt, dt);
	e = f.Bz(3, 4, 0) - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}



 