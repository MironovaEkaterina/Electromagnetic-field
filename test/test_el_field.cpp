#include "el_field.h"
#include <gtest.h>

const double ax = 0, bx = 5E+10, ay = 0, by = 5E+10,az=0,bz=0;

//Для тестов беру узел [3][4]

TEST(Field, right_Ex_axis_x)
{
	double dt = 0.001, t = 0.1, real_value,e,dx = 1E+9, dy = 1E+10,dz=1;
	real_value = 0; 
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Ex[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e)<0.1);
}

TEST(Field, right_Ey_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10,dz = 1;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Ey[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Ez_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Ez[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bx_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Bx[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_By_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.By[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bz_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10, dz = 1;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ey[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.Bz[i][j][0] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Bz[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}


TEST(Field, right_Ex_axis_y) //тоже беру узел [3][4]
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1;
	real_value = sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ex[i][j][0] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.Bz[i][j][0] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Ex[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Ey_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ex[i][j][0] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.Bz[i][j][0] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Ey[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Ez_axis_y) 
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ex[i][j][0] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.Bz[i][j][0] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Ez[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bx_axis_y) 
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ex[i][j][0] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.Bz[i][j][0] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Bx[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_By_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1;
	real_value = 0;
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ex[i][j][0] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.Bz[i][j][0] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.By[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bz_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9, dz = 1;
	real_value = -sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	Field f(dx, dy, dz, ax, ay, az, bx, by, bz);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.Ex[i][j][0] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.Bz[i][j][0] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.Bz[3][4][0] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}




 