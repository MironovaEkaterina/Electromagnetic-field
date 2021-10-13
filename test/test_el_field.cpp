#include "el_field.h"
#include <gtest.h>

const double ax = 0, bx = 5E+10, ay = 0, by = 5E+10;

//Для тестов беру узел [3][4]

TEST(Field, right_Ex_axis_x)
{
	double dt = 0.001, t = 0.1, real_value,e,dx = 1E+9, dy = 1E+10;
	real_value = 0; 
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[0][3][4] - real_value;
	EXPECT_EQ(true, abs(e)<0.1);
}

TEST(Field, right_Ey_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[1][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Ez_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[2][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bx_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.B[0][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_By_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.B[1][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bz_axis_x)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+9, dy = 1E+10;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.B[2][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}


TEST(Field, right_Ex_axis_y) //тоже беру узел [3][4]
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9;
	real_value = sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[0][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Ey_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[1][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Ez_axis_y) 
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[2][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bx_axis_y) 
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.B[0][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_By_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9;
	real_value = 0;
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.B[1][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, right_Bz_axis_y)
{
	double dt = 0.001, t = 0.1, real_value, e, dx = 1E+10, dy = 1E+9;
	real_value = -sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	Field f(dx, dy, ax, ay, bx, by);
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = -sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.B[2][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

TEST(Field, speed_comparison) 
{
	double dt = 0.001, dx = 1E+9, dy = 1E+10;
	double lambda = bx - ax;
	double k = 2 * M_PI / lambda;
	double v = (2 / (k * dt)) * asin((c * dt / dx) * sin(k * dx / 2));
	std::cout << v<<std::endl;
	EXPECT_EQ(true, v<c);
}





 