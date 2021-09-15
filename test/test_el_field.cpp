#include "el_field.h"
#include <gtest.h>

//Для тестов беру узел [3][4]

TEST(Field, right_Ex)
{
	Field f;
	double dt = 1E-11, t = 1E-10, real_value,e;
	real_value = 0; 
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

TEST(Field, right_Ey)
{
	Field f;
	double dt = 1E-10, t = 1E-10, real_value, e;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
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

TEST(Field, right_Ez)
{
	Field f;
	double dt = 1E-10, t = 1E-10, real_value, e;
	real_value = 0;
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

TEST(Field, right_Bx)
{
	Field f;
	double dt = 1E-10, t = 1E-10, real_value, e;
	real_value = 0;
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

TEST(Field, right_By)
{
	Field f;
	double dt = 1E-10, t = 1E-10, real_value, e;
	real_value = 0;
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

TEST(Field, right_Bz)
{
	Field f;
	double dt = 1E-10, t = 1E-10, real_value, e;
	real_value = sin(2 * M_PI * (3 * dx - c * t) / (bx - ax));
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



 