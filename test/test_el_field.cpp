#include "el_field.h"
#include <gtest.h>

//Для тестов беру узел [3][4]

TEST(Field, right_Ex)
{
	Field f;
	double dt = 0.001, t = 0.1, real_value,e;
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
	double dt = 0.001, t = 0.1, real_value, e;
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
	double dt = 0.001, t = 0.1, real_value, e;
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
	double dt = 0.001, t = 0.1, real_value, e;
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
	double dt = 0.001, t = 0.1, real_value, e;
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
	double dt = 0.001, t = 0.1, real_value, e;
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

TEST(Field, boundary_value) { //в узле [0][0]
	Field f1, f2;
	double dt = 0.001, t = 0.001, e;
	Field_characteristics res1,res2;
	for (int i = 0; i < f1.nx; i++)
		for (int j = 0; j < f1.ny; j++) {
			f1.v.E[0][i][j] = sin(2.0 * M_PI * (dy * j) / (by - ay));
			f1.v.B[2][i][j] = sin(2.0 * M_PI * (dy * j) / (by - ay));
		}
	//у меня квадратная сетка, в точке [0][0] значение при выбранном t и dt должно быть такое же, как получалось при распространении волны вдоль оси х, но оно получается такое же с другим знаком
	for (int i = 0; i < f2.nx; i++)
		for (int j = 0; j < f2.ny; j++) {
			f2.v.E[1][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
			f2.v.B[2][i][j] = sin(2 * M_PI * (dx * i) / (bx - ax));
		}

	for (double i = 0; i < t; i += dt) {
		res1 = f1.FDTD(dt);
		res2 = f2.FDTD(dt);
	}
	std::cout << res1.E[0][0][0] << ' ' << res2.E[1][0][0] << std::endl;
	e = res1.E[0][0][0] - res2.E[1][0][0];
	EXPECT_EQ(true, abs(e) < 0.001);

}

TEST(Field, right_Ex_axis_y) //беру узел [3][4]
{
	Field f;
	double dt = 0.001, t = 0.1, real_value, e;
	real_value = real_value = sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));;
	Field_characteristics res;
	for (int i = 0; i < f.nx; i++)
		for (int j = 0; j < f.ny; j++) {
			f.v.E[0][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
			f.v.B[2][i][j] = sin(2 * M_PI * (dy * j) / (by - ay));
		}
	for (double i = 0; i < t; i += dt)
		res = f.FDTD(dt);
	e = res.E[0][3][4] - real_value;
	EXPECT_EQ(true, abs(e) < 0.1);
}

//со значениями B аналогичная проблема





 