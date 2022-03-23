#include "el_field.h"
#include <gtest.h>

//const double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az=0,bz=0;

//Для тестов беру узел (30,4)

class Test_FDTD_axis_x : public ::testing::Test
{
public:
	double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 0;
	double dt = 0.001, t = 0.17, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	Field f;
	
	void SetUp()
	{
		f.Initialize(dx, dy, dz, ax, ay, az, bx, by, bz, (bx - ax) / dx + 1, (by - ay) / dy + 1, 1);
		for (int i = 0; i < f.nx; i++)
			for (int j = 0; j < f.ny; j++) {
				f.Ey(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
				f.Bz(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			}
	}

	void TearDown(){}
};

TEST_F(Test_FDTD_axis_x, FDTD_right_Ex_axis_x)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.Ex(3, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_right_Ey_axis_x)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = sin(2 * M_PI * (30 * dx - c * t) / (bx - ax));
	EXPECT_NEAR(f.Ey(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_right_Ez_axis_x)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.Ez(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_right_Bx_axis_x)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.Bx(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_right_By_axis_x)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.By(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_right_Bz_axis_x)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = sin(2 * M_PI * (30 * dx - c * t) / (bx - ax));
	EXPECT_NEAR(f.Bz(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_with_shift_right_Ey_axis_x)
{
	real_value = sin(2 * M_PI * (30 * dx - c * t) / (bx - ax));
	f.FDTD_with_shift(dt / 2, dt);
	for (double i = dt; i < t; i += dt)
		f.FDTD_with_shift(dt, dt);
	EXPECT_NEAR(f.Ey(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_x, FDTD_with_shift_right_Bz_axis_x)
{
	real_value = sin(2 * M_PI * (30 * dx - c * (t - dt / 2)) / (bx - ax));
	f.FDTD_with_shift(dt / 2, dt);
	for (double i = dt; i < t; i += dt)
		f.FDTD_with_shift(dt, dt);
	EXPECT_NEAR(f.Bz(30, 4, 0), real_value, 0.1);
}


class Test_FDTD_axis_y : public ::testing::Test
{
public:
	double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 0;
	double dt = 0.001, t = 0.17, real_value, dx = 1E+9, dy = 1E+10, dz = 1E+9;
	Field f;

	void SetUp()
	{
		f.Initialize(dx, dy, dz, ax, ay, az, bx, by, bz, (bx - ax) / dx + 1, (by - ay) / dy + 1, 1);
		for (int i = 0; i < f.nx; i++)
			for (int j = 0; j < f.ny; j++) {
				f.Ex(i, j, 0) = sin(2 * M_PI * (dy * j) / (by - ay));
				f.Bz(i, j, 0) = -sin(2 * M_PI * (dy * j) / (by - ay));
			}
		for (double i = 0; i < t; i += dt)
			f.FDTD(dt);
	}

	void TearDown() {}
};

TEST_F(Test_FDTD_axis_y, FDTD_right_Ex_axis_y) //тоже беру узел (30,4)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	EXPECT_NEAR(f.Ex(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_y, FDTD_right_Ey_axis_y)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.Ey(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_y, FDTD_right_Ez_axis_y)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.Ez(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_y, FDTD_right_Bx_axis_y)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.Bx(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_y, FDTD_right_By_axis_y)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = 0.0;
	EXPECT_NEAR(f.By(30, 4, 0), real_value, 0.1);
}

TEST_F(Test_FDTD_axis_y, FDTD_right_Bz_axis_y)
{
	for (double i = 0; i < t; i += dt)
		f.FDTD(dt);
	real_value = -sin(2.0 * M_PI * (4 * dy - c * t) / (by - ay));
	EXPECT_NEAR(f.Bz(30, 4, 0), real_value, 0.1);
}

class Test_PSTD_axis_x : public ::testing::Test
{
public:
	double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 1E+9;
	double dt = 0.001, t = 0.17, real_value, dx = 1E+9, dy = 1E+9, dz = 1E+9;
	Field f;

	void SetUp()
	{
		f.Initialize(dx, dy, dz, ax, ay, az, bx, by, bz, (bx - ax) / dx + 1, (by - ay) / dy + 1, 1);
		for (int i = 0; i < f.nx; i++)
			for (int j = 0; j < f.ny; j++) {
				f.Ey(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
				f.Bz(i, j, 0) = sin(2 * M_PI * (dx * i) / (bx - ax));
			}
		f.PSTD(dt, t);
	}

	void TearDown() {}
};

TEST_F(Test_PSTD_axis_x, PSTD_right_Ex_axis_x)
{
	real_value = 0.0;
	EXPECT_NEAR(f.Ex(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_x, PSTD_right_Ey_axis_x)
{
	real_value = sin(2 * M_PI * (30 * dx - c * t) / (bx - ax));
	EXPECT_NEAR(f.Ey(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_x, PSTD_right_Ez_axis_x)
{
	real_value = 0.0;
	EXPECT_NEAR(f.Ez(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_x, PSTD_right_Bx_axis_x)
{
	real_value = 0.0;
	EXPECT_NEAR(f.Bx(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_x, PSTD_right_By_axis_x)
{
	real_value = 0.0;
	EXPECT_NEAR(f.By(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_x, PSTD_right_Bz_axis_x)
{
	real_value = sin(2 * M_PI * (30 * dx - c * (t - dt / 2)) / (bx - ax));
	EXPECT_NEAR(f.Bz(30, 40, 0), real_value, 0.1);
}

class Test_PSTD_axis_y : public ::testing::Test
{
public:
	double ax = 0, bx = 5E+10, ay = 0, by = 5E+10, az = 0, bz = 1E+9;
	double dt = 0.001, t = 0.17, dx = 1E+9, dy = 1E+9, dz = 1E+9, real_value;
	Field f;

	void SetUp()
	{
		f.Initialize(dx, dy, dz, ax, ay, az, bx, by, bz, (bx - ax) / dx + 1, (by - ay) / dy + 1, 1);
		for (int i = 0; i < f.nx; i++)
			for (int j = 0; j < f.ny; j++) {
				f.Ex(i, j, 0) = sin(2 * M_PI * (dy * j) / (by - ay));
				f.Bz(i, j, 0) = -sin(2 * M_PI * (dy * j) / (by - ay));
			}
		f.PSTD(dt,t);
	}

	void TearDown() {}
};

TEST_F(Test_PSTD_axis_y, PSTD_right_Ex_axis_y) 
{
	real_value = sin(2.0 * M_PI * (40 * dy - c * t) / (by - ay));
	EXPECT_NEAR(f.Ex(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_y, PSTD_right_Ey_axis_y)
{
	real_value = 0.0;
	EXPECT_NEAR(f.Ey(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_y, PSTD_right_Ez_axis_y)
{
	real_value = 0.0;
	EXPECT_NEAR(f.Ez(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_y, PSTD_right_Bx_axis_y)
{
	real_value = 0.0;
	EXPECT_NEAR(f.Bx(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_y, PSTD_right_By_axis_y)
{
	real_value = 0.0;
	EXPECT_NEAR(f.By(30, 40, 0), real_value, 0.1);
}

TEST_F(Test_PSTD_axis_y, PSTD_right_Bz_axis_y)
{
	real_value = -sin(2.0 * M_PI * (40 * dy - c * (t - dt / 2)) / (by - ay));
	EXPECT_NEAR(f.Bz(30, 40, 0), real_value, 0.1);
}