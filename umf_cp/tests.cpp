/*tests.cpp*/
#include "tests.h"

namespace tests
{
	double Tests::Lambda(double x, double y, double t)
	{
		return 1;
	}
	double Tests::Sigma(double x, double y, double t)
	{
		return 1;
	}
	double Tests::Gamma(double x, double y, double t)
	{
		return 1;
	}
	double Tests::Ug(double x, double y, double t)
	{
		if (test == 1) return x+y+1;
		if (test == 2) return x*x + y*y + 1;
		if (test == 3) return x*x*x + y*y*y + 1;
		if (test == 4) return pow(x, 4) + pow(y, 4) + 1;

		if (test == 5) return t;
		if (test == 6) return t*t;
		if (test == 7) return t*t*t;

		if (test == 8) return x*x + y*y + 1;

		if (test == 9) return x*x*x + y*y*y + 1;

	}
	double Tests::Betta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//����� �����
			case 0: return 0;
				//������ �����
			case 1: return 0;
				//������ ����� 
			case 2: return 0;
				//������� �����
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//����� �����
			case 0: return 1;
				//������ �����
			case 1: return 1;
				//������ ����� 
			case 2: return 1;
				//������� �����
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//����� �����
			case 0: return 1;
				//������ �����
			case 1: return 1;
				//������ ����� 
			case 2: return 1;
				//������� �����
			case 3: return 1;
			}
	}
	double Tests::Ubetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//����� �����
			case 0: return 0;
				//������ �����
			case 1: return 0;
				//������ ����� 
			case 2: return 0;
				//������� �����
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//����� �����
			case 0: return x + y - 1;
				//������ �����
			case 1: return x + y + 1;
				//������ ����� 
			case 2: return x + y - 1;
				//������� �����
			case 3: return x + y + 1;
			}
		if (test == 3)
			switch (formNumber)
			{
				//����� �����
			case 0: return 0;
				//������ �����
			case 1: return 2*exp(x + y);
				//������ ����� 
			case 2: return 0;
				//������� �����
			case 3: return 2 * exp(x + y);
			}
	}
	double Tests::Tetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
			//����� �����
			case 0: return 0;
			//������ �����
			case 1: return 0;
			//������ ����� 
			case 2: return 0;
			//������� �����
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//����� �����
			case 0: return -1;
				//������ �����
			case 1: return 1;
				//������ ����� 
			case 2: return -1;
				//������� �����
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//����� �����
			case 0: return -exp(x+y);
				//������ �����
			case 1: return exp(x + y);
				//������ ����� 
			case 2: return -exp(x + y);
				//������� �����
			case 3: return exp(x + y);
			}
	}
	double Tests::Fi(double u, double x, double y, double t)
	{
		if (test == 1) return u - x - y - 1;
		if (test == 2) return u - (x*x + y*y + 1) - 4;
		if (test == 3) return u - (x*x*x + y*y*y + 1) - 6 * x - 6 * y;
		if (test == 4) return u - (pow(x,4) + pow(y,4) + 1) - 12 * x*x - 12 * y*y;

		if (test == 5) return u - t + 1;
		if (test == 6) return u - t*t + 2*t;
		if (test == 7) return u - t*t*t + 3*t*t;

		if (test == 8) return u*u - pow(x*x + y*y + 1,2) - 4;

		if (test == 9) return u*u - pow(x*x*x + y*y*y + 1, 2) - 6*x -6*y;

		return 0;
	}
	double Tests::dFdq(double u, double x, double y, double t)
	{
		if (test == 1) return 1;
		if (test == 2) return 1;
		if (test == 3) return 1;
		if (test == 4) return 1;

		if (test == 5) return 1;
		if (test == 6) return 1;
		if (test == 7) return 1;

		if (test == 8) return 2 * u;

		if (test == 9) return 2 * u*u;
		return 0;
	}
	Tests::Tests()
	{
		ifstream fo;
		fo.open("number_of_test.txt");
		fo >> test;
	}
}