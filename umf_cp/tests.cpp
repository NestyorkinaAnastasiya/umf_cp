/*tests.cpp*/
#include "tests.h"

namespace tests
{
	double Tests::Lambda(double x, double y, double t)
	{
		switch (test)
		{// По пространству
		case 1: return 1;
		case 2: return 2;
		case 3: return 1;
		case 4: return 1;
		// По времени
		case 5: return 1;
		case 6: return 1;
		case 7: return 1;
		}
		return 1;
	}
	double Tests::Sigma(double x, double y, double t)
	{
		switch (test)
		{
		case 1: return 1;
		case 2: return 1;
		case 3: return 1;
		case 4: return 1;
		// По времени
		case 5: return 1;
		case 6: return 1;
		case 7: return 1;
		}
		return 1;
	}
	double Tests::Gamma(double x, double y, double t)
	{
		switch (test)
		{// По пространству
		case 1: return 1;
		case 2: return 3;
		case 3: return 1;
		case 4: return 4;
		// По времени
		case 5: return 1;
		case 6: return 1;
		case 7: return 1;
		}
		return 1;
	}
	double Tests::Ug(double x, double y, double t)
	{
		switch (test)
		{// По пространству
		case 1: return Gamma(x, y, t)*(2 * x - 3 * y + 1);
		case 2: return 2 * pow(x, 2) - 3 * pow(y, 2) + 1;
		case 3: return 2 * pow(x, 3) - 3 * pow(y, 3) + 1;
		case 4: return 2 * pow(x, 4) - 3 * pow(y, 4) + 1;
		// По времени
		case 5: return 2*t;
		case 6: return 2*pow(t,2);
		case 7: return 2 * pow(t, 3);
		}
		return 0;
	}
	double Tests::Betta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 0;
				//правое ребро
			case 1: return 0;
				//нижнее ребро 
			case 2: return 0;
				//верхнее ребро
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 1;
				//правое ребро
			case 1: return 1;
				//нижнее ребро 
			case 2: return 1;
				//верхнее ребро
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 1;
				//правое ребро
			case 1: return 1;
				//нижнее ребро 
			case 2: return 1;
				//верхнее ребро
			case 3: return 1;
			}
	}
	double Tests::Ubetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 0;
				//правое ребро
			case 1: return 0;
				//нижнее ребро 
			case 2: return 0;
				//верхнее ребро
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//левое ребро
			case 0: return x + y - 1;
				//правое ребро
			case 1: return x + y + 1;
				//нижнее ребро 
			case 2: return x + y - 1;
				//верхнее ребро
			case 3: return x + y + 1;
			}
		if (test == 3)
			switch (formNumber)
			{
				//левое ребро
			case 0: return 0;
				//правое ребро
			case 1: return 2*exp(x + y);
				//нижнее ребро 
			case 2: return 0;
				//верхнее ребро
			case 3: return 2 * exp(x + y);
			}
	}
	double Tests::Tetta(int formNumber, double x, double y)
	{
		if (test == 1)
			switch (formNumber)
			{
			//левое ребро
			case 0: return 0;
			//правое ребро
			case 1: return 0;
			//нижнее ребро 
			case 2: return 0;
			//верхнее ребро
			case 3: return 0;
			}

		if (test == 2)
			switch (formNumber)
			{
				//левое ребро
			case 0: return -1;
				//правое ребро
			case 1: return 1;
				//нижнее ребро 
			case 2: return -1;
				//верхнее ребро
			case 3: return 1;
			}

		if (test == 3)
			switch (formNumber)
			{
				//левое ребро
			case 0: return -exp(x+y);
				//правое ребро
			case 1: return exp(x + y);
				//нижнее ребро 
			case 2: return -exp(x + y);
				//верхнее ребро
			case 3: return exp(x + y);
			}
	}
	double Tests::Fi(double x, double y, double t)
	{
		switch (test)
		{// По пространству
		case 1: return Gamma(x, y, t)*(2 * x - 3 * y + 1);
		case 2: return 2 * Lambda(x, y, t) + 
			Gamma(x, y, t)*(2 * pow(x, 2) - 3 * pow(y, 2) + 1);
		case 3: return -Lambda(x, y, t)*(12 * x - 18 * y) + 
			Gamma(x, y, t)*(2 * pow(x, 3) - 3 * pow(y, 3) + 1);
		case 4: return -Lambda(x, y, t)*(24 * pow(x, 2) - 36 * pow(y, 2)) +
			Gamma(x, y, t)*(2 * pow(x, 4) - 3 * pow(y, 4) + 1);
		// По времени
		case 5: return 2 * Sigma(x, y, t) + Gamma(x, y, t) * 2 * t;
		case 6: return 4 * Sigma(x, y, t)*t + Gamma(x, y, t) * 2 * pow(t, 2);
		case 7: return 6 * Sigma(x, y, t)*pow(t,2) + Gamma(x, y, t) * 2 * pow(t, 3);
		}
		return 0;
	}
	Tests::Tests()
	{
		ifstream fo;
		fo.open("number_of_test.txt");
		fo >> test;
	}
}