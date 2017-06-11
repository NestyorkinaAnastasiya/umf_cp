/*basis.h*/
#pragma once
#include "grid.h"

using namespace grid;
namespace basis
{
	const int nFunc = 9;
	const int nFunc1D = 3;

	struct Basis
	{
		//��������� �� ������� ���������� �������� ������� � �����
		array<function<double(double, double)>, nFunc> phi;
		//��������� �� ������� ���������� d/dksi �������� ������� � �����
		array<function<double(double, double)>, nFunc> dphiksi;
		//��������� �� ������� ���������� d/detta �������� ������� � �����
		array<function<double(double, double)>, nFunc> dphietta;
		Basis();
	};
}