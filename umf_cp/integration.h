/*integration.h*/
#pragma once
#include "matrix.h"
namespace integration
{
	const int n_ip = 25;
	const int n_ip1D = 5;
	struct GaussIntegration
	{
		//����� ������
		array<array<double, n_ip>, 2> gaussPoints;
		//���� ������
		array<double, n_ip> gaussWeights;
		//����� ������
		array<double, n_ip1D> gaussPoints1;
		//���� ������
		array<double, n_ip1D> gaussWeights1;

		GaussIntegration();
	};
}