/*slae.h*/
#pragma once
#include "tests.h"
using namespace matrix;
using namespace basis;
using namespace integration;
using namespace tests;

namespace slae
{
	class SLAE : private Basis, private GaussIntegration
	{
		//����������� ������
		int n;
		//������������ ���������� �������� � ��������
		int maxiter = 10000;
		//�������� ������� ����
		const double eps = 1e-10;
		//�����
		Grid grid;
		//��������� �������� �������
		Tests tests;
		//���������� �������
		Matrix A;
		Matrix globalM1;
		Matrix globalM2;
		//��������� �������
		//������� ��������
		array<array<double, 9>, 9> G;
		//������� �����
		array<array<double, 9>, 9> M;
		array<array<double, 9>, 9> Mg;
		//��������� ������ ������ �����
		array <double, 9> locF;
		//���������� ������ ������ �����
		vector <double> F;
		array <double, 9> newF;
		//������ ������������� ������� �� ���������� ��������
		//�� ������������
		vector <double> u_n;
		//������� ������������� ������� �� ���������� ���������
		//�� �������
		vector <double> u_t1;
		vector <double> u_t2;
		//������ ������������� �������
		vector <double> u;
		//����� ������� ������ �����
		double normF;

		//������ ��������� ������ ��������
		void CalculateG(int elementNumber);
		//������ ��������� ������ ����
		void CalculateM(int elementNumber);
		void CalculateMg(int elementNumber);
		//������ ��������� ������ ������
		void CalculateLocalF(int elementNumber);
		//������� ���������� �������� � ����������
		void AddElementToGlobalMatrix(Matrix &B, int i, int j, double element);
		//������ ��������� ������(��������) � ���������� � ����������
		void CalculateLocals(int elementNumber);


		//������ ������� ����� ��� ������� �������� 
		array<double, 3> g;
		//���������� ������ ����� ��� 1��� �������� �������
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//���������� 1��� �������� ������� ��� ������ ����
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//���� ������� �������� �������
		void CalculateBoundaries1(int number);

		//���������� ������� � �������������
		vector <double> L;
		vector <double> D;
		vector <double> U;

		//������� �������� �������
		double t;
		double deltat;
		double deltat0;
		double deltat1;
		//������ �������
		vector <double> r;
		//������ ������
		vector <double> z;

		//���������� ����� �������
		double Norm(const vector<double>& x);
		//��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);

		//��������� ���� �� i-�� �������� �� �������
		void GenerateSLAE();
		//LU-������������
		void LU();
		//��������������� ������� ��� ��������
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
		double Rel_Discrepancy();
		//�������� ��� � LU-�������������
		void LULOS();

		double StopIteration();

		void CalculateTimeVector(double time, vector<double> &result);

	public:
		SLAE();

		void TSolve();

		~SLAE() {};
	};
}
