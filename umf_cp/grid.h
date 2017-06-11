/*grid.h*/
#pragma once
#include <stdio.h>
#include <vector>
#include <fstream>
#include <functional>
#include <array>
using namespace std;

namespace grid
{
	//��������� ��� ��������� �������: 
	//�����,������, �����, ������
	//�� - �������� �������
	//�.�. - �������� �������

	//�������
	struct Area
	{
		//x-���������� ����� �������
		int leftX;
		//x-���������� ������ �������
		int rightX;
		//�-���������� ������ �������
		int lowY;
		//�-���������� ������� �������
		int upY;
		//1 - 1 ����
		//2 - 2 ����
		//3 - 3 ����
		//-1 - ��� ������� �������
		array<int, 4> ku;
		//C����� ������� ������� �������
		//-1 - ��� ������� �������
		array<int, 4> kuForm;
		void Input(FILE *fo);
	};

	//������������ ����� ��������
	struct AreasLines
	{
		//x,y - ������� �������� ���� ������������ 
		//� �������������� ������ ����������� Wi

		vector <double> x;
		vector <double> y;
		AreasLines();
	};

	//�������� �������
	struct Element
	{
		//����
		array<int, 4> nodes;
		//������� �������
		array<int, 9> dof;
		//����� �������
		int numberOfArea;
		//�������� ��������
		array<int, 4> neighbors;

		Element& operator=(Element element)
		{
			for (int i = 0; i < 4; i++)
				nodes[i] = element.nodes[i];
			for (int i = 0; i < 9; i++)
				dof[i] = element.dof[i];
			numberOfArea = element.numberOfArea;
			for (int i = 0; i < 4; i++)
				neighbors[i] = element.neighbors[i];
			return *this;
		}
		//����� ����������� ������ �.�. dofsNumber � ��������
		bool SearchDof(int dofsNumber);
	};

	//������� ������� ��� ��������� ��������
	struct BoundaryCondition
	{
		//����� ��������
		int elem;
		//������� �������(�� �����): 1 - ����, 0 - ���
		array<int, 4> edges;
		//����� ������� ��� �������
		array<int, 4> formNumber;
	};

	struct Point
	{
		double x;
		double y;
		Point() {};
		~Point() {};
		Point(double xx, double yy)
		{
			x = xx;
			y = yy;
		}

		bool operator==(Point point)
		{
			if (point.x == x && point.y == y)
				return true;
			else
				return false;
		}

	};

	//��������� ����� �������
	class Grid
	{
		//��������� ���������� � ������ ��������� �� ���� ����������
		void Partition�oordinate(vector <double> &ci, vector <double> areasLines,
			vector <double> coefficient, vector <int> nIntervals);
		void CreateGridTime();
		//���������� ����
		void PushNode(double x, double y);
		//���������� �����
		void BuildGrid();

		//��������� ����������� ������
		int GetGlobalNumber(int elementNumber, int localNumber);
		
		//���������� ������ �������, � ������� ����� �������� �����
		int FindArea(double x, double y);
		//���������� �������� �������� ���������
		void FindNeighbors(int elementNumber);
		//���������� �������� ���������
		void ComputeElements();

		//������������ �������� ������� �������
		void FormKU();

	public:
		Grid();
		//������ ��������
		vector <Area> areas;
		//������������ �����
		AreasLines areasLines;
		//������ �������� ���������
		vector <Element> elements;
		//������ �����
		vector <Point> nodes;
		vector <double> time;
		//������ ������� �������
		array<vector <BoundaryCondition>, 3> ku;
		//����� ���������� ���������� �� x � �� y
		//����� ���������� �����
		int nx, ny;
		//��������� ����������� ������ ��� �������� �������
		int GetGlobalFuncNumber(int elNum, int localFuncNum);

		//����������� ��������� ��
		void DoPartition();
		//������������ ������ ���������, ���������� ���������� ����� �.�.
		//������ dofsNumber
		void SearchElements(int dofsNumber, vector<int>& elList);

		~Grid();
	};

}