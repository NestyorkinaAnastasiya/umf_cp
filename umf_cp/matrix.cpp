/*matrix.cpp*/
#include "matrix.h"

namespace matrix
{
	//��������� �������� �������
	void Matrix::CreatePortret(int slaeSize, Grid grid)
	{
		vector <int> elList, unzeroNumbersList;

		vector<vector <int>> list;
		list.reserve(slaeSize);

		for (int i = 0; i < slaeSize; i++)
		{
			//�����,� ����� ��������� ������������ ������ i-� ������� �������
			grid.SearchElements(i, elList);
			//������� ������ ����� ���� ���������, ������� �� ������
			for (unsigned int j = 0; j < elList.size(); j++)
			{
				for (int k = 0; k < 9; k++)
				{
					//���� ������ ������ ��� ��� � �� ������ �� ������,�� ���������
					if (find(unzeroNumbersList.begin(), unzeroNumbersList.end(), grid.elements[elList[j]].dof[k])
						== unzeroNumbersList.end() && grid.elements[elList[j]].dof[k] < i)
						unzeroNumbersList.push_back(grid.elements[elList[j]].dof[k]);
				}
			}

			sort(unzeroNumbersList.begin(), unzeroNumbersList.end());
			list.push_back(unzeroNumbersList);
			unzeroNumbersList.clear();
			elList.clear();
		}

		//��������� ����������� ggl,ggu
		int gg_size = 0;
		for (int i = 0; i < slaeSize; i++)
		{
			if (!list[i].empty())
				gg_size += list[i].size();
		}

		//�������������� ������� � �������� �������
		Initialize(slaeSize, gg_size);

		ig[0] = 0;

		for (int i = 0; i < n; i++)
		{
			if (!list[i].empty())
				ig[i + 1] = ig[i] + list[i].size();
			else
				ig[i + 1] = ig[i];
		}

		int k = 0;
		for (int i = 0; i < n; i++)
		{
			if (!list[i].empty())
			{
				for (unsigned int j = 0; j < list[i].size(); j++)
				{
					jg[k] = list[i][j];
					k++;
				}
				list[i].clear();
			}
		}

		list.clear();
	}

	//������������� ������� ����� ��������� ��������
	void Matrix::Initialize(int size1, int size2)
	{
		n = size1; size = size2;

		ggl.resize(size);
		ggu.resize(size);
		di.resize(n);
		ig.resize(n + 1);
		jg.resize(size);
	}

	//��������� ������� �� ������
	void Matrix::MultiplyAx(const vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;

		for (i = 0; i < n; i++)
		{
			//������ i-�� ������(�������)
			l = ig[i];
			//������ (i+1)-�� ������(�������)
			iend = ig[i + 1];
			//���������� ��������� � i ������(�������)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//�������� �� ���� ��������� i ������ (�������)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggl[l] * a[j];
				result[j] += ggu[l] * a[i];
			}
		}
	}

	//��������� ����������������� ������� �� ������
	void Matrix::MultiplyATx(vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;
		for (i = 0; i < n; i++)
		{
			//������ i-�� ������(�������)
			l = ig[i];
			//������ (i+1)-�� ������(�������)
			iend = ig[i + 1];
			//���������� ��������� � i ������(�������)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//�������� �� ���� ��������� i ������ (�������)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggu[l] * a[j];
				result[j] += ggl[l] * a[i];
			}
		}
	}

}

