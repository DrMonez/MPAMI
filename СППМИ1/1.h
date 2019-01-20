#pragma once
#include <conio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <functional>
#include <math.h>
#include <set>

/* ���������� ����������� ������������ ���� */
using namespace std; 

struct point
{
	double x;
	double y;
};

struct node: point
{
	bool bound1 = false;
};

struct edge
{
	double length;
	int nIntervals;
	double coeff;
	vector <double> points;
};

/* �������� ��������� */
struct pointSources : point
{ 
	double P;
};

struct element
{
	/* ���������� �� � � � ��� �������� ������� */
	double hx, hy;
	/* ������� ���������� ������� ����� ��������� �������� � ������� ����� */
	array <int, 4> globalNodes;
	/* ����� �������� */
	int number;
};

struct pointReceivers : point
{
	double rA;
	double rB;
	double z;
	double signal;
	array <double, 2> values;
	array <element, 2> elements;
};

class basis
{
public:
	//������� ���������� �������� �� � ����� ������-���������, �� �����������
	//��������� �� ������� ���������� �������� ������� � �����
	array <function <double(double, double)>, 4> psi;
	//��������� �� ������� ���������� d/dksi �������� ������� � �����
	array <function <double(double, double)>, 4> dpsidx;
	//��������� �� ������� ���������� d/detta �������� ������� � �����
	array <function <double(double, double)>, 4> dpsidy;

	basis();
};

class mesh
{
protected:
	/* ����� � ������ ����� / ������ � ������� ����� */
	edge edgeY, edgeX;
	double sigma;
	array <pointSources, 2> sources; 
	array <pointReceivers, 3> receivers;
	vector <element> elements;
	/* ����� ������ �����: ������ - ���� ����� ���� */
	vector <node> nodes;

public:
	void Input();
	void BuildMesh();
	void FragmentationOfEdge(edge &edges);
	void GenerateParametersOfElement();
	void GetBound1Nodes();
	void GenerateNodes();
};

class gaussIntegration
{
public:
	array <array<double, 9>, 2> gaussPoints;
	array <double, 9> gaussWeights;
	array <double, 3> gaussPoints1;
	array <double, 3> gaussWeights1;

	gaussIntegration();
};

class slae : public mesh, protected basis, protected gaussIntegration
{
protected:
	vector <double> gg, ig, jg, di, f, x;
	array <array <double, 4>, 4> G;

	/* ���������� ������� � ������������� */
	vector <double> L, Di, U;

	/* ������� �������, ������ */
	vector <double> r, z, p;
	double pp, rr;

	void CreatePortrait();
	void CreateGlobalMatrix();
	void CreateLocalMatrix(int elementNumber);
	void AddElementToGlobalMatrix(int i, int g, double value);
	void AddSourcesToF();
	void AddBoundaries();
	void SolveDT();

	void LULOS();
	void LU();
	void LYF(const vector <double> &C, vector <double> &yl);
	void UXY(const vector <double> &C, vector <double> &yu);
	void Multiply(vector <double> &a, vector <double> &res);

	void Output();

	int GetElement(double x, double y);

	friend class sigma;
	friend class power;
	template<typename testParameter>
	friend class inverseTask;
public:

	slae();
	~slae();
};


