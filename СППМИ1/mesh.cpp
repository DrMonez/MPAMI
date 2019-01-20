#include "1.h"

void mesh::Input()
{
	ifstream file("in.txt");
	string s;
	vector <double> array(6);

	/* получаем длину, число интервалов и коэф разрядки */
	getline(file, s);
	istringstream iss(s);
	copy_n(istream_iterator < double >(iss), 3, array.begin());
	iss.clear();

	edgeX.length = array[0];
	edgeX.nIntervals = (int)array[1];
	edgeX.coeff = array[2];

	/* получаем глубину, число интервалов и коэф разрядки */
	getline(file, s);
	iss.str(s);
	copy_n(istream_iterator < double >(iss), 3, array.begin());
	iss.clear();

	edgeY.length = array[0];
	edgeY.nIntervals = (int)array[1];
	edgeY.coeff = array[2];

	/* получаем проводимость */
	getline(file, s);
	sigma = atof(s.c_str());

	/* получаем координаты точек A и B */
	getline(file, s);
	iss.str(s);
	copy_n(istream_iterator < double >(iss), 6, array.begin());
	iss.clear();

	for (unsigned int i = 0; i < 2; i++) {
		sources[i].x = array[3 * i];
		sources[i].y = array[3 * i + 1];
		sources[i].P = array[3 * i + 2];
	}

	/* получаем координаты приёмников */
	getline(file, s);
	iss.str(s);
	copy_n(istream_iterator < double >(iss), 6, array.begin());

	for (unsigned int i = 0; i < 3; i++) {
		receivers[i].x = array[2 * i];
		receivers[i].y = array[2 * i + 1];
	}
	file.close();
}

void mesh::FragmentationOfEdge(edge &edges)
{
	double h;

	if (edges.coeff == 1)
		h = edges.length / (edges.nIntervals);
	else if (edges.coeff == 0)
		h = edges.length;
	else
		if (edges.coeff > 1)
			h = edges.length*(edges.coeff - 1) / (edges.coeff*(pow(edges.coeff, edges.nIntervals - 1) - 1) + edges.coeff - 1);
		else
			h = edges.length*(1 - edges.coeff) / (1 - pow(edges.coeff, edges.nIntervals));

	double x = 0, y = 0;

	edges.points.push_back(x);

	for (int j = 1; j <= edges.nIntervals; j++) {
		x += h;
		h *= edges.coeff;
		edges.points.push_back(x);
	}
}

void mesh::GenerateNodes()
{
	unsigned int nNodes = edgeX.points.size() * edgeY.points.size();
	nodes.resize(nNodes);
	unsigned int k = 0;

	/* записываем координаты узлов */
	for (unsigned int j = 0; j < edgeY.points.size(); j++)
		for (unsigned int i = 0; i < edgeX.points.size(); i++) {
			nodes[k].x = edgeX.points[i];
			nodes[k].y = edgeY.points[j];
			k++;
		}
}

void mesh::GetBound1Nodes()
{
	/* правая граница */
	for (unsigned int j = 1; j <= edgeY.points.size(); j++) {
		unsigned int num = edgeX.points.size() * j - 1;
		nodes[num].bound1 = true;
	}
	/* нижняя граница */
	for (unsigned int i = 0; i < edgeX.points.size(); i++) {
		unsigned int num = i;
		nodes[num].bound1 = true;
	}
}

void mesh::GenerateParametersOfElement()
{
	unsigned int nElements = (edgeY.points.size() - 1) * (edgeX.points.size() - 1);
	elements.resize(nElements);
	unsigned int k = 0;
	for (unsigned int j = 0; j < edgeY.points.size() - 1; j++)
		for (unsigned int i = 0; i < edgeX.points.size() - 1; i++) {
			/* определяем глобальные номера узлов и hx, hy каждого элемента */
			elements[k].globalNodes[0] = edgeX.points.size() * j + i;
			elements[k].globalNodes[1] = edgeX.points.size() * j + i + 1;
			elements[k].globalNodes[2] = edgeX.points.size() * (j + 1) + i;
			elements[k].globalNodes[3] = edgeX.points.size() * (j + 1) + i + 1;

			elements[k].hx = nodes[elements[k].globalNodes[1]].x - nodes[elements[k].globalNodes[0]].x;
			elements[k].hy = nodes[elements[k].globalNodes[2]].y - nodes[elements[k].globalNodes[0]].y;
			k++;
		}
}

void mesh::BuildMesh()
{
	/* разбиваем грани */
	FragmentationOfEdge(edgeY);
	FragmentationOfEdge(edgeX);

	/* генерируем узлы, их координаты */
	GenerateNodes();

	/* устанавливаем первые краевые = 0, там где, они есть */
	GetBound1Nodes();

	/* генерируем массив элементов, хранящих свои узлы и hx, hy */
	GenerateParametersOfElement();

}