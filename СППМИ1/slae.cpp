#include "1.h"

void slae::CreatePortrait()
{
	int ggSize = 0;
	/* векторы для хранения связей узлов */
	vector<vector <int>> nodesLinks;
	nodesLinks.resize(nodes.size());

	/* идем по всем элементам */
	for (unsigned int i = 0; i < elements.size(); i++)
		/* идем по всем узлам элемента */
		for (int j = 0; j < 4; j++)
			/* для предыдущих узлов текущего элемента устанавливаем связь с текущим узлом */
			for (int k = 0; k < j; k++)
				nodesLinks[elements[i].globalNodes[k]].push_back(elements[i].globalNodes[j]);

	/* сортируем векторы связей и устраняем дублирование в них */
	for(vector <int> nodeLinks : nodesLinks)
	{
		sort(nodeLinks.begin(), nodeLinks.end());
		nodeLinks.erase(unique(nodeLinks.begin(), nodeLinks.end()), nodeLinks.end());
		ggSize += nodeLinks.size();
	}
	
	int n = nodes.size();
	ig.resize(n + 1);
	di.resize(n);
	f.resize(n);
	x.resize(n);
	jg.resize(ggSize);
	gg.resize(ggSize);

	unsigned int jgPosition = 0;
	for (unsigned int i = 0; i < n; i++)
	{
		unsigned int k = 0;
		for (unsigned int j = 0; j <= i; j++)
			/* узнаем, есть ли связь между i-ым и j-ым элементами */
			if (find(nodesLinks[j].begin(), nodesLinks[j].end(), i) != nodesLinks[j].end())
			{
			jg[jgPosition] = j;
			jgPosition++;
			k++;
			}
			/* количество ненулевых элементов в строке */
			ig[i + 1] = ig[i] + k;
	}
	nodesLinks.clear();

}

void slae::AddElementToGlobalMatrix(int i, int j, double value)
{
	if (i == j) di[i] += value;
	else
	{
		for (int id = ig[i]; id < ig[i + 1]; id++)
			if (jg[id] == j)
			{
				gg[id] += value;
				return;
			}
	}
}

void slae::CreateGlobalMatrix()
{
	int sizeGg = gg.size();
	int size = di.size();
	gg.clear();
	di.clear();
	f.clear();
	x.clear();
	gg.resize(sizeGg);
	di.resize(size);
	f.resize(size);
	x.resize(size);

	for (int i = 0; i < elements.size(); i++)
	{
		CreateLocalMatrix(i);
		for (int k = 0; k < 4; k++)
			for (int l = 0; l <= k; l++)
				AddElementToGlobalMatrix(elements[i].globalNodes[k], elements[i].globalNodes[l], G[k][l]);
	}
}

void slae::CreateLocalMatrix(int elementNumber)
{
	double hx = elements[elementNumber].hx;
	double hy = elements[elementNumber].hy;
	double jacobian = hx * hy / 4.0;

	for (int i = 0; i < 4; i++)
		for (int j = i; j < 4; j++)
		{
			double result = 0;
			for (int k = 0; k < 9; k++)
			{
				double ksi = 0.5 + 0.5 * gaussPoints[0][k],
					etta = 0.5 + 0.5 * gaussPoints[1][k];
				double x_ = nodes[elements[elementNumber].globalNodes[0]].x + ksi * elements[elementNumber].hx,
					y_ = nodes[elements[elementNumber].globalNodes[0]].y + etta * elements[elementNumber].hy;
				result += x_ * gaussWeights[k] * sigma * (dpsidx[i](etta, hx) * dpsidx[j](etta, hx) + dpsidy[i](ksi, hy) * dpsidy[j](ksi, hy));
			}

			G[i][j] = result * jacobian;
		}

	/* матрица симметричная, заполняем нижний треугольник */
	for (int i = 1; i < 4; i++)
		for (int j = 0; j < i; j++)
			G[i][j] = G[j][i];
}

int slae::GetElement(double x, double y)
{
	for (int i = 0; i < elements.size() || !exit; i++)
		if (x <= nodes[elements[i].globalNodes[1]].x && y <= nodes[elements[i].globalNodes[2]].y)
			if (x >= nodes[elements[i].globalNodes[0]].x && y >= nodes[elements[i].globalNodes[0]].y)
				return i;
	return -1;
}

void slae::AddSourcesToF()
{
	point tmp = { 1e-5,  edgeY.length };
	double ksi, etta;

	int i = GetElement(tmp.x, tmp.y);
	if (i == -1) {
		cout << "ERROR!!!";
		return;
	}
	ksi = (tmp.x - nodes[elements[i].globalNodes[0]].x) / elements[i].hx;
	etta = (tmp.y - nodes[elements[i].globalNodes[0]].y) / elements[i].hy;
	
	/* добавка к правой части */
	for (int j = 0; j < 4; j++)
		f[elements[i].globalNodes[j]] += sources[0].P * psi[j](ksi, etta);
}

void slae::AddBoundaries()
{
	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].bound1)
		{
			di[i] = 1.0;
			f[i] = 0.0;
			for (int k = ig[i]; k < ig[i + 1]; k++)
				gg[k] = 0.0;
			for (int l = i + 1; l < nodes.size(); l++)
				for (int p = ig[l]; p < ig[l + 1]; p++)
					if (jg[p] == i)
						gg[p] = 0.0;
		}
	}
}

void slae::SolveDT()
{
	CreateGlobalMatrix();
	AddSourcesToF();
	AddBoundaries();
	LULOS();
}

void slae::Output()
{
	ofstream file("out.txt");
	for (int i = 0; i < nodes.size(); i++)
		file << x[i] << endl;
	file.close();
}

/* генерируем эксперименттальные значения потенциала V, решая прямую задачу */
slae::slae()
{
	Input();
	BuildMesh();
	CreatePortrait();
	SolveDT();
}

slae::~slae()
{
	Output();
}



