#include "1.h"

/* скалярное произведение */
double Scalyar(vector <double> &a, vector <double> &b)
{
	double Res = 0.0;
	int o = a.size();
	for (int i = 0; i < o; i++)
		Res += a[i] * b[i];
	return Res;
}

/* умножение матрицы на вектор */
void slae::Multiply(vector <double> &a, vector <double> &res)
{
	int i0, i1, kol;
	for (int i = 0; i < nodes.size(); i++)
	{
		//начало i-ой строки(столбца)
		i0 = ig[i];
		//начало (i+1)-ой строки(столбца)
		i1 = ig[i + 1];
		//количество элементов в i строке(столбце)
		kol = i1 - i0;
		res[i] = di[i] * a[i];
		//проходим по всем элементам i строки (столбца)
		for (int k = 0; k < kol; k++, i0++)
		{
			int j = jg[i0];
			res[i] += gg[i0] * a[j];
			res[j] += gg[i0] * a[i];
		}
	}
}

void slae::LU()
{
	int i, i0, j0, iend, num, ki, kj;
	double suml, sumu, sumdg;
	L.resize(gg.size());
	L = gg;
	U.resize(gg.size());
	U = gg;
	Di.resize(di.size());
	Di = di;
	int n = nodes.size();
	for (i = 0; i < n; i++)
	{
		i0 = ig[i];
		iend = ig[i + 1];
		for (num = i0, sumdg = 0; num < iend; num++)
		{
			j0 = ig[jg[num]]; //в зависимости от номера фиксируем столбец,какой столбец l,такого столбца  ищем начальный эл у u 
			int jend = ig[jg[num] + 1];
			ki = i0; kj = j0;
			for (suml = 0, sumu = 0, ki = i0; ki < num; ki++) //для num учитываются все предыдущие элементы
				for (int m = kj; m < jend; m++)
					if (jg[ki] == jg[m]) //ищем соответствующие ненулевые элементы для умножения
					{
						suml += L[ki] * U[m];
						sumu += L[m] * U[ki];//для симметричного элемента из U
					}
			L[num] -= suml;
			U[num] = (U[num] - sumu) / Di[jg[num]];
			sumdg += L[num] * U[num];//умножаются симметричные элементы	
		}
		Di[i] -= sumdg;
	}
}

void slae::LYF(const vector <double> &C, vector <double> &yl)
{
	int i, i0, iend; //i0-адрес начала строки, iend-адрес конца строки
	double sum;
	for (i = 0; i < nodes.size(); i++)
	{
		i0 = ig[i]; iend = ig[i + 1];

		yl[i] = C[i];

		for (i0, sum = 0; i0 < iend; i0++)
			yl[i] -= yl[jg[i0]] * L[i0];
		yl[i] /= Di[i];
	}
}

void slae::UXY(const vector <double> &C, vector <double> &yu)
{
	int i, i0, iend;

	for (i = 0; i < nodes.size(); i++)
		yu[i] = 0.0;

	for (i = nodes.size() - 1; i >= 0; i--)//проход по столбцам с конца
	{
		yu[i] += C[i];

		i0 = ig[i]; iend = ig[i + 1]; iend--;

		for (; iend >= i0; iend--)//идём по столбцу с конца
			yu[jg[iend]] -= yu[i] * U[iend];
	}
}

void slae::LULOS()
{
	double a, b, pp, dis;
	int i;
	int n = nodes.size();
	vector <double> Ax(n);
	vector <double> C(n);
	r.resize(n);
	z.resize(n);
	p.resize(n);

	int maxiter = 10000;
	LU();
	Multiply(x, Ax);					//Ax0
	for (i = 0; i < n; i++)				//f-Ax0
		r[i] = f[i] - Ax[i];
	LYF(r, r);							//r0=L^(-1)(f-Ax0)
	UXY(r, z);							//z0=U^(-1)r0->r0=Uz0
	//p0=L^(-1)Az0
	Multiply(z, Ax);					//Az0
	LYF(Ax, p);

	rr = Scalyar(r, r);
	dis = Scalyar(r, r) / rr;
	dis = sqrt(dis);

	for (int iter = 1; dis > 1e-14 && iter <= maxiter; iter++)
	{
		pp = Scalyar(p, p);				//Аk
		a = Scalyar(p, r) / pp;

		for (i = 0; i < n; i++)
		{
			x[i] += a * z[i];			//Xk, Rk
			r[i] -= a * p[i];
		}
		UXY(r, C);						//UY=rk->Y=U^(-1)rk
		Multiply(C, Ax);				//AU^(-1)rk=Ax
		LYF(Ax, Ax);					//L^(-1)AU^(-1)rk=Y2->L^(-1)B=Y2->LY2=B->Y2=L^(-1)AU^(-1)rk
		b = -Scalyar(p, Ax) / pp;		//bk
		for (i = 0; i < n; i++)
		{
			z[i] = C[i] + b * z[i];		//zk=U^(-1)rk+bkz[k-1]				
			p[i] = Ax[i] + b * p[i];	//pk
		}
		dis = Scalyar(r, r) / rr;
		dis = sqrt(dis);
	}
}

