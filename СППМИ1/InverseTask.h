#pragma once
#include "1.h"

template <typename testParameter>
class inverseTask
{
protected:
	/* число параметров среды */
	const static int nParameters = 1;
	double delta = 0.05;
	double gamma = 1e-3;
	double eps = 1e-10;
	/* вектор параметров среды */
	vector<testParameter> parameters;
	slae slae;
	/* матрица СЛАУ для обратной задачи */
	array <array <double, nParameters>, nParameters> A;
	/* вектор правой части для обратной задачи */
	array <double, nParameters> f;
	/* направление поиска delta_u */
	array <double, nParameters> deltaU;
	/* нахождение позиции приёмников относительно каждого электрода */
	double CalculateR(double x, double y, double xAB, double yAB);
	/* создание приёмников и расчёт напряжённости в них */
	void CreateReceivers();
	/* расчёт компонентов напряженности E=-gradV */
	void CalculateE(int i, array<double, 2>& returnValue);
	/* установка средних значений параметров среды */
	void SetParameters();
	/* нахождение функционала Ju */
	double GetJu();
	/* генерация СЛАУ для обратной задачи */
	inline void GetInverseSLAE(array <array<double, 2>, 3> &prevE, array <array <array<double, 2>, 3>, nParameters> &newE);
public:
	inverseTask();
	~inverseTask();
};


struct sigma
{
	double initialSigma = 0.1;
	double GetValue(slae& slae) { return slae.sigma; };
	void SetValue(slae& slae, double sigma) { slae.sigma = sigma; };
	double GetParameter() { return initialSigma; };
};

struct power
{
	double initialPower = 2.5;
	double GetValue(slae& slae) { return slae.sources[0].P; };
	void SetValue(slae& slae, double power) { slae.sources[0].P = slae.sources[1].P = power; };
	double GetParameter() { return initialPower; };
};

template <typename testParameter>
double inverseTask<testParameter>::CalculateR(double x, double y, double xAB, double yAB)
{
	return sqrt(pow((x - xAB), 2) + pow((y - yAB), 2));
}

template <typename testParameter>
void inverseTask<testParameter>::CreateReceivers()
{
	for (int i = 0; i < slae.receivers.size(); i++) {
		/* позиция электрода A в r,z координатах */
		slae.receivers[i].rA = CalculateR(slae.receivers[i].x, slae.receivers[i].y, slae.sources[0].x, slae.sources[0].y);
		/* позиция электрода B в r,z координатах */
		slae.receivers[i].rB = CalculateR(slae.receivers[i].x, slae.receivers[i].y, slae.sources[1].x, slae.sources[1].y);
		slae.receivers[i].z = slae.edgeY.length;

		slae.receivers[i].signal = 1.0;

		/* элементы, в которых находятся координаты приемников, в зависимости от поля источника  */
		int searchIndex = slae.GetElement(slae.receivers[i].rA, slae.receivers[i].z);
		slae.receivers[i].elements[0] = slae.elements[searchIndex];

		searchIndex = slae.GetElement(slae.receivers[i].rB, slae.receivers[i].z);
		slae.receivers[i].elements[1] = slae.elements[searchIndex];

		slae.receivers[i].values[0] = 0;
		slae.receivers[i].values[1] = 0;
		/* находим напряженность Е = -gradV электрического поля */
		/* помним, что электрод A - положительный, а B - отрицательный */
		CalculateE(i, slae.receivers[i].values);
	}

	/* зашумляем данные */
	slae.receivers[0].signal = 1.0;
	double noise = 1.0 - slae.receivers[0].signal;
	slae.receivers[0].values[0] *= (1.0 + noise);
	slae.receivers[0].values[1] *= (1.0 + noise);

	slae.receivers[1].signal = slae.receivers[0].signal;
	slae.receivers[1].values[0] *= (1.0 + noise);
	slae.receivers[1].values[1] *= (1.0 + noise);
}

template <typename testParameter>
void inverseTask<testParameter>::CalculateE(int i, array<double, 2> &returnValue)
{
	int sign = -1;
	double x_ = slae.receivers[i].rA;

	/* -gradV поля от A - (-gradV поля от В) */
	for (int k = 0; k < 2; k++)
	{
		double ksi = (x_ - slae.nodes[slae.receivers[i].elements[k].globalNodes[0]].x) / slae.receivers[i].elements[k].hx;
		double etta = (slae.receivers[i].z - slae.nodes[slae.receivers[i].elements[k].globalNodes[0]].y) / slae.receivers[i].elements[k].hy;

		double g1 = 0, g2 = 0, q;
		for (int j = 0; j < 4; j++)
		{
			q = slae.x[slae.receivers[i].elements[k].globalNodes[j]];
			double a = slae.dpsidx[j](etta, slae.receivers[i].elements[k].hx);
			double b = slae.dpsidy[j](ksi, slae.receivers[i].elements[k].hy);
			g1 += q * slae.dpsidx[j](etta, slae.receivers[i].elements[k].hx);
			g2 += q * slae.dpsidy[j](ksi, slae.receivers[i].elements[k].hy);
		}
		x_ = slae.receivers[i].rB;
		returnValue[0] += sign * g1;
		returnValue[1] += sign * g2;
		sign = 1;
	}
}

template<>
void inverseTask<sigma>::SetParameters()
{
	sigma s;
	/* задаем начальное значение сигма */
	parameters[0].SetValue(slae, s.initialSigma);
}

template<>
void inverseTask<power>::SetParameters()
{
	power p;
	/* задаем начальное значение P */
	parameters[0].SetValue(slae, p.initialPower);
}

template<typename testParameter>
double inverseTask<testParameter>::GetJu()
{
	/* решаем прямую задачу */
	slae.SolveDT();

	/* значение функционала J и веса с целью его минимизации */
	double Ju = 0.0, w;
	array <double, 2> sol;

	for (unsigned int i = 0; i < slae.receivers.size(); i++) {

		sol = { 0, 0 };

		CalculateE(i, sol);

		/* по r */
		w = slae.receivers[i].signal / slae.receivers[i].values[0];
		Ju += pow(w * (slae.receivers[i].values[0] - sol[0]), 2);

		/* по z */
		w = slae.receivers[i].signal / slae.receivers[i].values[1];
		Ju += pow(w * (slae.receivers[i].values[1] - sol[1]), 2);
	}
	return Ju;
}

/* расчёт элементов матрицы и правой части СЛАУ для параметров (без альфы) */
template<typename testParameter>
inline void inverseTask<testParameter>::GetInverseSLAE(array <array<double, 2>, 3> &prevE, array <array <array<double, 2>, 3>, nParameters> &newE)
{
	for (unsigned int i = 0; i < nParameters; i++) {
		f[i] = 0;
		for (unsigned int j = 0; j < nParameters; j++)
			A[i][j] = 0;
	}
	/* w2 - уровень погрешности при приёме сигнала в k-ом источнике */
	double w2 = 0;
	/* prevE - значения E=-gradV на предыдущей итерации */
	/* newE - значения E=-gradV при приращении параметра(ов) */
	for (unsigned int i = 0; i < nParameters; i++) {
		for (unsigned int k = 0; k < slae.receivers.size(); k++) {
			for (unsigned int j = 0; j < nParameters; j++)  {
				w2 = (slae.receivers[k].signal * slae.receivers[k].signal) / (slae.receivers[k].values[0] * slae.receivers[k].values[0]);
				A[i][j] += w2 *
					(prevE[k][0] - newE[i][k][0]) / (parameters[i].GetValue(slae) * delta) *
					(prevE[k][0] - newE[j][k][0]) / (parameters[j].GetValue(slae) * delta);
				w2 = (slae.receivers[k].signal * slae.receivers[k].signal) / (slae.receivers[k].values[1] * slae.receivers[k].values[1]);
				double a = prevE[k][1];
				double b = newE[i][k][1];
				double c = newE[j][k][1];
				double d = parameters[i].GetValue(slae) * delta;
				double e = parameters[j].GetValue(slae) * delta;
				A[i][j] += w2 *
					(prevE[k][1] - newE[i][k][1]) / (parameters[i].GetValue(slae) * delta) *
					(prevE[k][1] - newE[j][k][1]) / (parameters[j].GetValue(slae) * delta);
			}
			w2 = (slae.receivers[k].signal * slae.receivers[k].signal) / (slae.receivers[k].values[0] * slae.receivers[k].values[0]);///ПОМЕНЯТЬ МЕСТАМИ
			f[i] -= w2 *
				(slae.receivers[k].values[0] - prevE[k][0]) *
				(prevE[k][0] - newE[i][k][0]) / (parameters[i].GetValue(slae) * delta);
			w2 = (slae.receivers[k].signal * slae.receivers[k].signal) / (slae.receivers[k].values[1] * slae.receivers[k].values[1]);
			f[i] -= w2 *
				(slae.receivers[k].values[1] - prevE[k][1]) *
				(prevE[k][1] - newE[i][k][1]) / (parameters[i].GetValue(slae) * delta);
		}
	}
}


template<typename testParameter>
inverseTask<testParameter>::inverseTask()
{
	cout << "Count of elements: " << slae.elements.size() << endl;
	cout << "Count of nodes: " << slae.x.size() << endl;

	parameters.resize(nParameters);
	/* создаем приёмники и считаем снятые в них данные E = -gradV */
	CreateReceivers();
	/* устанавливем начальные значения искомых параметров (сигма/мощность) */
	SetParameters();
	/* решаем с ними прямую задачу */
	slae.SolveDT();

	cout.precision(17);

	/* искомые параметры на предыдущей итерации */
	array <double, nParameters> prevParameters;
	for (unsigned int i = 0; i < nParameters; i++)
		prevParameters[i] = parameters[i].GetValue(slae);

	/* значения E=-gradV на предыдущей итерации */
	int nReceivers = slae.receivers.size();
	array <array<double, 2>, 3> prevE;
	for (unsigned int i = 0; i < nReceivers; i++) {
		for (int k = 0; k < 2; k++)
			prevE[i][k] = 0;
		CalculateE(i, prevE[i]);
	}

	/* значения E = -gradV при приращении параметров */
	array <array <array<double, 2>, 3>, nParameters> newE;
	for (int i = 0; i < nParameters; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 2; k++)
				newE[i][j][k] = 0;

	/* получаем J(u) */
	double Ju = GetJu();

	bool flag = true;

	for (unsigned int iter = 1; iter<2; iter++) {
		/* расчет новых значений при приращении параметров */
		for (unsigned int i = 0; i < nParameters; i++) {
			double u = parameters[i].GetValue(slae);
			parameters[i].SetValue(slae, u + u * delta);
			slae.SolveDT();
			for (unsigned int j = 0; j < slae.receivers.size(); j++)
			{
				for (int k = 0; k < 2; k++)
					newE[i][j][k] = 0;
				CalculateE(j, newE[i][j]);
			}
			parameters[i].SetValue(slae, u);
		}
		/* расчёт элементов матрицы и правой части СЛАУ для параметров (без альфы) */
		GetInverseSLAE(prevE, newE);

		/* выбираем альфа */
		double alpha = 0.0;
		double sum = 0.0;
		for (unsigned int i = 0; i < nParameters; i++) 
			sum += pow((prevParameters[i] - parameters[0].GetParameter()), 2);
		
		if (sum) alpha = gamma * Ju / sum;

		/* учитываем регуляризацию в СЛАУ */
		for (unsigned int i = 0; i < nParameters; i++) {
			A[i][i] += alpha; // alpha * I
			f[i] -= alpha * (prevParameters[i] - parameters[0].GetParameter());
		}

		/* находим направление поиска дельта_u (с) */
		deltaU[0] = f[0] / A[0][0];


		/* выбираем параметр релаксации бетта и новый вектор параметров среды */
		double betta = 1.0, JuNew;

		while (betta > eps) {
			for (unsigned int i = 0; i < nParameters; i++)
				parameters[i].SetValue(slae, prevParameters[i] + betta * deltaU[i]);

			JuNew = GetJu();
			if (JuNew > Ju)
				betta /= 3;
			else
				break;
		}
		Ju = JuNew;

		
		/* если не удалось подобрать бету, фиксируем старые параметры и выходим из цикла */
		if (betta <= eps) {
			cout << "Stagnation detected." << endl;
			cout << "Number of iterations: " << iter << endl;
			for (unsigned int i = 0; i < nParameters; i++)
				cout << "param[" << i << "] = " << prevParameters[i] << endl;
			system("pause");
			flag = false;
			break;
		}

		unsigned int count = 0;
		for (unsigned int i = 0; i < nParameters; i++)
			if (fabs((prevParameters[i] - parameters[i].GetValue(slae)) / prevParameters[i]) > eps)
				count++;

		if (count == 0) {
			cout << "Stagnation detected." << endl;
			cout << "Number of iterations: " << iter << endl;
			for (unsigned int i = 0; i < nParameters; i++)
				cout << "param[" << i << "] = " << prevParameters[i] << endl;
			system("pause");
			flag = false;
			break;
		}

		/* обновляем старые значения */
		slae.SolveDT();
		for (unsigned int i = 0; i < nReceivers; i++)
			CalculateE(i, prevE[i]);
		for (unsigned int i = 0; i < nParameters; i++)
			prevParameters[i] = parameters[i].GetValue(slae);

		for (unsigned int i = 0; i < nParameters; i++)
			cout << "<" << iter << "> param[" << i << "] = " << prevParameters[i] << endl;

		if (Ju <= eps) {
			cout << "Solution found." << endl;
			cout << "Number of iterations: " << iter << endl;
			for (unsigned int i = 0; i < nParameters; i++)
				cout << "param[" << i << "] = " << prevParameters[i] << endl;
			system("pause");
			flag = false;
		}

	}
}

template<typename testParameter>
inverseTask<testParameter>::~inverseTask()
{

}

