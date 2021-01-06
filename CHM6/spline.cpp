#include "spline.h"

bool read_data(std::vector<double>& Y_var, std::vector<double>& X_var)
{
	std::ifstream datafile("data.txt");
	int n = 0;

	if (datafile.is_open())
	{
		datafile >> n;
		Y_var.resize(n);
		X_var.resize(n);

		for (int i = 0; i < n; i++)
		{
			datafile >> X_var[i];
			datafile >> Y_var[i];
		}
	}
	else
		return false;

	datafile.close();

	return true;
}

bool type_data(std::vector<double>& X_var, spline_data& spl_data)
{
	std::ofstream datafile("spline.txt");
	int n = X_var.size();

	if (datafile.is_open())
	{
		for (int i = 0; i < n; i++)
		{
			datafile << X_var[i] << " ";
			datafile << spl_data.s_k[i] << std::endl;
		}
	}
	else
		return false;

	datafile.close();

	return true;
}

void spline_coeff(std::vector<double>& Y_var, std::vector<double>& X_var, spline_data& spl_data)
{
	int n = Y_var.size();

	spl_data.h_k.resize(n - 1);
	spl_data.b_k.resize(n);
	spl_data.AL.resize(n - 1);
	spl_data.AD.resize(n);
	spl_data.AU.resize(n - 1);
	spl_data.s_k.resize(n);

	for (int i = 1; i < n; i++)
		spl_data.h_k[i - 1] = X_var[i] - X_var[i - 1];

	spl_data.b_k[0] = (-3 * spl_data.h_k[0] - 2 * spl_data.h_k[1]) / (spl_data.h_k[0] * (spl_data.h_k[0] + spl_data.h_k[1])) * Y_var[0];
	spl_data.b_k[0] += (spl_data.h_k[0] + spl_data.h_k[1]) / (spl_data.h_k[0] * spl_data.h_k[1]) * Y_var[1];
	spl_data.b_k[0] -= spl_data.h_k[0] / ((spl_data.h_k[0] + spl_data.h_k[1]) * spl_data.h_k[1]) * Y_var[2];
	spl_data.b_k[0] /= 2;

	spl_data.b_k[n - 1] = spl_data.h_k[n - 2] / ((spl_data.h_k[n - 2] + spl_data.h_k[n - 3]) * spl_data.h_k[n - 3]) * Y_var[n - 3];
	spl_data.b_k[n - 1] -= (spl_data.h_k[n - 2] + 2 * spl_data.h_k[n - 3]) * Y_var[n - 2] / (spl_data.h_k[n - 2] * spl_data.h_k[n - 3]);
	spl_data.b_k[n - 1] += (3 * spl_data.h_k[n - 2] + 2 * spl_data.h_k[n - 3]) * Y_var[n - 1] / (spl_data.h_k[n - 2] * (spl_data.h_k[n - 2] + spl_data.h_k[n - 3]));
	spl_data.b_k[n - 1] /= 2;

	spl_data.AD[0] = 1;
	spl_data.AD[n - 1] = 1;
	spl_data.AU[0] = 0;

	for (int i = 1; i < n - 1; i++)
	{
		spl_data.AL[i - 1] = 2 / spl_data.h_k[i - 1];
		spl_data.AD[i] = 4 * (1 / spl_data.h_k[i - 1] + 1 / spl_data.h_k[i]);
		spl_data.AU[i] = 2 / spl_data.h_k[i];

		spl_data.b_k[i] = -6 * Y_var[i - 1] / (spl_data.h_k[i - 1] * spl_data.h_k[i - 1]);
		spl_data.b_k[i] += 6 * Y_var[i] * (1 / (spl_data.h_k[i - 1] * spl_data.h_k[i - 1]) - 6 / (spl_data.h_k[i] * spl_data.h_k[i]));
		spl_data.b_k[i] += 6 * Y_var[i + 1] / (spl_data.h_k[i] * spl_data.h_k[i]);
	}
}

void lingauss(spline_data& spl_data)
{
	unsigned long long n = spl_data.b_k.size();
	std::vector<double> x;
	x.resize(n);

	double temp;
	int k = 0;

	for (unsigned int i = 0; i < n - 1; i++)
	{
		temp = spl_data.AL[i] / spl_data.AD[i];
		spl_data.AL[i] -= temp * spl_data.AD[i];
		spl_data.AD[i + 1] -= temp * spl_data.AU[i];
		spl_data.b_k[i + 1] -= temp * spl_data.b_k[i];
	}

	spl_data.s_k[n - 1] = spl_data.b_k[n - 1] / spl_data.AD[n - 1];

	for (int i = n - 2; i >= 0; i--)
		spl_data.s_k[i] = (spl_data.b_k[i] - spl_data.s_k[i + 1] * spl_data.AU[i]) / spl_data.AD[i];
}

void calc_spline()
{
	spline_data spl_data;
	std::vector<double> Y_var;
	std::vector<double> X_var;

	read_data(Y_var, X_var);
	spline_coeff(Y_var, X_var, spl_data);
	lingauss(spl_data);
	type_data(X_var, spl_data);
}