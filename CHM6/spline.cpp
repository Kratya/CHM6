#include "spline.h"
#include <iostream>


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

bool type_data(std::vector<double>& X_var, std::vector<double>& Y_var, spline_data& spl_data)
{
	std::ofstream datafile("spline.txt");
	int n = X_var.size();
	double h = 0.25;
	double ksi, ksi_2, ksi_3;
	int i = 0;
	int k;

	if (datafile.is_open())
	{
		for (double x = X_var[0]; x <= X_var[n - 1]; x = X_var[0] + (++i) * h)
		{
			for (k = 0; k < n - 1; k++)
				if (X_var[k] <= x && x <= X_var[k + 1])
					break;
			if (k == n - 1)
			{
				std::cout << "FAIL";
				continue;
			}

			ksi = (x - X_var[k]) / spl_data.h_k[k];
			ksi_2 = ksi * ksi;
			ksi_3 = ksi_2 * ksi;

			datafile << x << " ";
			datafile << Y_var[k] * (1.0 - 3.0 * ksi_2 + 2.0 * ksi_3) +
				spl_data.b_k[k] * spl_data.h_k[k] * (ksi - 2.0 * ksi_2 + ksi_3) +
				Y_var[k + 1] * (3.0 * ksi_2 - 2.0 * ksi_3) +
				spl_data.b_k[k + 1] * spl_data.h_k[k] * (-ksi_2 + ksi_3) << "\n";
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
	//spl_data.AL.resize(n - 1);
	//spl_data.AD.resize(n);
	//spl_data.AU.resize(n - 1);
	spl_data.s_k.resize(n);

	spl_data.A.resize(n);
	for (int i = 0; i < n; i++)
		spl_data.A[i].resize(n);

	for (int i = 1; i < n; i++)
		spl_data.h_k[i - 1] = X_var[i] - X_var[i - 1];

	spl_data.b_k[0] = (-3.0 * spl_data.h_k[0] - 2.0 * spl_data.h_k[1]) / (spl_data.h_k[0] * (spl_data.h_k[0] + spl_data.h_k[1])) * Y_var[0];
	spl_data.b_k[0] += (spl_data.h_k[0] + 2.0 * spl_data.h_k[1]) / (spl_data.h_k[0] * spl_data.h_k[1]) * Y_var[1];
	spl_data.b_k[0] -= spl_data.h_k[0] / ((spl_data.h_k[0] + spl_data.h_k[1]) * spl_data.h_k[1]) * Y_var[2];
	spl_data.b_k[0] /= 2.0;

	spl_data.b_k[n - 1] = spl_data.h_k[n - 2] / ((spl_data.h_k[n - 2] + spl_data.h_k[n - 3]) * spl_data.h_k[n - 3]) * Y_var[n - 3];
	spl_data.b_k[n - 1] -= (spl_data.h_k[n - 2] + 2.0 * spl_data.h_k[n - 3]) * Y_var[n - 2] / (spl_data.h_k[n - 2] * spl_data.h_k[n - 3]);
	spl_data.b_k[n - 1] += (3.0 * spl_data.h_k[n - 2] + 2.0 * spl_data.h_k[n - 3]) * Y_var[n - 1] / (spl_data.h_k[n - 2] * (spl_data.h_k[n - 2] + spl_data.h_k[n - 3]));
	spl_data.b_k[n - 1] /= 2.0;

	spl_data.A[0][0] = 1.0;
	spl_data.A[n - 1][n - 1] = 1.0;

	for (int i = 1; i < n - 1; i++)
	{
		spl_data.A[i][i - 1] = 2.0 / spl_data.h_k[i - 1];
		spl_data.A[i][i] = 4.0 * (1.0 / spl_data.h_k[i - 1] + 1.0 / spl_data.h_k[i]);
		spl_data.A[i][i + 1] = 2.0 / spl_data.h_k[i];

		spl_data.b_k[i] = -6.0 * Y_var[i - 1] / (spl_data.h_k[i - 1] * spl_data.h_k[i - 1]);
		spl_data.b_k[i] += 6.0 * Y_var[i] * (1.0 / (spl_data.h_k[i - 1] * spl_data.h_k[i - 1]) - 1.0 / (spl_data.h_k[i] * spl_data.h_k[i]));
		spl_data.b_k[i] += 6.0 * Y_var[i + 1] / (spl_data.h_k[i] * spl_data.h_k[i]);
	}
	int i = 0;
}

void lingauss(spline_data& spl_data)
{
	unsigned long long n = spl_data.b_k.size();
	std::vector<double> x;
	x.resize(n);

	double temp;
	int k = 0;

	//for (unsigned int i = 0; i < n - 1; i++)
	//{
	//	temp = spl_data.AL[i] / spl_data.AD[i];
	//	spl_data.AL[i] -= temp * spl_data.AD[i];
	//	spl_data.AD[i + 1] -= temp * spl_data.AU[i];
	//	spl_data.b_k[i + 1] -= temp * spl_data.b_k[i];
	//}
	//
	//spl_data.s_k[n - 1] = spl_data.b_k[n - 1] / spl_data.AD[n - 1];
	//
	//for (int i = n - 2; i >= 0; i--)
	//	spl_data.s_k[i] = (spl_data.b_k[i] - spl_data.s_k[i + 1] * spl_data.AU[i]) / spl_data.AD[i];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j != i)
			{
				temp = spl_data.A[j][i] / spl_data.A[i][i];
				for (int k = i; k < n; k++)
					spl_data.A[j][k] -= temp * spl_data.A[i][k];
				spl_data.b_k[j] -= temp * spl_data.b_k[i];
			}
		}
	}

	for (int i = n - 1; i >= 0; i--)
		spl_data.b_k[i] /= spl_data.A[i][i];
}

void calc_spline()
{
	spline_data spl_data;
	std::vector<double> Y_var;
	std::vector<double> X_var;

	read_data(Y_var, X_var);
	spline_coeff(Y_var, X_var, spl_data);
	lingauss(spl_data);
	type_data(X_var, Y_var, spl_data);
}