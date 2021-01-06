#pragma once
#include <fstream>
#include <vector>

struct _spline_data
{
	std::vector<double> b_k;
	std::vector<double> h_k;
	std::vector<double> s_k;
	std::vector<double> AL;
	std::vector<double> AU;
	std::vector<double> AD;
} typedef  spline_data;

void calc_spline();