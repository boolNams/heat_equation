#include <iostream>
#include "HEAT.h"
using namespace std;

double f(double x, double t)
{
	return pow(t,2) * sqrt(x);
}

double psi0(double t)
{
	return 1.0 - pow(t,3);
}

double psi1(double t)
{
	return t;
}

double phi (double x)
{
	return 1.0 - x;
}

int main(void)
{
	int N = 16, M = 16;
	bool plot = true;

	ImplicitScheme(N, M, plot);

	cout << "STOP" << endl;
	return 0;
}
