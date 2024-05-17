#include "HEAT.h"
#include <math.h>
#include <fstream>
using namespace std;

struct sign_out
{
   double x;
   sign_out() {this->x = 0.0;}
   sign_out(double x){this->x = x;}
};

template<class T> T& operator<<(T &out,sign_out x)
{
   if(x.x > 0.0 || abs(x.x) < 2.2e-308) out<<"+"<<x.x;
   else out<<x.x;
   return out;
}

void save_txt(const char *name, double *x, double *t, double **u, int N, int M)
{
	ofstream out; out.open(name); out.precision(6);
	if(out.is_open())
	{
		out << scientific;
		out << N + 1 << " " << M + 1 << endl;
		for (int i = 0; i < N + 1; ++i) for (int j = 0; j < M + 1; ++j)
		{
			out << sign_out(x[i]) << " " << sign_out(t[j]) << " " << sign_out(u[i][j]);
			out << endl;
		}
	}
}
