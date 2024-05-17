#include "HEAT.h"
#include <math.h>

const char *graph = "python3 graph.py";

void ExplicitScheme(int N, int M, bool G)	// M should be >= 4*N*N for stability
{
	/*  ExplicitScheme
	-N	number of segments along the x axis
	-M	number of segments along the y axis
	-G	plotting parameter
	*/

	double   h  = 1.0 / N, tau = 1.0 / M;
	double  *x  = new double[N + 1];
	double  *t  = new double[M + 1];
	double **u  = new double*[N + 1];

	for (int i = 0; i < N + 1; ++i)
	{
		x[i] = i*h;
		u[i] = new double[M + 1];
		u[i][0] = phi(x[i]);
	}
	for (int j = 0; j < M + 1; ++j)
	{
		t[j] = j*tau;
		u[0][j] = psi0(t[j]);
		u[N][j] = psi1(t[j]);
	}

	double k = tau / pow(h,2);
	for (int j = 0; j < M; ++j)
	{
		for (int i = 1; i < N; ++i)
		{
			u[i][j + 1] = u[i][j] + k * (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) + tau * f (x[i], t[j]);
		}
	}

	if(G)
	{
		save_txt("data.txt", x, t, u, N, M);
		system(graph);
	}

	delete [] x, t, u;
}

void ImplicitScheme(int N, int M, bool G)
{
	/*  ExplicitScheme
	-N	number of segments along the x axis
	-M	number of segments along the y axis
	-G	plotting parameter
	*/

	double   h   = 1.0 / N, tau = 1.0 / M;
	double  *x  = new double[N + 1];
	double  *t  = new double[M + 1];
	double **u = new double*[N + 1];
	for (int i = 0; i < N + 1; ++i) {x[i] = i*h; u[i] = new double[M + 1]; u[i][0] = phi(x[i]);}
	for (int j = 0; j < M + 1; ++j) {t[j] = j*tau; u[0][j] = psi0(t[j]); u[N][j] = psi1(t[j]);}

	double sigma  = 0.5; int n = N - 1;
	double *A = new double[3*n]; double *B = new double[n];

	double k1 = -sigma / pow(h, 2), k2 = 1.0 / tau + 2.0 * sigma / pow(h, 2);
	double q1 = (1 - sigma) / pow(h, 2), q2 = (1.0 / tau) - 2.0 * (1.0 - sigma) / pow(h, 2);
	A[1      ] = k2; A[2      ] = k1; A[0      ] = 0.0;
	A[3*n - 3] = k1; A[3*n - 2] = k2; A[3*n - 1] = 0.0;
	for (int i = 1; i < n - 1; ++i)
	{
		A[3*i + 0] = k1;
		A[3*i + 1] = k2;
		A[3*i + 2] = k1;
	}

	double *y = new double[n + 1], *ksi = new double[n + 1], *eta = new double[n + 1];
	for (int j = 0; j < M; ++j)
	{

		for (int i = 0; i < n; ++i)
		{
			B[i] = (1 - sigma) * f(x[i + 1], t[j]) + sigma * f(x[i + 1], t[j + 1]) + q1*(u[i][j] + u[i + 2][j]) + q2 *u[i + 1][j];
		}
		B[0]     -= k1 * u[0][j + 1];
		B[n - 1] -= k1 * u[N][j + 1];

		//========================================================================== Tridiagonal matrix algorithm

		ksi[0] = 0.0; eta[0] = 0.0; y[n] = 0.0;
		for (int i = 0; i < n; ++i)
		{
			ksi[i + 1] = A[3*i + 2] / (-A[3*i + 1] - A[3*i] * ksi[i]);
			eta[i + 1] = (A[3*i] * eta[i] - B[i]) / (-A[3*i + 1] - A[3*i] * ksi[i]);
		}

		for (int i = n - 1; i >= 0; --i) y[i] = ksi[i + 1] * y[i + 1] + eta[i + 1];
		for (int i = 0; i < n; ++i) u[i + 1][j + 1] = y[i];
	}
	delete [] y, ksi, eta;

	if(G)
	{
		save_txt("data.txt", x, t, u, N, M);
		system(graph);
	}

	delete [] x, t, A, B, u;
}
