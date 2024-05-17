#pragma once
#include <math.h>

double f(double x, double t);
double psi0(double t);
double psi1(double t);
double phi (double x);

void ExplicitScheme(int N, int M, bool G);
void ImplicitScheme(int N, int M, bool G);

void save_txt(const char *name, double *x, double *t, double **u, int N, int M);
