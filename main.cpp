#include <iostream>
using namespace std;

const int n = 4;
double A[n][n] = {10, 6, 2, 0,
		  5, 1, -2, 4,
		  3, 5, 1, -1,
		  0, 6, -2, 2};
double B[n] = {25,
	       14,
	       10,
	       8};
double L[n][n] = {};
double U[n][n] = {};

double T1[n][n] = {};
double T2[n] = {};

void LU(double (&A)[n][n]){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			U[i][j] = A[i][j];

	for(int k = 1; k < n; k++)
	{
		for(int i = k - 1; i < n; i++)
			for(int j = i; j < n; j++)
				L[j][i] = U[j][i] / U[i][i];

		for(int i = k; i < n; i++)
			for(int j = k - 1; j < n; j++)
				U[i][j] = U[i][j] - L[i][k-1] * U[k-1][j];
	}
}

void Modif(double (&A)[n][n], double (&B)[n]){
	for (int i = 0; i < n; i++){
		T2[i] = B[i];
		for (int j = 0; j < n; j++){
			T1[i][j] = A[i][j];
		}
	}

	for (int j = 0; j < n - 1; j++){
		double imax = A[j][j];
		int nmax = j;
		for (int i = 0; i < n; i++){
			if (A[i][j] > imax) nmax = i;
		}

		double temp[n] = {};
		double btemp;

		for (int k = 0; k < n; k++){
			temp[k] = A[nmax][k];
		}
		btemp = B[nmax];

		for (int k = 0; k < n; k++){
			A[nmax][k] = A[j][k];
		}
		B[nmax] = B[j];

		for (int k = 0; k < n; k++){
			A[j][k] = temp[k];
		}
		B[j] = btemp;
	}
}

void Mprod(double (&A)[n][n], double (&B)[n][n]){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			T1[i][j] = 0;

	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			for(int k = 0; k < n; k++)
				T1[i][j] += A[i][k] * B[k][j];
}

void Mprod(double (&A)[n][n], double (&B)[n]){
	for (int i = 0; i < n; i++)
		T2[i] = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			T2[i] += A[i][j] * B[j];
}

void Show(double (&A)[n][n]){
	for (int i = 0; i < n; i++){
		cout << "|";
		for (int j = 0; j < n; j++){
			cout << A[i][j] << "\t";
		}
		cout << "|" << endl;
	}
}

void Show(double (&A)[n]){
	for (int i = 0; i < n; i++)
		cout << "|" << A[i] << "\t|" << endl;
}

void BiShow(double (&A)[n][n], double (&B)[n]){
	for (int i = 0; i < n; i++){
		cout << "|";
		for (int j = 0; j < n; j++){
			cout << A[i][j] << "\t";
		}
		cout << "|" << "x" << i << "|\t" << B[i] << endl;
	}
}

void Invert(double (&A)[n][n]){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			T1[i][j] = A[i][j];

	double temp;
	double E[n][n];
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			E[i][j] = 0;
			if (i == j) E[i][j] = 1;
		}
	}

	for (int k = 0; k < n; k++){
		temp = T1[k][k];
		for (int j = 0; j < n; j++){
			T1[k][j] /= temp;
			E[k][j] /= temp;
		}
		for (int i = k + 1; i < n; i++){
			temp = T1[i][k];
			for (int j = 0; j < n; j++){
				T1[i][j] -= T1[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = n - 1; k > 0; k--){
		for (int i = k - 1; i >= 0; i--){
			temp = T1[i][k];
			for (int j = 0; j < n; j++){
				T1[i][j] -= T1[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			T1[i][j] = E[i][j];
}

void main(){
	setlocale(LC_ALL, "ru");
	cout.setf(ios::fixed);
	cout.fill('0');
	cout.width(3);
	cout.precision(2);

	cout << "???????????????? ?????????????? ??????????????????" << endl;
	cout << "[A] * [x] = [B]" << endl;
	BiShow(A, B);
	cout << "???? ???????????? ???????????? ???????????????????????? ?????????????? ?? ??????, ?????????? ???? ?????????????????? ???? ???????? ??????????:" << endl;
	Modif(A, B);
	BiShow(A, B);
	cout << "???????????????????? ?????????????? ?? ????????: [L] * [U] * [x] = [B]" << endl;
	LU(A);
	cout << "?????? ?????????? ???????????? ?????????????? L:" << endl;
	Show(L);
	cout << "?? ?????????????? U:" << endl;
	Show(U);
	cout << "?????????? [y] = [U] * [x], ?????????? [L] * [y] = [B]" << endl;
	cout << "???????????? ???????????? [y] = [L-1] * [B]" << endl;
	cout << "[L-1] =" << endl;
	Invert(L);
	Show(T1);
	cout << "[y] =" << endl;
	Mprod(T1, B);
	Show(T2);
	cout << "???????????? ?????????? [U] * [x] = [y]; [x] = [U-1] * [y]" << endl;
	cout << "[U-1] =" << endl;
	Invert(U);
	Show(T1);
	cout << "?????????? [x] =" << endl;

	double temp[n];
	for (int i = 0; i < n; i++)
		temp[i] = T2[i];

	Mprod(T1, temp);
	Show(T2);

	system("pause");
}
