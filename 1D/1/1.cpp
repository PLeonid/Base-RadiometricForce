#include <iostream> 
#include <fstream> 
#include <cmath> 
#include <stdlib.h> 

using namespace std; 


double const R = 8.3145;
double const u = 4 * 0.001;

int const N = 1001;			//Количество координат 
int const K = 1000;			//Количество проекций скоростей
int const M = 6100;			//Количество временных итераций
int const w = 10;			//Частота вывода в файл и в терминал

double const T1 = 800;		//Температура левой стенки
double const T2 = 200;		//Температура правой стенки
double const T  = 100;		//Температура газа начальная

double pop = 1;				//Поправочный коэфф-т

double Tmax; 
	
	
//Функция проекции скорости от номера
double Vx(int m)
{
	double V;
	
	Tmax = max(T1,T2);
	Tmax = max(Tmax,T);

	
	V = 5 * pow ((R * Tmax / u), 0.5) / (K/2) * (m -(K/2 - 0.5) + 0.5 * abs(m - (K/2 - 0.5))/(m - (K/2 - 0.5)));
	
	return V;
}

int main () 
{ 
	ofstream fout("file.txt"); 

///Зачем?
	Tmax = max(T1,T2);
	Tmax = max(Tmax,T);

	
	double N0 = 1;				//Начальная концентрация  
	double N1;					//Концентрация

	double const L = 0.5;		//Число Куранта максимальное

	int i;						//Координатный счётчик
	int k;						//Скоростной счётчик
	int j;						//Временной счётчик

	double t;					//Временной шаг
	t = abs( L / Vx(0));


	//Массив данных
	double **A = new double* [N+1];
	for (i = 0; i < N+1; i++)
	{
		A[i] = new double [K];
	}


	//Массив данных вспомогательный
	double **B = new double* [N+1];
	for (i = 0; i < N+1; i++)
	{
		B[i] = new double [K];
	}

  
	//Массив температур
	double *Tem = new double [N+1];          


	//Вспомогательные переменные
	double SumA = 0;
	double  C;
	double S1,S2,S3;


	cout << "Start of data generation"<< endl;


	//Задание начальных данных
	for (i=1; i<N; i=i+1)
	{
		SumA = 0;
		//Поиск плотности вероятности из уравнения Максвелла без нормировочного коэффициента
		for (k=0; k<K; k++)
		{
			A[i][k]= exp( - Vx(k) * Vx(k) * u / (2 * R * T));
			SumA += A[i][k];
		}

		//Поиск нормировочного коэффициента
		C = N0/SumA;

		//Подстановка нормировочного коэфф-та
		for (k=0; k<K; k++)
		{
			A[i][k] = A[i][k] * C;
		}

	}
	
/*	
		for (k=0; k<K; k++)
		{ 
			for (i=0; i<N+1; i++)
			{  
				fout << A[i][k] << " ";
			}
			fout << endl;
		}
		fout << endl << endl << endl;
*/

	//Временной цикл
	for (j=1; j<M; j++)
	{

		//Создание 0-ой вспомогательной ячейки
		//Обнуление счётчиков
		SumA = 0;
		N1 = 0;

		//Поиск суммы приходящей волны
		for (k = 0; k < K / 2; k++)
		{
			N1 += A[1][k];
		}

		//Поиск суммы уходящей волны
		for (k = K / 2; k < K; k++)
		{
			A[0][k]= exp( - Vx(k) * Vx(k) * u /(2 * R * T1));		
			SumA += A[0][k];
		}

		//Определение норм. коэфф-та и умножение на него
		C = N1 / SumA;
		for (k = K/2; k<K; k++)
		{
			A[0][k] = C * A[0][k];
		}



		//Создание N-ной вспомогательной ячейки
		//Обнуление счётчиков
		SumA = 0;
		N1 = 0;

		//Поиск суммы приходящей волны
		for (k = K / 2; k < K; k++)
		{
			N1 += A[N-1][k];
		}

		//Поиск суммы уходящей волны
		for (k=0; k<K/2; k++)
		{
			A[N][k]= exp(- Vx(k) * Vx(k) * u/(2 * R * T2));	
			SumA += A[N][k];
		}
	
		//Определение норм. коэфф-та и умножение на него
		C=N1/SumA;
		for (k=0; k<K/2; k++)
		{
			A[N][k]= C * A[N][k];
		}
 

		//Определение следующих значений функции
		for (k=0; k<K; k++)
		{ 
			for (i=1; i<N; i++)
			{  

				//B[i][k]=A[i][k][j+1]=f(A[i][k][j],A[i-1][k][j],A[i+1][k][j])
				B[i][k] = A[i][k] - Vx(k) * t * (A[i+1][k] - A[i-1][k]) + abs(Vx(k)) * t * (A[i+1][k] - 2 * A[i][k] + A[i-1][k]);

			}
		}


		//Возврат к традиционным обозначениям
		for (k=0; k<K; k++)
		{ 
			for (i=1; i<N; i=i+1)
			{  
				A[i][k] = B[i][k];
			}
		}

		//Определение температуры
		for (i=1; i<N; i++)
		{ 

			//Обнуление сумм
			S1=0;
			S2=0;
			S3=0;   
    
			//Подсчёт сумм 
			for (k=0; k<K; k++)
			{ 
				S1 += A[i][k];
				S2 += A[i][k] * Vx(k) * Vx(k);
				S3 += A[i][k] * Vx(k);
			}

			//Подсчёт температуры по формуле T=m/3*(<e^2> - v^2)
			Tem[i] = ((S2 / S1) - (S3 * S3)/(S1 * S1))* u / R/ pop;   ///Проверка формулы
///			Tem[i] = (S2 / S1) * u / R/ pop;

		}
		
		//Граничные температуры
		Tem[0]=T1;
		Tem[N]=T2;

		
		//Вывод данных
		if ((j%w) == 0 ) 
		{
			fout << endl << endl << endl;
			for (i=0; i<N+1; i++)
			{
				fout << i << " "<< Tem[i] << endl;
			}
		}

/*
		for (k=0; k<K; k++)
		{ 
			for (i=0; i<N+1; i++)
			{  
				fout << A[i][k] << " ";
			}
			fout << endl;
		}
		fout << endl << endl << endl;
		
*/		
		//Вывод этапа исполнения в терминал
		if ((j%w) == 0)
		{
			cout << j << " / " << M << endl;

		}

	}   

	cout << "End of data generation"<< endl;

	fout.close();
	return 0; 
}
