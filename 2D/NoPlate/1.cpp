#include <iostream> 
#include <fstream> 
#include <cmath> 
#include <stdlib.h> 

using namespace std; 


double const R = 8.3145;		  	//Универсальная газовая постоянная
double const u = 4 * 0.001;			//Молярная масса газа

int const N = 100;					//Количество координат 
int const P = 100;					//Количество проекций скоростей
int const K = 510;					//Количество временных итераций
int const w = 1;					//Частота вывода в файл и в терминал

double const T1 = 100;				//Температура левой   стенки
double const T2 = 100;				//Температура правой  стенки
double const T3 = 50;				//Температура нижней  стенки
double const T4 = 50;				//Температура верхней стенки

double const Tl = 100;				//Температура левой  стенки пластинки
double const Tp = 100;				//Температура правой стенки пластинки

double const T  = 10;				//Температура газа начальная

double Tmax; 						//Максимальная температура системы

double const pop = 1.04154;				//Поправочный коэффициент
	
//Функция проекции скорости от номера
double V(int m)
{
	double v;
	
	Tmax = max(T1,T2);
	Tmax = max(Tmax,T3);
	Tmax = max(Tmax,T4);
	Tmax = max(Tmax,Tl);
	Tmax = max(Tmax,Tp);
	Tmax = max(Tmax,T);
	
	v = 5 * pow ((R * Tmax / u), 0.5) / (P/2) * (m -(P/2 - 0.5) + 0.5 * abs(m - (P/2 - 0.5))/(m - (P/2 - 0.5)));
	
	return v;
}

int main () 
{ 
	ofstream fout("file.txt"); 

///Зачем?
	Tmax = max(T1,T2);
	Tmax = max(Tmax,T3);
	Tmax = max(Tmax,T4);
	Tmax = max(Tmax,Tl);
	Tmax = max(Tmax,Tp);
	Tmax = max(Tmax,T);

	
	double N0 = 1;					//Начальная концентрация  
	double N1;						//Концентрация

	double const L = 0.5;			//Число Куранта максимальное

	int i,j;						//Кооргинатные счётчики
	int vi,vj;						//Скоростные счётчики
	int k;							//Временной счётчик

	double t;						//Временной шаг
	t = abs( L / 2 / V(0));


	//Массив данных
	//double A[N+1][N+1][P][P];
	
	double ****A = new double*** [N+1];
	for (i = 0; i < N+1; i++)
	{
		A[i] = new double** [N+1];
		for (j = 0; j < N+1; j++)
		{
			A[i][j] = new double* [P];
			for (vi = 0; vi < P; vi++)
			{
				A[i][j][vi] = new double [P];
			}
		}
	}
	
	
	//double B[N+1][N+1][P][P];
	
	//Массив данных вспомогательный
	double ****B = new double*** [N+1];
	for (i = 0; i < N+1; i++)
	{
		B[i] = new double** [N+1];
		for (j = 0; j < N+1; j++)
		{
			B[i][j] = new double* [P];
			for (vi = 0; vi < P; vi++)
			{
				B[i][j][vi] = new double [P];
			}
		}
	}
	
	
//	double Tem[N+1][N+1];
	
	//Массив температур
	double **Tem = new double* [N+1];          
	for (i = 0; i < N+1; i++)
	{
		Tem[i] = new double [N+1];
	}
	
	
	//Вспомогательные переменные
	double SumA;
	double  C;
	double S1,S2,S3;


	cout << "Start of data generation"<< endl;


	//Задание начальных данных
	for (i=1; i<N; i++)
	{
		for (j=1; j<N; j++)
		{
			SumA = 0;
			
			//Поиск плотности вероятности из уравнения Максвелла без нормировочного коэффициента
			for (vi=0; vi<P; vi++)
			{
				for (vj=0; vj<P; vj++)
				{
					A[i][j][vi][vj] = exp( -( (V(vi) * V(vi)) + (V(vj) * V(vj)) ) * u / (2 * R * T));
					SumA = SumA + A[i][j][vi][vj];
				}
			}

			//Поиск нормировочного коэффициента
			C = N0 / SumA;

			//Подстановка нормировочного коэфф-та
			for (vi=0; vi<P; vi++)
			{
				for (vj=0; vj<P; vj++)
				{
					A[i][j][vi][vj] = C * A[i][j][vi][vj];
					B[i][j][vi][vj] = A[i][j][vi][vj];
				}
			}
		}
	}


	//Временной цикл
	for (k=0; k < K; k++)
	{

		//Создание левой вспомогательной строки
		//Обнуление счётчиков
		for (j=0; j < N+1; j++)
		{
			for (vj = 0; vj < P; vj++)
			{	
				SumA = 0;
				N1 = 0;

				//Поиск суммы приходящей волны
				for (vi = 0; vi < P/2; vi++)
				{
					N1 += A[1][j][vi][vj];
				}

				//Поиск суммы уходящей волны
				for (vi = P/2; vi < P; vi++)
				{
					A[0][j][vi][vj]= exp( - ( (V(vi) * V(vi)) + (V(vj) * V(vj)) ) * u /(2 * R * T1));		
					SumA += A[0][j][vi][vj];
				}
		
				//Определение норм. коэфф-та и умножение на него
				C = N1 / SumA;
				for (vi = P/2; vi < P; vi++)
				{
					A[0][j][vi][vj] = C * A[0][j][vi][vj];
				}
			}
		}

		//Создание правой вспомогательной строки
		//Обнуление счётчиков
		for (j=0; j<N+1; j++)
		{
			for (vj = 0; vj < P; vj++)
			{	
				SumA = 0;
				N1 = 0;

				//Поиск суммы приходящей волны
				for (vi = P/2; vi < P; vi++)
				{
					N1 += A[N-1][j][vi][vj];
				}

				//Поиск суммы уходящей волны
				for (vi = 0; vi < P/2; vi++)
				{
					A[N][j][vi][vj]= exp( - ( (V(vi) * V(vi)) + (V(vj) * V(vj)) ) * u /(2 * R * T2));		
					SumA += A[N][j][vi][vj];
				}
		
				//Определение норм. коэфф-та и умножение на него
				C = N1 / SumA;
				for (vi = 0; vi < P/2; vi++)
				{
					A[N][j][vi][vj] = C * A[N][j][vi][vj];
				}
			}
		}

		
		
		//Создание нижней вспомогательной строки
		//Обнуление счётчиков
		for (i=0; i < N+1; i++)
		{
			for (vi = 0; vi < P; vi++)
			{	
				SumA = 0;
				N1 = 0;

				//Поиск суммы приходящей волны
				for (vj = 0; vj < P/2; vj++)
				{
					N1 += A[i][1][vi][vj];
				}

				//Поиск суммы уходящей волны
				for (vj = P/2; vj < P; vj++)
				{
					A[i][0][vi][vj]= exp( - ( (V(vi) * V(vi)) + (V(vj) * V(vj)) ) * u /(2 * R * T3));		
					SumA += A[i][0][vi][vj];
				}
		
				//Определение норм. коэфф-та и умножение на него
				C = N1 / SumA;
				for (vj = P/2; vj < P; vj++)
				{
					A[i][0][vi][vj] = C * A[i][0][vi][vj];
				}
			}
		}

		//Создание верхней вспомогательной строки
		//Обнуление счётчиков
		for (i=0; i<N+1; i++)
		{
			for (vi = 0; vi < P; vi++)
			{	
				SumA = 0;
				N1 = 0;

				//Поиск суммы приходящей волны
				for (vj = P/2; vj < P; vj++)
				{
					N1 += A[i][N-1][vi][vj];
				}

				//Поиск суммы уходящей волны
				for (vj = 0; vj < P/2; vj++)
				{
					A[i][N][vi][vj]= exp( - ( (V(vi) * V(vi)) + (V(vj) * V(vj)) ) * u /(2 * R * T4));		
					SumA += A[i][N][vi][vj];
				}
		
				//Определение норм. коэфф-та и умножение на него
				C = N1 / SumA;
				for (vj = 0; vj < P/2; vj++)
				{
					A[i][N][vi][vj] = C * A[i][N][vi][vj];
				}
			}
		}


		//Определение следующих значений функции
		for (vi=0; vi<P; vi++)
		{
			for (vj=0; vj<P; vj++)
			{
				for (i=1; i<N; i++)
				{
					for (j=0; j<N+1; j++)
					{ 
						//B[i][k]= A[i][k][j+1] =f(A[i][k][j],A[i-1][k][j],A[i+1][k][j])
						B[i][j][vi][vj] = A[i][j][vi][vj] - V(vi) * t * (A[i+1][j][vi][vj] - A[i-1][j][vi][vj]) + abs(V(vi)) * t * (A[i+1][j][vi][vj] - 2 * A[i][j][vi][vj] + A[i-1][j][vi][vj]);
					}
				}
			}
		}
		for (vj=0; vj<P; vj++)
		{
			for (vi=0; vi<P; vi++)
			{
				for (j=1; j<N; j++)
				{
					for (i=1; i<N; i++)
					{ 
						//B[i][k]= A[i][k][j+1] =f(A[i][k][j],A[i-1][k][j],A[i+1][k][j])
						A[i][j][vi][vj] = B[i][j][vi][vj] - V(vj) * t * (B[i][j+1][vi][vj] - B[i][j-1][vi][vj]) + abs(V(vj)) * t * (B[i][j+1][vi][vj] - 2 * B[i][j][vi][vj] + B[i][j-1][vi][vj]);
					}
				}
			}
		}
/*
		//Возврат к традиционным обозначениям
		for (i=1; i<N; i++)
		{
			for (j=1; j<N; j++)
			{
				for (vi=0; vi<P; vi++)
				{
					for (vj=0; vj<P; vj++)
					{
						A[i][j][vi][vj] = B[i][j][vi][vj];
					}
				}
			}
		}
*/		
		
		//Определение температуры
		for (i=1; i<N; i++)
		{ 
			for (j=1; j<N; j++)
			{
				//Обнуление сумм
				S1=0;
				S2=0;
				S3=0;   
    
				//Подсчёт сумм 
				for (vi=0; vi<P; vi++)
				{
					for (vj=0; vj<P; vj++)
					{
						S1 += A[i][j][vi][vj];
						S2 += A[i][j][vi][vj] * ((V(vi) * V(vi)) + (V(vj) * V(vj)));
//						S3 += A[i][j][vi][vj] * Vx(k);
					}
				}
				
///Проверка формул
				//Подсчёт температуры по формуле T=m/3*(<e^2> - v^2)
//				Tem[i][j] = ((S2 / S1) - (S3 * S3)/(S1 * S1)) * u / R / pop;  
				Tem[i][j] = (S2 / S1) * u / R/ 2 / pop;
			}
		}
		
		//Граничные температуры
		for (j=0; j<N+1; j++)
		{	
			Tem[0][j] = T1;
		}
		
		for (j=0; j<N+1; j++)
		{	
			Tem[N][j] = T2;
		}
		
		for (i=0; i<N+1; i++)
		{	
			Tem[i][0] = T3;
		}
		
		for (i=0; i<N+1; i++)
		{	
			Tem[i][N] = T4;
		}
		
		
		//Вывод данных
		if ((k%w) == 0) 
		{
			fout << endl << endl << endl;
			for (i=0; i<N+1; i++)
			{
				for (j=0; j<N+1; j++)
				{
					fout << i << " " << j << " " << Tem[i][j] << endl;
				}
				fout << endl;
			}
		}

		//Вывод этапа исполнения в терминал
		if ((k%w) == 0)
		{
			cout << k << " / " << K << endl;
		}

	}   

	cout << "End of data generation" << endl;

	fout.close();
	return 0; 
}
