#include "pch.h"
#include <fstream>   /* Работа с файлами */
#include <iostream>  /* Поток ввода / вывода */
#include <iomanip>
#include <string>    /* Строки */
#include <cmath>     /* Математика */

#include "LineAlgebra.h" /* Моя библиотека линейной алгебры */

using std::string;
using std::cout; /* Для вывода в консоль */
using std::endl; /* Для перехода на другую строку */

using std::ofstream; /* Для файла */

typedef double type;

/*typedef type (*Name_function) (type argument); /* Здесь Name_function - имя функции, *Name_function - указатель на функцию,
											   (type argument) - аргументы которые принимает функция на которую мы указываем,
											   другими словами это можно назвать */

const static int n = 127; /* Число отрезков разбиения */
const static int between_nodeGrid = 1; /* Число точек между узлами */
type all_point = n * pow(2, between_nodeGrid);


static type a = 0; /* Левая граница отрезка */
static type b = 0; /* Правая граница отрезка */

string fileNameUnGrid = "Uniform_Grid.txt";  /* Файл с равномерной сеткой */


                                     /* Если надо использовать функцию std::remove, 
									  то использовать только тип "const char" 
									  std::remove - удаляет файл */

string fileNameChebyshevGrid = "Chebyshev_Grid.txt";  /* Файл с чебышевской сеткой */

string fileNameInterpolated = "Interpolated_function_L.txt"; /* Файл с координатами интерполянта методом лагранжа */

string fileNameInterpolateCube = "Interpolated_function_C.txt";

string VARIBLE_NAME = "NULL"; /* Глобальная переменная в которой
							  будет храниться имя переменной*/

#define SAVE_VARIBLE_NAME(varible) VARIBLE_NAME = string(#varible); /* Создано для того, чтобы 
																	 записывать имя самой переменной
																	 в VARIBLE_NAME и выводить на экран */

/* Функции см. в методичке */
 /* -1 <= x <= 1 */
type y_1(type x) {
	return x * x;
}

/* -1 <= x <= 1 */
type y_2(type x) {
	return 1 / (1 + x * x);
}

/* -3 <= x <= 3 */
type y_3(type x) {
	return 1 / atan(1 + 10 * x*x);
}

 /* -3 <= x <= 3 */
type y_4(type x) {

	if (4 * x*x*x + 2 * x* x - 4 * x + 2 <= 0) {
		return -1000;
	} else {
       	return  pow(4*x*x*x + 2*x*x - 4*x + 2, sqrt(2)) + asin(1 / (5 + x - x * x)) - 5;
	}

}

type y_5(type x) {
	return 1;
}

 /* -1 <= x <= 1 */
type y_6(type x) {
	return 1 / (1 + 25 * x * x);
}

/* -1 <= x <= 1 */
type y_7(type x) {
	return x;
}

 /* -1 <= x <= 1 */
type y_8(type x) {
	return sin(x);
}


/*-------------------------------------------------------------------------------------------------------------------------------*/

/* Запись данных в файл */
/* Здесь записываются значения в 2 столбца : x   y(x) */
void WritingDataToFile(string fileName, type* l_grid_x, type* l_grid_y, int N) {
	                                         /* type* l_grid_x - аргумент в сетки
											 type* l_grid_y - значение в сетке
											 string fileName - имя файла, куда производится запись значений интерполяции */


	ofstream File(fileName);                 /* Создаем файл */
	File.close();                            /* Закрываем файл */

	File.open(fileName, std::ios_base::out); /* ios_base::out - открытие файла для записи */

	if (!File.is_open()) {
		cout << "ERROR! FILE IS NOT OPEN! \n";
	} else {
		for (size_t i = 0; i < N + 1; i++) {
			if (i == N) {
				File << l_grid_x[i] << "\t" << l_grid_y[i];
			} else {
				File << l_grid_x[i] << "\t" << l_grid_y[i] << "\n";
			}
		}
	}

	File.close();
}

type error(type(*function) (type x), type* arg, type* val, int N) {

	type result = 0;
	type tmp = 0;

	for (size_t i = 0; i < N; i++) {

		tmp = abs(function(arg[i]) - val[i]);

		if (tmp > result) {
			result = tmp;
		}

	}

	return result;
}

/* Создание равномерной сетки */
void CreateUniformGrid( type(*function) (type x),  type* l_grid_x, type* l_grid_y, type& a, type b, int N, string fileName) {
	/* 									 type (*function) (type x) - функция интерполирования 
	                                     type* l_grid_x - аргумент в сетки 
										 type* l_grid_y - значение в сетке 
										 type a - левая часть отрезка, где проводится интерполяция 
										 type b - правая часть отрезка, где проводится интерполяция
										 int N - число отрезков разбиения
										 string fileName - имя файла, куда производится запись значений интерполяции */
	
	type h = (b - a) / N; /* Шаг сетки */

	for (size_t i = 0; i < N + 1; i++) {
		l_grid_x[i] = a + h * i;
		l_grid_y[i] = function(l_grid_x[i]);
	}


	WritingDataToFile(fileName, l_grid_x, l_grid_y, N);
}

/* Создание Чебышевской сетки */
void CreateChebyshevGrid( type(*function) (type x), type* l_grid_x,  type* l_grid_y, type a, type b, int N, string fileName) {
	/* 								 type (*function) (type x) - функция интерполирования
									 type* l_grid_x - аргумент в сетки
									 type* l_grid_y - значение в сетке
									 type a - левая часть отрезка, где проводится интерполяция
									 type b - правая часть отрезка, где проводится интерполяция
									 int N - число отрезков разбиения
									 string fileName - имя файла, куда производится запись значений интерполяции */

	type h = (b - a) / N; /* Шаг сетки */
    constexpr type PI = 3.14; /* constexpr рассчитывает объект на этапе компиляции */

	for (size_t i = 0; i < N + 1 ; i++) {

		if (i == 0) {

			l_grid_x[i] = a;
			l_grid_y[i] = function(l_grid_x[i]);

		} else if (i == N) {
			l_grid_x[i] = b;
			l_grid_y[i] = function(l_grid_x[i]);

		} else {
			l_grid_x[i] = (a + b) / 2.0 + (b - a) / 2.0 * cos((2.0*(i - 1) + 1.0)*PI / (2.0*(N - 1)));
			l_grid_y[i] = function(l_grid_x[i]);
		}
	}

	WritingDataToFile(fileName, l_grid_x, l_grid_y, N);
}

/* Сортировка по возрастанияю сетки */
void OrderingGrid(type* l_grid_x, type* l_grid_y, int N) {
	/* 								 type* l_grid_x - аргумент в сетки
									 type* l_grid_y - значение в сетке
									 int N - количество отрезков разбиения
									 */

	type tmp = 0;

	for (size_t i = 0; i < N - 1; i++) {
		for (size_t j = 0; j < N - i - 1; j++) {

			if (l_grid_x[j] > l_grid_x[j + 1]) {
				tmp = l_grid_x[j];
				l_grid_x[j] = l_grid_x[j + 1];
				l_grid_x[j + 1] = tmp;

				tmp = l_grid_y[j];
				l_grid_y[j] = l_grid_y[j + 1];
				l_grid_y[j + 1] = tmp;
			}
		}
	}

	for (size_t i = 0; i < N; i++) {
		cout << l_grid_x[i] << "\t" << l_grid_y[i] << "\n";
	}
}

/* Создание полинома Ланграндж с подстановкой в него числа */
type Polinom_Lagrange(const type* l_grid_x, const type* l_grid_y, type argument, type a, type b, int N) {
	/*								 type* l_grid_x - аргумент в сетки
									 type* l_grid_y - значение в сетке
									 type argument - точка в которой считается поном 
									 type a - левая часть отрезка, где проводится интерполяция
									 type b - правая часть отрезка, где проводится интерполяция
									 int N - число отрезков разбиения
									 */
	type h = (b - a) / N;

	type c_k = 0;  /*   _                                            ^
				       | | (x - x_i)        i != k                   |
			     c_k = | |-----------               x = argument     |
				       | |(x_k - x_i) */

	type y_i = 0; /*
				       ___
					   \    c_k(x_i) * y_k, где y_k - узел сетки
				y_i  = /
					   ---      */
	
	for (size_t k = 0; k < N + 1; k++) {
		c_k = 1;
		for (size_t i = 0; i < N + 1; i++) {
			if (i != k) {
				c_k *= (argument - l_grid_x[i]) / (l_grid_x[k] - l_grid_x[i]);
			}
		}

		y_i += c_k * l_grid_y[k];
	}

	return y_i;
}

/* Интерполяция методом Лагранджа */
void Interpolation_Lagrange(type(*function) (type x), type* l_grid_x, type* l_grid_y, type a, type b, int N, int numberPointInterpol, string fileName) {                                       
	                                              /* type (*function) (type x) - функция интерполирования 
												     type* l_grid_x - аргумент в сетки 
													 type* l_grid_y - значение в сетке 
													 type a - левая часть отрезка, где проводится интерполяция 
													 type b - правая часть отрезка, где проводится интерполяция
													 int N - число отрезков разбиения
													 int numberPointInterpol - всего точек на отрезке
													 string fileName - имя файла, куда производится запись значений интерполяции */


	int k = 0;

	type h = (b - a) / numberPointInterpol;

	type* val_build_f = nullptr;
	val_build_f = new type[numberPointInterpol + 1];   /* Значение интреполированной функции, т.е.
												      той которую постороим по нашей сетке */

	type* arg_by_build = nullptr;
	arg_by_build = new type[numberPointInterpol + 1];  /* Аргумент интреполированной функции, т.е.
												      той которую постороим по нашей сетке */

	type* step = nullptr;
	step = new type[numberPointInterpol + 1];

	type tmp = 0;

	step[k] = a;

	for (size_t i = 0; i < N; i++) {

		tmp = l_grid_x[i + 1] - l_grid_x[i];
		tmp = tmp / (numberPointInterpol / n);

		for (size_t j = 0; j < (numberPointInterpol / n); j++) {
			k++;
			step[k] = step[k - 1] + tmp;
		}
	}

	 /* Добавляем + 1 точку, чтобы вошла правая  граница отрезка интерполяции */
	for (size_t i = 0; i < numberPointInterpol + 1; i++) {
		arg_by_build[i] = step[i];
		val_build_f[i] = Polinom_Lagrange( l_grid_x, l_grid_y, arg_by_build[i], a, b, N);
	}

	WritingDataToFile(fileName, arg_by_build, val_build_f, numberPointInterpol);

	cout << "error of interpolation : \n" << error(function, arg_by_build, val_build_f, numberPointInterpol) << "\n";

	delete[] val_build_f;
	delete[] arg_by_build;
	delete[] step;
}

/* Решение СЛАУ методом прогонки */
type* Direct_sweep(type* a, type* b, type* c, type* d, int N, type* x) {
		                        /* 
							      a - нижняя диагональ в трехдиагональной матрице 
							      b - диагональ в трехдиагональной матрице 
							      c - верхняя диагональ в трехдиагональной матрице 
							      d - правая часть системы */

	type* alpha = nullptr;
	type* beta = nullptr;

	alpha = new type[N];
	beta = new type[N];

	/* Алгоритм метода смотри в книге Дубинского стр. 162 */

	for (size_t i = 1; i < N; i++) {
		if (i == 1) {
			alpha[i] = -c[i] / b[i];
			beta[i] = d[i] / b[i];

		} else if (i == N - 1) {
			beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] + a[i] * alpha[i - 1]);
		}
		else {
			alpha[i] = -c[i] / (b[i] + a[i] * alpha[i - 1]);
			beta[i] = (d[i] - a[i] * beta[i - 1]) / (b[i] + a[i] * alpha[i - 1]);
		}
	}

	for (int i = N - 1; i >= 1; i--) {

		if (i == N - 1) {
			x[i] = beta[i];

		} else {
			x[i] = alpha[i] * x[i + 1] + beta[i];
		}
	}


	delete[] beta;
	delete[] alpha;

	return x;
}

/* Интерполяция методом кубического сплайна */
void Cube_Spline(type* l_grid_x, type* l_grid_y, type* a, type* b, type* c, type* d, int N, int numberPointInterpol, type(*function) (type x))
{
	                        /* h - шаг сетки 
							   under_diag - нижняя диагональ в трехдиагональной матрице 
							   diag - диагональ в трехдиагональной матрице 
							   over_diag - верхняя диагональ в трехдиагональной матрице 
							   rPart - правая часть системы 
							   point - все точки на отрезке [a, b]
							   s - интерполянт */
	type* h = nullptr; /* Шаг сетки */

	int k = 0;
	type tmp = 0;
	type null = 0.0;

	type* under_diag = nullptr;  /* h_i-1 */
	type* diag = nullptr;        /* h_i-1 + h_i */
	type* over_diag = nullptr;   /* h_i */
	type* rPart = nullptr;       /* g_i - g_i-1  */
	type* point = nullptr;        /* Все точки на отрезке */
	type* s = nullptr;           /* Интерполянт */
	type* g = nullptr;


	under_diag = new type[N + 1];
	diag = new type[N + 1];
	over_diag = new type[N + 1];
	rPart = new type[N + 1];
	point = new type[numberPointInterpol + 1];
	s = new type[numberPointInterpol + 2];
	h = new type[N + 1];
	g = new type[N + 1];

	for (size_t i = 0; i < N + 1 ; i++) {
		h[i] = l_grid_x[i + 1] - l_grid_x[i];
		g[i] = (l_grid_y[i + 1] - l_grid_y[i]) / h[i];
	}

	for (size_t i = 1; i < N ; i++) {
		under_diag[i] = h[i];
		diag[i] = 2 * (h[i - 1] + h[i]);
		over_diag[i] = h[i];

		rPart[i] = 3 * (g[i] - g[i - 1]) ; /* g_i */
	}

	under_diag[N] = null;
	under_diag[0] = null;

	diag[0] = 1;
	diag[N] = 1;

	over_diag[0] = null;
	over_diag[N] = null;

	rPart[0] = 1;
	rPart[N] = 1;

	type* x = nullptr;
	x = new type[N + 1];

	c = Direct_sweep(under_diag, diag, over_diag, rPart, N + 1, x);

	delete[] x;

	point[k] = l_grid_x[k];

	for (size_t i = 0; i < N   ; i++) {
		tmp = l_grid_x[i + 1] - l_grid_x[i];
		tmp = tmp / (numberPointInterpol / n);

		for (size_t j = 0; j < (numberPointInterpol / n); j++) {
			k++;
			point[k] = point[k - 1] + tmp;
		}
	}

	point[k] = l_grid_x[N];

	c[0] = null;
	c[N] = null;

	for (size_t i = 0; i < N  ; i++) {
		b[i] =  g[i] - (c[i + 1] + 2* c[i]) * h[i] / 3;
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
		a[i] = l_grid_y[i];
	}

	a[N] = l_grid_y[N];

	Print(a, N + 1);
	Print(b, N + 1);

	delete[] under_diag;
	delete[] diag;
	delete[] over_diag;
	delete[] rPart;

	k = 0;
	for (size_t i = 0; i < N + 1 ; i++) {

		/* Делаем столько шагов, сколько точек между двумя значениями шага + 1*/
		for (size_t j = 0; j < numberPointInterpol / n ; j++) {
			s[k] = a[i] + b[i] * (point[k] - l_grid_x[i]) + c[i] * pow((point[k] - l_grid_x[i]), 2) + d[i] * pow((point[k] - l_grid_x[i]), 3);
			k++;
		}
		
	}

	WritingDataToFile("Interpolation_CubeSpline.txt", point, s, numberPointInterpol);

	cout << "error of interpolation : \n" << error(function, point, s, numberPointInterpol) << "\n";

	delete[] s;
	delete[] h;
	delete[] g;
	delete[] point;
}



int main()
{
	setlocale(LC_ALL, "Russian");

	type(*ptrFunction) (type argument) = NULL;   /* Здесь ptrFunction - имя функции, *ptrFunction - указатель на функцию,
											   (type argument) - аргументы которые принимает функция на которую мы указываем,
											   другими словами это можно назвать функциональным типом. 
											   Сделано для того чтобы передавать функцию в процедуру/ функцию */

	int number = 0;              /* Номер функции */
	bool SelectFunction = false; /* Для того, чтобы пользователь выбрал функцию, независимо от его уровня интелекта */

	cout << " 1) f(x) = x^2 \n";
	cout << " 2) f(x) = 1/(1 + x^2) \n";
	cout << " 3) f(x) = 1/(atan[1 + 10*x^2]) \n";
	cout << " 4) f(x) = pow(4 * x* x*x + 2 * x* x - 4 * x + 2, sqrt(2)) + asin(1 / [5 + x - x * x]) \n";
	cout << " 5) f(x) = 1 \n";
	cout << " 6) f(x) = 1/(1 + 25x^2) \n";
	cout << " 7) f(x) = x \n";
	cout << " 8) f(x) = sin(п*x) \n";
	cout << " U must select function (input number) : ";
	std::cin >> number;

	while (!SelectFunction)
	{
		switch (number)
		{
		case 1:
			ptrFunction = &y_1;  /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1;
			b = -a;
			SelectFunction = true;
			cout << "Ur function f(x) = x^2 \n";
			break;
		case 2:
			ptrFunction = &y_2; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1;
			b = -a;
			SelectFunction = true;
			cout << "Ur function f(x) = 1/(1 + x^2) \n";
			break;
		case 3:
			ptrFunction = &y_3; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -3;
			b = -a;
			SelectFunction = true;
			cout << "Ur function f(x) = 1/(atan[1 + 10*x^2]) \n";
			break;
		case 4:
			ptrFunction = &y_4; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1.43;
			b = 3;
			SelectFunction = true;
			cout << "Ur function : f(x) = pow(4 * x* x*x + 2 * x* x - 4 * x + 2, sqrt(2)) + asin(1 / [5 + x - x * x]) \n";
			break;
		case 5:
			ptrFunction = &y_5; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1;
			b = 1;
			SelectFunction = true;
			cout << "Ur function : f(x) = 1 \n";
			break;
		case 6:
			ptrFunction = &y_6; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1;
			b = 1;
			SelectFunction = true;
			cout << "Ur function : f(x) = 1/(1 + 25x^2) \n";
			break;
		case 7:
			ptrFunction = &y_7; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1;
			b = 1;
			SelectFunction = true;
			cout << "Ur function : f(x) = x \n";
			break;
		case 8:
			ptrFunction = &y_8; /* В наш указатель на функцию записываю необходимый нам указатель */
			a = -1.25;
			b = 1.25;
			SelectFunction = true;
			cout << "Ur function : f(x) = sin(п*х) \n";
			break;
		default:
			cout << "U didn't select function, try again!\n";
			cout << "Number :";
			std::cin >> number;
			break;
		}
	}

	type* grid_x = nullptr; /* Координата х сетки */
	type* grid_y = nullptr; /* Коoрдината у сетки */

	grid_x = new type[n + 1];
	grid_y = new type[n + 1];

	//CreateChebyshevGrid(ptrFunction, grid_x, grid_y, a, b, n, fileNameChebyshevGrid);  /* Создание Чебышевской сетки */
	//OrderingGrid(grid_x, grid_y, n + 1);                                               /* Сортировка по возрастанию */

	CreateUniformGrid(ptrFunction, grid_x, grid_y, a, b, n, fileNameUnGrid);           /* Создание равномерной сетки */
	//Interpolation_Lagrange(ptrFunction, grid_x, grid_y, a, b, n, all_point, fileNameInterpolated);

	/* Коэф создания интерполянта кубическим сплайном */
	type* a = nullptr; 
	type* b = nullptr; 
	type* c = nullptr; 
	type* d = nullptr; 

	a = new type[n + 1];
	b = new type[n + 1];
	c = new type[n + 1];
	d = new type[n + 1];

	/* Процедура кубического сплайна */
	Cube_Spline(grid_x, grid_y, a, b, c, d, n, all_point, ptrFunction);

	delete[] grid_x;
	delete[] grid_y;

	delete[] a;
	delete[] b;
	delete[] d;
	delete[] c;

		
	system("pause");
}
