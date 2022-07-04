#ifndef LineAlgebra
#define LineAlgebra
#include "pch.h"

using std::string;

typedef double type;

//const type ZeroEps = 0.1e-5; // null for float
const type ZeroEps = 1e-7; // null for double
const type Eps = 1e-4; 
const type EpsForTri = 1e-4;


static bool sys_denerate = true;  /* static говороит о том, что данная 
								  переменная видна только в этом файле */

void Print(type** inPtrMatr, int inRow); /* Вывод матрицы в консоль */

void Print(const type* inPtrRight, int inRow, string var); /* Вывод вектора в консоль */

void Print(const type* inPtrRight, int inRow); /* Вывод массива в строку */

void Print(type** inPtrMatr, const type* inPtrRight, int inRow, int inColumn); /* Вывод СЛАУ в консоль */

void CopyMatrix(type* b, type* b_2, int N); /* Копирование вектора */

void CopyMatrix(type** ptrMatr, type** ptrMatr_2, int N); /* Копирование Матрицы */

void CopyMatrix(type** ptrMatr, type** ptrMatr_2, type* b, type* b_2, int N); /* Копирование всей СЛАУ */

type Norma_Euclidean(type* b, int N); /* Евклидова норма 
									     для вектора */

type Norma_1(type* b, int N); /* Норма вида ||.||_1 
							  Смотрим max_sum по столбцам (абсолютное значение) 
							  в случае вектора сумма элементов столбца */

type Norma_infinity(type* b, int N);  /* Норма вида ||.||_inf 
							           Смотрим max_sum по строкам (абсолютное значение) 
							           в случае вектора max элемент в векторе */

type Norma_infinity(type** ptrMatr, int N);   /* Норма вида ||.||_inf 
							                Смотрим max_sum по строкам (абсолютное значение) */
							             
type Norma_1(type** ptrMatr, int N);   /* Норма вида ||.||_1 
							           Смотрим max_sum по столбцам (абсолютное значение) */

type** CreateOneMatrix(type** ptrMatr, int N); /* Заполнение единичной матрицы */

void TransposeMatrix(type** ptrMatr, int N); /* Транспонирование матрицы */

type** MultMatrix(type** ptrMatr, type** ptrMatrSecond, int N); /* Умножение матриц */

type* MultMatrix(type** ptrMatr, type* vector, int N); /* Умножение матрицы на вектор */

type** MultMatrix(type** Matrix, type numeric, int N); /* Умножение матрицы на число */

void StaggeredMatrix(type** ptrMatr, type* ptrRight, int N);   /* Приведение марицы 
															   к ступенчатому виду */

void SolutionVector_or_returnStroke(type** ptrMatr, type* ptrRight, type* x, int N); /* Поиск решения для 
																					   метода Гаусса  */

type* MGausses(type** ptrMatr, type* ptrRight, type* x, int N); /* Реализация метода Гаусса */

type** SearchInverseMatrix(type** ptrMatr, type** InverseMatr, int N); /* Поиск обратной матрицы */

type* DifferVector(type* vector_1, type* vector_2, int N); /* Разность двух векторов */

type* SumMatrix(type* vector_1, type* vector_2, int N); /* Сумма векторов */

type** SumMatrix(type** Matrix_1, type** Matrix_2, int N); /* Сумма матриц */

void VectorNevyazki(type** ptrMatr, type* ptrRight, type* x, int N);

void DLU(type** A, type** D, type** L, type** U, int N); /* Разложение матрицы на верхнетреугольную - U, 
														 нижнетреугольную - L и диагональную - D */

void FillTheMatrix(type* x, int N); /* Заполнение массива нулями */

void FillTheMatrix(type** Matrix, int N); /* Заполнение матрицы нулями */

void UpperTriangleMatrix(type** Matrix, type** UpTri, int N, bool WriteScreen );   /* Выделение верхнетреугольной 
																						матрицы из исходной */

void LowerTriangleMatrix(type** Matrix, type** DownTri, int N, bool WriteScreen );  /* Выделение нижнетреугольной
																						матрицы из исходной */

void CorrectionMatrix(type** Matrix, type* vector, int N);


#endif