#include "pch.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "LineAlgebra.h"


void Print(type** inPtrMatr, int inRow)
{
	for (size_t i = 0; i < inRow; i++)
	{
		for (size_t j = 0; j < inRow; j++)
		{
			std::cout << std::setw(10) << inPtrMatr[i][j] << '\t';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Print(const type* inPtrRight, int inRow, string var)
{
	for (size_t i = 0; i < inRow; i++)
	{
		std::cout << var <<"[" << i << "] = " << inPtrRight[i] << std::endl;
	}
}

void Print(const type* inPtrRight, int inRow)
{
	for (size_t i = 0; i < inRow; i++)
	{
		std::cout << std::setw(10) << inPtrRight[i] << '\t';
	}
	
	std::cout << "\n";
}

void Print(type** inPtrMatr, const type* inPtrRight, int inRow, int inColumn)
{
	bool indicator_first;

	std::cout << std::endl;
	for (size_t i = 0; i < inRow; i++)
	{
		indicator_first = true;
		for (size_t j = 0; j < inColumn; j++)
		{
			if (inPtrMatr[i][j] != 0)
			{
				if (indicator_first)              // Условие на то, что множитель = 0 или нет
				{
					if (inRow - j == 0)                                              // Условие того, что элемент является правй частью уравнения
					{
						std::cout << "= " << inPtrRight[i] << std::endl;    // т.к. храню массив указателей и массив элементов, где
					}
					else                                                             // последний элемент в каждом массиве - правый элемент ур-ия
					{
						if (inPtrMatr[i][j] == 1)
						{
							std::cout << "x" << j + 1 << " ";        // Вывод элементов с множителем 1
						}
						else
						{
							std::cout << inPtrMatr[i][j] << "*x" << j + 1 << " ";
						}

					}
				}
				else
				{
					if (inPtrMatr[i][j] != 0)              // Условие на то, что множитель = 0 или нет
					{
						if (inRow - j == 0)                                              // Условие того, что элемент является правй частью уравнения
						{
							std::cout << "= " << inPtrRight[i] << std::endl;    // т.к. храню массив указателей и массив элементов, где
						}
						else                                                             // последний элемент в каждом массиве - правый элемент ур-ия
						{
							if (inPtrMatr[i][j] == 1)
							{
								std::cout << "+ x" << j + 1 << " ";        // Вывод элементов с множителем 1
							}
							else
							{
								if (inPtrMatr[i][j] > 0)
								{
									std::cout << "+ " << inPtrMatr[i][j] << "*x" << j + 1 << " ";
								}
								else
								{
									std::cout << inPtrMatr[i][j] << "*x" << j + 1 << " ";
								}

							}

						}
					}
				}
				indicator_first = false;
			}

		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void CopyMatrix(type* b, type* b_2, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		b_2[i] = b[i];
	}
}

void CopyMatrix(type** ptrMatr, type** ptrMatr_2, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			ptrMatr_2[i][j] = ptrMatr[i][j];
		}
	}
}

void CopyMatrix(type** ptrMatr, type** ptrMatr_2, type* b, type* b_2, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		b_2[i] = b[i];
		for (size_t j = 0; j < N; j++)
		{
			ptrMatr_2[i][j] = ptrMatr[i][j];
		}
	}
}

type Norma_Euclidean(type* b, int N)
{
	type result = 0;

	for (size_t i = 0; i < N; i++)
	{
		result += pow(b[i], 2);
	}

	return result;
}

type Norma_1(type* b, int N)
{
	type sum = 0;

	for (size_t i = 0; i < N; i++)
	{
		sum += abs(b[i]);
	}
	return sum;
}

type Norma_infinity(type* vector, int N)
{
	type max = 0;

	for (size_t i = 0; i < N; i++)
	{
		if (max < abs(vector[i]))
		{
			max = abs(vector[i]);
		}
	}
	return max;
}

type Norma_infinity(type** Matrix, int N)
{
	type result = 0;
	type max = 0;

	for (size_t i = 0; i < N; i++)
	{
		result = 0;
		for (size_t j = 0; j < N; j++)
		{
			result = result + abs(Matrix[i][j]);
		}
		if (result > max)
		{
			max = result;
		}
	}
	return max;
}

type Norma_1(type** Matrix, int N)
{
	type result = 0;
	type max = 0;
	for (size_t i = 0; i < N; i++)
	{
		result = 0;
		for (size_t j = 0; j < N; j++)
		{
			result = result + abs(Matrix[j][i]);
		}
		if (result > max)
		{
			max = result;
		}
	}
	return max;
}

type** CreateOneMatrix(type** Matrix, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				Matrix[i][j] = 1;
			}
			else
			{
				Matrix[i][j] = 0;
			}
		}
	}
	return Matrix;
}

void TransposeMatrix(type** Matrix, int N)
{
	type temp = 0;
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = i + 1; j < N; j++)
		{
			if (i != j)
			{
				temp = Matrix[i][j];
				Matrix[i][j] = Matrix[j][i];
				Matrix[j][i] = temp;

			}
		}
	}
}

void StaggeredMatrix(type** ptrMatr, type* ptrRight, int N)    // N - количество неизвестных или размер матрицы
{
	type coef;       // Вспомогательные элементы для записи, чтобы не обращаться к памяти

	int i_max = 0;   // Номер максимального элемента в столбце i

	bool indicator_denerate;  // Индикатор вырожденности матрицы

	type *x = nullptr;    // Неизвестный вектор
	type *HelpVector = nullptr;  // Вспомогательный вектор, для смены мест строк

	type HelpVar = 0;           // Вспомогательная переменная для смены мест правых частей 

	x = new type[N];

	for (size_t k = 0; k < N; k++)
	{
		i_max = k;
		// ------------------------------------------------ Ищем ведущий элемент в столбце ------------------------------------------
		for (size_t m = k + 1; m < N; m++)
		{
			if (fabs(ptrMatr[m][k]) > fabs(ptrMatr[i_max][k]))          // Сравнение элементов в столбце i, для поиска максимального
			{
				i_max = m;
			}
		}

		if (fabs(ptrMatr[i_max][k]) < ZeroEps)              // Если максимальный элемент нуль? Все - выход.
		{
			indicator_denerate = false;
			std::cout << "Matrix is degenerate!!!";     // Матрица вырождена
			std::cout << "cond(A) = infinity. " << std::endl;
			sys_denerate = false;
			return;
		}
		else
		{                                                    // 2
			indicator_denerate = true;    // Матрица не вырождена
			if (i_max != k)
			{
				HelpVector = ptrMatr[i_max];       // Меняем места указателей в массиве указателей
				ptrMatr[i_max] = ptrMatr[k];       //<=>
				ptrMatr[k] = HelpVector;           //меням местами строки

				HelpVar = ptrRight[i_max];       // Меняем местами правые части
				ptrRight[i_max] = ptrRight[k];
				ptrRight[k] = HelpVar;

			}

		}                                                 // 2

		if (indicator_denerate = true)
		{
			for (size_t i = k + 1; i < N; i++)
			{
				double coef = ptrMatr[i][k] / ptrMatr[k][k];     // Чтобы использовать меньше памяти
				for (size_t j = k; j < N; j++)
				{
					ptrMatr[i][j] -= ptrMatr[k][j] * coef;
				}

				ptrRight[i] = ptrRight[i] - coef* ptrRight[k];
			}
		}
	}
}

void SolutionVector_or_returnStroke(type** Metrix, type* vector, type* x, int N)
{
	type temp = 0;             //Экономим память,записывая промежуточные вычисления в переменную

	for (int i = N - 1; i >= 0; i--) // Производим поиск неизвестного вектора обратным ходом Гаусса
	{
		temp = vector[i];                        // Экономим память
		for (int j = i + 1; j < N; j++)
		{
			temp = temp - Metrix[i][j] * x[j];
		}
		x[i] = temp / Metrix[i][i];
		if (fabs(x[i]) < ZeroEps)
		{
			x[i] = 0.0;
		}
	}
}

type* MGausses(type** Matrix, type* vector, type* x, int N)
{
	StaggeredMatrix(Matrix, vector, N);

	SolutionVector_or_returnStroke(Matrix, vector, x, N);

	return x;
}

type** SearchInverseMatrix(type** ptrMatr, type** InverseMatr, int N)
{
	type* ptrSolution = nullptr;
	ptrSolution = new type[N];

	type** ptrOneMatrix = nullptr;
	ptrOneMatrix = new type*[N];

	type** ptrFM = nullptr;
	ptrFM = new type*[N];

	for (size_t i = 0; i < N; i++)
	{
		ptrOneMatrix[i] = new type[N];
		ptrFM[i] = new type[N];
	}

	CopyMatrix(ptrMatr, ptrFM, N);
	ptrOneMatrix = CreateOneMatrix(ptrOneMatrix, N);

	for (size_t i = 0; i < N; i++)
	{
		CopyMatrix(ptrFM, ptrMatr, N);  // По скольку нам нужна всегда исходная матрица
		ptrSolution = MGausses(ptrMatr, ptrOneMatrix[i], ptrSolution, N);
		for (size_t j = 0; j < N; j++)
		{
			InverseMatr[j][i] = ptrSolution[j];
		}
	}

	CopyMatrix(ptrFM, ptrMatr, N);  // Чтобы оставить исходную матрицу без изменений
	for (size_t i = 0; i < N; i++)
	{
		delete[] ptrOneMatrix[i];
		delete[] ptrFM[i];
	}
	delete[] ptrOneMatrix;
	delete[] ptrSolution;
	delete[] ptrFM;

	return InverseMatr;
}

type* DifferVector(type* vector_1, type* vector_2, int N)
{
	type* result;
	result = new type[N];
	for (size_t i = 0; i < N; i++)
	{
		result[i] = vector_1[i] - vector_2[i];
	}

	return result;

	delete[] result;
}

type** MultMatrix(type** Matrix_1, type** Matrix_2, int N)
{
	type **Matrix_res;
	Matrix_res = new type*[N];

	type null = 0.0;

	for (int i = 0; i < N; i++)
	{
		Matrix_res[i] = new type[N];
		for (int j = 0; j < N; j++)
		{
			type t = null;
			for (int k = 0; k < N; k++)
				if (fabs(Matrix_1[i][k]) < ZeroEps || fabs(Matrix_2[k][j]) < ZeroEps)
				{
					t += null;
				}
				else
				{
					t += Matrix_1[i][k] * Matrix_2[k][j];
				}
			Matrix_res[i][j] = t;
		}
	}

	return Matrix_res;

	delete[] Matrix_res;
	for (size_t i = 0; i < N; i++)
	{
		delete[] Matrix_res[i];
	}
}

type* MultMatrix(type** Matrix, type* vector, int N)
{
	type temp = 0;
	type* result = nullptr;

	result = new type[N];

	for (int i = 0; i < N; i++)
	{
	    temp = 0;
		for (int j = 0; j < N; j++)
		{
			temp += Matrix[i][j] * vector[j];
		}
		result[i] = temp;
	}

	return result;

	delete[] result;
} 

type** MultMatrix(type** Matrix, type numeric, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			Matrix[i][j] = numeric * Matrix[i][j];
		}
	}

	return Matrix;
}

type* SumMatrix(type* vector_1, type* vector_2, int N)
{
	type* result;
	result = new type[N];
	for (size_t i = 0; i < N; i++)
	{
		result[i] = vector_1[i] + vector_2[i];
	}

	return result;

	delete[] result;
}

type** SumMatrix(type** Matrix_1, type** Matrix_2, int N)
{
	type** result = nullptr;

	result = new type*[N];
	for (size_t i = 0; i < N; i++)
	{
		result[i] = new type[N];
		for (size_t j = 0; j < N; j++)
		{
			result[i][j] = Matrix_1[i][j] + Matrix_2[i][j];
		}
	}

	return result;

	for (size_t i = 0; i < N; i++)
	{
		delete[] result[i];
	}

	delete[] result;
}

void VectorNevyazki(type** Matrix, type* vector, type* x, int N)
{
	type temp = 0;  // Сумма левой части умноженной на вектор решения
	type norm = 0;
	type norm_1 = 0;
	type norm_infinity = 0;

	type* Nevyazka;

	Nevyazka = new type[N];

	for (size_t i = 0; i < N; i++)
	{
		temp = 0;
		for (size_t j = 0; j < N; j++)
		{
			temp += Matrix[i][j] * x[j];   // Считаем что получается при найденном решении с правой стороны
		}
		Nevyazka[i] = abs(temp - vector[i]);
		norm += pow(Nevyazka[i], 2);
	}
	norm_1 = Norma_1(Nevyazka, N);
	norm_infinity = Norma_infinity(Nevyazka, N);
	norm = sqrt(norm);

	std::cout << std::endl;
	std::cout << "Норма вектора невязки ||b1 - b||_2 = " << norm << std::endl;
	std::cout << "Норма вектора невязки ||b1 - b||_1 = " << norm_1 << std::endl;
	std::cout << "Норма вектора невязки ||b1 - b||_inf = " << norm_infinity << std::endl;
	std::cout << std::endl << std::endl;

	delete[] Nevyazka;

}

void DLU(type** A, type** D, type** L, type** U, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				D[i][j] = A[i][j];
			}
			else
			{
				D[i][j] = 0;
			}

			if (i < j)
			{
				U[i][j] = A[i][j];
			}
			else
			{
				U[i][j] = 0;
			}

			if (i > j)
			{
				L[i][j] = A[i][j];
			}
			else
			{
				L[i][j] = 0;
			}
		}
	}


}

void UpperTriangleMatrix(type** Matrix, type** UpTri, int N, bool WriteScreen = false)
{
	type null = 0.0;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i < j)
			{
				UpTri[i][j] = Matrix[i][j];
			}
			else
			{
				UpTri[i][j] = null;
			}
		}
	}

	if (WriteScreen)
	{
		type norm_inf = 0;
		norm_inf = Norma_infinity(UpTri, N);

		std::cout << " ||C_U||_inf = " << norm_inf << std::endl;
		std::cout << "\n";
	}
}

void LowerTriangleMatrix(type** Matrix, type** DownTri, int N, bool WriteScreen = false)
{
	type null = 0.0;

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i > j)
			{
				DownTri[i][j] = Matrix[i][j];
			}
			else
			{
				DownTri[i][j] = null;
			}
		}
	}

	if (WriteScreen)
	{
		type norm_inf = 0;
		norm_inf = Norma_infinity(DownTri, N);

		std::cout << " ||C_L||_inf = " << norm_inf << std::endl;
		std::cout << "\n";
	}
}

void FillTheMatrix(type* vector, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		vector[i] = 0;
	}
}

void FillTheMatrix(type** Matrix, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			Matrix[i][j] = 0;
		}
	}
}

void CorrectionMatrix(type** Matrix, type* vector,  int N) 
{
	type null = 0;

	type* sum;
	sum = new type[N];
	for (int i = 0; i < N; i++)
	{
		sum[i] = null;
		for (int j = 0; j < N; j++)
		{
			if (i != j)
				sum[i] += abs(Matrix[i][j]);
		}
	}

	for (int i = 0; i < N; i++)
	{
		if (Matrix[i][i] < null && abs(Matrix[i][i]) > sum[i])
		{
			vector[i] = -vector[i];
			for (int j = 0; j < N; j++)
				Matrix[i][j] = -Matrix[i][j];
		}
	}

	delete[] sum;
}