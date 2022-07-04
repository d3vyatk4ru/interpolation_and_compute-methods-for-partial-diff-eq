#ifndef LineAlgebra
#define LineAlgebra
#include "pch.h"

using std::string;

typedef double type;

//const type ZeroEps = 0.1e-5; // null for float
const type ZeroEps = 1e-7; // null for double
const type Eps = 1e-4; 
const type EpsForTri = 1e-4;


static bool sys_denerate = true;  /* static �������� � ���, ��� ������ 
								  ���������� ����� ������ � ���� ����� */

void Print(type** inPtrMatr, int inRow); /* ����� ������� � ������� */

void Print(const type* inPtrRight, int inRow, string var); /* ����� ������� � ������� */

void Print(const type* inPtrRight, int inRow); /* ����� ������� � ������ */

void Print(type** inPtrMatr, const type* inPtrRight, int inRow, int inColumn); /* ����� ���� � ������� */

void CopyMatrix(type* b, type* b_2, int N); /* ����������� ������� */

void CopyMatrix(type** ptrMatr, type** ptrMatr_2, int N); /* ����������� ������� */

void CopyMatrix(type** ptrMatr, type** ptrMatr_2, type* b, type* b_2, int N); /* ����������� ���� ���� */

type Norma_Euclidean(type* b, int N); /* ��������� ����� 
									     ��� ������� */

type Norma_1(type* b, int N); /* ����� ���� ||.||_1 
							  ������� max_sum �� �������� (���������� ��������) 
							  � ������ ������� ����� ��������� ������� */

type Norma_infinity(type* b, int N);  /* ����� ���� ||.||_inf 
							           ������� max_sum �� ������� (���������� ��������) 
							           � ������ ������� max ������� � ������� */

type Norma_infinity(type** ptrMatr, int N);   /* ����� ���� ||.||_inf 
							                ������� max_sum �� ������� (���������� ��������) */
							             
type Norma_1(type** ptrMatr, int N);   /* ����� ���� ||.||_1 
							           ������� max_sum �� �������� (���������� ��������) */

type** CreateOneMatrix(type** ptrMatr, int N); /* ���������� ��������� ������� */

void TransposeMatrix(type** ptrMatr, int N); /* ���������������� ������� */

type** MultMatrix(type** ptrMatr, type** ptrMatrSecond, int N); /* ��������� ������ */

type* MultMatrix(type** ptrMatr, type* vector, int N); /* ��������� ������� �� ������ */

type** MultMatrix(type** Matrix, type numeric, int N); /* ��������� ������� �� ����� */

void StaggeredMatrix(type** ptrMatr, type* ptrRight, int N);   /* ���������� ������ 
															   � ������������ ���� */

void SolutionVector_or_returnStroke(type** ptrMatr, type* ptrRight, type* x, int N); /* ����� ������� ��� 
																					   ������ ������  */

type* MGausses(type** ptrMatr, type* ptrRight, type* x, int N); /* ���������� ������ ������ */

type** SearchInverseMatrix(type** ptrMatr, type** InverseMatr, int N); /* ����� �������� ������� */

type* DifferVector(type* vector_1, type* vector_2, int N); /* �������� ���� �������� */

type* SumMatrix(type* vector_1, type* vector_2, int N); /* ����� �������� */

type** SumMatrix(type** Matrix_1, type** Matrix_2, int N); /* ����� ������ */

void VectorNevyazki(type** ptrMatr, type* ptrRight, type* x, int N);

void DLU(type** A, type** D, type** L, type** U, int N); /* ���������� ������� �� ����������������� - U, 
														 ���������������� - L � ������������ - D */

void FillTheMatrix(type* x, int N); /* ���������� ������� ������ */

void FillTheMatrix(type** Matrix, int N); /* ���������� ������� ������ */

void UpperTriangleMatrix(type** Matrix, type** UpTri, int N, bool WriteScreen );   /* ��������� ����������������� 
																						������� �� �������� */

void LowerTriangleMatrix(type** Matrix, type** DownTri, int N, bool WriteScreen );  /* ��������� ����������������
																						������� �� �������� */

void CorrectionMatrix(type** Matrix, type* vector, int N);


#endif