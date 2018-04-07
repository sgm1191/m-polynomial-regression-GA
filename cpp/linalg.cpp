#include <stdio.h>
#include <stdlib.h>

typedef struct{
	double* V;
	int size;
}Vector;

typedef struct{
	double** M;
	int rows;
	int cols;
}Matrix;

Vector* new_vector(double *v_pt, int size)
{
	Vector *v = (Vector*)malloc(sizeof(Vector));
	v->V = v_pt;
	v->size = size;
	return v;
}

Vector* new_vector(int size)
{
	Vector *v = (Vector*)malloc(sizeof(Vector));
	v->V = (double*)malloc(sizeof(double)*size);//new double[size];
	v->size = size;
	return v;

}

Matrix* new_matrix(int rows, int cols)
{
	Matrix *matrix = (Matrix*)malloc(sizeof(Matrix));
	double **A = (double**) malloc(sizeof(double*)*rows);
	for (int i = 0; i < rows; ++i) A[i] =  (double*) malloc(sizeof(double)*cols);
	matrix->M = A;
	matrix->rows = rows;
	matrix->cols = cols;
	return matrix;
}

Matrix* new_I_matrix(int size)
{
	Matrix *M = new_matrix(size,size);
	for (int i = 0; i < size; ++i) M->M[i][i] = 1;
	return M;
}

Matrix* transpose(Matrix *M)
{
	Matrix *T = new_matrix(M->cols,M->rows);
	printf("check1\n");
	for (int i = 0; i < T->rows; ++i)
		for (int j = 0; j < T->cols; ++j) T->M[i][j] = M->M[j][i];
	return T;
}

/*gets a column as a 1-D array (useful for operations with columns)*/
Vector* get_col(Matrix *matrix, int col_idx)
{
	if(col_idx >= matrix->cols)
	{
		printf("Error in get_col(M,i): index of column is out of bounds.\n");
		return NULL;
	}
	double* col_pt = new double[matrix->rows];
	for (int i = 0; i < matrix->rows; ++i) col_pt[i] = *(*(matrix->M+i)+col_idx);
	return new_vector(col_pt,matrix->rows);
}

/* sum of two Vectors, Vectors must be same size */
Vector* sum(Vector *u, Vector *v)
{
	if(u->size != v->size)
	{
		printf("Error in sum(u,v): Vector size must be the same\n");
		return NULL;
	}
	double *w_pt = new double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] += (*(u->V+i))+(*(v->V+i));
	return new_vector(w_pt,u->size);
}

/* sum of the elements of a Vector */
double sum(Vector *u)
{
	double sum = 0;
	for (int i = 0; i < u->size; ++i) sum += *(u->V+i);
	return sum;
}

/* sum of a Vector and a scalar */
Vector* sum(Vector* u, double scalar)
{
	double *u_pt = new double[u->size];
	for (int i = 0; i < u->size; ++i) u_pt[i] = u->V[i] + scalar;
	return new_vector(u_pt,u->size);
}

/* mul of a Vector and a scalar */
Vector* mul(Vector* u, double scalar)
{
	double *u_pt = new double[u->size];
	for (int i = 0; i < u->size; ++i) u_pt[i] = u->V[i] * scalar;
	return new_vector(u_pt,u->size);
}

/* elemtn wise product of two Vectors u,v */
Vector* mul(Vector* u, Vector* v)
{
	if(u->size != v->size)
	{
		printf("Error in mul(u,v): Vector sizes must be the same\n");
		return NULL;
	}
	double *w_pt = new double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] = (*(u->V+i)) * (*(v->V+i));
	return new_vector(w_pt,u->size);
}

/* product the elements of a Vector */
double mul(Vector *u)
{
	double mul = 1;
	for (int i = 0; i < u->size; ++i) mul *= *(u->V+i);
	return mul;
}

/* dot product of two Vectors u,v*/
double dot(Vector *u, Vector *v)
{
	if(u->size != v->size)
	{
		printf("Error in dot(u,v): Vector sizes must be the same\n");
		return -.10101;
	}
	return sum(mul(u,v));
}

/* product of matrices C(n,p) = A(n,m)*B(m,p) */
Matrix* mul(Matrix *A, Matrix *B)
{
	if (A->cols != B->rows)
	{ 
		printf("Error in mul(A,B): Matrix A cols must be same size as B rows\n"); 
		return NULL; 
	}
	Matrix *C = new_matrix(A->rows, B->cols);
	for (int i = 0; i < A->rows; ++i)
		for (int j = 0; j < B->cols; ++j)
			C->M[i][j] = dot(new_vector(*(A->M+i),A->cols),get_col(B,j));
	return C;
}

/* product of matrix M and a Vector u*/
Vector* mul(Matrix *A, Vector *u)
{
	if (A->cols != u->size)
	{ 
		printf("Error in mul(A,v): Matrix A cols must be same size as u rows\n"); 
		return NULL; 
	}
	double* v_pt = new double[A->rows];
	for (int i = 0; i < A->rows; ++i)
		v_pt[i] = dot(new_vector(*(A->M+i),A->cols),u);
	return new_vector(v_pt,A->rows);
}
Vector* mul(Vector *u, Matrix *A)
{
	if (A->rows != u->size)
	{ 
		printf("Error in mul(v,A): Matrix A rows must be same size as u rows\n"); 
		return NULL; 
	}
	return mul(transpose(A),u);
}

/* prints a Vector */
void print(Vector *u)
{
	printf("[ ");
	for (int i = 0; i < u->size; ++i) printf("%f ", *(u->V+i));
	printf(" ]\n");
}

/* prints a matrix */
void print(Matrix *A)
{
	for (int i = 0; i < A->rows; ++i) print(new_vector(*(A->M+i),A->cols));
}


/*
reads a csv file and return data matrix **M
*/
double** read_csv(char *filename, char separator)
{

}


int main(int argc, char const *argv[])
{
	int v_size = 3;
	int rows = 3;
	int cols = 3;
	Vector *u = new_vector(v_size-1);
	Vector *v = new_vector(v_size);
	Matrix *A = new_matrix(2,3);
	Matrix *B = new_matrix(3,2);
	Matrix *I = new_I_matrix(3);
	A->M[0][0] = 2;
	A->M[0][1] = 2;
	A->M[0][2] = 2;
	A->M[1][0] = 1;
	A->M[1][1] = 0;
	A->M[1][2] = 1;
	
	B->M[0][0] = 2;
	B->M[0][1] = 2;
	B->M[1][0] = 1;
	B->M[1][1] = 0;
	B->M[2][0] = 1;
	B->M[2][1] = 1;

	u->V[0] = 3;
	u->V[1] = 2;
	v->V[0] = 1;
	v->V[1] = 2;
	v->V[2] = 3;
	printf("\nVector v:\n");
	print(v);
	printf("\nV:Vector u:\n");
	print(u);
	printf("\nMatrix A:\n");
	print(A);
	printf("\nmul(A,v):\n");
	print(mul(A,v));
	printf("\nmul(B,u):\n");
	print(mul(B,u));
	printf("\nmul(u,A):\n");
	print(mul(u,A));
	printf("\nmul(v,B):\n");
	print(mul(v,B));


	
	return 0;
}

