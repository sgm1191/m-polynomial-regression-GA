#include <stdio.h>
#include <stdlib.h>

typedef struct vect{
	double* v;
	int size;
} vect;

typedef struct matrix{
	double** M;
	int rows;
	int cols;
} matrix;

double** new_matrix(int rows, int cols)
{
	double **M = (double**) malloc(sizeof(double*)*rows);
	for (int i = 0; i < rows; ++i) M[i] =  (double*) malloc(sizeof(double)*cols);
	return M;
}

double** new_I_matrix(int size)
{
	double **M = new_matrix(size,size);
	for (int i = 0; i < size; ++i) M[i][i] = 1;
	return M;
}

/*gets a column as a 1-D array (useful for operations with columns)*/
double* get_col(double** M, int col_idx, int size)
{
	double* col = new double[size];
	for (int i = 0; i < size; ++i) col[i] = *(*(M+i)+col_idx);
	return col;
}

double* sum(double *u, double *v, int size)
{
	double *w = new double[size];
	for (int i = 0; i < size; ++i,u++,v++) w[i] = *u + *v;
	return w;
}

double sum(double *u, int size)
{
	double sum = 0;
	for (int i = 0; i < size; ++i, u++) sum += *u;
	return sum;
}

double* mul(double *u, double *v, int size)
{
	double *w = new double[size];
	for (int i = 0; i < size; ++i,u++,v++) w[i] = (*u)*(*v);
	return w;
}

double mul(double *u, int size)
{
	double mul = 1;
	for (int i = 0; i < size; ++i, u++) mul *= *u;
	return mul;
}

double dot(double *u, double *v, int size)
{
	return sum(mul(u,v,size), size);
}

/* multiplication of matrices C(n,p) = A(n,m)*B(m,p) */
double** mul(double **A, int rows_a, int cols_a, double **B, int rows_b, int cols_b)
{
	if (cols_a != rows_b)
	{ 
		printf("Error in mul(A,B)!: Matrix A cols must be same size as B rows\n"); 
		return NULL; 
	}
	double** C = new_matrix(rows_a, cols_b);
	for (int i = 0; i < rows_a; ++i,A++)
		for (int j = 0; j < cols_b; ++j)
			C[i][j] = dot(*A,get_col(B,j,rows_b),rows_b);
	return C;
}

/* prints a vector */
void print(double *u, int size)
{
	printf("[ ");
	for (int i = 0; i < size; ++i, u++) printf("%f ", *u);
	printf(" ]\n");
}

/* prints a matrix */
void print(double **M, int r, int c)
{
	for (int i = 0; i < r; ++i, M++)
	{
		print(*M,c);
	}
}


/*
reads a csv file and return data matrix **M
*/
double** read_csv(char *filename, char separator)
{

}


int main(int argc, char const *argv[])
{
	int v_size = 5;
	int rows = 3;
	int cols = 3;

	vect a;
	a.size = v_size;
	double *u = (double*)malloc(sizeof(double)*v_size);
	double *v = new double[v_size];
	double **A = new_matrix(2,3);
	double **B = new_matrix(6,2);
	double **Ix = new_I_matrix(3);
	
	A[0][0] = 3;
	A[0][1] = 1;
	A[0][2] = 2;
	A[1][0] = 2.4;
	A[1][1] = 3.1;
	A[1][2] = 6.8;

	B[0][0] = 2;
	B[0][1] = 3;
	B[1][0] = 1;
	B[1][1] = 2;
	B[2][0] = 5;
	B[2][1] = 4;
	B[3][0] = 1;
	B[3][1] = 1;
	B[4][0] = 1;
	B[4][1] = 0;
	B[5][0] = 0;
	B[5][1] = 1;

	u[0] = 1;
	u[1] = 2;
	u[2] = 3;
	u[3] = 4;
	u[4] = 5;

	v[0] = 2;
	v[1] = 2;
	v[2] = 2;
	v[3] = 2;
	v[4] = 2;

	a.v = v;

	printf("vector B(6,2) * A(2,3):\n");
	double** C = mul(B,6,2,A,2,3);
	C = mul(C,6,3,Ix,3,3);
	printf("matrix Ix:\n"); 
	print(Ix,3,3);
	printf("matrix B:\n"); 
	print(B,6,2);
	printf("matrix C:\n"); 
	print(C,6,3);
	double *w = sum(a.v,u,v_size);
	printf("\nsuma de v + u:\n");
	print(w,v_size);
	printf("suma interna de u:\n%f",sum(u,v_size));
	printf("\nproducto punto de u dot v \n%f\n",dot(v,u,v_size));




	return 0;
}
