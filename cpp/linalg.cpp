#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*########################### GLOBAL VARIABLES ############################*/

srand(time(NULL));

/*########################### OTHER OPERATORS ############################*/

/* char array must end with '\0' special character */
int str2int(char *num)
{
	int res = 0;
	int pw = strlen(num);
	while(pw--) res += (int(*(num++))-48)*pow(10,pw);
	return res;
}

/* char array must end with '\0' special character */
long double str2double(char *num)
{
	long double res = 0;
	int pw = 0;
	int sign = 1;
	if(strchr(num, '.') == NULL) return (long double) str2int(num);
	if (*num == '-'){ sign = -1; num++; }
	for (int i = 0; *(num+i) != '.'; ++i, pw++);
	while(*(num)!='.') res+= (int(*(num++))-48)*pow(10,--pw);
	num++;
	while(*(num)>='0' && *(num)<='9') res += (int(*(num++))-48)*pow(10,--pw);
	return sign*res;
}

/* get sign of a number */
int sign(long double number)
{
	if(number < 0) return -1;
	return 1;
}

/* calculates logarithm in base 2 of a number */
double log2( double n )
{  
    return log( n ) / log( 2 );  
} 

/* generates a random double in (0,1] */
double rand01()
{
	return (double)rand() / (double)RAND_MAX;
}

/* generates a random integer in [a,b] */
int randint(int min, int max)
{
	return rand() % (max + 1 - min) + min;
}

/*#######################  LINEAR ALGEBRA MODULE  ############################*/

typedef struct{
	long double* V;
	int size;
}Vector;

typedef struct{
	long double** M;
	int rows;
	int cols;
}Matrix;

Vector* new_vector(long double *v_pt, int size)
{
	Vector *v = (Vector*)malloc(sizeof(Vector));
	v->V = v_pt;
	v->size = size;
	return v;
}
/* gets new vector from another one, it may be a subsection */
Vector* new_vector(long double *v_pt, int ini, int fin)
{
	Vector *v = (Vector*)malloc(sizeof(Vector));
	v->V = (long double*)malloc(sizeof(long double)*(fin-ini));
	int j = 0;
	for (int i = ini; i < fin; i++,j++) *(v->V+j) = *(v_pt+i);
	v->size = fin-ini;
	return v;
}

Vector* new_vector(int size)
{
	Vector *v = (Vector*)malloc(sizeof(Vector));
	v->V = (long double*)malloc(sizeof(long double)*size);//new long double[size];
	v->size = size;
	return v;

}

Matrix* new_matrix(int rows, int cols)
{
	Matrix *matrix = (Matrix*)malloc(sizeof(Matrix));
	long double **A = (long double**) malloc(sizeof(long double*)*rows);
	for (int i = 0; i < rows; ++i) A[i] =  (long double*) malloc(sizeof(long double)*cols);
	matrix->M = A;
	matrix->rows = rows;
	matrix->cols = cols;
	return matrix;
}

int** new_Ix(int size)//Ix is squared
{
	int **A = (int**) malloc(sizeof(int*)*size);
	for (int i = 0; i < rows; ++i) *(A+i) =  (int*) calloc(size,sizeof(int));
	for (int i = 0; i < size; ++i) *(*(A+i)+i) = 1;
	return A;
}

long double** new_array_f(int rows, int cols)
{
	long double **A = (long double**) malloc(sizeof(long double*)*rows);
	for (int i = 0; i < rows; ++i) *(A+i) =  (long double*) calloc(sizeof(long double)*cols);
	return A;
}

int** new_array_d(int rows, int cols)
{
	int **A = (int**) malloc(sizeof(int*)*rows);
	for (int i = 0; i < rows; ++i) *(A+i) =  (int*) calloc(sizeof(int)*cols);
	return A;
}

long double* new_array_f(int size)
{
	long double *u = (long double*) calloc(sizeof(long double*)*size);
	return u;
}


// void destroy(Vector *u)
// {
// 	free(u->V);
// }

/* print shape */
void print_shape(Matrix *A, char* label)
{
	printf("shape of %s= (%d,%d)\n", label, A->rows, A->cols);
}

void print_shape(Vector *v, char* label)
{
	printf("shape of %s = (%d, )\n",label, v->size);
}

Matrix* copy(Matrix *A)
{
	Matrix *M = new_matrix(A->rows,A->cols);
	for (int i = 0; i < A->rows; ++i)
	{
		for (int j = 0; j < A->cols; ++j)
		{
			*(*(M->M+i)+j) = *(*(A->M+i)+j);
		}
	}
	return M;
}

long double** transpose(long double **M,int mr,int mc)
{
	for (int i = 0; i < mr; ++i)
		for (int j = 0; j < mc; ++j) 
		{
			long double temp = *(*(M+i)+j);
			*(*(M+i)+j) = *(*(M+j)+i);
			*(*(M+j)+i) = temp;
		}
	return M;
}
/* gets a subsection of the vector */
Vector* sub(Vector *v, int ini, int fin)
{
	if (ini == 0 && fin == v->size) return v;
	if ( ini < 0 || fin > v->size)
	{
		printf("ERROR in sub(Vector, int, int): indices must be inside of vector bounds\n");
		return NULL;
	}
	if (ini >= fin)
	{
		printf("ERROR in sub(Vector, int, int): ini index must be lower than fin index \n");
		return NULL;
	}
	int len = fin-ini;
	long double *v_pt = new long double[len];
	for (int i = 0; i < len; i++, ini++) *(v_pt+i) = *(v->V+ini);
	return new_vector(v_pt,len);
}

/* gets a subsection of a matrix
** starting at indices _x to _y-1 */
Matrix* sub(Matrix *A, int i_x, int f_x, int i_y, int f_y)
{
	if(i_x == 0 && i_y == 0 && f_x == A->rows && f_y == A->cols) return A;
	if ( i_x < 0 || f_x > A->rows || i_y < 0 || f_y > A->cols )
	{
		printf("ERROR in sub(Matrix, int, int, int, int): indices must be inside of matrix bounds\n");
		return NULL;
	}
	if (i_x >= f_x || i_y >= f_y)
	{
		printf("ERROR in sub(Matrix, int, int, int, int): init indices must be lower than end indices \n");
		return NULL;
	}
	Matrix *B = new_matrix((f_x - i_x),(f_y - i_y));
	int ib = 0;
	for (int i = i_x; i < f_x; i++,ib++)
	{
		int jb = 0;
		for (int j = i_y; j < f_y; j++,jb++)
			*(*(B->M+ib)+jb) = *(*(A->M+i)+j);
	}
	return B;
}

Matrix* append_c(Matrix *A, Vector *b)
{
	if (A->rows != b->size)
	{
		printf("Error in append_c(A,b): rows of A and size of b must be the same. \n");
	}
	Matrix *B = new_matrix(A->rows,A->cols+1);
	for (int r_i = 0; r_i < A->rows; ++r_i)
	{
		for (int c_i = 0; c_i < A->cols; ++c_i)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i);
		*(*(B->M+r_i)+B->cols-1) = *(b->V+r_i);
	}
	return B;
}

Matrix* append_r(Matrix *A, Vector *b)
{
	if (A->cols != b->size)
	{
		printf("Error in append_r(A,b): cols of A and size of b must be the same. \n");
	}
	Matrix *B = new_matrix(A->rows+1,A->cols);
	for (int c_i = 0; c_i < A->cols; ++c_i)
	{
		for (int r_i = 0; r_i < A->rows; ++r_i)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i);
		*(*(B->M + B->rows-1) + c_i) = *(b->V+c_i);
	}
	return B;
}

/* appends a vector at the end of a matrix 
   in the axis given(0 for rows, 1 for cols)*/
Matrix* append(Matrix *A, Vector *b, int axis)
{
	if (axis==0) return append_r(A,b);
	return append_c(A,b);
}

/* appends at the end of the vector the value given */
Vector* append(Vector *u, long double value)
{
	Vector *v = new_vector(u->size+1);
	int i = 0;
	for (i; i < u->size; i++) *(v->V+i) = *(u->V+i);
	*(v->V+i) = value;
	return v;
}

Vector* insert(Vector *u, int index, long double value)
{
	Vector *v = new_vector(u->size+1);
	int i = 0;
	for (i; i < index; i++) *(v->V+i) = *(u->V+i);
	*(v->V+i) = value;
	i++;
	for (i; i < v->size; i++) *(v->V+i) = *(u->V+i-1);
	return v;
}

Matrix* insert_r(Matrix *A, int index, Vector *b)
{
	if(index > A->rows)
	{
		printf("Error in insert_r(A,index,b): index must be lower than number of rows of matrix A\n");
		return NULL;
	}
	Matrix *B = new_matrix(A->rows+1,A->cols);
	int r_i = 0;
	for (r_i; r_i < index; r_i++)
		for (int c_i = 0; c_i < A->cols; ++c_i)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i);

	for (int c_i = 0; c_i < A->cols; ++c_i)
		*(*(B->M + r_i) + c_i) = *(b->V+c_i);
	r_i++;
	for (r_i; r_i-1 < A->rows; r_i++)
		for (int c_i = 0; c_i < A->cols; ++c_i)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i-1)+c_i);
	return B;
}

Matrix* insert_c(Matrix *A, int index, Vector *b)
{
	if(index > A->cols)
	{
		printf("Error in insert_c(A,index,b,axis=1): index must be lower than number of cols of matrix A\n");
		return NULL;
	}
	Matrix *B = new_matrix(A->rows,A->cols+1);
	int c_i = 0;
	for (c_i; c_i < index; c_i++)
		for (int r_i = 0; r_i < A->rows; r_i++)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i);

	for (int r_i = 0; r_i < A->rows; r_i++)
		*(*(B->M + r_i) + c_i) = *(b->V+r_i);
	c_i++;

	for (c_i; c_i < B->cols; c_i++)
		for (int r_i = 0; r_i < A->rows; r_i++)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i-1);
	return B;
}
/* inserts a vector in the axis given at the index given */
Matrix* insert(Matrix *A, int index, Vector *b, int axis)
{
	if (axis==0)
	{
		return insert_r(A,index,b);
	}
	return insert_c(A,index,b);
}

Vector* remove(Vector *u, int index)
{
	if(index >= u->size)
	{
		printf("Error in remove(u,index): index must be lower than size of vector u.\n");
		return NULL;
	}
	Vector *v = new_vector(u->size-1);
	int i = 0;
	for (i; i < index; ++i) *(v->V+i) = *(u->V+i);
	i++;
	for (i; i-1 < v->size; ++i) *(v->V+i-1) = *(u->V+i);
	return v;
}

Matrix* remove_r(Matrix *A, int index)
{
	if(index >= A->rows)
	{
		printf("Error in remove_r(A,index,b): index must be lower than number of rows of matrix A\n");
		return NULL;
	}
	Matrix *B = new_matrix(A->rows-1,A->cols);
	int r_i = 0;
	for (r_i; r_i < index; r_i++)
		for (int c_i = 0; c_i < A->cols; ++c_i)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i);
	r_i++;
	for (r_i; r_i-1 < B->rows; r_i++)
		for (int c_i = 0; c_i < A->cols; ++c_i)
			*(*(B->M+r_i-1)+c_i) = *(*(A->M+r_i)+c_i);
	return B;
}

Matrix* remove_c(Matrix *A, int index)
{
	if(index > A->cols)
	{
		printf("Error in remove_c(A,index,b,axis=1): index must be lower than number of cols of matrix A\n");
		return NULL;
	}
	Matrix *B = new_matrix(A->rows,A->cols-1);
	int c_i = 0;
	for (c_i; c_i < index; c_i++)
		for (int r_i = 0; r_i < A->rows; r_i++)
			*(*(B->M+r_i)+c_i) = *(*(A->M+r_i)+c_i);
	c_i++;
	for (c_i; c_i-1 < B->cols; c_i++)
		for (int r_i = 0; r_i < A->rows; r_i++)
			*(*(B->M+r_i)+c_i-1) = *(*(A->M+r_i)+c_i);
	return B;
}
/* removes a column or row in the axis given at the index given */
Matrix* remove(Matrix *A, int index, int axis)
{
	if (axis==0)
	{
		return remove_r(A,index);
	}
	return remove_c(A,index);
}

/*gets a column as a 1-D array (useful for operations with columns)*/
Vector* get_col(Matrix *matrix, int col_idx)
{
	if(col_idx >= matrix->cols)
	{
		printf("Error in get_col(M,i): index of column is out of bounds.\n");
		return NULL;
	}
	long double* col_pt = new long double[matrix->rows];
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
	long double *w_pt = new long double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] += (*(u->V+i))+(*(v->V+i));
	return new_vector(w_pt,u->size);
}

/* sum of the elements of a Vector */
long double sum(Vector *u)
{
	long double sum = 0;
	for (int i = 0; i < u->size; ++i) sum += *(u->V+i);
	return sum;
}

/* sum of a Vector and a scalar */
Vector* sum(Vector* u, long double scalar)
{
	long double *u_pt = new long double[u->size];
	for (int i = 0; i < u->size; ++i) u_pt[i] = u->V[i] + scalar;
	return new_vector(u_pt,u->size);
}

/* mul of a Vector and a scalar */
Vector* mul(Vector* u, long double scalar)
{
	long double *u_pt = new long double[u->size];
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
	long double *w_pt = new long double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] = (*(u->V+i)) * (*(v->V+i));
	return new_vector(w_pt,u->size);
}

/* elemtn wise division of two Vectors u,v */
Vector* div(Vector* u, Vector* v)
{
	if(u->size != v->size)
	{
		printf("Error in div(u,v): Vector sizes must be the same\n");
		return NULL;
	}
	long double *w_pt = new long double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] = (*(u->V+i)) / (*(v->V+i));
	return new_vector(w_pt,u->size);
}

/* element-wise avsolute of a vector */
Vector* absolute(Vector *u) 
{
	long double *w_pt = new long double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] = abs(*(u->V+i));
	return new_vector(w_pt,u->size);
}

/* get only the signs of the elemesnt of a vector */
Vector* sign(Vector *u) 
{
	long double *w_pt = new long double[u->size];
	for (int i = 0; i < u->size; ++i) w_pt[i] = sign(*(u->V+i));
	return new_vector(w_pt,u->size);
}

/* product the elements of a Vector */
long double mul(Vector *u)
{
	long double mul = 1;
	for (int i = 0; i < u->size; ++i) mul *= *(u->V+i);
	return mul;
}

/* dot product of two Vectors u,v*/
long double dot(long double *u, int ul, long double *v, int vl)
{
	if(ul != vl)
	{
		printf("Error E05 in dot(u,v): Vector sizes must be the same\n");
		return -.10101;
	}
	long double sum = 0;
	for (int i = 0; i < ul; i++)
		sum += (*(u+i)) * (*(v+i));
	return sum;
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
long double* mul(long double **A, int Ar, int Ac, long double *u, int ul)
{
	if (Ac != ul)
	{ 
		printf("Error E04 in mul(A,v): Matrix A cols must be same size as u rows\n"); 
		return NULL; 
	}
	long double* v_pt = new_array_f(Ar);
	for (int i = 0; i < Ar; ++i)
		*(v_pt+i) = dot(*(A+i),Ac,u,ul);
	return v_pt;
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

/* swap rows of a matrix */
void swap_r(Matrix *&M, int row_a, int row_b)
{
	for (int col_i = 0; col_i < M->cols; ++col_i)
	{
		long double temp = *(*(M->M+row_a)+col_i);
		*(*(M->M+row_a)+col_i) = *(*(M->M+row_b)+col_i);
		*(*(M->M+row_b)+col_i) = temp;
	}
}
/* swap columns of a matrix */
void swap_c(Matrix *&M, int col_a, int col_b)
{
	for (int row_i = 0; row_i < M->rows; ++row_i)
	{
		long double temp = *(*(M->M+row_i)+col_a);
		*(*(M->M+row_i)+col_a) = *(*(M->M+row_i)+col_b);
		*(*(M->M+row_i)+col_b) = temp;
	}
}
/* swap two columns or rows of a vector */
void swap(Matrix *&M, int a, int b, int axis)
{
	if (axis==0) swap_r(M, a, b);
	else swap_c(M, a, b);
}
/* swap  2 elements of a vector */
void swap(Vector* u, int a, int b)
{
	long double temp = *(u->V+a);
	*(u->V+a) = *(u->V+b);
	*(u->V+b) = temp;
}
/* swap 2 rows from different matrices */
void swap(Matrix *&A, Matrix *&B, int row_a, int row_b)
{
	if (A->cols != B->cols)
	{
		printf("Error in swap(A,B,row_a,row_b): matrices A and B must have same number of columns \n");
		return;
	}
	for (int col = 0; col < A->cols; ++col)
	{
		long double temp = *(*(A->M+row_a)+col);
		*(*(A->M+row_a)+col) = *(*(B->M+row_b)+col);
		*(*(B->M+row_b)+col) = temp;
	}
}

/* gets a random sample of a matrix A.
   inplace = substracts the sample from the original */
Matrix* sample(Matrix *&A, int samp_size, bool inplace)
{
	//srand(90418); // semilla random
	srand(time(NULL));
	Matrix *S = new_matrix(samp_size, A->cols);
	Matrix *new_A = new_matrix(A->rows-samp_size, A->cols);
	int iA = 0; // index for original data matrix
	int iS = 0; // index for sample matrix
	int inA = 0; // index for remaining data matrix
	int steps = int(A->rows/samp_size);
	while (samp_size)
	{
		int samp_v = iA + rand() % steps;
		while(iA < samp_v) *(new_A->M + inA++) = *(A->M + iA++);
		samp_size--;
		*(S->M + iS++) = *(A->M + iA++);
	}
	while(iA < A->rows) *(new_A->M + inA++) = *(A->M + iA++);;
	if (inplace) A = new_A;
	return S;
}

/* get first n rows of a matrix substracting it from the original
** The array A is the original one and S is the sample array where the first elements
** are going to be substracted */
void first(long double **&A,int &ar,int ac, long double **&S, int n)
{
	for (int i = 0; i < n; i++)
			*(*(S+ i)+ j) = *(*(A+ i)+ j);
	ar -= n;
	long double **residual = new_array_f(ar, ac);
	for (int i = 0; i < ar; i++)
		for (int j = 0; j < ac; j++)
			*(*(residual+i)+j) = *(*(A+n+i)+j);
	A = residual;
}

long double min(Vector* u, int &index)
{
	long double min = *(u->V);
	for (int i = 1; i < u->size; ++i)
	{
		if (*(u->V+i) > min)
		{
			min = *(u->V+i);
			index = i;
		}
	}
	return min;
}

long double max(Vector* u, int &index)
{
	long double max = *(u->V);
	for (int i = 1; i < u->size; ++i)
	{
		if (*(u->V+i) > max)
		{
			max = *(u->V+i);
			index = i;
		}
	}
	return max;
}

/* prints a Vector */
void print(long double *u, int len)
{
	printf("[ ");
	for (int i = 0; i < len; ++i) printf("%10.12Lf ", *(u+i));
	printf("]\n");
}

/* prints a matrix */
void print(long double **A, int ar, int ac)
{
	for (int i = 0; i < ar; i++) print(*(A+i),ac);
}

/* returns Vector c of coefficients */
void solve(long double **XY,int xyr,int xyc, long double *&c)
{
	int m = xyr;

	/* loop for the generation of upper triangular matrix*/
	for(int j=0; j<m; j++) 
	{
		for(int i=j+1; i<m; i++)
		{
			long double c= *(*(XY+i)+j)/(*(*(XY+j)+j));
			for(int k=0; k<m+1; k++)
				*(*(XY+i)+k) = *(*(XY+i)+k)-c*(*(*(XY+j)+k));
		}
	}
	
	*(x->V+m-1) = *(*(XY+m-1)+m)/(*(*(XY+m-1)+m-1));

	/* this loop is for backward substitution*/
	for(int i=m-2; i>=0; i--)
	{
		long double sum=0;
		for(int j=i+1; j<m; j++) sum += *(*(XY+i)+j) * (*(x->V+j));
		*(x+i) = (*(*(XY+i)+m)-sum)/(*(*(XY+i)+i));
	}
}

long double** inverse(long double **A, int A_r, int A_c)
{
	int dim = A_c;
	long double **Ix = (long double**)new_Ix(dim);

	// scale matrices
	for (int i = 0; i < dim; i++)
	{
		long double r_max = fabsl(*(*(A+i)));
		for (int j = 1; j < dim; j++)
			r_max = fmaxl(fabsl(*(*(A+i)+j)), r_max);
		if (r_max == 0)
			printf("Error E03 in inverse(A): can't get inverse, unstable matrix.\n");
		long double scale = 1/r_max;
		for (int j = 0; j < dim; j++)
		{
			*(*(A+i)+j) *= scale;
			if (i==j) *(*(Ix+i)+j) = scale;
		}
	}

	// put largest element in pivot position
	int ipiv;
	for (int k = 0; k < dim-1; k++)
	{
		long double temp;
		long double big = 0;
		for (int i = k; i < dim; i++)
		{
			temp = fabsl(*(*(A+ i)+ k));
			if (big < temp)
			{
				big = temp;
				ipiv = i;
			}
		}
		if (big == 0)
			printf("Error E03 in inverse(A): can't get inverse, unstable matrix.\n");
		if (ipiv != k)
		{
			swap(X,ipiv,k,0);
			swap(Ix,ipiv,k,0);
		}

		// eliminate X(k) from equations k+1, k+2, ..., k+dim
		long double quot;
		for (int i = k+1; i < dim; i++)
		{
			quot = *(*(A+ i)+ k) / (*(*(A+ k)+ k));

			for (int j = k+1; j < dim; j++)
				*(*(A+ i)+ j) -= (*(*(A+ k)+ j)) * quot;

			for (int j = 0; j < dim; j++)
				*(*(Ix+ i)+ j) -= (*(*(Ix+ k)+ j)) * quot;
		}
	}
	if (*(*(A+dim-1)+dim-1) == 0)
			printf("Error E03 in inverse(A): can't get inverse, unstable matrix.\n");
	// back substitution
	for (int l = 0; l < dim; l++)
	{
		*(*(Ix+dim-1)+l) /= (*(*(A+dim-1)+dim-1));
		for (int i = dim-2; i >= 0; i--)
		{
			long double sum = 0;
			for (int j = i+1; j < dim; j++)
				sum += *(*(A+i)+j) * (*(*(Ix+j)+l));
			*(*(Ix+i)+l) = (*(*(Ix+i)+l) - sum)/(*(*(A+i)+i));
		}
	}
	return Ix;

}



/* replicates the matrix x times in x axis, and y times in y axis
** NOTE: it shouldn't enter a 0 in any of the inputs */
Matrix* tile(Matrix *A, int x, int y)
{
	int c = A->cols;
	int r = A->rows;
	Matrix *B = new_matrix(A->rows*x, A->cols*y);
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			// loop for matrix A
			int stp_x = i*r;
			int stp_y = j*c;
			for (int ir = 0; ir < r; ir++)
			{
				for (int ic = 0; ic < c; ic++)
				{
					*(*(B->M+ (ir+stp_x))+ (ic+stp_y)) = *(*(A->M+ ir)+ ic);
				}
			}
		}
	}
	return B;
}
/* replicates the matrix t times
** NOTE: it shouldn't enter a 0 in any of the inputs */
Vector* repeat(Vector *u, int t)
{
	Vector *v = new_vector(u->size*t);
	for (int i = 0; i < t; i++)
		for (int j = 0; j < u->size; j++)
			*(v->V +j+(u->size*i)) = *(u->V+j);
	return v;
}

/*#######################  FILE I/O  #########################*/

/*
reads a csv file and return data matrix M
note: all numbers must be in float format
for example, instead of 1 it must be writen as 1.0
*/
long double** read_csv(char *filename, char separator, int rows, int fields)
{
	FILE* stream = fopen(filename, "r");
	if(stream == NULL)
	{
		printf("Error E00 in read_csv(filename,sep,rows,fields): file does not exists.\n");
	}
	long double **data = new_array_f(rows, fields);
	int c = fgetc(stream);
	int row = 0;
	int col = 0;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < fields; ++j)
		{
			char *num = (char*)malloc(30);
			int i_2 = 0;
			while(c != separator && c != EOF && c != '\n')
			{
				*(num + i_2++) = c;
				c = fgetc(stream);
			}
			*(num+i_2) = '\0';
			*(*(data+i)+j) = str2double(num);
			if (c == EOF) break;
			c = fgetc(stream);
			free(num);
		}
		if (c == EOF) break;
	}
	fclose(stream);
	return data;
}

/*###################### ASCEND ALGORITHM ########################*/

/* solve system of linear equations for finding the minimax signs the first time 
** note: not implemented, instead solve is used, it is not implemented yet because
** it has bugs. (°-°)" 
*/
Vector* lassol(Matrix *XY_cap)
{
	Matrix *XY = copy(XY_cap);
	Vector *lambda = new_vector(XY->rows);
	int dim = XY->rows;
	int dimp1 = XY->cols;
	// scale each row to its max elements
	for (int i = 0; i < dim; i++)
	{
		long double r_max = fabsl(*(*(XY->M+i)));
		for (int j = 1; j < dim; j++)
			r_max = fmaxl(fabsl(*(*(XY->M+i)+j)), r_max);
		if (r_max == 0)
		{
			printf("Error in lassol(A): can't solve, unstable system.\n");
			return NULL;
		}
		long double scale = 1/r_max;
		for (int j = 0; j < dim; j++)
			*(*(XY->M+ i)+ j) *= scale;
	}

	// put largest element in column i
	int ipiv;
	for (int k = 0; k < dim-1; k++)
	{
		long double temp;
		long double big = 0;
		for (int i = k; i < dim; i++)
		{
			temp = fabsl(*(*(XY->M+ i)+ k));
			if (big < temp)
			{
				big = temp;
				ipiv = i;
			}
		}
		if (big == 0)
		{
			printf("Error in lassol(A): can't solve, unstable system.\n");
			return NULL;
		}
		// exchange column with the largest element
		if (ipiv != k) swap(XY,ipiv,k,0);

		// eliminate all in column except first
		long double quot;
		for (int i = k+1; i < dim; i++)
		{
			quot = *(*(XY->M+ i)+ k) / (*(*(XY->M+ k)+ k));
			for (int j = k+1; j < dim; j++)
				*(*(XY->M+ i)+ j) -= (*(*(XY->M+ k)+ j)) * quot;
		}
	}
	
	if (*(*(XY->M+dim-1)+dim-1) == 0)
	{
		printf("Error in lassol(A): can't solve, unstable system.\n");
		return NULL;
	}
	// back substitution
	*(lambda->V+dim-1) = (*(*(XY->M+dim-1)+dim))/(*(*(XY->M+dim-1)+dim-1));
	for (int i = dim-2; i >= 0; i--)
	{
		long double sum = 0;
		for (int j = i+1; j < dim; j++)
			sum += *(*(XY->M+i)+j) * (*(lambda->V+j));
		*(lambda->V+i) = (*(*(XY->M+i)+dim) - sum)/(*(*(XY->M+i)+i));
	}
	return lambda;
}

/* get minimax signs by te 4th method of the minimax theory document 
** and appends the signs found to the inner set and asigns it to inner_*/
void get_signs(long double **&inner,int ir,int ic, long double **&inner_)
{
	long double **innT = transpose(inner);
	long double *signs = new_array_f(ir);
	*(signs+ir-1) = -1;
	signs = solve(innT,ic,ir);
	free(innT);
	
	// append signs to inner_ first column
	for (int i = 0; i < ir; i++) *(*(inner_+i)) = sign(*(signs+i));
	for (int i = 0; i < ir; i++)
		for (int j = 0; j < ic j++)
			*(*(inner_+i)+j+1) = *(*(inner+i)+j);
	// printf("Minimax signs:\n");
	// print(signs);
}

/* calculates all the combinations of the degrees of the variables */
Matrix* get_terms(Vector *deg)
{
	int comb = 1;
	int temp = 1;
	for(int i = 0; i < deg->size; i++) *(deg->V+i) = *(deg->V+i) + 1;
	for(int i = 0; i < deg->size; i++) comb *= *(deg->V+i);
	Matrix *terms = new_matrix(comb, deg->size);
	for (int ci = terms->cols-1; ci >=0 ; ci--)
	{
		for (int ri = 0; ri < terms->rows; ri++)
		{
			int max = *(deg->V+ci);
			for (int di=0; di < max; di++,ri++)
			{
				for (int k = 0; k < temp; k++,ri++)
				{
					*(*(terms->M+ri)+ci) = di;
				}
				ri--;
			}
			ri--;
		}
		temp *= *(deg->V+ci);
	}
	return terms;
}

/* maps the data into new dataset with the degrees of the variables given 
** input:
**    data: array of data of dimensions (dr,dc)
**    dr: data rows
**    dc: data columns
**    terms: matrix with the following format.
** 		if there are originally 5 variables
** 		each array is a term which 
** 		has the exponent for each variable.
** 		example:
**
** 		       v1, v2, v3, v4, v5=vc
** 		t1	[ [ 0,  1,  2,  1,  3],
** 		t2	  [ 1,  2,  3,  1,  0],
** 		...	  ...
** 		tr	  [ 1,  2,  3,  4,  1] ]
** 	  tr: terms rows
** 	  tc: terms columns   
*/
void map(long double **&data,int dr,int dc, int **terms, int tr, int tc)
{
	// checkMM for memmory liberation
	if (dc-1 != tc)
	{
		printf("Error E01 in map(data,terms): data fields and terms columns must have the same size.\n");
		return NULL;
	}
	
	long double **m_data = new_array_f(dr,tr+1); // +1 for the dependent variable

	for (int i = 0; i < m_dr; i++)
	{
		for (int j = 0; j < m_dc-1; j++)
		{
			long double value = 1;
			for (int k = 0; k < tc; k++)
				value *= powl((*(*(data+i)+k)),(*(*(terms+j)+k)));

			*(*(m_data+i)+j) = value;
		}
		*(*(m_data+i) + m_dc-1) = *(*(data+i) + dc-1);
	}
	data = m_data;
}

void get_coeff(long double  **B,int Br,int Bc, long double *f, int fl, long double *c, int cl long double &eps_th)
{
	long double *C = mul(B,Br,Bc, f,fl);
	eps_th = *C;
	for (int i = 0; i < cl; ++i)
		*(c+i) = *(C+i+1);
	free(C);
	// printf("coefficients:\n");
	// print(C);
	// printf("\n");
}

/* test coefficients in the outter set.
   input
      -outter: outter set
      -coefficients: vector of coefficients(with the error removed)
      -eps_ph: refernce variable where error is stored 
      -sgn: reference variable where sign of the error is stored
      -idx: index of the vector with the maximum error */
long double test_coeff(Matrix *outter, Vector *coefficients, int &sgn, int &idx)
{	
	long double error, abs_err;
	long double eps_ph;
	eps_ph = -100000;
	for (int i = 0; i < outter->rows; ++i)
	{
		long double y_cap = 0;
		for (int j = 0; j < coefficients->size; j++)
			y_cap += (*(coefficients->V+j)) * (*(*(outter->M+i)+j));

		error = *(*(outter->M+i)+outter->cols-1) - y_cap;
		// printf("Error at index %d = %Lf\n", i, error);
		abs_err = fabsl(error);
		if(abs_err > eps_ph)
		{
			eps_ph = abs_err;
			idx = i;
			sgn = sign(error);
		}
	}
	return eps_ph;
}

/* gets new inverse with the lambdas and the index  of the maximum value
   at the inner set.
   input
   	   -B: inverse matrix of A(inner set)
   	   -betha: index of the maximum internal error
   	   -lambdas: vector of lambdas calculated in the swapping step */
void get_new_inverse(Matrix *&B, int betha, Vector *lambdas)
{
	for (int i = 0; i < B->rows; ++i)
		*(*(B->M+i)+betha) =  *(*(B->M+i)+betha) / *(lambdas->V+betha);

	for (int i = 0; i < B->rows; ++i)
		for (int j = 0; j < B->rows; ++j)
			if (i!=betha) *(*(B->M+j)+i) = *(*(B->M+j)+i) - (*(lambdas->V+i)) * (*(*(B->M+j)+betha));
}

/* stabilizes data by adding neglectable(<=10e-6) random value */
void stabilize(long double **&A, int ar, int ac, long double factor)
{
	double random;
	for (int i = 0; i < ar; i++)
	{
		for (int j = 0; j < ac; j++)
		{
			random = rand01();
			if(*(*(A+i)+j) == 0) 
				*(*(A+i)+j) = *(*(A+i)+j) + random*factor;
			else *(*(A+i)+j) = (*(*(A+i)+j)) * (1+random*factor);
		}
	}
}

/* swaps a vector from the inner set for the one in the outter set
   whit the maximum external error.
   input
   	   -outter: outter set
   	   -inner: inner set(matrix A)
   	   -B: inverse matrix of A
   	   -mu: sign of the external error 
   	   -IE: index of the maximum external error
   	output - None, the Matrices are passed by reference */
void swap_vector(Matrix *&outter, Matrix *&inner, Matrix *&B, long double mu, int IE)
{

	Vector *amp1 = insert(new_vector(*(remove(outter,outter->cols-1,1)->M+IE),outter->cols-1),0,mu);
	Vector *lambdas = mul(amp1,B);
	long double betha_max = -10000;
	long double betha;
	int bmi = -1; // betha max index
	for (int i = 0; i < B->cols; ++i)
	{
		betha = mu * (*(lambdas->V+i)/(*(*(B->M)+i)));
		if (betha > betha_max)
		{
			betha_max = betha;
			bmi = i;
		}
	}

	//printf(" Swapping %d for %d\n", bmi+1, IE+inner->cols);
	
	*(*(inner->M+bmi)) = mu;
	for (int i = 0; i < outter->cols; ++i)
	{
		long double temp = *(*(inner->M+bmi)+i+1);
		*(*(inner->M+bmi)+i+1) = *(*(outter->M+IE)+i);
		*(*(outter->M+IE)+i) = temp;
	}
	// printf("\nLambda Vector:\n");
	// print(lambdas);
	// printf("\n");

	// calculate new inverse
	get_new_inverse(B,bmi,lambdas);
}

/* ascend algorithm returns vector of coefficients C 
   input:
   		terms: matrix with the following format.
   		if there are originally 5 variables
   		each array is a term which 
   		has the exponent for each variable.
   		example:

   		       v1, v2, v3, v4, v5
   		t1	[ [ 0,  1,  2,  1,  3],
   		t2	  [ 1,  2,  3,  1,  0],
   		...	  ...
   		tn	  [ 1,  2,  3,  4,  1] ]
   	output:
   		Vector of coefficients for each term of the polynomial,
   		the first element is the minimax internal error.
   		[epsilon_theta, c1, c2, ..., cn] */
Vector* ascend(Matrix *terms, long double &eps_th, long double &eps_ph)
{
	int IE, mu;
	int m = terms->rows;
	int M = m+1;
	long double epsilon_th,epsilon_ph;
	Vector *c; // coefficient vector

	long double stab_fac = 1e-6;
	long double quasi = 0.05;
	bool Q_F = true;

	char* filename = (char*) malloc(100);
	filename = (char*)"z3Vars.dat";
	// dimensions of the original dataset
	int rows = 300; 
	int fields = 4;
	// dimensions of the mapped dataset
	int dr,dc;
	
	long double *c = new_array_f(m);


	// printf("Training filename: \n");
	// scanf("%s",filename);
	// printf("data shape(rows fields):");
	// scanf("%d %d",&rows,&fields);
	// read data
	long double **outter = read_csv(filename,'\t', rows, fields);
	free(filename);
	
	// map data
	map(outter,rows,fields, terms,tr,tc);
	dr = rows; // rows of the mapped dataset
	dc = tr+1;

	// stabilize data
	stabilize(outter,dr,dc, stab_fac);
	
	// split data
	long double **inner = new_array_f(M,dc);
	first(outter,dr,dc, inner_,M);
	
	// printf("Inner set:\n");
	// print(inner);
	// printf("\nOutter set:\n");
	// print(outter);
	
	// get minimax signs
	int Ar = M;
	int Ac = dc+1;
	long double **A = new_array_f(Ar,Ac);
	get_signs(inner,M,dc, A);
	free(inner);

	// get matrix A
	// get 1st inverse
	long double **B = inverse(A,A_r,A_c-1);
	// printf("Identity:\n");
	// print(mul(B,remove(inner,inner->cols-1,1)));
	int iteration = 1;
	char cont_flag;

	while(true)
	{
		//printf("check1\n");
		get_coeff(B,Ar,Ac, get_col(inner,inner->cols-1), epsilon_th);
		//printf("check2\n");
		// printf("coefficients:\n");
		// print(c);
		epsilon_ph = test_coeff(outter, c, mu, IE);
		//printf("\n IT[%d]: eps_th = %4.10Lf eps_ph = %4.10Lf ",iteration, epsilon_th, epsilon_ph);
		if ((epsilon_th >= epsilon_ph) || (Q_F && fabsl(epsilon_th - epsilon_ph) <= quasi))
		{
			//if (Q_F)
			//{
				//printf("quasi minimax achieved\n");
			//}
			break;
		}
		else
		{
			swap_vector(outter, inner, B, mu, IE);
			// printf(" \n");
			// print(mul(remove(inner,inner->cols-1,1), B));
		}
		// printf("Next? (y/n): ");
		// scanf("%c",&cont_flag);
		// if(cont_flag=='n') break;
		iteration++;
	}
	eps_ph = epsilon_ph;
	eps_th = epsilon_th;
	free(B);
	free(outter);
	free(inner);
	return c;
}

/*###################### GENETIC ALGORITHM(EGA) ########################*/

/* gen population N= number of individuals, L= Length of each individual,
** nt= number of terms, nv= number of variables, bid= bits in digit,
** max_deg= maximum degree */
Matrix *gen_population(int N, int L, int nt, int nv, int bid, int max_deg)
{
	Matrix *pop = new_matrix(N,L);

	for (int ii = 0; ii < N; ii++)
	{
		// naive generation of individuals(fully random)
		for (int bi = 0; bi < L; bi++)
			*(*(pop->M+ii)+bi) = round(rand01());
		// making a single individual
		//1 chose degree form a normal distribution
		//2 find the combination of the powers wich adds to chosen degree
		//3 shufle powers
		//4 if ok insert ocurrences separated by ','
		//5 concatenate valid combinations without ','
	}
	return pop;
}

/* decodes the individual into a matrix of degrees of variables 
** input:
**       individual: pointer to array of 1s and 0s 
**       bid: bytes in digit
**       nt: number of terms in polynomial
**       nv: number of variables */
Matrix* decode(long double *ind, int nt, int nv, int bid)
{
	Matrix *terms = new_matrix(nt,nv);
	int step_t = nv*bid;
	for (int t = 0; t < nt; t++) // for each term
	{
		for (int v = 0; v < nv; v++) // for each variable
		{
			// binary to decimal
			int power = 1;
			int decimal = 0;
			for (int b = (t*step_t)+(bid*(v+1))-1; b >= (t*step_t)+(bid*v); b--) // for each bit in digit
			{
				decimal += *(ind+b)*power;
				power<<=1;
			}
			*(*(terms->M+t)+v) = decimal;	
		}
	}
	return terms;
}

/* returns a vector of fitness values starting at individual ini-th
** to the (fin-1)-th individual */
void evaluate(Matrix *pop, Vector *&fit, int ini, int fin, int nt, int nv, int bid)
{
	long double eps_th,eps_ph;
	for (int ii = ini; ii < fin; ii++)
	{
		Matrix *terms = decode(*(pop->M+ii),nt,nv,bid);
		Vector *c = ascend(terms,eps_th,eps_ph);

		*(fit->V+ii) = eps_ph; // minimizing by minimax error
		// test_RMS(test_filename,coefficients) - minimizing by rms error
	}
}

/* annular crossover */
void annular_cross(Matrix *&pop, int N, int L, double Pc)
{
	int n = N-1;
	int L_2 = (int)(L/2);
	int N_2 = (int)(N/2);
	for (int i = 0; i < N_2; i++,n--)
	{
		double p = rand01();
		if (p <= Pc)
		{
			int xp = randint(0,L_2);
			// swap the middle chunk starting at xp and finishing at L_2
			for (int j = xp; j <= L_2; j++)
			{
				// crossing ind i-th with ind n-th
				long double temp = *(*(pop->M+i)+j);
				*(*(pop->M+i)+j) = *(*(pop->M+n)+j);
				*(*(pop->M+n)+j) = temp;
			}

		}

	}
}

/* uniform mutation */
void mutate(Matrix *&pop, int N, int L, int b2m)
{
	while(b2m--)
	{
		int p1 = round(rand01()*L)-1;
		int p2 = round(rand01()*N)-1;
		if(*(*(pop->M+p2)+p1) == 0) *(*(pop->M+p2)+p1) = 1;
		else *(*(pop->M+p2)+p1) = 0;
	}
}

void sort(Matrix *&pop, Vector *&fit)
{
	int n = fit->size;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (*(fit->V+j) > *(fit->V+i))
			{
				// swap fitness values
				double tmp = *(fit->V+i);
				*(fit->V+i) = *(fit->V+j);
				*(fit->V+j) = tmp;
				// swap individuals
				long double *ind = *(pop->M+i);
				*(pop->M+i) = *(pop->M+j);
				*(pop->M+j) = ind;
			}  
		}
	}
}

void print_scores(Vector *fit, int n)
{
	printf("Best %d scores:\n",n);
	for (int i = 0; i < n; i++)
	{
		printf( "%d: %Lf\t",i,*(fit->V+i));
	}
}

void run_ega(/*Vector *&best_ind, long double &best_fit*/)
{
	// parameters of the EGA

	// file definition parameters
	char* f_train = (char*)"DB24-glass/TRAIN.TXT";
	char* f_test = (char*)"DB24-glass/TEST.TXT";
	
	double Pc = 1; // crossover probability
	double Pm = .05; // mutation probability
	int gen = 100; // number of generations
	int N = 50; // number of individuals
	int max_deg = 8; // maximum degree of the variables

	// parameters inferred from the dataset
	int NV = 3; // number of independent variables
	int NT = 8; // number of terms

	// derivated parameters
	int BID = ceil(log2(max_deg));
	int L = NT*NV*BID;
	int b2m = round(L*N*Pm); // parameter for later calculations

	// generate random population
	Matrix *pop = gen_population(N,L, NT, NV, BID, max_deg);
	// evaluate population
	Vector *fitness = new_vector(N); 
	evaluate(pop,fitness, 0, N, NT, NV, BID);

	for (int g = 0; g < gen; g++)
	{
		// duplicate population and fitness
		//pop = tile(sub(pop,0,N,0,L),2,1);
		//fitness = repeat(sub(fitness,0,N),2);
		print_shape(pop,(char*)"population");
		printf("Size of fitness: %d\n",fitness->size );
		// cross
		//annular_cross(pop, N, L, Pc);
		// mutation
		//mutate(pop, N, L, Pm);
		// evaluate new population
		evaluate(pop,fitness, 0, N, NT, NV, BID);
		// sort population and fitness
		//sort(pop, fitness);
		// print results
		system("clear");
		printf("Generation %d:\n",g);
		print_scores(fitness, 20);

	}
	printf("Best individual: \n");
	Matrix *terms = decode(*(pop->M),NT,NV,BID);
	long double eps_th,eps_ph;
	Vector *c = ascend(terms,eps_th,eps_ph);
	printf("Minimax error: %Lf\n", *(fitness->V));
	printf("coefficients:\n");
	print(c);

}

/*###################### Main entrance of the program ########################*/

int main(int argc, char const *argv[])
{
	// long double epsilon_th;
	// long double epsilon_ph;
	// Vector *deg = new_vector(3);
	// deg->V[0] = 3;
	// deg->V[1] = 3;
	// deg->V[2] = 3;
	
	// Matrix *terms = get_terms(deg);
	// print_shape(terms,(char*)"terms");
	// print(terms);
	// Vector *c = ascend(terms, epsilon_th, epsilon_ph);
	// printf("\nEpsilon Theta: %Lf\nEpsilon Phi: %Lf\n", epsilon_th, epsilon_ph);
	// printf("\ncoefficients found\n");
	 print(c);

	run_ega();
	
	// printf("Sub Matrix A_cap(1,3,2,5):\n");
	// print(sub(A,0,5,0,5));


	return 0;
}



