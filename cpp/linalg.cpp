#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* ########################### OTHER OPERATORS ############################*/

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

Matrix* new_I_matrix(int size)
{
	Matrix *M = new_matrix(size,size);
	for (int i = 0; i < size; ++i) M->M[i][i] = 1;
	return M;
}

Matrix* transpose(Matrix *M)
{
	Matrix *T = new_matrix(M->cols,M->rows);
	for (int i = 0; i < T->rows; ++i)
		for (int j = 0; j < T->cols; ++j) T->M[i][j] = M->M[j][i];
	return T;
}
/* gets a subsection of the vector */
Vector* sub(Vector *v, int ini, int fin)
{
	int len = fin-ini;
	long double *v_pt = new long double[len];
	for (int i = 0; i < len; i++, ini++) *(v_pt+i) = *(v->V+ini);
	return new_vector(v_pt,len);
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
long double dot(Vector *u, Vector *v)
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
	long double* v_pt = new long double[A->rows];
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

/* swap rows of a matrix */
void swap_r(Matrix *M, int row_a, int row_b)
{
	for (int col_i = 0; col_i < M->cols; ++col_i)
	{
		long double temp = *(*(M->M+row_a)+col_i);
		*(*(M->M+row_a)+col_i) = *(*(M->M+row_b)+col_i);
		*(*(M->M+row_b)+col_i) = temp;
	}
}
/* swap columns of a matrix */
void swap_c(Matrix *M, int col_a, int col_b)
{
	for (int row_i = 0; row_i < M->rows; ++row_i)
	{
		long double temp = *(*(M->M+row_i)+col_a);
		*(*(M->M+row_i)+col_a) = *(*(M->M+row_i)+col_b);
		*(*(M->M+row_i)+col_b) = temp;
	}
}
/* swap two columns or rows of a vector */
void swap(Matrix *M, int a, int b, int axis)
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
void print(Vector *u)
{
	printf("[ ");
	for (int i = 0; i < u->size; ++i) printf("%10.12Lf ", *(u->V+i));
	printf("]\n");
}

/* prints a matrix */
void print(Matrix *A)
{
	for (int i = 0; i < A->rows; ++i) print(new_vector(*(A->M+i),A->cols));
}

/* print shape */
void print_shape(Matrix *A, char* label)
{
	printf("shape of %s= (%d,%d)\n", label, A->rows, A->cols);
}
void print_shape(Vector *v, char* label)
{
	printf("shape of %s = (%d, )\n",label, v->size);
}

Matrix* getCofactor(Matrix *M, int p, int q)
{
	int i = 0, j = 0;
	Matrix *temp = new_matrix(M->rows-1,M->cols-1);
 
	// Looping for each element of the matrix
	for (int row = 0; row < M->rows; row++) 
	{
		for (int col = 0; col < M->cols; col++) 
		{
			if (row != p && col != q) 
			{
				*(*(temp->M+i)+(j++)) = *(*(M->M+row)+col);
				if (j == M->cols-1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
	return temp;
}

long double det(Matrix *A)
{
	if (A->rows == 1) return *(*(A->M));
	long double D = 0;
	int sign = 1; 
 
	for (int f = 0; f < A->cols; f++) {
		D += sign * (*(*(A->M)+f)) * det(getCofactor(A, 0, f));
		sign = -sign;
	}
	return D;
}

/* returns Vector x of coefficients */
Vector* solve(Matrix *A, Vector *b)
{
	if (A->rows != A->cols)
	{
		printf("Error in solve(A,b): Matrix A must be squared.\n");
		return NULL;
	}
	if (A->rows != b->size)
	{
		printf("Error in solve(A,b): Vector b must be same size as Matrix A order.\n");
		return NULL;
	}
	if(det(A) == 0)
	{
		printf("Error in solve(A,b): Matrix A is singular and can't be solved.\n");
		return NULL;
	}

	int m = A->rows;
	Matrix *A_temp = append(A,b,1);
	Vector *x = new_vector(m);

	/* loop for the generation of upper triangular matrix*/
	for(int j=0; j<m; j++) 
	{
		for(int i=j+1; i<m; i++)
		{
			long double c= *(*(A_temp->M+i)+j)/(*(*(A_temp->M+j)+j));
			for(int k=0; k<m+1; k++)
				*(*(A_temp->M+i)+k) = *(*(A_temp->M+i)+k)-c*(*(*(A_temp->M+j)+k));
		}
	}
	
	*(x->V+m-1) = *(*(A_temp->M+m-1)+m)/(*(*(A_temp->M+m-1)+m-1));

	/* this loop is for backward substitution*/
	for(int i=m-2; i>=0; i--)
	{
		long double sum=0;
		for(int j=i+1; j<m; j++) sum += *(*(A_temp->M+i)+j) * (*(x->V+j));
		*(x->V+i) = (*(*(A_temp->M+i)+m)-sum)/(*(*(A_temp->M+i)+i));
	}

	return x;
}

/* adjoint of a matrix */
Matrix* adjoint(Matrix *A)
{
	if (A->cols == 1)
	{
		Matrix *adj = new_matrix(1,1);
		*(*adj->M) = 1;
		return adj;
	}
	Matrix *adj = new_matrix(A->cols,A->cols);
	int sign = 1;
	for (int i = 0; i < A->cols; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			sign = ((i+j)%2==0)? 1: -1;
			*(*(adj->M+j)+i) = sign*det(getCofactor(A,i,j));
		}
	}
	return adj;
}

/* inverse of a matrix */
Matrix* inverse(Matrix *A)
{
	if (A->rows != A->cols)
	{
		printf("Error in inverse(A): Matrix A must be squared.\n");
		return NULL;
	}
	long double D = det(A);
	if (D == 0)
	{
		printf("Error in inverse(A): can't get inverse of singular  matrix.\n");
	}
	Matrix *inverse = new_matrix(A->cols, A->cols);
	A = adjoint(A);
	for (int i = 0; i < A->cols; ++i)
		for (int j = 0; j < A->cols; ++j)
			*(*(inverse->M+i)+j) = (*(*(A->M+i)+j))/D;
	return inverse;
}

/*#######################  FILE READING  #########################*/

/*
reads a csv file and return data matrix M
*/
Matrix* read_csv(char *filename, char separator, int rows, int fields)
{
	FILE* stream = fopen(filename, "r");
	if(stream == NULL)
	{
		printf("Error in read_csv(filename,sep,rows,fields): file does not exists.\n");
	}
	Matrix *data = new_matrix(rows, fields);
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
			*((*(data->M+i))+j) = str2double(num);
			if (c == EOF) break;
			c = fgetc(stream);
		}
		if (c == EOF) break;
	}
	fclose(stream);
	return data;
}

/*###################### ASCEND ALGORITHM ########################*/

/* get minimax signs by te 4th method of the minimax theory document */
Vector* get_signs(Matrix *inner)
{
	inner = transpose(inner);
	Matrix *B = remove(inner,inner->cols-1,1);
	Vector *f = get_col(inner,inner->cols-1);
	Vector *signs = solve(B,f);
	signs = append(signs,-1);
	signs = sign(signs);
	printf("Minimax signs:\n");
	print(signs);
	return signs;
}

/* calculates all the combinations of the degrees of the variables */
Matrix* get_terms(Vector *deg)
{
	int comb = 1;
	int temp = 1;
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

/* maps the data into new dataset with the degrees of the variables given */
Matrix* map(Matrix *data, Matrix *terms)
{
	if (data->cols-1 != terms->cols)
	{
		printf("Error in map(data,terms): data fields and terms columns must have the same size.\n");
		return NULL;
	}
	
	Matrix *m_data = new_matrix(data->rows,terms->rows+1); // +1 for the dependent variable

	for (int i = 0; i < m_data->rows; i++)
	{
		for (int j = 0; j < m_data->cols-1; j++)
		{
			long double value = 1;
			for (int k = 0; k < terms->cols; k++)
			{
				//value = data[i][k]^terms[j][k]
				value *= powl((*(*(data->M+i)+k)),(*(*(terms->M+j)+k)));
			}
			*(*(m_data->M+i)+j) = value;
		}
		*(*(m_data->M+i) + m_data->cols-1) = *(*(data->M+i) + data->cols-1);
	}
	return m_data;
}

Vector* get_coeff(Matrix *B, Vector *f, long double &eps_th)
{
	Vector *C = mul(B,f);
	eps_th = *C->V;
	return remove(C,0);
}

/* test coefficients in the outter set.
   input
      -outter: outter set
      -coefficients: vector of coefficients(with the error removed)
      -eps_ph: refernce variable where error is stored 
      -sgn: reference variable where sign of the error is stored
      -idx: index of the vector with the maximum error */
void test_coeff(Matrix *outter, Vector *coefficients, long double &eps_ph, int &sgn, int &idx)
{	
	long double error, abs_err;
	eps_ph = -10000;
	for (int i = 0; i < outter->rows; ++i)
	{
		Vector *out_x = new_vector(*(outter->M+i),0,outter->cols-1);

		long double y = *(*(outter->M+i)+outter->cols-1);
		error = y - dot(coefficients,out_x);
		abs_err = abs(error);
		if(abs_err > eps_ph)
		{
			eps_ph = abs_err;
			idx = i;
			sgn = sign(error);
		}
	}
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
		*(*(B->M+i)+betha) /= *(lambdas->V+betha);
	for (int i = 0; i < B->rows; ++i)
		for (int j = 0; j < B->rows; ++j)
			if (i!=betha) *(*(B->M+j)+i) -= (*(lambdas->V+i))*(*(*(B->M+j)+betha));
}

/* stabilizes data by adding neglectable random value */
Matrix* stabilize(Matrix *A, long double factor)
{
	Matrix *S = new_matrix(A->rows,A->cols);
	srand(time(NULL));
	double random;
	for (int i = 0; i < S->rows; ++i)
	{
		for (int j = 0; j < S->cols; ++j)
		{
			random = (double)rand() / (double)RAND_MAX ;
			if(*(*(A->M+i)+j) == 0) 
				*(*(S->M+i)+j) = *(*(A->M+i)+j) + random*factor;
			else *(*(S->M+i)+j) = (*(*(A->M+i)+j)) * (1+random*factor);
		}
	}
	return S;
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
	print_shape(inner,(char*)"inner");
	print_shape(outter,(char*)"outter");
	
	*(*(inner->M+bmi)) = mu;
	for (int i = 0; i < outter->cols; ++i)
	{
		long double temp = *(*(inner->M+bmi)+i+1);
		*(*(inner->M+bmi)+i+1) = *(*(outter->M+IE)+i);
		*(*(outter->M+IE)+i) = temp;
	}

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
Vector *ascend(Matrix *terms, long double &eps_th, long double &eps_ph)
{
	int IE, mu;
	int m = terms->rows;
	int M = m+1;
	// long double eps_th,eps_ph;
	Vector *c; // coefficient vector

	long double stab_fac = 1e-6;
	long double quasi = 0.05;
	bool Q_F = true;

	char* filename = (char*) malloc(100);
	int rows = 0;
	int fields = 0;
	
	printf("Training filename: \n");
	scanf("%s",filename);
	printf("data shape(rows fields):");
	scanf("%d %d",&rows,&fields);

	// read data
	Matrix *outter = read_csv(filename,'\t', rows, fields);
	
	// map data
	outter = map(outter, terms);

	// stabilize data
	outter = stabilize(outter, stab_fac);

	// split data
	Matrix *inner = sample(outter,M,true);
	
	// get minimax signs
	Vector *signs = get_signs(remove(inner,inner->cols-1,1));

	// get matrix A
	inner = insert(inner,0,signs,1); 

	// get 1st inverse
	Matrix *B = inverse(remove(inner,inner->cols-1,1)); 

	int iteration = 1;

	while(true)
	{
		//printf("check1\n");
		c = get_coeff(B, get_col(inner,inner->cols-1), eps_th);
		//printf("check2\n");
		test_coeff(outter, c, eps_ph, mu, IE);
		printf("IT[%d]: eps_th = %2.6Lf eps_ph = %2.6Lf\n",iteration, eps_th, eps_th);
			if ((eps_th >= eps_ph) || (Q_F && abs(eps_th - eps_ph) <= quasi))
		{
			if (Q_F)
			{
				printf("quasi minimax achieved\n");
			}
			break;
		}
		else
		{
			swap_vector(outter, inner, B, mu, IE);
		}
		iteration++;
	}
	return NULL;

}


int main(int argc, char const *argv[])
{
	long double epsilon_th;
	long double epsilon_ph;
	Vector *deg = new_vector(3);
	deg->V[0] = 2;
	deg->V[1] = 2;
	deg->V[2] = 2;
	
	Matrix *terms = get_terms(deg);
	print_shape(terms,(char*)"terms");
	print(terms);
	Vector *c = ascend(terms, epsilon_th, epsilon_ph);
	printf("Epsilon Theta: %Lf\nEpsilon Phi: %Lf\n", epsilon_th, epsilon_ph);
	printf("coefficients found\n");
	print(c);
	// char *filename;
	// scanf("%s",filename);
	// Matrix *terms = read_csv(filename,'\t',3,3);
	// print(algo);
	// printf("\n");
	// print(insert(algo,5,999.9));

	return 0;
}



