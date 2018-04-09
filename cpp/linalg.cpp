#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>


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
double str2double(char *num)
{
	double res = 0;
	int pw = 0;
	int sign = 1;
	if(strchr(num, '.') == NULL) return (double) str2int(num);
	if (*num == '-'){ sign = -1; num++; }
	for (int i = 0; *(num+i) != '.'; ++i, pw++);
	while(*(num)!='.') res+= (int(*(num++))-48)*pow(10,--pw);
	num++;
	while(*(num)>='0' && *(num)<='9') res += (int(*(num++))-48)*pow(10,--pw);
	return sign*res;
}

/*#######################  LINEAR ALGEBRA MODULE  ############################*/

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

Matrix* append(Matrix *A, Vector *b, int axis)
{
	if (axis==0) return append_r(A,b);
	return append_c(A,b);
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

/* swap rows of a matrix */
void swap_r(Matrix *M, int row_a, int row_b)
{
	for (int col_i = 0; col_i < M->cols; ++col_i)
	{
		double temp = *(*(M->M+row_a)+col_i);
		*(*(M->M+row_a)+col_i) = *(*(M->M+row_b)+col_i);
		*(*(M->M+row_b)+col_i) = temp;
	}
}
/* swap columns of a matrix */
void swap_c(Matrix *M, int col_a, int col_b)
{
	for (int row_i = 0; row_i < M->rows; ++row_i)
	{
		double temp = *(*(M->M+row_i)+col_a);
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
	double temp = *(u->V+a);
	*(u->V+a) = *(u->V+b);
	*(u->V+b) = temp;
}

double min(Vector* u, int &index)
{
	double min = *(u->V);
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

double max(Vector* u, int &index)
{
	double max = *(u->V);
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
	for (int i = 0; i < u->size; ++i) printf("%f ", *(u->V+i));
	printf(" ]\n");
}

/* prints a matrix */
void print(Matrix *A)
{
	for (int i = 0; i < A->rows; ++i) print(new_vector(*(A->M+i),A->cols));
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

bool isSingular(Matrix *A)
{
	if (A->rows == 1) return *(*(A->M));
	double D = 0;
	int sign = 1; 
 
	for (int f = 0; f < A->cols; f++) {
		D += sign * (*(*(A->M)+f)) * isSingular(getCofactor(A, 0, f));
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
	if(isSingular(A) == 0)
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
			double c= *(*(A_temp->M+i)+j)/(*(*(A_temp->M+j)+j));
			for(int k=0; k<m+1; k++)
				*(*(A_temp->M+i)+k) = *(*(A_temp->M+i)+k)-c*(*(*(A_temp->M+j)+k));
		}
	}
	
	*(x->V+m-1) = *(*(A_temp->M+m-1)+m)/(*(*(A_temp->M+m-1)+m-1));

	/* this loop is for backward substitution*/
	for(int i=m-2; i>=0; i--)
	{
		double sum=0;
		for(int j=i+1; j<m; j++) sum += *(*(A_temp->M+i)+j) * (*(x->V+j));
		*(x->V+i) = (*(*(A_temp->M+i)+m)-sum)/(*(*(A_temp->M+i)+i));
	}

	return x;
}

/*#######################  FILE READING  #########################*/

/*
reads a csv file and return data matrix **M
*/
Matrix* read_csv(char *filename, char separator, int rows, int fields)
{
	FILE* stream = fopen(filename, "r");
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



int main(int argc, char const *argv[])
{
	// int v_size = 3;
	// int rows = 3;
	// int cols = 3;
	// Vector *u = new_vector(v_size-1);
	// Vector *v = new_vector(v_size);
	// Vector *b = new_vector(3);
	// Matrix *A = new_matrix(3,3);
	// Matrix *B = new_matrix(3,2);
	// Matrix *I = new_I_matrix(3);
	// A->M[0][0] = 1;
	// A->M[0][1] = 2;
	// A->M[0][2] = 3;
	// A->M[1][0] = 4;
	// A->M[1][1] = 5;
	// A->M[1][2] = 6;
	// A->M[2][0] = 7;
	// A->M[2][1] = 8;
	// A->M[2][2] = 9;
	
	// B->M[0][0] = 2;
	// B->M[0][1] = 2;
	// B->M[1][0] = 1;
	// B->M[1][1] = 0;
	// B->M[2][0] = 1;
	// B->M[2][1] = 1;

	// u->V[0] = 3;
	// u->V[1] = 2;
	// v->V[0] = 1;
	// v->V[1] = 2;
	// v->V[2] = 3;
	// b->V[0] = 4;
	// b->V[1] = 3;
	// b->V[2] = 5;
	char *filename_ = (char*)malloc(sizeof(char)*1000);
	printf("\nname of the file: ");
	scanf("%s",filename_);
	//printf("number parsed: %f\n", str2double(filename_));
	Matrix *data = read_csv(filename_,'\t',1000,23);
	printf("data read: \n");
	print(data);
	// Matrix *data = read_csv(filename_,',');

	return 0;
}


