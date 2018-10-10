#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*#########################################################################*/
/*########################### GLOBAL VARIABLES ############################*/
/*#########################################################################*/

// constants
const int powers[] = {0, 1, 3, 5, 7, 9,
					  11,15,21,25,27,
					  33,35,45,49,55,
					  63,77,81,99,121};
const int MMX = 1;
const int RMS = 2;
const int REC_TST = 1;
const int REC_TRN = 0;
const double rseed = 500221.000;

// file parameters
char* train_file;
int tr_rows;
int tr_cols;
char* test_file;
int ts_rows;
int ts_cols;
int recall_data = REC_TST;
bool HAS_TEST = 1;
bool HAS_RECALL = 0;

// parameters of the ascend algorithm
long double stab_fac;
long double quasi;
bool Q_F;
int minimize;

// parameters of the EGA
double Pc; // crossover probability
double Pm; // mutation probability
int gen; // number of generations
int N; // number of individuals
int max_deg; // maximum degree of the variables
double mu; // standar deviation for the degrees in the terms of the individuals
double sigma;
long double **coeffs;


// parameters inferred from the dataset
int NV; // number of independent variables
int NT; // number of terms




/*########################################################################*/
/*########################### OTHER OPERATORS ############################*/
/*########################################################################*/

int get_degree(int max_d){
	for(int i = 0; i < 21; i++) {
		if (max_d <= powers[i])
		{
			return powers[i];
		}
	}
	return powers[20];
}

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
	return ((long double)rand())/((long double)RAND_MAX);
}

/* generates a random integer in [a,b] */
int randint(int min, int max)
{
	return (rand() % (max + 1 - min)) + min;
}

/**/
int get_Off()
{
	double Y = 0;
	for (int i = 1; i <= 12; i++) Y += rand01();
	Y = round(fabs(sigma*(Y-6)+mu)+.5);
	return Y;
}

/* shuffles an array */
void shuffle(long double *&array, int n)
{
    if (n > 1) 
    {
        int i;
        for (i = 0; i < n - 1; i++) 
        {
          int j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = *(array+j);
          *(array+j) = *(array+i);
          *(array+i) = t;
        }
    }
}

/*##################################################################################*/
/*########################### LINEAL ALGEBRA MODULE ################################*/
/*##################################################################################*/

long double** new_Ix(int size)//Ix is squared
{
	long double **A = (long double**) malloc(sizeof(long double*)*size);
	for (int i = 0; i < size; ++i) *(A+i) =  (long double*) calloc(size,sizeof(long double));
	for (int i = 0; i < size; ++i) *(*(A+i)+i) = 1;
	return A;
}

long double** new_array_f(int rows, int cols)
{
	//long double **A = (long double**) calloc(rows,sizeof(long double));
	long double **A = (long double**) calloc(rows,sizeof(long double*)*rows);
	for (int i = 0; i < rows; ++i) *(A+i) =  (long double*) calloc(cols,sizeof(long double));
	return A;
}

long double* new_array_f(int size)
{
	long double *u = (long double*)calloc(size,sizeof(long double));
	return u;
}

long double dot(long double *u, int ul, long double *v, int vl)
{
	if(ul != vl)
	{
		printf("Error E05 in dot(u,v): Vector sizes must be the same\n");
		return 0;
	}
	long double sum = 0;
	for (int i = 0; i < ul; i++)
		sum += (*(u+i)) * (*(v+i));
	return sum;
}

/* product of matrix M and a Vector u*/
long double* mul(long double **A, int ar, int ac, long double *u, int ul)
{
	if (ac != ul)
	{ 
		printf("Error E04 in mul(A,v): Matrix A cols must be same size as u rows\n"); 
		return NULL; 
	}
	long double* v_pt = new_array_f(ar);
	for (int i = 0; i < ar; ++i)
		*(v_pt+i) = dot(*(A+i),ac,u,ul);
	return v_pt;
}

long double* mul(long double *u, int ul, long double **A, int ar, int ac)
{
	if (ac != ul)
	{ 
		printf("Error E04 in mul(A,v): Matrix A cols must be same size as u rows\n"); 
		return NULL; 
	}
	long double* v_pt = new_array_f(ar);
	long double sum;
	for (int i = 0; i < ac; i++)
	{
		sum = 0;
		for (int j = 0; j < ar; j++)
			sum += *(*(A+j)+i) * *(u+j);
		*(v_pt+i) = sum;
	}
		
	return v_pt;
}

long double** transpose(long double **M,int mr,int mc)
{
	long double** MT = new_array_f(mc,mr);
	for (int i = 0; i < mc; ++i)
		for (int j = 0; j < mr; ++j) 
			*(*(MT+i)+j) = *(*(M+j)+i);
	return MT;
}

/* solves a system of equations and stores the result in c
***/
bool lassol(long double **XY,int xyr,int xyc, long double *&c)
{
	int m = xyr;
	int mp1 = xyc;
	long double rowmax, scale, big, temp, quot, sum;
	int ipiv,kp1,ii;
	long double *temp_a;
	// scale each row to its max element
	for (int i = 0; i < m; i++)
	{
		rowmax = fabsl(*(*(XY)+i));
		for (int j = 1; j < m; j++)
			rowmax = fmaxl(fabsl(*(*(XY+i)+j)), rowmax);
		if (rowmax == 0)
		{
			printf("Error E03 in solve(XY): unstable system.\n");
			return false;
		}
		scale = 1/rowmax;
		for (int j = 0; j < mp1; j++)
			*(*(XY+i)+j) *= scale;
	}

	/* Largest element in column i */
	for (int k = 0; k < m; k++)
	{
		big = 0;
		for (int i = k; i < m; i++)
		{
			temp = fabsl(*(*(XY+i)+k));
			if (big < temp)
			{
				big = temp;
				ipiv = i;
			}
		}
		if (big == 0)
		{
			printf("Error E03 in solve(XY): unstable system.\n");
			return false;
		}

		if (ipiv != k) // swap all in column except first
		{
			for (int i = k; i < mp1; ++i)
			{
				temp = *(*(XY+k)+i);
				*(*(XY+k)+i) = *(*(XY+ipiv)+i);
				*(*(XY+ipiv)+i) = temp;
			}
		}

		// eliminate all in column except first
		kp1 = k+1;
		for (int i = kp1; i < m; i++)
		{
			quot = (*(*(XY+i)+k))/(*(*(XY+k)+k));
			for (int j = kp1; j < mp1; j++)
				*(*(XY+i)+j) -= (*(*(XY+k)+j)) * quot;
		}
	}
	if (*(*(XY+m-1)+m-1)==0)
	{
		printf("Error E03 in solve(XY): unstable system.\n");
		return false;
	}
	

	// back substitution
	*(c+m-1) = (*(*(XY+m-1)+m))/(*(*(XY+m-1)+m-1));
	
	for (int ib = 1; ib <= m; ib++)
	{
		ii 	= m-ib;
		sum = 0;
		for (int j = ii+1; j < m; j++)
			sum += (*(*(XY+ii)+j)) * (*(c+j));
		////printf("check2 ii=%d, m=%d\n",ii,m);
		*(c+ii) = (*(*(XY+ii)+m)-sum)/(*(*(XY+ii)+ii));
	}
	return true;
}

long double** inverse(long double **A_, int A_r, int A_c)
{

	int dim = A_c;
	long double ** A = new_array_f(A_r,A_c);
	for (int i = 0; i < A_r; i++)
		for (int j = 0; j < A_c; j++)
			*(*(A+i)+j) = *(*(A_+i)+j);
	long double **Ix = (long double**)new_Ix(dim);

	// scale matrices
	for (int i = 0; i < dim; i++)
	{
		long double r_max = fabsl(*(*(A+i)));
		for (int j = 1; j < dim; j++)
			r_max = fmaxl(fabsl(*(*(A+i)+j)), r_max);
		if (r_max == 0)
			printf("Error E03.1 in inverse(A): can't get inverse, unstable matrix.\n");
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
			printf("Error E03.2 in inverse(A): can't get inverse, unstable matrix.\n");
		if (ipiv != k)
		{// swap rows
			long double *temp_a;
			temp_a = *(A+ipiv);
			*(A+ipiv) = *(A+k);
			*(A+k) = temp_a;

			temp_a = *(Ix+ipiv);
			*(Ix+ipiv) = *(Ix+k);
			*(Ix+k) = temp_a;
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
			printf("Error E03.3 in inverse(A): can't get inverse, unstable matrix.\n");
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
	for(int i = 0; i < A_r; i++) free(*(A+i));
	free(A);
	return Ix;

}


/*###############################################################################*/
/*############################### I/O MODULE ####################################*/
/*###############################################################################*/

void print(long double *u, int len)
{
	printf("[ ");
	for (int i = 0; i < len; i++) printf("%10.12Lf ", *(u+i));
	printf("]\n");
}

void print(long double **A, int ar, int ac)
{
	for (int i = 0; i < ar; i++) print(*(A+i),ac);
}

void printInt(long double *u, int len)
{
	printf("[ ");
	for (int i = 0; i < len; i++) printf("%3.0Lf ", *(u+i));
	printf("]\n");
}

void printInt(long double **A, int ar, int ac)
{
	for (int i = 0; i < ar; i++) printInt(*(A+i),ac);
}


void print_params()
{
	printf("--- FILE DEFINITION\n\n");
	printf("TRAIN FILE = %s\n",train_file);
	printf("TEST FILE  = %s\n",test_file);
	if (HAS_RECALL) printf("RECALL DATA  = %s\n",((recall_data == REC_TST)? "TEST":"TRAIN"));
	else printf("RECALL DATA  = %s\n","none");
	printf("NUMBER OF FIELDS = %d\n",tr_cols);
	printf("ROWS IN TRAIN FILE = %d\n",tr_rows);
	printf("ROWS IN TEST FILE = %d\n",(HAS_TEST)? ts_rows : 0);
	printf("--------------------------------\n\n");
	printf("STABILIZATION FACTOR = %Lf\n",stab_fac);
	printf("QUASI MINIMAX CONDITION = %s\n",(char*)((Q_F == 1)? "true":"false"));
	printf("QUASI MINIMAX VALUE = %Lf\n",quasi);
	printf("OPTIMIZE = %s\n",(char*)((minimize == MMX)? "minimax":"rms"));
	printf("--------------------------------\n\n");
	printf("CROSS PROBABILITY = %f\n",Pc);
	printf("MUTATION PROBABILITY = %f\n",Pm);
	printf("GENERATIONS = %d\n",gen);
	printf("INDIVIDUALS = %d\n",N);
	printf("MAXIMUM DEGREE = %d\n",max_deg);
	printf("--------------------------------\n\n");
	printf("NUMBER OF TERMS = %d\n",NT);
}

/*
** reads a csv file and return data matrix M
** WARNING!!!: the file must end with a number, 
** NOT with a '\n' or a blank space, 
** this may lead to wrong values reading.
*/
long double** read_csv(char *filename, char separator, int rows, int fields)
{
	FILE* stream = fopen(filename, "r");
	if(stream == NULL)
	{
		printf("Error E00.1 in read_csv(filename,sep,rows,fields): file does not exists.\n");
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

void to_text(char *filename, long double *u, int len)
{
	FILE *csv = fopen(filename, "w");
	for(int j=0;j<len;j++)
	    fprintf(csv,"%10.12Lf\t",*(u+j));
	fprintf(csv,"\n");
	fclose(csv);
}

void to_text(char *filename, long double **A, int ar,int ac)
{
	FILE *csv = fopen(filename, "w");
	for(int i=0;i<ar;i++) 
	{
		for(int j=0;j<ac;j++)
		    fprintf(csv,"%10.12Lf\t",*(*(A+i)+j));
		fprintf(csv,"\n");
	}
	fclose(csv);
}

void get_rows_and_cols(char *filename,char sep, int &rows, int &cols)
{
	FILE* stream = fopen(filename, "r");
	if(stream == NULL)
	{
		printf("Error E00.2 in get_rows_and_cols(filename,sep): file %s does not exists.\n",filename);
		exit(0);
	}
	int r = 0;
	int c = 1;
	char ch = fgetc(stream);
	while(ch != '\n')
	{
		if(ch==sep) c++;
		ch = fgetc(stream);
	}
	r++;
	ch = fgetc(stream);
	while(ch != EOF)
	{
		////printf("check2\n");
		while(ch != '\n') ch = fgetc(stream);
		r++;
		ch = fgetc(stream);
	}
	fclose(stream);
	rows = r;
	cols = c;

}

void read_input()
{
	FILE* stream = fopen((char*)"configuration.conf", "r");
	if(stream == NULL)
	{
		printf("Error E09 in read_input(): file %s does not exists.\n", "configuration.conf");
	}
	char c;
	for (int i = 0; i < 13; ++i) // 12 inputs, check configuration file
	{
		int j = 0;
		char* input = (char*)malloc(sizeof(char)*30);
		char* value = (char*)malloc(sizeof(char)*30);
		c = fgetc(stream);
		while(c != ' ') { input[j++] = c; c = fgetc(stream); }
		input[j] = '\0';
		while(c == ' ' || c == '=') c = fgetc(stream);
		j = 0;
		while(c != '\n' && c != EOF ) { value[j++] = c; c = fgetc(stream); }
		value[j] = '\0';
		char none[50]; strcpy(none, "none");
		char test[50]; strcpy(test, "test");
		char train[50]; strcpy(train, "train");
		if (i == 0) train_file = value;
		else if (i == 1) { test_file = value; if (strcmp(none, test_file) == 0) HAS_TEST = 0; else HAS_TEST = 1; }
		else if (i == 2) { 
			if (strcmp(value,none) == 0) HAS_RECALL = 0; 
			else if (strcmp(value,train) == 0) { HAS_RECALL = 1; recall_data = REC_TRN; }
			else if (strcmp(value,test) == 0 && HAS_TEST) { HAS_RECALL = 1; recall_data = REC_TST; }
			else { printf("\n\nRecall test data must be defined, edit configuration file and try again.\n"); exit(0);}
		}
		else if (i == 3) stab_fac = str2double(value);
		else if (i == 4) quasi = str2double(value);
		else if (i == 5) Q_F = ((str2int(value) == 1));
		else if (i == 6) minimize = (1 == str2int(value))? MMX : RMS;
		else if (i == 7) Pc = str2double(value);
		else if (i == 8) Pm = str2double(value);
		else if (i == 9) gen = str2int(value);
		else if (i == 10) N = str2int(value);
		else if (i == 11) max_deg = get_degree(str2int(value));
		else if (i == 12) NT = str2int(value);
		free(input);
	}
	//if(stream == NULL) printf("mamadas\n");
	//print_params();
	fclose(stream);
	printf("checkpoint2\n");
	sigma = max_deg/sqrt(10);
	mu = 0;
	get_rows_and_cols(train_file,'\t',tr_rows,tr_cols);
	printf("checkpoint1\n");
	if(HAS_TEST) get_rows_and_cols(test_file,'\t',ts_rows,ts_cols);
	NV = tr_cols - 1;

}

/*################################################################*/
/*###################### ASCEND ALGORITHM ########################*/
/*################################################################*/

/* get minimax signs by te 4th method of the minimax theory document 
** and appends the signs found to the inner set and asigns it to inner_*/
void get_signs(long double **inner,int ir,int ic, long double **&inner_)
{
	long double **inner_cap = new_array_f(ir,ic-1); // matrix to store inner - f column of inner
	for (int i = 0; i < ir; i++)
		for (int j = 0; j < ic-1; j++)
			*(*(inner_cap+i)+j) = *(*(inner+i)+j);

	long double **innT = transpose(inner_cap,ir,ic);
	long double *signs = new_array_f(ir);
	*(signs+ir-1) = -1;
	lassol(innT,ic-1,ir, signs);

	for(int i = 0; i < ic; i++) free(*(innT+i));
	for(int i = 0; i < ir; i++) free(*(inner_cap+i));
	free(innT);
	free(inner_cap);
	
	
	// append signs to inner_ first column
	for (int i = 0; i < ir; i++) *(*(inner_+i)) = sign(*(signs+i));
	for (int i = 0; i < ir; i++)
		for (int j = 0; j < ic; j++)
			*(*(inner_+i)+j+1) = *(*(inner+i)+j);
}

/* calculates all the combinations of the degrees of the variables given 
** input.
**	deg: array of degrees for each variable
**	dl: length of deg array
**  combinations: variable where the number of rows(or combinations) of 
**  the terms are gonna be stored*/
long double** get_terms(long double *deg,int dl, int &combinations)
{
	int comb = 1;
	int temp = 1;
	for(int i = 0; i < dl; i++) *(deg+i) = *(deg+i) + 1;
	for(int i = 0; i < dl; i++) comb *= *(deg+i);
	combinations = comb;
	long double **terms = new_array_f(comb, dl);
	for (int ci = dl-1; ci >=0 ; ci--)
	{
		for (int ri = 0; ri < comb; ri++)
		{
			int max = *(deg+ci);
			for (int di=0; di < max; di++,ri++)
			{
				for (int k = 0; k < temp; k++,ri++)
				{
					*(*(terms+ri)+ci) = di;
				}
				ri--;
			}
			ri--;
		}
		temp *= *(deg+ci);
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
long double** map(long double **data,int dr,int dc, long double **terms, int tr, int tc)
{
	// checkMM for memmory liberation
	if (dc-1 != tc)
	{
		printf("Error E01 in map(data,terms): data fields and terms columns must have the same size.\n");
	}
	int m_dr = dr;
	int m_dc = tr+1;
	long double **m_data = new_array_f(m_dr,m_dc); // +1 for the dependent variable

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
	return m_data;
}

void get_coeff(long double **B,int Br,int Bc, long double **A, long double *c, int cl, long double &eps_th)
{
	long double *f = new_array_f(Br);
	for (int i = 0; i < Br; i++)
		*(f+i) = *(*(A+i)+Bc);
	long double *C = mul(B,Br,Bc, f,Br);
	eps_th = *C;
	for (int i = 0; i < cl; ++i)
		*(c+i) = *(C+i+1);
	free(C);
	free(f);
}

/* test coefficients in the outter set.
   input
      -outter: outter set
      -coefficients: vector of coefficients(with the error removed)
      -eps_ph: refernce variable where error is stored 
      -sgn: reference variable where sign of the error is stored
      -idx: index of the vector with the maximum error */
long double test_coeff(long double **outter,int o_r,int o_c, long double *coefficients, int cl, int &sgn, int &idx)
{	
	long double error, abs_err;
	long double eps_ph;
	eps_ph = -100000;
	for (int i = 0; i < o_r; i++)
	{
		long double y_cap = 0;
		for (int j = 0; j < cl; j++)
			y_cap += (*(coefficients+j)) * (*(*(outter+i)+j));

		error = *(*(outter+i)+o_c-1) - y_cap;
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
void get_new_inverse(long double **&B,int br, int betha, long double *lambdas)
{
	for (int i = 0; i < br; i++)
		*(*(B+i)+betha) =  *(*(B+i)+betha) / *(lambdas+betha);

	for (int i = 0; i < br; ++i)
		for (int j = 0; j < br; ++j)
			if (i!=betha) *(*(B+j)+i) = *(*(B+j)+i) - (*(lambdas+i)) * (*(*(B+j)+betha));
}

/* stabilizes data by adding neglectable(<=10e-6) random value */
void stabilize(long double **&A, int ar, int ac, long double factor)
{
	double random;
	for (int i = 0; i < ar; i++)
	{
		for (int j = 0; j < ac-1; j++) // -1 to ommit dependent variable column
		{
			random = rand01();
			if(*(*(A+i)+j) == 0 || *(*(A+i)+j) < 1e-8) 
				*(*(A+i)+j) += random*factor;
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
void swap_vector(long double **&outter,int o_r,int o_c, long double **&inner,int ir,int ic, long double **&B,int br,int bc, long double mu, int IE)
{
	long double *amp1 = new_array_f(bc);
	*(amp1) = mu;
	for (int i = 1; i < bc; ++i) *(amp1+i) = *(*(outter+IE)+i-1);
	long double *lambdas = mul(amp1,bc, B,br,bc);
	long double betha_max = -1e38;
	long double betha;
	int bmi = -1; // betha max index
	for (int i = 0; i < bc; ++i)
	{
		betha = mu * (*(lambdas+i) / (*(*(B)+i)));
		if (betha > betha_max)
		{
			betha_max = betha;
			bmi = i;
		}
	}

	// printf("%d\t%d",bmi+1, IE+bc+1);
	
	*(*(inner+bmi)) = mu;
	for (int i = 0; i < o_c; ++i)
	{
		long double temp = *(*(inner+bmi)+i+1);
		*(*(inner+bmi)+i+1) = *(*(outter+IE)+i);
		*(*(outter+IE)+i) = temp;
	}
	// printf("\nLambda Vector:\n");
	// print(lambdas, bc);
	// printf("\n");

	// calculate new inverse
	get_new_inverse(B,br,bmi,lambdas);
	free(lambdas);
	free(amp1);
}

/* get first n rows of a matrix substracting it from the original
** The array A is the original one and S is the sample array where the first elements
** are going to be substracted */
void first(long double **&A,int &ar,int ac, long double **&S, int n, long double **&outter)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < ac; j++)
			*(*(S+ i)+ j) = *(*(A+ i)+ j);
	ar -= n;
	for (int i = 0; i < ar; i++)
		for (int j = 0; j < ac; j++)
			*(*(outter+i)+j) = *(*(A+n+i)+j);
}

/*
** obtains test RMS and test Minimax errors
** input:
**   terms: matrix of terms to map the test data
**   nt,nv: number of terms and number of variables in each term
**   coef: 1D-array of coefficients of length nt
**   rms,mmx: variables where the errors are going to be stored
*/
void get_tst_errors(long double **terms, int nt,int nv, long double *coef, long double &rms, long double &mmx)
{
	long double **data = read_csv(test_file,'\t', ts_rows, ts_cols);
	long double **mapped = map(data,ts_rows, ts_cols, terms,nt,nv);
	for(int i = 0; i < ts_rows; i++) free(*(data+i));
	free(data);
	long double rms_error = 0;
	long double mmx_error, abs_err, error;
	mmx_error = 0;

	for (int i = 0; i < ts_rows; i++)
	{
		long double y_cap = 0;
		for (int j = 0; j < nt; j++)
			y_cap += (*(coef+j)) * (*(*(mapped+i)+j));
		error = *(*(mapped+ i)+ nt) - y_cap;
		rms_error += error*error;
		
		abs_err = fabsl(error);
		if(abs_err > mmx_error)
			mmx_error = abs_err;
	}
	for(int i = 0; i < ts_rows; i++) free(*(mapped+i));
	free(mapped);
	rms_error /= ts_rows;
	rms = sqrtl(rms_error);
	mmx = mmx_error;
}

/*
** returns the train RMS error
** input:
**   MAP: mapped train data
**   mr,mc: dimensions of MAP, mr=MAP rows, mc=MAP columns
**   c: coefficients
**   cl: coefficients' length
*/
long double get_trn_RMS(long double **MAP,int mr,int mc, long double *c, int cl)
{
	long double rms_error = 0;
	long double error;

	for (int i = 0; i < mr; i++)
	{
		long double y_cap = 0;
		for (int j = 0; j < cl; j++)
			y_cap += (*(c+j)) * (*(*(MAP+i)+j));
		error = *(*(MAP+ i)+ mc-1) - y_cap;
		rms_error += error*error;
		
	}
	// printf("rms_error = %Lf\n",rms_error );
	rms_error /= mr;
	rms_error = sqrtl(rms_error);
	return rms_error;
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
long double* ascend(long double **terms, int tr, int tc, long double &eps_th, long double &eps_ph, 
	long double &trn_rms, long double &tst_rms, long double &tst_mmx)
{
	int IE, mu;
	int m = tr;
	int M = m+1;
	long double epsilon_th,epsilon_ph;
	long double *c = new_array_f(m); // coefficient vector
	char z;

	char* filename = train_file;
	// dimensions of the original dataset
	int rows = tr_rows;
	int fields = tr_cols;
	// dimensions of the mapped dataset
	int dr,dc;
	
	// read data
	long double **data = read_csv(filename,'\t', rows, fields);
	// free(filename);
	
	// map data
	long double **MAP = map(data,rows,fields, terms,tr,tc);
	for(int i = 0; i < rows; i++) free(*(data+i));
	free(data);
	dr = rows; // rows of the mapped dataset
	dc = tr+1; // cols of mapped dataset = number of terms plus one

	// stabilize data
	stabilize(MAP,dr,dc, stab_fac);
	// to_text((char*)"MAPPED.xls",MAP,dr,dc);

	// split data
	long double **inner_ = new_array_f(M,dc);
	long double **outter = new_array_f(dr-M,dc);
	first(MAP,dr,dc, inner_,M, outter);

	// to_text((char*)"inner.xls",inner_,M,dc);
	// to_text((char*)"outter.xls",outter,dr-M,dc);
	
	// get minimax signs
	int ar = M;
	int ac = dc+1; // +1 because of the signs column
	long double **A = new_array_f(ar,ac); // A = Matrix augmented with column of signs 


	// printf("M = %d, dc = %d\n", M,dc);
	// print(inner_, M,dc);
	get_signs(inner_,M,dc, A);
	for(int i = 0; i < M; i++) free(*(inner_+i));
	free(inner_);

	// to_text((char*)"A_plus_signs.xls",A,ar,ac);

	// get matrix A
	// get 1st inverse
	// print(A, ar,ac-1);
	long double **B = inverse(A,ar,ac-1);
	int iteration = 1;
	char cont_flag;

	// to_text((char*)"1stB.xls",B,ar,ac-1);
	long double prev_th; // para guardar los errores anteriores y checar convergencia
	prev_th = 10;
	while(true)
	{
		get_coeff(B,ar,ac-1, A, c, m, epsilon_th);
		epsilon_ph = test_coeff(outter,dr,dc, c, m, mu, IE);
		// printf("\n IT[%d]: %4.10Lf\t%4.10Lf\t",iteration, epsilon_th, epsilon_ph);
		
		if (epsilon_th >= epsilon_ph) { /*printf("\nABSOLUTE\n");*/ break; } // Absolute convvergence
		if (Q_F && fabsl(epsilon_ph/epsilon_th-1) <= quasi) { /*printf("\nQUASI-MMX\n");*/ break; } // quasi convergence
		if (fabsl(epsilon_th-prev_th) < 0.0000001) { /*printf("\nDELTA %1.12Lf,%1.12Lf\n",epsilon_th,prev_th);*/ break; } // delta convergence
		
		swap_vector(outter,dr,dc, A,ar,ac, B,ar,ac-1, mu, IE);
		prev_th = epsilon_th;
		// printf("Next? (y/n): ");
		// scanf("%c",&cont_flag);
		// if(cont_flag=='n') break;
		iteration++;
	}
	eps_ph = epsilon_ph;
	eps_th = epsilon_th;

	trn_rms = get_trn_RMS(MAP,dr+M,dc, c,m);
	if(HAS_TEST) get_tst_errors(terms,tr,tc, c, tst_rms, tst_mmx);
	else { tst_rms = -1; tst_mmx = -1; }
	

	for(int i = 0; i < ar; i++) {free(*(A+i)); free(*(B+i));}
	for(int i = 0; i < rows; i++) free(*(MAP+i));
	for(int i = 0; i < dr; i++) free(*(outter+i));
	free(B);
	free(A);
	free(outter);
	free(MAP);
	return c;
}

/*######################################################################*/
/*###################### GENETIC ALGORITHM(EGA) ########################*/
/*######################################################################*/

/*
** true if term is already in individual, false otherwise
*/
bool term_in_ind(long double *ind,int il, long double *T, int nv)
{
	int total_terms = il/nv;
	for (int i = 0; i < total_terms; i++)
	{
		bool equal = true;
		for (int j = 0; j < nv; j++)
		{
			if(*(T+j) != *(ind+(i*nv)+j)) {equal = false; break;}
		}
		if(equal) return true;
	}
	return false;
}

/*
** creates a valid term
*/

long double* gen_valid_term(int nv)
{
	//1 chose degree form a normal distribution
	int off_i = get_Off()-1;
	int deg_i = powers[off_i];
	int sum_i = 0;
	int nv_i = 0;
	//2 find the combination of the powers wich adds to chosen degree
	long double *term = new_array_f(nv);
	if (deg_i == 0)
	{
		for (int i = 0; i < nv; i++)
		{
			*(term+i) = 0;
		}
		return term;
	}
	while(nv_i < nv)
	{
		int r = randint(1,deg_i);
		sum_i += r;
		if(sum_i > deg_i)
		{
			sum_i -= r;
			break;
		}
		*(term+nv_i) = r;
		nv_i++;
	}
	if(nv_i <= nv-1) *(term+nv_i) = deg_i-sum_i;
	if(nv_i == nv && sum_i < deg_i) *(term+nv_i-1) += deg_i-sum_i;
	//3 shufle powers
	shuffle(term,nv);
	return term;
}

/* gen population N= number of individuals, L= Length of each individual,
** nt= number of terms, nv= number of variables, bid= bits in digit,
** max_deg= maximum degree */
long double** gen_population(int N, int L, int nt, int nv/*, int bid, int max_deg*/)
{
	long double **pop = new_array_f(2*N,L);
	for (int ii = 0; ii < N; ii++)
	{
		// naive generation of individuals(fully random)
		// for (int bi = 0; bi < L; bi++)
		// 	*(*(pop+ii)+bi) = round(rand01());
		// making a single individual
		for (int ti = 0; ti < nt; ti++) // for each term
		{
			//generate valid term
			long double *term = gen_valid_term(nv);
			if (term_in_ind(*(pop+ii),L,term,nv)) // if term is already there get another one
			{
				ti--;
				continue;
			}
			for (int t = 0; t < nv; t++)
				*(*(pop+ii)+ t +(ti*nv)) = *(term+t);
			free(term);
		}

	}
	return pop;
}

void repair(long double **&pop,int N,int L,int nt,int nv)
{
	for (int ii = 0; ii < N; ii++)
	{
		for (int ti = 0; ti < nt; ti++)
		{
			//printf("check.1\n");
			long double *term = new_array_f(nv);
			//printf("check.2\n");
			for (int vi = 0; vi < nv; vi++) *(term+vi) = *(*(pop+ii)+(ti*nv)+vi);
			//printf("check.3\n");
			if(term_in_ind(*(pop+ii),L,term,nv))
			{
				//printf("check.4\n");
				long double *temp = gen_valid_term(nv);
				//printf("check.5\n");
				while(term_in_ind(*(pop+ii),L,temp,nv)) {free(temp); temp = gen_valid_term(nv);}
				for (int vi = 0; vi < nv; vi++) *(*(pop+ii)+(ti*nv)+vi) = *(temp+vi);
					//printf("check.6\n");
				free(temp);
			}
			free(term);
		}
	}
}

/* decodes the individual into a matrix of degrees of variables 
** input:
**       individual: pointer to array of 1s and 0s 
**       bid: bytes in digit
**       nt: number of terms in polynomial
**       nv: number of variables */
long double** decode_bin(long double *ind, int nt, int nv, int bid)
{
	long double **terms = new_array_f(nt,nv);
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
			*(*(terms+t)+v) = decimal;	
		}
	}
	return terms;
}

long double** decode_dec(long double *ind, int nt, int nv)
{
	int bid = 0;
	long double **terms = new_array_f(nt,nv);
	for (int t = 0; t < nt; t++) // for each term
		for (int v = 0; v < nv; v++) // for each variable
			*(*(terms+t)+v) = *(ind+(t*nv)+v);	
	return terms;
}

/* returns a vector of fitness values starting at individual ini-th
** to the (fin-1)-th individual */
void evaluate(long double **pop, long double *&fit_tr_mmx, long double *&fit_tr_rms, 
	long double *&fit_ts_mmx, long double *&fit_ts_rms, 
	int ini, int fin, int nt, int nv)
{
	long double eps_th,trn_mmx,trn_rms,tst_mmx,tst_rms;
	for (int ii = ini; ii < fin; ii++)
	{
		long double **terms = decode_dec(*(pop+ii),nt,nv);
		
		// printf("\nEvaluating:\n");
		// print(terms, nt,nv);
		*(coeffs+ii) = ascend(terms,nt,nv,eps_th,trn_mmx,trn_rms,tst_rms,tst_mmx);
		// long double *c = ascend(terms,nt,nv,eps_th,trn_mmx,trn_rms,tst_rms,tst_mmx);
		free(terms);
		*(fit_ts_mmx+ii) = tst_mmx; // test minimax error
		*(fit_ts_rms+ii) = tst_rms; // test RMS error
		*(fit_tr_mmx+ii) = trn_mmx; // train minimax error
		*(fit_tr_rms+ii) = trn_rms; // train RMS error
	}
}

/* annular crossover */
void annular_cross(long double **&pop, int N, int L, double Pc)
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
				long double temp = *(*(pop+i)+j);
				*(*(pop+i)+j) = *(*(pop+n)+j);
				*(*(pop+n)+j) = temp;
			}

		}

	}
}

/* uniform mutation */
void mutate(long double **&pop, int N, int L, int b2m)
{
	while(b2m--)
	{
		int p1 = round(rand01()*L)-1;
		int p2 = round(rand01()*N)-1;

		// if(*(*(pop+p2)+p1) == 0) *(*(pop+p2)+p1) = 1;
		// else *(*(pop+p2)+p1) = 0;

		int temp = randint(0,max_deg);
		while (temp == *(*(pop+p2)+p1)) temp = randint(0,max_deg);
		*(*(pop+p2)+p1) = temp;
	}
}

/*
** sorts pop, fit2, fit3, fit4 according to values in fit
*/
void sort(long double **&pop, long double *&fit, long double *&fit2,long double *&fit3,long double *&fit4, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (*(fit+j) > *(fit+i))
			{
				// swap fitness1 values
				double tmp = *(fit+i);
				*(fit+i) = *(fit+j);
				*(fit+j) = tmp;
				// swap fitness2 values
				tmp = *(fit2+i);
				*(fit2+i) = *(fit2+j);
				*(fit2+j) = tmp;
				// swap fitness3 values
				tmp = *(fit3+i);
				*(fit3+i) = *(fit3+j);
				*(fit3+j) = tmp;
				// swap fitness4 values
				tmp = *(fit4+i);
				*(fit4+i) = *(fit4+j);
				*(fit4+j) = tmp;
				// swap individuals
				long double *ind = *(pop+i);
				*(pop+i) = *(pop+j);
				*(pop+j) = ind;
				//free(ind);
				// swap coefficients
				long double *c = *(coeffs+i);
				*(coeffs+i) = *(coeffs+j);
				*(coeffs+j) = c;
				//free(c);
			}  
		}
	}
}
/*
** Memmory of pop must be reserved for double the size(2N)
** gen_population already does this.
*/
void duplicate(long double **&pop,int N,int L, long double *&fit, long double *&fit2, long double *&fit3, long double *&fit4)
{

	for (int i = 0; i < N; i++)
	{
		*(fit+i+N) = *(fit+i);
		*(fit2+i+N) = *(fit2+i);
		*(fit3+i+N) = *(fit3+i);
		*(fit4+i+N) = *(fit4+i);		
	    *(coeffs+i+N) = *(coeffs+i);
		for (int j = 0; j < L; j++)
		{
			*(*(pop+i+N)+j) = *(*(pop+i)+j);
		}

	}
}

void print_scores(long double *fit, int n)
{
	printf("Best %d scores:\n",n);
	for (int i = 0; i < n; i++)
	{
		printf( "%d: %Lf\t",i,*(fit+i));
	}
}

void recall(char* prefix, long double *c, long double **terms, int NT, int NV)
{
	if (!HAS_RECALL) return;
	char filename[50];
	int dr,dc;
	dc = tr_cols;
	if (recall_data == REC_TRN) { strcpy(filename, train_file); dr = tr_rows; }
	else if (recall_data == REC_TST) { strcpy(filename, test_file); dr = ts_rows; }
	long double **data = read_csv(filename,'\t',dr,dc);
	long double **MAP = map(data,dr,dc, terms,NT,NV);
	int mc = NT+1; // numero de columnas de los datos mapeados, NT + columna de variable dependiente
	to_text((char*)"data.xls",data,dr,dc);
	for (int i = 0; i < dr; i++)
	{
		long double sum = 0;
		for (int j = 0; j < NT; j++) sum += *(*(MAP+i)+j) * (*(c+j));
		*(*(data+i)+dc-1) = sum;
	}

	char recall_file[100];
	strcpy(recall_file, prefix);
	strcat(recall_file,"_recall.xls");
	//strcat(recall_file, (recall_data == REC_TRN)? "train.xls":"test.xls");
	to_text(recall_file,data,dr,dc);
	for (int i = 0; i < dr; ++i){ free(*(MAP+i)); free(*(data+i)); }
	free(MAP);
	free(data);
}
/* return index of the best fitness */
int get_best(long double *fitness, int n)
{
	int min = *(fitness);
	int index = 0;
	for (int i = 1; i < n; i++)
		if (*(fitness+i) < min) { min = *(fitness+i); index = i; }
	return index;
}

void run_ega(/*Vector *&best_ind, long double &best_fit*/)
{
	// parameters of the EGA

	// file definition parameters
	char* f_train = train_file; // (char*)"DB24-glass/TRAIN.TXT";
	char* f_test = test_file; // (char*)"DB24-glass/TEST.TXT";

	// derivated parameters
	int L = NT*NV;
	int b2m = round(L*N*Pm); // parameter for later calculations

	// generate random population
	long double **pop = gen_population(N,L, NT, NV);//, BID, max_deg);
	// arrays to store fitnesses and coefficients of each individual
	long double *fit_tr_mmx = new_array_f(N*2);
	long double *fit_tr_rms = new_array_f(N*2);
	long double *fit_ts_mmx = new_array_f(N*2);
	long double *fit_ts_rms = new_array_f(N*2);
	// evaluate population
	evaluate(pop,fit_tr_mmx,fit_tr_rms,fit_ts_mmx,fit_ts_rms, 0, N, NT, NV);
	for (int g = 0; g < gen; g++)
	{
		printf("\nGeneration %d:\n",g);

		// duplicate population and fitness
		//printf("check1\n");
		duplicate(pop,N,L,fit_tr_mmx,fit_tr_rms,fit_ts_mmx,fit_ts_rms);
		//printf("check2\n");
		// cross
		annular_cross(pop, N, L, Pc);
		//printf("check3\n");
		// mutation
		//printf("check4\n");
		mutate(pop, N, L, Pm);
		// repair population
		//printf("check5\n");
		repair(pop, N, L, NT, NV);
		// evaluate new population
		//printf("check6\n");
		evaluate(pop,fit_tr_mmx,fit_tr_rms,fit_ts_mmx,fit_ts_rms, 0, N, NT, NV);
		// sort population and fitness
		// and print results
		switch(minimize)
		{
			case RMS:
				if ( HAS_TEST ) {
					//printf("check7\n");
					sort(pop, fit_ts_rms,fit_tr_mmx,fit_tr_rms,fit_ts_mmx, 2*N);
					//printf("check8\n");
					print_scores(fit_ts_rms, 20);
				} else {
					//printf("check9\n");
					sort(pop, fit_tr_rms,fit_ts_rms,fit_tr_mmx,fit_ts_mmx, 2*N);
					//printf("check10\n");
					print_scores(fit_tr_rms, 20);
				}
				break;
			case MMX:
				if ( HAS_TEST ) {
					//printf("check11\n");
					sort(pop, fit_ts_mmx,fit_ts_rms,fit_tr_mmx,fit_tr_rms, 2*N);
					//printf("check12\n");
					print_scores(fit_ts_mmx, 20);
				} else {
					//printf("check13\n");
					sort(pop, fit_tr_mmx,fit_ts_mmx,fit_ts_rms,fit_tr_rms, 2*N);
					//printf("check14\n");
					print_scores(fit_tr_mmx, 20);
				}
				break;
		}

		system("clear");
		printf("\nbest train rms = %Lf\tbest test rms = %Lf\nbest train mmx = %Lf\tbest test mmx = %Lf\n",
			*(fit_tr_rms),*(fit_ts_rms),*(fit_tr_mmx),*(fit_ts_mmx));

	}
	// printf("Best individual: \n");
	// long double **terms = decode_dec(*(pop),NT,NV);
	// printf("Terms:\n");
	// printInt(terms,NT,NV);
	// long double eps_th,trn_mmx,trn_rms,tst_mmx,tst_rms;
	// long double *c = ascend(terms,NT,NV,eps_th,trn_mmx,trn_rms,tst_rms,tst_mmx);
	// printf("coefficients:\n");
	// print(c,NT);
	// printf("\nbest train rms = %Lf\tbest test rms = %Lf\nbest train mmx = %Lf\tbest test mmx = %Lf\n",
	// 		trn_rms,tst_rms,trn_mmx,tst_mmx);
	
	// recall((char*)"TST_",c,terms,NT,NV);

	int best;
	long double **terms;

	best = get_best(fit_tr_mmx,N);
	terms = decode_dec(*(pop+best),NT,NV);
	printf("\n\nTRN MMX error: %2.10Lf\n", *(fit_tr_mmx+best));
	printf("terms: \n");
	printInt(terms,NT,NV);
	printf("\t coefficients:\n");
	print(*(coeffs+best),NT);
	recall((char*)"TRN_MMX",*(coeffs+best),terms,NT,NV);

	best = get_best(fit_ts_mmx,N);
	terms = decode_dec(*(pop+best),NT,NV);
	printf("\n\nTST MMX error: %2.10Lf\n", *(fit_ts_mmx+best));
	printf("terms: \n");
	printInt(terms,NT,NV);
	printf("\t coefficients:\n");
	print(*(coeffs+best),NT);
	recall((char*)"TST_MMX",*(coeffs+best),terms,NT,NV);

	best = get_best(fit_tr_rms,N);
	terms = decode_dec(*(pop+best),NT,NV);
	printf("\n\nTRN RMS error: %2.10Lf\n", *(fit_tr_rms+best));
	printf("terms: \n");
	printInt(terms,NT,NV);
	printf("\t coefficients:\n");
	print(*(coeffs+best),NT);
	recall((char*)"TRN_RMS",*(coeffs+best),terms,NT,NV);

	best = get_best(fit_ts_rms,N);
	terms = decode_dec(*(pop+best),NT,NV);
	printf("\n\nTST RMS error: %2.10Lf\n", *(fit_ts_rms+best));
	printf("terms: \n");
	printInt(terms,NT,NV);
	printf("\t coefficients:\n");
	print(*(coeffs+best),NT);
	recall((char*)"TST_RMS",*(coeffs+best),terms,NT,NV);

}

/*############################################################################*/
/*###################### Main entrance of the program ########################*/
/*############################################################################*/


int main(int argc, char const *argv[])
{
	srand(time(NULL));
	//srand(rseed);
	read_input();
	// array of coefficients
	coeffs = new_array_f(N*2, NT);
	// run_ega();

	// long double eps_th,trn_mmx,trn_rms,tst_mmx,tst_rms; // errores interno y externo
	// long double **terms = read_csv((char*)"terms.conf",'\t',NT,NV);
	// long double *vars = new_array_f(NV);
	// vars[0] = 2;
	// vars[1] = 2;
	// vars[2] = 2;
	// //vars[3] = 1;
	// // long double **terms = get_terms(vars,NV,NT);
	// // printInt(terms,NT,NV);	
	// print_params();
	// printf("Continue? ");
	// char r;
	// scanf("%c",&r);
	// if (r != 'y') return 0;
	// long double *coef = ascend(terms,NT,NV, eps_th,trn_mmx,trn_rms,tst_rms,tst_mmx);
	// printf("coefficients found:\n");
	// print(coef,NT);
	// printf("train minimax=%1.12Lf test minimax=%1.12Lf\ntrain rms=%1.12Lf test rms=%1.12Lf\n", trn_mmx, tst_mmx, trn_rms, tst_rms);
	

	system("clear");
	print_params();
	char r;
	printf("\nProceed? [y/n]: ");
	scanf("%c",&r);
	if (r == 'y')
	{
		run_ega();
	}
	else if (r == 'n')
	{
		printf("\nProcess aborted. Edit configuration.conf file if you wish and try again.\n\n");
	}
	else
	{
		printf("\nFATAL ERROR!!! Invalid Option. Program has crashed.\n\n");
	}


	// long double *term = gen_valid_term(20);
	// printf("generated term:\n");
	// print(term,20);
	// int suma = 0;
	// for (int i = 0; i < 20; ++i) suma += term[i];
	// printf("suma %d\n", suma);

	// long double **population = gen_population(10,30, 5, 6, 0, 20);
	// printf("poblacion:\n");
	// print(population,10,30);
	// // printf("decoded:\n");
	// print(decode_dec(population))

	// long double *ind = new_array_f(12); // 4 terms 3 vars
	// ind[0] = 1;
	// ind[1] = 2;
	// ind[2] = 3;
	// ind[3] = 4;
	// ind[4] = 5;
	// ind[5] = 6;
	// ind[6] = 7;
	// ind[7] = 8;
	// ind[8] = 9;
	// ind[9] = 10;
	// ind[10] = 11;
	// ind[11] = 12;
	// long double *term = new_array_f(3);
	// term[0] = 4;
	// term[1] = 5;
	// term[2] = 6;
	// if (term_in_ind(ind,12,term,3))
	// {
	// 	printf("term is in individual\n");
	// }
	// else printf("term is not in individual\n");

	// long double **Ix = new_Ix(lx);
	// long double **A = new_array_f(lx,ly);
	// *(*(A+0)+0) = 1;
	// *(*(A+0)+1) = 1;
	// *(*(A+0)+2) = 1;
	// *(*(A+1)+0) = 13;
	// *(*(A+1)+1) = 16;
	// *(*(A+1)+2) = 2;
	// *(*(A+2)+0) = 7;
	// *(*(A+2)+1) = 8;
	// *(*(A+2)+2) = 9;
	// printf("Matrix A:\n");
	// print(A,lx,ly);
	// printf("inverse of A:\n");
	// print(inverse(A,lx,ly),lx,ly);

	// *(b) = 2;
	// *(b+1) = 2;
	// *(b+2) = 2;
	// printf("Vector b:\n");
	// print(b,lx);
	// long double **terms = get_terms(b,lx,comb);
	// printf("terms(%d,%d):\n",comb);
	// print(terms,comb,lx);

	// long double *c = new_array_f(lx);
	// *(*(A+0)+0) = 1.9;
	// *(*(A+0)+1) = 2.8;
	// *(*(A+0)+2) = 3.7;
	// *(*(A+1)+0) = 4.6;
	// *(*(A+1)+1) = 5.5;
	// *(*(A+1)+2) = 6.4;
	// *(*(A+2)+0) = 7.3;
	// *(*(A+2)+1) = 8.2;
	// *(*(A+2)+2) = 9.1;
	// *(b) = 3;
	// *(b+1) = 2;
	// *(b+2) = 1;
	// printf("Matrix A:\n");
	// print(A,lx,ly);
	// printf("Vector b:\n");
	// print(b,lx);
	// printf("Product b*A:\n");
	// print(mul(b,lx, A,lx,ly),lx);
	// printf("Product A*b:\n");
	// print(mul(A,lx,ly, b,lx),lx);
	
	// to_text((char*)"prueba1.csv",Ix,lx,lx);

	// printf("Matrix A from z3Vars.dat:\n");
	// print(A,lx,ly);
	
	// long double **a = new_array_f(4,5);
	// long double *x = new_array_f(4);
	// *(*(a+0)+0) = 0.1234;
	// *(*(a+0)+1) = 0.5836;
	// *(*(a+0)+2) = 0.3461;
	// *(*(a+0)+3) = 0.0924;
	// *(*(a+0)+4) = 0.0183;
	// *(*(a+1)+0) = 0.4215;
	// *(*(a+1)+1) = 48.0;
	// *(*(a+1)+2) = 17.0;
	// *(*(a+1)+3) = 0.0;
	// *(*(a+1)+4) = -2.51;
	// *(*(a+2)+0) = 0.2341;
	// *(*(a+2)+1) = 0.9124;
	// *(*(a+2)+2) = -0.5424;
	// *(*(a+2)+3) = 0.1246;
	// *(*(a+2)+4) = 0.1307;
	// *(*(a+3)+0) = -2.0;
	// *(*(a+3)+1) = -11.0;
	// *(*(a+3)+2) = 0.13;
	// *(*(a+3)+3) = 0.65;
	// *(*(a+3)+4) = 0.9916;
	// lassol(a,4,5,x);
	// printf("A:\n");
	// print(a, 4,5);
	// printf("\nx:\n");
	// print(x,4);

	return 0;
}


