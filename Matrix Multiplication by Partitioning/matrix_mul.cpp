/******************************************************************************
FileName: matrix_mul.cpp
Assignment: Assignment-1
Author: M DILEEPKUMAR(CS18M031)
Description: This file multiplies matrix A and B using naive method, single Partition Method and Strassen Matrix Multiplication method and determines
time taken by those methods and computes forbenious norm between naive and strssen resulted matrices
*******************************************************************************/

#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
using namespace std;



void DivideMatrix(float **temp, float **temp11, float **temp12, float **temp21, float **temp22, long int max1)
{
	
	for(int i=0; i<max1; i++)
	{
		for(int j=0; j<max1; j++)
		{

			temp11[i][j]=temp[i][j];
			temp12[i][j]=temp[i][j+max1];
			temp21[i][j]=temp[i+max1][j];
			temp22[i][j]=temp[i+max1][j+max1];
		}
	}
}


void MergeMatrix(float **temp, float **temp11, float **temp12, float **temp21, float **temp22, long int max1)
{
	for(int i=0; i<max1; i++)
	{
		for(int j=0; j<max1; j++)
		{
			temp[i][j]=temp11[i][j];
			temp[i][j+max1]=temp12[i][j];
			temp[i+max1][j]=temp21[i][j];
			temp[i+max1][j+max1]=temp22[i][j];
		}
	}
}

void AddMatrix(float **temp, float **temp1, float **temp2, long int max1)
{
	for(int i=0; i<max1; i++)
	{
		for(int j=0; j<max1; j++)
		{
			temp[i][j]=temp1[i][j]+temp2[i][j];
		}
	}
}

void SubMatrix(float **temp, float **temp1, float **temp2, long int max1)
{
	for(int i=0; i<max1; i++)
	{
		for(int j=0; j<max1; j++)
		{
			temp[i][j]=temp1[i][j]-temp2[i][j];
		}
	}
}





void Multiply(float **C, float **A, float **B, long int m_A, long int n_B, long int p_B)
{
	for(int i=0; i< m_A; i++)
	{
		for(int j=0; j< p_B; j++)
		{
			double temp_result=0.0;
			for(int k=0; k<n_B; k++)
			{
				temp_result= temp_result+ A[i][k]*B[k][j];	
			}
			C[i][j]=temp_result;
		}
	}
}

void WriteOutput(std::ofstream& fp_matC, float **C, long int m_A, long int p_B)
{
	fp_matC << m_A << " " << p_B <<"\n";
	for(int i=0; i<m_A; i++)
	{
		for(int j=0; j<p_B; j++)
		{
			fp_matC << C[i][j] <<" ";
		}
		fp_matC<<"\n";
	}

}


void Strassen(float **A, float **B, float **C, long int max1)
{
	
	if(max1==1)
	{
		//Multiply(A,B,C, 4);
		//cout<<C[0][0]<<" ";
		C[0][0]=A[0][0]*B[0][0];
	}
	else
	{
		float **A11 = new float *[max1/2];
		float **A12 = new float *[max1/2];
		float **A21 = new float *[max1/2];
		float **A22 = new float *[max1/2];
		float **B11 = new float *[max1/2];
		float **B12 = new float *[max1/2];
		float **B21 = new float *[max1/2];
		float **B22 = new float *[max1/2];
		float **C11 = new float *[max1/2];
		float **C12 = new float *[max1/2];
		float **C21 = new float *[max1/2];
		float **C22 = new float *[max1/2];
		float **S1 = new float *[max1/2];
		float **S2 = new float *[max1/2];
		float **S3 = new float *[max1/2];
		float **S4 = new float *[max1/2];
		float **S5 = new float *[max1/2];
		float **S6 = new float *[max1/2];
		float **S7 = new float *[max1/2];
		float **S8 = new float *[max1/2];
		float **S9 = new float *[max1/2];
		float **S10 = new float *[max1/2];
		float **P1 = new float *[max1/2];
		float **P2 = new float *[max1/2];
		float **P3 = new float *[max1/2];
		float **P4 = new float *[max1/2];
		float **P5 = new float *[max1/2];
		float **P6 = new float *[max1/2];
		float **P7 = new float *[max1/2];



	

		for(int i=0; i<max1/2; i++)
		{
			A11[i] =new float[max1/2];
			A12[i] =new float[max1/2];
			A21[i] =new float[max1/2];
			A22[i] =new float[max1/2];
			B11[i] =new float[max1/2];
			B12[i] =new float[max1/2];
			B21[i] =new float[max1/2];
			B22[i] =new float[max1/2];
			C11[i] =new float[max1/2];
			C12[i] =new float[max1/2];
			C21[i] =new float[max1/2];
			C22[i] =new float[max1/2];
			S1[i] =new float[max1/2];
			S2[i] =new float[max1/2];
			S3[i] =new float[max1/2];
			S4[i] =new float[max1/2];
			S5[i] =new float[max1/2];
			S6[i] =new float[max1/2];
			S7[i] =new float[max1/2];
			S8[i] =new float[max1/2];
			S9[i] =new float[max1/2];
			S10[i] =new float[max1/2];
			P1[i] =new float[max1/2];
			P2[i] =new float[max1/2];
			P3[i] =new float[max1/2];
			P4[i] =new float[max1/2];
			P5[i] =new float[max1/2];
			P6[i] =new float[max1/2];
			P7[i] =new float[max1/2];
		}
		DivideMatrix(A, A11, A12, A21, A22, max1/2);
		DivideMatrix(B, B11, B12, B21, B22, max1/2);
		
		
		
		SubMatrix(S1, B12, B22, max1/2);
		AddMatrix(S2, A11, A12, max1/2);
		AddMatrix(S3, A21, A22, max1/2);
		SubMatrix(S4, B21, B11, max1/2);
		AddMatrix(S5, A11, A22, max1/2);
		AddMatrix(S6, B11, B22, max1/2);
		SubMatrix(S7, A12, A22, max1/2);
		AddMatrix(S8, B21, B22, max1/2);
		SubMatrix(S9, A11, A21, max1/2);
		AddMatrix(S10, B11, B12, max1/2);
		
		
		Strassen(A11, S1,P1, max1/2);
		
		Strassen(S2, B22,P2, max1/2);
		
		Strassen(S3, B11,P3, max1/2);
		Strassen(A22, S4,P4, max1/2);
		Strassen(S5, S6,P5, max1/2);
		Strassen( S7, S8,P6, max1/2);
		Strassen( S9, S10,P7, max1/2);
	
		
		AddMatrix(C11, P5, P4, max1/2);
		SubMatrix(C11, C11, P2, max1/2);
		AddMatrix(C11, C11, P6, max1/2);
		AddMatrix(C12, P1, P2, max1/2);
		AddMatrix(C21, P3, P4, max1/2);
		AddMatrix(C22, P5, P1, max1/2);
		SubMatrix(C22, C22, P3, max1/2);
		SubMatrix(C22, C22, P7, max1/2);
	
	
		MergeMatrix(C, C11, C12, C21, C22, max1/2);
		delete(A11);delete(A12);delete(A21);delete(A22);
		delete(B11);delete(B12);delete(B21);delete(B22);
		delete(C11);delete(C12);delete(C21);delete(C22);
		delete(S1);delete(S2);delete(S3);delete(S4);delete(S5);
		delete(S6);delete(S7);delete(S8);delete(S9);delete(S10);
		delete(P1);delete(P2);delete(P3);delete(P4);delete(P5);
		delete(P6);delete(P7);
		

	}
}



int main()
{
   	ifstream fp_matA, fp_matB; 
	ofstream fp_matC_N, fp_matC_P,fp_matC_RP, fp_output;

	fp_matA.open("matA.txt");
	fp_matB.open("matB.txt");
	fp_matC_N.open("matC_N.txt");
	fp_matC_P.open("matC_P.txt");
	fp_matC_RP.open("matC_RP.txt");
	
	fp_output.open("output.txt");
	
    	if (!fp_matA && !fp_matB) {
        	cout << "Unable to open file";
        	return 1; // terminate with error
   	 }
	
	double temp,temp1;
	long int m_A, n_A, n_B, p_B, maxA, maxB, max1;

	
	
	
	
	fp_matA >> m_A;
	fp_matA >> n_A;
	fp_matB >> n_B;
	fp_matB >> p_B;
	if(n_A!=n_B)
	{
		cout<<"Matrix multiplication not compatable";
		return 1;
	}
	
	maxA=m_A>n_A?m_A:n_A;
	maxB=n_B>p_B?n_B:p_B;
	max1=maxA>maxB?maxA:maxB;
	if(max1!=0 && max1-1!=0)
	{
		max1=ceil(log2(max1));
		max1=pow(2, max1);
	}
	


	float **A = new float *[max1];
	float **B = new float *[max1];
	float **C_N = new float *[max1];
	float **C_P = new float *[max1];
	float **C_RP = new float *[max1];
	float **C_F = new float *[max1];
	
	for(int i=0; i<max1; i++)
	{
		A[i] =new float[max1];
		B[i] =new float[max1];
		C_N[i] =new float[max1];
		C_P[i] =new float[max1];
		C_RP[i] =new float[max1];
		C_F[i] =new float[max1];
		for(int j=0; j<max1; j++)
		{
			A[i][j]=0;
			B[i][j]=0;
			C_N[i][j]=0;
			C_P[i][j]=0;
			C_RP[i][j]=0;
			C_F[i][j]=0;
		}
		
	}
	
	
	
	
	for(int i=0; i<m_A; i++)
	{
		for(int j=0; j<n_A; j++)
		{
			fp_matA >> A[i][j];
		}
	}

	for(int i=0; i<n_B; i++)
	{
		for(int j=0; j<p_B; j++)
		{
			fp_matB >> B[i][j];
		}
	}

	clock_t start_N, end_N, start_P, end_P, start_RP, end_RP;
	double cpu_time_used_N, cpu_time_used_P, cpu_time_used_RP;
	
	//Naive Method

	
	start_N = clock();

	Multiply(C_N, A, B, m_A, n_B, p_B);
	

	end_N = clock();
	cpu_time_used_N = ((double) (end_N - start_N)) / CLOCKS_PER_SEC;
	fp_output<<cpu_time_used_N<<"\n";

	WriteOutput(fp_matC_N, C_N, m_A, p_B);
	
	//Single Level Partition
	
	start_P = clock();

	float **A11 = new float *[max1/2];
	float **A12 = new float *[max1/2];
	float **A21 = new float *[max1/2];
	float **A22 = new float *[max1/2];
	float **B11 = new float *[max1/2];
	float **B12 = new float *[max1/2];
	float **B21 = new float *[max1/2];
	float **B22 = new float *[max1/2];
	float **C11 = new float *[max1/2];
	float **C12 = new float *[max1/2];
	float **C21 = new float *[max1/2];
	float **C22 = new float *[max1/2];
	float **S1 = new float *[max1/2];
	float **S2 = new float *[max1/2];
	float **S3 = new float *[max1/2];
	float **S4 = new float *[max1/2];
	float **S5 = new float *[max1/2];
	float **S6 = new float *[max1/2];
	float **S7 = new float *[max1/2];
	float **S8 = new float *[max1/2];
	float **S9 = new float *[max1/2];
	float **S10 = new float *[max1/2];
	float **P1 = new float *[max1/2];
	float **P2 = new float *[max1/2];
	float **P3 = new float *[max1/2];
	float **P4 = new float *[max1/2];
	float **P5 = new float *[max1/2];
	float **P6 = new float *[max1/2];
	float **P7 = new float *[max1/2];
	
	for(int i=0; i<max1/2; i++)
	{
		A11[i] =new float[max1/2];
		A12[i] =new float[max1/2];
		A21[i] =new float[max1/2];
		A22[i] =new float[max1/2];
		B11[i] =new float[max1/2];
		B12[i] =new float[max1/2];
		B21[i] =new float[max1/2];
		B22[i] =new float[max1/2];
		C11[i] =new float[max1/2];
		C12[i] =new float[max1/2];
		C21[i] =new float[max1/2];
		C22[i] =new float[max1/2];
		S1[i] =new float[max1/2];
		S2[i] =new float[max1/2];
		S3[i] =new float[max1/2];
		S4[i] =new float[max1/2];
		S5[i] =new float[max1/2];
		S6[i] =new float[max1/2];
		S7[i] =new float[max1/2];
		S8[i] =new float[max1/2];
		S9[i] =new float[max1/2];
		S10[i] =new float[max1/2];
		P1[i] =new float[max1/2];
		P2[i] =new float[max1/2];
		P3[i] =new float[max1/2];
		P4[i] =new float[max1/2];
		P5[i] =new float[max1/2];
		P6[i] =new float[max1/2];
		P7[i] =new float[max1/2];
	}

	DivideMatrix(A, A11, A12, A21, A22, max1/2);
	DivideMatrix(B, B11, B12, B21, B22, max1/2);

	SubMatrix(S1, B12, B22, max1/2);
	AddMatrix(S2, A11, A12, max1/2);
	AddMatrix(S3, A21, A22, max1/2);
	SubMatrix(S4, B21, B11, max1/2);
	AddMatrix(S5, A11, A22, max1/2);
	AddMatrix(S6, B11, B22, max1/2);
	SubMatrix(S7, A12, A22, max1/2);
	AddMatrix(S8, B21, B22, max1/2);
	SubMatrix(S9, A11, A21, max1/2);
	AddMatrix(S10, B11, B12, max1/2);

	Multiply(C_N, A, B, m_A, n_B, p_B);
	Multiply(P1, A11, S1, max1/2, max1/2, max1/2);
	Multiply(P2, S2, B22, max1/2, max1/2, max1/2);
	Multiply(P3, S3, B11, max1/2, max1/2, max1/2);
	Multiply(P4, A22, S4, max1/2, max1/2, max1/2);
	Multiply(P5, S5, S6, max1/2, max1/2, max1/2);
	Multiply(P6, S7, S8, max1/2, max1/2, max1/2);
	Multiply(P7, S9, S10, max1/2, max1/2, max1/2);
	

	AddMatrix(C11, P5, P4, max1/2);
	SubMatrix(C11, C11, P2, max1/2);
	AddMatrix(C11, C11, P6, max1/2);
	AddMatrix(C12, P1, P2, max1/2);
	AddMatrix(C21, P3, P4, max1/2);
	AddMatrix(C22, P5, P1, max1/2);
	SubMatrix(C22, C22, P3, max1/2);
	SubMatrix(C22, C22, P7, max1/2);
	
	
	MergeMatrix(C_P, C11, C12, C21, C22, max1/2);	
	
		delete(A11);delete(A12);delete(A21);delete(A22);
		delete(B11);delete(B12);delete(B21);delete(B22);
		delete(C11);delete(C12);delete(C21);delete(C22);
		delete(S1);delete(S2);delete(S3);delete(S4);delete(S5);
		delete(S6);delete(S7);delete(S8);delete(S9);delete(S10);
		delete(P1);delete(P2);delete(P3);delete(P4);delete(P5);
		delete(P6);delete(P7);
	

	end_P = clock();
	cpu_time_used_P = ((double) (end_P - start_P)) / CLOCKS_PER_SEC;
	
	fp_output<<cpu_time_used_P<<"\n";

	WriteOutput(fp_matC_P, C_P, m_A, p_B);

	//Stratssen method

	start_RP = clock();		

	Strassen(A, B, C_RP, max1);
	

	end_RP = clock();
	cpu_time_used_RP = ((double) (end_RP - start_RP)) / CLOCKS_PER_SEC;
	
	fp_output<<cpu_time_used_RP<<"\n";
	//cout<<cpu_time_used_N<<"\n"<<cpu_time_used_P<<"\n"<<cpu_time_used_RP;
	WriteOutput(fp_matC_RP, C_RP, m_A, p_B);

	SubMatrix(C_F, C_N, C_RP, max1);
	double forb_norm=0.0;
	for(int i=0; i<m_A;i++)
	{
		for(int j=0; j<p_B;j++)
		{
			forb_norm+=(C_F[i][j] * C_F[i][j]);
		}
	}
	
	forb_norm= sqrt(forb_norm);
	fp_output<<forb_norm<<"\n";
	//cout<<"\n"<<forb_norm;
	
   
    return 1;
}




