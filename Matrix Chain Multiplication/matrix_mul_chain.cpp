/******************************************************************************
FileName: matrix_chain_mul.cpp
Assignment: Assignment-2
Author: M DILEEPKUMAR(CS18M031)
Description: This file finds no of scalar multiplications by Optimal, Left to right sequence of multiplications of matrices and Non optimal pairing
*******************************************************************************/


#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <bits/stdc++.h>

using namespace std;

//optpaimal stores the optimum sequence of pairing for multiplication of matrices
int optpar[100000];
void PrintOptimalParenthesis(std::ofstream& fp_output , int start, int end, int &matNo, int **optimalbreak)
{
    static int ind=0;
    if(start==end)
    {
        int a=matNo++;
        fp_output<<a;
        optpar[ind++]=a;
        return;
    }
    else
    {
        fp_output<<"(";
        optpar[ind++]=-1;
        PrintOptimalParenthesis(fp_output, start, optimalbreak[start][end], matNo, optimalbreak);
        fp_output<<",";
        optpar[ind++]=-3;
        PrintOptimalParenthesis(fp_output, optimalbreak[start][end]+1, end, matNo, optimalbreak);
        fp_output<<")";
        optpar[ind++]=-2;
    }
}

//findOptimum finds the no of scalar multiplication for he optimum sequence of pairing for multiplication of matrices
int findOptimum(int *d, int **optimalbreak, int noofmatrices)
{
    int mul[noofmatrices+1][noofmatrices+1];

    for(int i=1; i<=noofmatrices; i++)
    {
        mul[i][i]=0;
    }
    for(int length=2; length<=noofmatrices; length++)
    {
        for(int i=1; i<=noofmatrices-length+1; i++)
        {
            int j=i+length-1;
            mul[i][j]=INT_MAX;
            for(int k=i; k<j; k++)
            {
                int cost=mul[i][k]+mul[k+1][j] + d[i-1]*d[k]*d[j];
                if(mul[i][j]>cost)
                {
                    mul[i][j]=cost;
                    optimalbreak[i][j]=k;
                }
            }
        }
    }
    return mul[1][noofmatrices];
}

//findNonOptimumPairing finds sequence of matrix multiplication by performing a sequence of non-optimal pairing from left to right
int matno=1;
int nonoptpar[100000];
static int ind1=0;
int findNonOptimumPairing(int *d,int start, int endm, int N)
{

    int mul=0;
    if(start+1==endm && endm<=N)
    {
        nonoptpar[ind1++]=matno;
        matno+=1;
        nonoptpar[ind1++]=-3;
        nonoptpar[ind1++]=matno;
        matno+=1;
        return d[start-1]*d[start]*d[endm];
    }
    if(start==endm)
    {
        nonoptpar[ind1++]=matno;
        matno+=1;
        return 0;
    }
    int mid=(start+endm)/2;
        nonoptpar[ind1++]=-1;
        mul+=findNonOptimumPairing( d, start, mid, N);
        nonoptpar[ind1++]=-2;
         nonoptpar[ind1++]=-1;
        if(endm<=N)
        {
            mul+=findNonOptimumPairing( d, mid+1, endm, N)+d[start-1]*d[mid]*d[endm];
        }
        else
        {
            mul+=findNonOptimumPairing( d, mid+1, N, N)+d[start-1]*d[mid]*d[N];
        }
        nonoptpar[ind1++]=-2;
    return mul;
}

//findSequenceScalar finds the no of scalar multiplication for  sequence of pairing for multiplication of matrices from left to right

int sequence[100000];

int findSequenceScalar(int *d, int noofmatrices)
{
    int cost=0;
    int ind=0;
    for(int i=2; i<=noofmatrices; i++)
    {
        cost=cost+d[0]*d[i-1]*d[i];
        sequence[ind++]=-1;
    }
    int matno=1;
    if(noofmatrices>=2)
    {
        int t=matno++;
        sequence[ind++]=t;
        sequence[ind++]=-3;
        t=matno++;
        sequence[ind++]=t;
         sequence[ind++]=-2;
        for(;matno<=noofmatrices; matno++)
        {
            sequence[ind++]=matno;
            sequence[ind++]=-2;
        }
    }

    return cost;
}


//for multiplying two matrices

void Multiply(int **A, int **B, long int m_A, long int n_B, long int p_B)
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
		}
	}
}

int top;
int temopstack[100000];

void push(int m, int n)
{
    temopstack[++top]=m;
    temopstack[++top]=n;
}

int pop()
{
    if(top!=-1)
    {
        return temopstack[top--];
    }
    else
    {
        return 0;
    }
}

void multiply(int *d, int *seq)
{
     top=-1;
     int i=0;
     while(seq[i]!='\0')
     {
          if(seq[i]==-2)
         {
             int rn=pop();
             int rm=pop();
             int ln=pop();
             int lm=pop();
            int **A=new int *[lm];
            for(int i=0; i<lm; i++)
            {
                A[i]=new int[ln];
                for(int j=0; j<ln; j++)
                {
                    A[i][j]=100+rand()%50;
                }
            }

            int **B=new int *[rm];
            for(int i=0; i<rm; i++)
            {
                B[i]=new int[rn];
                for(int j=0; j<rn; j++)
                {
                    B[i][j]=100+rand()%50;
                }
            }
             Multiply(A,B, lm, ln, rn);
             push(lm, rn);
        }
        else if(seq[i]!=-1 && seq[i]!=-3)
         {
             push(d[seq[i]-1], d[seq[i]]);
         }
        i++;
     }
}

int main()
{
    int N;
    //fp_input pointer to file for reading input file
    ifstream fp_input;
    //fp_ouptut pointer to file for writing to output file
    ofstream fp_output;
    fp_input.open("input.txt");
    fp_output.open("output.txt");
    if (!fp_input)
    {
        cout << "Unable to open file";
        return 1;
    }

    fp_input>>N;
    int d[N+1];

    for(int i=0; i<N+1; i++)
    {
        fp_input>>d[i];
    }

    //finding optimal scalar multiplications
    int **optimalbreak =new int *[N+1];
    for(int i=0; i<N+1; i++)
    {
        optimalbreak[i]=new int[N+1];
    }
    int optimal=findOptimum(d, optimalbreak, N);
    int matNo=1;

    //optimum pairing sequence
    PrintOptimalParenthesis(fp_output,  1, N, matNo, optimalbreak);
    fp_output<<"\n";

    clock_t start_opt, end_opt, start_seq, end_seq, start_non, end_non;
    float cpu_time_used_opt, cpu_time_used_seq, cpu_time_used_non;

    //finding optimal matrix multiplication time
    start_opt = clock();
    multiply(d, optpar);
    end_opt = clock();
    cpu_time_used_opt = ((float) (end_opt - start_opt)) / CLOCKS_PER_SEC;

    //sequence

    //sequencemul is no of scalare multiplications taken by sequence of matrix multiplication from left to right
    long int sequencemul=findSequenceScalar(d,N);;
    fp_output<<sequencemul;
    fp_output<<"\n";

    //finding sequence matrix multiplication time
    start_seq = clock();
    multiply(d, sequence);
    end_seq = clock();
    cpu_time_used_seq = ((float) (end_seq - start_seq)) / CLOCKS_PER_SEC;


    //non optimal

    int twopow=ceil(log2(N));
    twopow=pow(2, twopow);
    nonoptpar[ind1++]=-1;
    //nonoptimal is no of scalare multiplications taken by sequence of non optimal pairing matrix multiplication
    long int nonoptimal=findNonOptimumPairing(d,1, twopow, N);
    fp_output<<nonoptimal;
    nonoptpar[ind1++]=-2;
     fp_output<<"\n";

    //finding nonoptimal matrix multiplication time
     start_non= clock();
   multiply(d, nonoptpar);
    end_non = clock();
    cpu_time_used_non = ((float) (end_non - start_non)) / CLOCKS_PER_SEC;

    //optimal value
    fp_output<<optimal;

    fp_output<<"\n";
    fp_output<<setprecision(4)<<cpu_time_used_seq;
    fp_output<<"\n";
    fp_output<<setprecision(4)<<cpu_time_used_non;
    fp_output<<"\n";
    fp_output<<setprecision(4)<<cpu_time_used_opt;

    return 0;
}
