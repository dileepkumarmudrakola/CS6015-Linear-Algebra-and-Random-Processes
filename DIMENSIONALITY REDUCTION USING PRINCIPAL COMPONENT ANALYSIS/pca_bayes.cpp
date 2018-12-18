/******************************************************************************
FileName: pca_bayes.cpp
Assignment: 4
Author: M DILEEPKUMXR(CS18M031)
Date: 27-11-18
*******************************************************************************/

#include <iostream>
#include <fstream>
#include<cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace std;

//Finds no of lines in text file
int Find_No_Lines(std::ifstream& fp_temp)
{
    int count_lines=0;
    string line="";

    while (getline(fp_temp, line))
        count_lines++;
    fp_temp.close();
    return count_lines;
}

//finds dimension
int Find_No_Words(std::ifstream& fp_temp, int lines)
{
    int count_words=0;
    string words="";

    while (!fp_temp.eof())
    {
        fp_temp>>words;
        count_words++;
    }
    count_words=count_words/lines;
    fp_temp.close();
    return count_words;
}

//Finds the mean vector
void computeMean(double **temp_class, double *mean, int D, int N)
{

    for(int i=0; i<D; i++)
    {
        double sum=0.0;
        for(int j=0; j<N; j++)
        {
           sum+=temp_class[j][i];
        }
        mean[i]=sum/N;
    }
}

//Exchange row when needed in finding determinant
void exchangeRow(double **covar, int D, int i, int j)
{
    double temp=0.0;
    for(int l=0; l<D; l++)
    {
        temp=covar[i][l];
        covar[i][l]=covar[j][l];
        covar[j][l]=temp;
    }
}

//Finds determinant of given covariance matrix
double FindDeterminant(double **covar, int D)
{
    int No_of_exchanges=0;
    for(int i=0; i<D; i++)
    {
        //checking for max value in pivot
        for(int l=i+1; l<D; l++)
        {
            if(abs(covar[i][i]) < abs(covar[l][i]))
            {
                No_of_exchanges++;
                exchangeRow(covar, D, i, l);
            }
        }

        //checking pivot 0;
        if(covar[i][i]==0)
        {
            return 0;
        }
        else
        {
            for(int k=i+1; k<D; k++)
            {
                double temp=covar[k][i]/covar[i][i];
                for(int k1=0; k1<D; k1++)
                {
                    covar[k][k1]=covar[k][k1]-temp*covar[i][k1];
                }
            }
        }
    }

    double determinant=1.0;
    for(int i=0; i<D; i++)
    {
        determinant=determinant*covar[i][i];
    }
    if(No_of_exchanges%2!=0)
    {
        determinant=(-1)*determinant;
    }

    return determinant;
}

//Computes the covariance of given matrix
void ComputeCovariance(double **covar, double **tempclass, double *mean, int D, int N)
{
    for(int i=0; i<D; i++)
    {
        for(int j=0; j<D; j++)
        {
            if(i>j)
            {
                covar[i][j]=covar[j][i];
            }
            else
            {
                double temp1, temp2, sum=0.0;
                for(int k=0; k<N; k++)
                {
                    temp1=tempclass[k][i]-mean[i];
                    temp2=tempclass[k][j]-mean[j];
                    sum+=temp1*temp2;
                }
                covar[i][j]=sum/N;
            }

        }
    }
}

//Finds the inverse of covariance matrix
void FindInverse(double **InvCov, double **covar, int D, double det)
{
    double **tempInv=new double*[D];
    for(int i=0; i<D; i++)
    {
        tempInv[i]=new double[D];
    }

    for(int i=0; i<D; i++)
    {
        for(int j=0; j<D; j++)
        {
            double ** tempcovar=new double*[D-1];
            for(int k=0; k<D-1; k++)
            {
                tempcovar[k]=new double[D-1];
            }
            int p=0,l=0,m=0;
            while(p<D)
            {
                    m=0;
                    if(p!=i)
                    {
                        for(int k=0; k<D; k++)
                        {
                            if(k!=j)
                            {
                                tempcovar[l][m]=covar[p][k];
                                m++;
                            }
                        }
                        l++;
                    }

                p++;
            }

            tempInv[i][j]=FindDeterminant(tempcovar, D-1) * pow(-1, i+j);

        }
    }

    for(int i=0; i<D; i++)
    {
        for(int j=0; j<D; j++)
        {
            InvCov[i][j]=tempInv[j][i]/det;
        }
    }
}

//Finds the class of given data
int FindClass(double det1, double det2, double **inv1, double **inv2, double **totol_x, double *mean1, double *mean2, int D, int N, std::ofstream& fp_output)
{
    int count=0;
    for(int y=0; y<N; y++)
    {
        int x[D+1];
        for(int z=0; z<=D; z++)
        {
            x[z]=totol_x[y][z];
        }
        double x1[D], x2[D];
        for(int i=0; i<D; i++)
        {
            x1[i]=x[i]-mean1[i];
            x2[i]=x[i]-mean2[i];
        }
        double tempreult1[D]={0.0}, tempreult2[D]={0.0};
        for(int i=0; i<D; i++)
        {
            for(int j=0; j<D; j++)
            {
                tempreult1[i]+=x1[i]*inv1[i][j];
                tempreult2[i]+=x2[i]*inv2[i][j];
            }
        }
        double temp1=0.0, temp2=0.0;
        for(int i=0; i<D; i++)
        {
            temp1+=tempreult1[i]*x1[i];
            temp2+=tempreult2[i]*x2[i];
        }

        temp1=(-0.5)*temp1;
        temp2=(-0.5)*temp2;

        temp1=pow(2.71, temp1);
        temp2=pow(2.71, temp2);

        temp1=temp1/(pow(det1,0.5));
        temp2=temp2/(pow(det2,0.5));

        if(temp1>temp2)
        {
            fp_output<< "1 ";
            if(x[D]==1)
                count++;
        }
        else
        {
            fp_output<<"2 ";
            if(x[D]==2)
                count++;
        }

    }
    return count;
}


//Converts N*D to N*2 Dimension data
void Multiply(double **C, double **A, double **B, long int m_A, long int n_B, long int p_B)
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


int main(int argc, char **argv)
{
    ifstream fp_input1, fp_input2, fp_temp, fp_test;
    fp_temp.open("class1_train.txt");
    int N1=Find_No_Lines(fp_temp);
    fp_temp.open("class2_train.txt");
    int N2=Find_No_Lines(fp_temp);
    fp_temp.open("class1_train.txt");
    int D=Find_No_Words(fp_temp, N1);

    //TRAINING SET
    // int N1, N2, N1_Test,N2_Test;
    // N1=N1_total*80/100;
    // N2=N2_total*80/100;

    // //Test
    // N1_Test=N1_total*20/100;
    // N2_Test=N2_total*20/100;


    //Declaration of Training 
    double **tclass=new double*[N1+N2];
    for(int i=0; i<N1+N2; i++)
    {
        tclass[i]=new double[D];
    }

    double **tempclass1=new double*[N1];
    for(int i=0; i<N1; i++)
    {
        tempclass1[i]=new double[3];
    }

    double **tempclass2=new double*[N2];
    for(int i=0; i<N2; i++)
    {
        tempclass2[i]=new double[3];
    }

    double **tempclass=new double*[N1+N2];
    for(int i=0; i<N1+N2; i++)
    {
        tempclass[i]=new double[2];
    }

    //Declaration of Test
    // double **tclasstest=new double*[N1_Test+N2_Test];
    // for(int i=0; i<N1_Test+N2_Test; i++)
    // {
    //     tclasstest[i]=new double[D];
    // }
    // double ** class1_Test=new double*[N1_Test];
    // for(int i=0; i<N1_Test; i++)
    // {
    //     class1_Test[i]=new double[3];
    // }
    // double **class2_Test=new double*[N2_Test];
    // for(int i=0; i<N2_Test; i++)
    // {
    //     class2_Test[i]=new double[3];
    // }
    // double **tempclasstest=new double*[N1_Test+N2_Test];
    // for(int i=0; i<N1_Test+N2_Test; i++)
    // {
    //     tempclasstest[i]=new double[2];
    // }


    //Declaration of EigenVector Matrix
    double **eigvec=new double*[D];
    for(int i=0; i<D; i++)
    {
        eigvec[i]=new double[2];
    }

    //taking input from file for training

    fp_input1.open("class1_train.txt");
    for(int i=0; i<N1; i++)
    {
        for(int j=0; j<D; j++)
        {
            fp_input1>>tclass[i][j];
        }
    }

    fp_input2.open("class2_train.txt");
    for(int i=0; i<N2; i++)
    {
        for(int j=0; j<D; j++)
        {
            fp_input2>>tclass[i+N1][j];
        }
    }

    //taking input from file for test
    // for(int i=0; i<N1_Test; i++)
    // {
    //     for(int j=0; j<D; j++)
    //     {
    //         fp_input1>>tclasstest[i][j];
    //     }
    // }
    fp_input1.close();


    // for(int i=0; i<N2_Test; i++)
    // {
    //     for(int j=0; j<D; j++)
    //     {
    //         fp_input2>>tclasstest[i+N1_Test][j];
    //     }
    // }
    fp_input2.close();


    double mean1[2], mean2[2], tmean[D];
    //computing Mean
    computeMean(tclass, tmean, D, N1+N2);

    //covariance Matrix declaration
    double ** tcovar=new double*[D];
    for(int i=0; i<D; i++)
    {
        tcovar[i]=new double[D];
    }

    //Computing Covariance
    ComputeCovariance(tcovar, tclass, tmean, D, N1+N2);

    gsl_matrix *tcovmatrix = gsl_matrix_alloc (D, D);
    for(int i=0; i<D; i++)
        {
            for(int j=0; j<D; j++)
            {
               gsl_matrix_set(tcovmatrix, i, j, tcovar[i][j]);
            }
        }

  gsl_vector *eigenvalues = gsl_vector_alloc (D);
  gsl_matrix *eigenvectors = gsl_matrix_alloc (D, D);

  gsl_eigen_symmv_workspace * w = 
    gsl_eigen_symmv_alloc (D);
  
  gsl_eigen_symmv (tcovmatrix, eigenvalues, eigenvectors, w);

  gsl_eigen_symmv_free (w);

//sorts eigen values according to most significant eigen vectors
  gsl_eigen_symmv_sort (eigenvalues, eigenvectors, 
                        GSL_EIGEN_SORT_VAL_DESC);

    for(int i=0; i<D; i++)
        {
            for(int j=0; j<2; j++)
            {
                eigvec[i][j]=gsl_matrix_get(eigenvectors, i, j);
            }
        }

	
    //Dimensionality reduction tempclass consist of training data of dimension2 adn tempclasstess consist of testing data of dimension 2
    Multiply(tempclass, tclass, eigvec, N1+N2, D, 2);
    //Multiply(tempclasstest, tclasstest, eigvec, N1_Test+N2_Test, D, 2);

    //copying reduced dimension data to individual matrices
    //Training
    for(int i=0; i<N1; i++)
    {
        tempclass1[i][0]=tempclass[i][0];
        tempclass1[i][1]=tempclass[i][1];
        tempclass1[i][2]=1;
    }

    for(int i=0; i<N2; i++)
    {
        tempclass2[i][0]=tempclass[i+N1][0];
        tempclass2[i][1]=tempclass[i+N1][1];
        tempclass2[i][2]=2;
    }

    //Testing
    // for(int i=0; i<N1_Test; i++)
    // {
    //     class1_Test[i][0]=tempclasstest[i][0];
    //     class1_Test[i][1]=tempclasstest[i][1];
    //     class1_Test[i][2]=1;
    // }

    // for(int i=0; i<N2_Test; i++)
    // {
    //     class2_Test[i][0]=tempclasstest[i+N1_Test][0];
    //     class2_Test[i][1]=tempclass[i+N1_Test][1];
    //     class2_Test[i][2]=2;
    // }



    computeMean(tempclass1, mean1, 2, N1);
    computeMean(tempclass2, mean2, 2, N2);

    double ** covar1=new double*[2];
    double ** covar2=new double*[2];
   
    double ** Invcovar1=new double*[2];
    double ** Invcovar2=new double*[2];
    double ** tInvcovar=new double*[2];
    double ** tempcovar=new double*[2];
    for(int i=0; i<2; i++)
    {
        covar1[i]=new double[2];
        covar2[i]=new double[2];
        //tcovar[i]=new double[D];
        Invcovar1[i]=new double[2];
        Invcovar2[i]=new double[2];
        tInvcovar[i]=new double[2];
        tempcovar[i]=new double[2];
    }

    //Computing Covariance
    ComputeCovariance(covar1, tempclass1, mean1, 2, N1);
    ComputeCovariance(covar2, tempclass2, mean2, 2, N2);

    // computing Determinants
    double Det1, Det2;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
           tempcovar[i][j]=covar1[i][j];
        }
    }
    Det1=FindDeterminant(tempcovar, 2);

    for(int i=0; i<2; i++)
    {
        for(int j=0; j<2; j++)
        {
           tempcovar[i][j]=covar2[i][j];
        }

    }
    Det2=FindDeterminant(tempcovar, 2);




    FindInverse(Invcovar1, covar1, 2, Det1);
    FindInverse(Invcovar2, covar2, 2, Det2);

    fp_test.open(argv[1]);
    int N_Test=Find_No_Lines(fp_test);
    fp_test.close();

    //Declaration of Test
    double ** class_Test=new double*[N_Test];
    for(int i=0; i<N_Test; i++)
    {
        class_Test[i]=new double[D];
    }

    
    fp_test.open(argv[1]);
    //taking input from file for test
    for(int i=0; i<N_Test; i++)
    {
        for(int j=0; j<D; j++)
        {
            fp_test>>class_Test[i][j];
        }
        class_Test[i][D]=1;
    }
    fp_test.close();



    ofstream fp_output;
    fp_output.open("output.txt");
    
    FindClass(Det1, Det2, Invcovar1, Invcovar2, class_Test, mean1, mean2, 2, N_Test, fp_output); 

   
    // ofstream fp_output;
    // fp_output.open("output.txt");
    // double count=0;
    // count+=FindClass(Det1, Det2, Invcovar1, Invcovar2, tempclass1, mean1, mean2, 2, N1, fp_output);
    // fp_output<<" ";
    // count+=FindClass(Det1, Det2, Invcovar1, Invcovar2, tempclass2, mean1, mean2, 2, N2, fp_output);
    // double acc=((count)/(N1+N2))*100;
    // fp_output<<"\n"<<acc<<"\n";

    // count=0;
    // count+=FindClass(Det1, Det2, Invcovar1, Invcovar2, class1_Test, mean1, mean2, 2, N1_Test, fp_output);
    // fp_output<<" ";
    // count+=FindClass(Det1, Det2, Invcovar1, Invcovar2, class2_Test, mean1, mean2, 2, N2_Test, fp_output);
    // acc=((count)/(N1_Test+N2_Test))*100;
    // fp_output<<"\n"<<acc<<"\n";
    return 0;
}

