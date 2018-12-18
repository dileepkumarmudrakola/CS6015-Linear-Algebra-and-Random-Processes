#include <iostream>
#include <fstream>
#include<cmath>
using namespace std;

//Finds the no of data inputs in a given file
int Find_No_Lines(std::ifstream& fp_temp)
{
    int count_lines=0;
    string line="";

    while (getline(fp_temp, line))
        count_lines++;
    fp_temp.close();
    return count_lines;
}

//Finds the Dimension
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

//Finds the Mean of the features
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

//Swaps the rows in gaussian elimination
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

//Finds the determinant of covariance Matrix
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

//Compute the covariance feature Matrix
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

//Finds the inverse of covariance Matrix
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

//Finds the class of given Test file Data inputs and prints to the output file
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


int main(int argc,char **argv)
{
    ifstream fp_input1, fp_input2, fp_temp,fp_test;
    fp_temp.open("./class1_train.txt");
    int N1_total=Find_No_Lines(fp_temp);
    fp_temp.open("./class2_train.txt");
    int N2_total=Find_No_Lines(fp_temp);
    fp_temp.open("./class1_train.txt");
    int D=Find_No_Words(fp_temp, N1_total);

    //TRAINING SET
    int N1, N2;
    N1=N1_total;
    N2=N2_total;

    //Declaration of Training
    double ** class1=new double*[N1];
    for(int i=0; i<N1; i++)
    {
        class1[i]=new double[D+1];
    }
    double **class2=new double*[N2];
    for(int i=0; i<N2; i++)
    {
        class2[i]=new double[D+1];
    }

    //taking input from file for training
    fp_input1.open("./class1_train.txt");
    for(int i=0; i<N1; i++)
    {
        for(int j=0; j<D; j++)
        {
            fp_input1>>class1[i][j];
        }
        class1[i][D]=1;
    }


    fp_input2.open("./class2_train.txt");
    for(int i=0; i<N2; i++)
    {
        for(int j=0; j<D; j++)
        {
            fp_input2>>class2[i][j];
        }
        class2[i][D]=2;
    }

    
    fp_test.open(argv[1]);
    int N1_Test=Find_No_Lines(fp_test);
    fp_test.close();
    fp_test.open(argv[1]);

    //Declaration of Test
    double ** class1_Test=new double*[N1_Test];
    for(int i=0; i<N1_Test; i++)
    {
        class1_Test[i]=new double[D+1];
    }

    //taking input from file for test
    for(int i=0; i<N1_Test; i++)
    {
        for(int j=0; j<D; j++)
        {
            fp_test>>class1_Test[i][j];
        }
        class1_Test[i][D]=1;
    }
    fp_test.close();

    double mean1[D], mean2[D];

    //computing Mean
    computeMean(class1, mean1, D, N1);
    computeMean(class2, mean2, D, N2);

    //covariance Matrix declaration

    double ** covar1=new double*[D];
    double ** covar2=new double*[D];
    double ** Invcovar1=new double*[D];
    double ** Invcovar2=new double*[D];
    double ** tempcovar=new double*[D];
    for(int i=0; i<D; i++)
    {
        covar1[i]=new double[D];
        covar2[i]=new double[D];
        Invcovar1[i]=new double[D];
        Invcovar2[i]=new double[D];
        tempcovar[i]=new double[D];
    }

    //Computing Covariance
    ComputeCovariance(covar1, class1, mean1, D, N1);
    ComputeCovariance(covar2, class2, mean2, D, N2);

        //computing Determinants
    double Det1, Det2;
    for(int i=0; i<D; i++)
    {
        for(int j=0; j<D; j++)
        {
           tempcovar[i][j]=covar1[i][j];
        }
    }
    Det1=FindDeterminant(tempcovar, D);

    for(int i=0; i<D; i++)
    {
        for(int j=0; j<D; j++)
        {
           tempcovar[i][j]=covar2[i][j];
        }

    }
    Det2=FindDeterminant(tempcovar, D);

    //computing Inverses
    FindInverse(Invcovar1, covar1, D, Det1);
    FindInverse(Invcovar2, covar2, D, Det2);

    ofstream fp_output;
    fp_output.open("./output.txt");
    double count=0;
    
    //Finding class for test file
    count+=FindClass(Det1, Det2, Invcovar1, Invcovar2, class1_Test, mean1, mean2, D, N1_Test, fp_output);
    fp_output.close();
    return 0;
}

