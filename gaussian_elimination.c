#include<stdio.h>
#include<stdlib.h>
#include<time.h>

void computeGaussianElimination(int n,float arr[][n]);
int determinant(int n,float matrix[][n]);
void getCoFactorMatrix(int n,float matrix[][n],int index,float cofactor[][n-1]);

int main()
{
	int n;
	printf("Enter the number of unknowns : ");
	scanf("%d",&n);
	printf("Enter the augmented matrix : \n");
	float arr[n][n];
	for(int i=0 ; i<n ; i++)
	{
		for(int j=0 ; j<n+1 ; j++)
		{
			printf("Enter arr[%d][%d] : ",i,j);
			scanf("%f",&arr[i][j]);
		}
	}
	
	
	
	int det = determinant(n,arr);
	printf("%d\n",det);
	if(det==0)
	{
		printf("\nThe matrix is singular!\n");
	}
	
	else
	{
		computeGaussianElimination(n,arr);
	}
	
	
	
	
	return 0;
}



int determinant(int n,float matrix[][n])
{
	int D=0 ;
	
	if(n==1)
	{
		return matrix[0][0];
	}
	
	int sign=1 ;
	float cofactor[n][n];
	
	for(int i=0 ; i<n ; i++)
	{
		getCoFactorMatrix(n,matrix,i,cofactor);
		D+=sign*matrix[0][i]*determinant(n-1,cofactor);
		
		sign*=-1;
	}
	
	return D;
	
}

void getCoFactorMatrix(int n,float matrix[][n],int index,float cofactor[][n-1])
{
	int i=0 ;
	int j=0 ;
	for(int row=0 ; row<n ; row++)
	{
		for(int col=0 ; col<n ; col++)
		{
			if(row!=0 && col!=index)
			{
				cofactor[i][j++]=matrix[row][col];
				
				if(j==n-1)
				{
					++i;
					j=0;
				}
				
			}
		}
	}	
}

void computeGaussianElimination(int n,float arr[][n])
{
	clock_t start ;
	clock_t end ;
	double cpu_time_used=0 ;
	float solution[n];
	float lMatrix[n][n];
	for(int i=0 ; i<n ; i++)
	{
		for(int j=0 ; j<n ; j++)
		{
			lMatrix[i][j]=0;
		}
		lMatrix[i][i]=1;
	}
	//Gaussian elimination
	//Generate Upper Triangular Matrix
	start = clock() ;
	for(int i=0 ; i<n ; i++)
	{
		float pivot = arr[i][i] ;
		for(int j=i+1 ; j<n ; j++)
		{
			float multiplier = arr[j][i]/pivot ;
			lMatrix[j][i]=multiplier;

			for(int k=0 ; k<=n ; k++)
			{
				arr[j][k]= arr[j][k] - multiplier*arr[i][k];
			}
		}
	}
	solution[n-1] = arr[n-1][n] / arr[n-1][n-1] ; 

	

	//Back substitution
	float sum;
	for(int i=n-2 ; i>=0 ; --i)
	{
		sum=0;
		for(int j=n-1 ; j>i ; j--)
		{
			sum+= arr[i][j]*solution[j] ;
		}
		solution[i] = (arr[i][n]-sum)/arr[i][i] ;
	}
	
	
	
	end=clock() ;
	cpu_time_used = (double) (end-start)/CLOCKS_PER_SEC ;

	printf("\nThe L matrix : \n") ;	
	for(int i=0 ; i<n ; i++)
	{
		for(int j=0 ; j<n ; j++)
		{
			printf("%f\t",lMatrix[i][j]);
		}
		printf("\n");
	}
	
	printf("\nThe U matrix : \n ");
	for(int i=0 ; i<n ; i++)
	{
		for(int j=0 ; j<n ; j++)
		{
			printf("%f\t",arr[i][j]);
		}
		printf("\n");
	}
	
	
	for(int i=0 ; i<n ; i++)
	{
		printf("X[%d] : %f\n",i,solution[i]);
	}
	cpu_time_used = (double) (end-start)/CLOCKS_PER_SEC ;
	
	printf("The time taken for the computation is : %f \n",cpu_time_used);
}
