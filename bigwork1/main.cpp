#include <iostream>
#include <fstream>
#include "math.h"
using namespace std;

const double pi =3.141592;
const double e =1e-7 ;   /*收敛误差限*/

const int m =101; //在圆弧上取101个点
const int n_in=601; //内部比较密的点 6m
const int n_out=100 ;//外部比较稀疏的点 10m

const double de=0.01;
const double dn=de;

//翼型 NACA0012
const double t=0.122 ;
const double c=1;

const int n =701;
const double q=1.94;




int main()
{
	int i,j,s=0,s1=0;
   	double x[m][n],y[m][n],dx[m][n],dy[m][n],a[m][n],b[m][n],maxa,maxb;
	double alph[m][n],belta[m][n],gamma[m][n],bw[m][n],be[m][n],bs[m][n],bn[m][n],bp[m][n],cpx[m][n],cpy[m][n];
    /*-----------网格生成-----------*/

	cout<<"正在计算中******"<<endl;





   	for (i=0;i<m;i++)/*-----------翼型坐标-----------*/
	{
		x[i][0]=0.499*(1+cos(2.0*i*pi/(m-1))); //用正弦分布给出x点
		if(i<0.5*(m-1))
		{
				y[i][0]=5*t*c*(-0.1015*pow(x[i][0]/c,4)+0.2843*pow(x[i][0]/c,3)-0.3516*pow(x[i][0]/c,2)-0.1221*x[i][0]/c+0.2969*sqrt(x[i][0]/c));
		}
		else
		{  y[i][0]=-5*t*c*(-0.1015*pow(x[i][0]/c,4)+0.2843*pow(x[i][0]/c,3)-0.3516*pow(x[i][0]/c,2)-0.1221*x[i][0]/c+0.2969*sqrt(x[i][0]/c));
		}

	}

cout<<"ok"<<endl;

for(i=0;i<m;i++)/*-----------远场坐标-----------*/
	{
		  x[i][n-1]=16.5*cos(2.0*i*pi/(m-1))+0.5;
		  //f[i][n-1]=Vecioty*x[i][n-1];
		 if(i<0.5*(m-1)) y[i][n-1]=sqrt(16.5*16.5- (x[i][n-1]-0.5)*(x[i][n-1]-0.5) );
			    else
					y[i][n-1]=-sqrt(16.5*16.5- (x[i][n-1]-0.5)*(x[i][n-1]-0.5) );

	}

	for(i=0;i<m;i++)/*-----------分划坐标-----------*/
	{
		x[i][n_in-1]=6.5*cos(2.0*i*pi/(m-1))+0.5;
		//f[i][n-1]=Vecioty*x[i][n-1];
	 if(i<0.5*(m-1))  y[i][n_in-1]=sqrt(6.5*6.5- (x[i][n_in-1]-0.5)*(x[i][n_in-1]-0.5) );
				else
				y[i][n_in-1]=-sqrt(6.5*6.5- (x[i][n_in-1]-0.5)*(x[i][n_in-1]-0.5) );

	}
	//内场分为两个部分
	for(i=0;i<m;i++)

	{
		for(j=1;j<n_in;j++)
		{
			x[i][j]=x[i][0]+j*(x[i][n_in]-x[i][0])/(n_in);
	        y[i][j]=y[i][0]+j*(y[i][n_in]-y[i][0])/(n_in);

		}
	}

	for(i=0;i<m;i++)

	{
		for(j=n_in+1;j<n;j++)
		{

		x[i][j]=x[i][n_in]+j*(x[i][n-1]-x[i][n_in])/(n-n_in);
		y[i][j]=y[i][n_in]+j*(y[i][n-1]-y[i][n_in])/(n-n_in);

	}
}

	cout<<"正在生成网格数据中******"<<endl;
	cout<<"......"<<endl;

do        /*-----------内场点坐标计算-----------*/

    {


		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
			a[i][j]=x[i][j];
			b[i][j]=y[i][j];

			alph[i][j]=((x[i][j+1]-x[i][j-1])*(x[i][j+1]-x[i][j-1])+(y[i][j+1]-y[i][j-1])*(y[i][j+1]-y[i][j-1]))/(4*dn*dn);
			belta[i][j]=((x[i+1][j]-x[i-1][j])*(x[i][j+1]-x[i][j-1])+(y[i+1][j]-y[i-1][j])*(y[i][j+1]-y[i][j-1]))/(4*de*dn);
			gamma[i][j]=((x[i+1][j]-x[i-1][j])*(x[i+1][j]-x[i-1][j])+(y[i+1][j]-y[i-1][j])*(y[i+1][j]-y[i-1][j]))/(4*de*de);

			bw[i][j]=be[i][j]=alph[i][j]/(de*dn);
			bs[i][j]=bn[i][j]=gamma[i][j]/(de*dn);
			bp[i][j]=bw[i][j]+bs[i][j]+be[i][j]+bn[i][j];

			cpx[i][j]=-belta[i][j]*(x[i+1][j+1]-x[i+1][j-1]-x[i-1][j+1]+x[i-1][j-1])/(2*de*dn);
			cpy[i][j]=-belta[i][j]*(y[i+1][j+1]-y[i+1][j-1]-y[i-1][j+1]+y[i-1][j-1])/(2*de*dn);

			}
		}


		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				x[i][j]=(bw[i][j]*a[i-1][j]+be[i][j]*a[i+1][j]+bs[i][j]*a[i][j-1]+bn[i][j]*a[i][j+1]+cpx[i][j])/bp[i][j];
				y[i][j]=(bw[i][j]*b[i-1][j]+be[i][j]*b[i+1][j]+bs[i][j]*b[i][j-1]+bn[i][j]*b[i][j+1]+cpy[i][j])/bp[i][j];




			dx[i][j]=fabs(a[i][j]-x[i][j]);
	        dy[i][j]=fabs(b[i][j]-y[i][j]);

			}
		}

			 maxa=dx[1][1],maxb=dy[1][1];

		for(i=1;i<m-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				if(dx[i][j]>=maxa)
				{
					maxa=dx[i][j];
				}
			    if(dy[i][j]>=maxb)
				{
					maxb=dy[i][j];
				}
			}
		}
		s+=1;
	}while(maxa>e||maxb>e); /*-----------判断收敛-----------*/

	cout<<"网格收敛迭代次数: "<<s<<endl;




ofstream outfile;

	outfile.open("网格.dat");


        for(i=0;i<m;i++)
		{
			for(j=0;j<n;j++)
			{
			  outfile<<x[i][j]<<"  "<<y[i][j]<<endl;
            }
		}

       outfile.close();

       cout<<"计算结束!"<<endl;
	   return 0;
	}
