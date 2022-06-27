#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
// int n1,mx1,my1; //number of latttice nodes
char setupfile[]="setupfile.txt";
double u0;
int mstep; // The total number of time steps
double alpha;
// const int n=n1,mx=mx1,my=my1; //number of latttice nodes
int n,mx,my; //number of latttice nodes
double ResU=1.0,ResV=1.0,Resrho=1.0;//residuals of different variable
int i,j;
int dx=1,dy=1; //space and time step
double const rho0=5.00;
int static const norm=2;//different norm selection
double Re;
// double Re=u0*my/alpha;
double omega;
// double omega=1.0/(3.*alpha+0.5);
int t=0;
//double *cx,*cy,*w,*x,*y,*feq,**rho,**rhoo,**uo,**vo,**u,**v,***f;
void Readinitialfile(char* filename,int *n,int *mx,int *my,double *u0,double *alpha,int *mstep)
{ 
    FILE *file=NULL;
    file=fopen(filename,"r");
    if(!file)
    {
        std::cout<<"cannot find file"<<std::endl;
    }
    fscanf(file,"%d",n);
    fscanf(file,"%d",mx);
    fscanf(file,"%d",my);
    fscanf(file,"%lf",u0);
    fscanf(file,"%lf",alpha);
    fscanf(file,"%d",mstep);
    fclose(file);
}

void initialize(double *cx,double *cy,double *w,double *x,double *y,double **rho,double **u,double **v,double ***f)
{
Re=u0*my/alpha;
omega=1.0/(3.*alpha+0.5);
// std::cout << "rho, before" << std::endl;
// std::cout << rho[1][1] << std::endl;
// std::cout << "rho, after" << std::endl;

std::cout<<"Re= "<<Re<<std::endl;
std::cout<<"omega= "<<omega<<std::endl;
std::string filename = std::string("animation0") +std::string(".dat");
/*---------streaming vector--------*/
cx[0]=0,cx[1]=1,cx[2]=0,cx[3]=-1,cx[4]=0,cx[5]=1,cx[6]=-1,cx[7]=-1,cx[8]=1;
cy[0]=0,cy[1]=0,cy[2]=1,cy[3]=0,cy[4]=-1,cy[5]=1,cy[6]=1,cy[7]=-1,cy[8]=-1;
/*---------streaming vector--------*/
/*-------weight factor ------*/
w[0]=4./9;
    for(i=1;i<5;i++)
    {
        w[i]=1./9;
    }
    for(i=5;i<9;i++)
    {
        w[i]=1./36;
    }
/*---------weight factor ----------*/
//coordinate
x[0] =0.0; 
y[0] =0.0;
for (i=1;i<mx;i++)
{
    x[i]=x[i-1]+dx;
}
for (j=1;j<my;j++)
{
    y[j]=y[j-1]+dy;
}
/*initial condition--------------*/

    for(j=0;j<my;j++)
    {
        for (i=0;i<mx;i++)
        {
            rho[i][j]=rho0; // Initial value
            u[i][j]=0.0;
            v[i][j]=0.0;
        }
    }

    for (i=1;i<mx;i++)//top moving wall.
    {
        u[i][my-1]=u0;
        v[i][my-1]=0.0;
    }
//initial condition of distribution function 
    double strf[mx][my];
    for(int j=0;j<my;j++)
    {
       for(int i=0;i<mx;i++)
       {
           for(int k=0;k<9;k++)
           {
                f[k][i][j]=w[k]*rho[i][j];
                strf[i][j]=0;
           }
       }
    }
/*initial condition--------------*/
std::cout<<"n "<<n<<std::endl;
    std::cout<<"mx "<<mx<<std::endl;
    std::cout<<"my "<<my<<std::endl;
    std::cout<<"u0 "<<u0<<std::endl;
    std::cout<<"alpha "<<alpha<<std::endl;
    std::cout<<"mstep "<<mstep<<std::endl;
}
void collesion(double **u,double **v,double ***f,double *feq,double **rho,double omega,
double *w,double *cx,double *cy)
{
    double t1,t2;
    for (int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            t1=u[i][j]*u[i][j]+v[i][j]*v[i][j];
                for(int k=0;k<9;k++)
                {
                    t2=u[i][j]*cx[k]+v[i][j]*cy[k];
                    feq[k]=rho[i][j]*w[k]*(1.0+3.0*t2+4.50*t2*t2-1.50*t1);//ck=1=dt/dx
                    f[k][i][j]=omega*feq[k]+(1.-omega)*f[k][i][j];
                }
        }
    }
}

void streaming(double ***f)
{
    for (int j=0;j<my;j++)
    {
        for(int i=mx-1;i>0;i--)
        {
            f[1][i][j]=f[1][i-1][j];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[3][i][j]=f[3][i+1][j];
        }
    }

    for(int j=my-1;j>0;j--)
    {
        for(int i=0;i<mx;i++)
        {
            f[2][i][j]=f[2][i][j-1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f[5][i][j]=f[5][i-1][j-1];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[6][i][j]=f[6][i+1][j-1];
        }
    }

    for(int j=0;j<my;j++)
    {
        for(int i=0;i<mx;i++)
        {
            f[4][i][j]=f[4][i][j+1];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[7][i][j]=f[7][i+1][j+1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f[8][i][j]=f[8][i-1][j+1];
        }
    }
}

void bound(double ***f,double u0)
{
    double rhon;
    for(int j=0;j<my;j++)
    {
        f[1][0][j]=f[3][0][j];//left
        f[5][0][j]=f[7][0][j];
        f[8][0][j]=f[6][0][j];

        f[3][mx-1][j]=f[1][mx-1][j];//right
        f[7][mx-1][j]=f[5][mx-1][j];
        f[6][mx-1][j]=f[8][mx-1][j];
    }

    for (int i=0;i<mx;i++)
    {
        f[2][i][0]=f[4][i][0];//bottom
        f[5][i][0]=f[7][i][0];
        f[6][i][0]=f[8][i][0];
    }

    // TODO: check bounds of loop (should be mx-1?)
    for (int i=1;i<mx-1;i++)//top wall 
    {
        rhon=(f[0][i][my-1]+f[1][i][my-1]+f[3][i][my-1]+2.*(f[2][i][my-1]+f[6][i][my-1]+f[5][i][my-1]))/(1+u0);
        //zou-he
        f[4][i][my-1]=f[2][i][my-1]-2/3*rhon*u0;
        f[8][i][my-1]=f[6][i][my-1]+rhon*u0/6.0+1/2*(f[1][i][my-1]-f[3][i][my-1])-1/2*rhon*u0;//zou-he
        f[7][i][my-1]=f[5][i][my-1]-rhon*u0/6.0+1/2*(f[3][i][my-1]-f[1][i][my-1])+1/2*rhon*u0;//zou-he
    }

    //corner
    // for(int k=0;k<9;k++)
    // {
    //     f[k][0][0]=0;
    //     f[k][mx-1][0]=0;
    //     f[k][0][my-1]=0;
    //     f[k][mx-1][my-1]=0; 
    // }
    
   
}

void rhouv(double ***f,double **rho,double **u,double **v,double *cx,double *cy)
{
    double ssum,usum,vsum;
    for(int j=0;j<my;j++)
    {
       for(int i=0;i<mx;i++)
       {
           ssum=0;
           for(int k=0;k<9;k++)
           {
               ssum=ssum+f[k][i][j];
           }
           rho[i][j]=ssum;
       }
    }

    for(int i=1;i<mx;i++)//top density 
    {
        rho[i][my-1]=f[0][i][my-1]+f[1][i][my-1]+
        f[3][i][my-1]+2.*(f[2][i][my-1]+f[6][i][my-1]+f[5][i][my-1]);
    }

    for(int i=1;i<mx;i++)
    {
        for(int j=1;j<my;j++)
        {
            usum=0.0;
            vsum=0.0;
            for(int k=0;k<9;k++)
            {
                usum=usum+f[k][i][j]*cx[k];
                vsum=vsum+f[k][i][j]*cy[k];
            }
            u[i][j]=usum/rho[i][j];
            v[i][j]=vsum/rho[i][j];
        }
    }

    // TODO: override macroscopic boundary conditions here
    for(int j=0;j<my;j++)
    {
        u[0][j]  = 0.0;
        v[0][j]  = 0.0;
        u[mx-1][j] = 0.0;
        v[mx-1][j] = 0.0; 
    }

    for(int i=0;i<mx;i++)
    {
        u[i][0]  = 0.0;
        v[i][0]  = 0.0;
        u[i][my-1] = 0.1;
        v[i][my-1] = 0.0; 
    }

}

// TODO: implement equillibrium function writing
void result(double **u,double **v,double **rho,double *x, double *y,std::string filename,double ***f,double *feq,int time)
{
double rhoav,rhom;

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"U\",\"V\",\"density\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\",\"feq0\",\"feq1\",\"feq2\",\"feq3\",\"feq4\",\"feq5\",\"feq6\",\"feq7\",\"feq8\"," << std::endl;
  out << "ZONE T = \"fluid\", I=" << mx << ", J=" << my << ", F=POINT" << std::endl;
  out << "SOLUTIONTIME="<< time << std::endl;
  for (unsigned j = 0; j < my; ++j)
    for (unsigned i = 0; i < mx; ++i)
    {
      out << std::scientific << std::setprecision(5) << std::setw(15) << x[i];
      out << std::scientific << std::setprecision(5) << std::setw(15) << y[j];
      out << std::scientific << std::setprecision(5) << std::setw(15) << u[i][j];
      out << std::scientific << std::setprecision(5) << std::setw(15) << v[i][j];
      out << std::scientific << std::setprecision(5) << std::setw(15) << rho[i][j];
      //out << std::scientific << std::setprecision(5) << std::setw(15) << strf[i][j];
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << f[k][i][j];
      }
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << feq[k];
      }
      out << std::endl;
    }
  out.close();

}
//residual calculation function 
void calculationResiduals(int norm, double &ResU,double &ResV,double &Resrho, double **rhoo,double **uo,double **vo,double **u,double **v,double **rho)
{
    
    // using an appropriate vector norm (L0, L1, L2)
switch (norm) 
    {
        case 0:
        {
            // L0: Look for largest absolute difference
            double maxU,maxV,maxrho = -999999.9;
            double tempResidualU,tempResidualV,tempResidualrho=0;
            for (int i=0;i<mx;i++)
            {
                for (int j=0;j<my;j++)
                {
                    // repeat below for each variable
                    tempResidualU = std::fabs(u[i][j] - uo[i][j]);
                    if (tempResidualU > maxU)
                    maxU = tempResidualU;
                    tempResidualV = std::fabs(v[i][j] - vo[i][j]);
                    if (tempResidualV > maxV)
                    maxV = tempResidualV;
                    tempResidualrho = std::fabs(rho[i][j] - rhoo[i][j]);
                    if (tempResidualrho > maxrho)
                    maxrho = tempResidualrho;
                } 
            }
            ResU=maxU;
            ResV=maxV;
            Resrho=maxrho;
            break;
        }
        case 1:
        {
            // L1: look for average difference
            for (int i=0;i<mx;i++)
            {
                for (int j=0;j<my;j++)
                {
                    ResU = ResU + u[i][j];
                    ResV = ResV + v[i][j];
                    Resrho = Resrho + rho[i][j];                                       
                }
            }
            ResU = ResU / (mx * my);
            ResV = ResV / (mx * my);
            Resrho = Resrho / (mx * my);
            break;
        }
        case 2:
        {
            // L2: look for root-mean square difference
            for (int i=0;i<mx;i++)
            {
                for (int j=0;j<my;j++)
                {
                    //ResU = ResU + std::sqrt(std::pow(u[i][j], 2)) // this line is equivalent to the below
                    ResU = ResU + std::fabs(u[i][j]);
                    ResV = ResV + std::fabs(u[i][j]);
                    Resrho = Resrho + std::fabs(u[i][j]);
                }
            }
            ResU = ResU / (mx * my);
            ResV = ResV / (mx * my);
            Resrho = Resrho / (mx * my);
            break;
        }
        
    }
    
}

bool IfcalculationEnd(double ResU,double ResV,double Resrho)
{
    double epsilon = 1e-5;
    bool calculationFinished = false;
    if (ResU < epsilon && ResV < epsilon && Resrho < epsilon)
        calculationFinished = true;
    return calculationFinished;    
}

// void Readinitialfile(char* filename,int &n,int &mx,int &my,double &u0,double &alpha,int &mstep)
// {
//     FILE *file=NULL;
//     file=fopen(filename,"r");
//     if(!file)
//     {
//         std::cout<<"cannot find file"<<std::endl;
//     }
//     fscanf(file,"%d","%d","%d","%f","%f","%d",&n1,&mx1,&my1,&u0,&alpha,&mstep);
//     fclose(file);
// }

int main()
{
Readinitialfile(setupfile,&n,&mx,&my,&u0,&alpha,&mstep);
double *cx = new double[n]; 
double *cy = new double[n];
double *w = new double[n]; 
double *x = new double[mx]; 
double *y = new double[my];
double *feq = new double[n];

double **rho = new double *[mx];
double **rhoo = new double *[mx];
double **uo = new double *[mx];
double **vo = new double *[mx];
double **u = new double *[mx];
double **v = new double *[mx];
for(int s=0;s<mx;s++)
{
    rho[s] = new double [my];
    rhoo[s] = new double [my];
    uo[s] = new double [my];
    vo[s] = new double [my];
    u[s] = new double [my];
    v[s] = new double [my];
}
double ***f = new double**[n];
for(int s=0;s<n;s++)
{
    f[s]= new double *[mx];
    for(int ss=0;ss<mx;ss++)
    {
        f[s][ss]=new double [my];
    }
}
initialize(cx, cy, w, x, y, rho, u,v, f);
//main loop
int time;
for (time=1;time<mstep+1;++time)
{
//---copy data from previous time step uo,vo, rhoo----//
    for (i=0;i<mx;i++)
    {
        for (j=0;j<my;j++)
        {
            uo[i][j]   = u[i][j];  
            vo[i][j]   = v[i][j];  
            rhoo[i][j] = rho[i][j];
        }
    }
//----copy data from previous time step uo,vo, rhoo---//
collesion(u,v,f,feq,rho,omega,w,cx,cy);
streaming(f);
bound(f,u0);
rhouv(f,rho,u,v,cx,cy);
calculationResiduals(norm,ResU,ResV,Resrho,rhoo,uo,vo,u,v,rho);
//---------judge whether the residual is below the threshold--//
bool calculationFinsihed = IfcalculationEnd(ResU,ResV,Resrho);
if (calculationFinsihed)
    break;
//---------judge whether the residual is below the threshold--//
std::string filename = std::string("animation") + std::to_string(time) + std::string(".dat");
result(u,v, rho, x, y,filename, f, feq,time);

//print residuals
std::cout << std::fixed << std::setw(8) << time;
std::cout << std::scientific << std::setprecision(5) << std::setw(15) << ResU;
std::cout << std::scientific << std::setprecision(5) << std::setw(15) << ResV;
std::cout << std::scientific << std::setprecision(5) << std::setw(15) << Resrho;
std::cout << std::endl;
}//main loop end
//result(u,v,rho,x,y);
return 0;
}



        
        