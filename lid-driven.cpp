#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
int const static n=9,mx=40,my=40; //number of latttice nodes
double f[n][mx][my],feq[9],rho[mx][my],cx[n],cy[n],w[n],u[mx][my],v[mx][my],x[mx],y[my];
double rhoo[mx][my],uo[mx][my],vo[mx][my];//previous variable for resdiual calculation 
int main()
{
double ResU=1.0,ResV=1.0,Resrho=1.0;//residuals of different variable
int i,j;
int dx=1,dy=1; //space and time step
x[0] =0.0; 
y[0] =0.0;
double u0=0.1; 
double const rho0=5.00;
int static const norm=2;//different norm selection
//coordinate
for (i=1;i<mx;i++)
{
    x[i]=x[i-1]+dx;
}
for (j=1;j<my;j++)
{
    y[j]=y[j-1]+dy;
}
double const alpha=0.02;
double Re=u0*my/alpha;
std::cout<<"Re= "<<Re<<std::endl;
double omega=1.0/(3.*alpha+0.5);
std::cout<<"omega= "<<omega<<std::endl;
int mstep=4000; // The total number of time steps
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
//function declaration 
void collesion(double u[mx][my],double v[mx][my],double f[9][mx][my],double feq[9],double rho[mx][my],double omega,
double w[9],double cx[9],double cy[9]);
void streaming(double f[9][mx][my]);
void bound(double f[9][mx][my],double u0);
void rhouv(double f[9][mx][my],double rho[mx][my],double u[mx][my],double v[mx][my],double cx[9],double cy[9]);
void result(double u[mx][my],double v[mx][my],double rho[mx][my],double x[mx], double y[my],std::string filename,double f[9][mx][my],double feq[9],int time);
void calculationResiduals(int norm, double &ResU,double &ResV,double &Resrho, double rhoo[mx][my],double uo[mx][my],double vo[mx][my],double u[mx][my],double v[mx][my],double rho[mx][my]);
bool IfcalculationEnd(double ResU,double ResV,double Resrho);

/*---------streaming vector--------*/
cx[0]=0,cx[1]=1,cx[2]=0,cx[3]=-1,cx[4]=0,cx[5]=1,cx[6]=-1,cx[7]=-1,cx[8]=1;
cy[0]=0,cy[1]=0,cy[2]=1,cy[3]=0,cy[4]=-1,cy[5]=1,cy[6]=1,cy[7]=-1,cy[8]=-1;
/*---------streaming vector--------*/

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
std::string filename = std::string("animation0") +std::string(".dat");
int t=0;
result(u,v,rho,x,y,filename,f,feq,t);
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
result(u,v,rho,x,y,filename,f,feq,time);

// print residuals
std::cout << std::fixed << std::setw(8) << time;
std::cout << std::scientific << std::setprecision(5) << std::setw(15) << ResU;
std::cout << std::scientific << std::setprecision(5) << std::setw(15) << ResV;
std::cout << std::scientific << std::setprecision(5) << std::setw(15) << Resrho;
std::cout << std::endl;
}//main loop end
//result(u,v,rho,x,y);
return 0;
}

void collesion(double u[mx][my],double v[mx][my],double f[9][mx][my],double feq[9],double rho[mx][my],double omega,
double w[9],double cx[9],double cy[9])
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

void streaming(double f[9][mx][my])
{
    for (int j=0;j<my;j++)
    {
        for(int i=mx-1;i>0;i--)
        {
            f[1][i][j]=f[1][i-1][j];
        }
        for(int i=0;i<mx;i++)
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
        for(int i=0;i<mx;i++)
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
        for(int i=0;i<mx;i++)
        {
            f[7][i][j]=f[7][i+1][j+1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f[8][i][j]=f[8][i-1][j+1];
        }
    }
}

void bound(double f[9][mx][my],double u0)
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

void 
rhouv(double f[9][mx][my],double rho[mx][my],double u[mx][my],double v[mx][my],double cx[9],double cy[9])
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
void result(/*double f[n][mx][my], */double u[mx][my],double v[mx][my],double rho[mx][my],double x[mx], double y[my],std::string filename,double f[9][mx][my],double feq[9],int time)
{
double strf[mx][my],rhoav,rhom;
/*-----stream function calculation--------*/
strf[0][0]=0.0;
for(int i=0;i<mx-1;i++)
{
    rhoav=0.5*(rho[i-1][0]+rho[i][0]);
    if(i!=0) strf[i][0]=strf[i-1][0]-rhoav*0.5*(v[i-1][0]+v[i][0]);

    for(int j=1;j<my-1;j++)
    {
        rhom=0.5*(rho[i][j]+rho[i][j-1]);
        strf[i][j]=strf[i][j-1]+rhom*0.5*(u[i][j-1]+u[i][j]);
    }
}
/*-----stream function calculation--------*/

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
void calculationResiduals(int norm, double &ResU,double &ResV,double &Resrho, double rhoo[mx][my],double uo[mx][my],double vo[mx][my],double u[mx][my],double v[mx][my],double rho[mx][my])
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
        
        