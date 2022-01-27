#include <cstdlib>
#include "common.h"
/***************************************************************
 *               SIMPLE MATRIX MULTIPLICATION                  *
 ***************************************************************/
void linalg::multiply_Ax ( const double mat[3][3], double & ax, double & ay, 
		double & az, double & bx, double & by, double & bz )
{
	bx = mat[0][0]*ax + mat[0][1]*ay + mat[0][2]*az;
	by = mat[1][0]*ax + mat[1][1]*ay + mat[1][2]*az;
	bz = mat[2][0]*ax + mat[2][1]*ay + mat[2][2]*az;
}

// jacobi matrix inversion: stolen from vmd code - heh,heh!
/* jacobi.C, taken from Numerical Recipes and specialized to 3x3 case */

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

int linalg::jacobi(double a[3][3], double d[3], double v[3][3])
{
  int n=3;
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  b=new double[n];
  z=new double[n];
  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      delete [] z;
      delete [] b;
      return 0; // Exit normally
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((double)(fabs(h)+g) == (double)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=0;j<n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  delete [] b;
  delete [] z;
  return 1; // Failed to converge
}

#undef ROTATE

void linalg::normalise ( double x[3] )
{
	double modx = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
	x[0] /= modx;
	x[1] /= modx;
	x[2] /= modx;
}

void linalg::normalise ( double & x, double & y, double & z )
{
	double modx = sqrt( x*x + y*y + z*z );
	x /= modx;
	y /= modx;
	z /= modx;
}

double linalg::modx ( double x[3] )
{
	return sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
}

double linalg::dotp (double x[3], double y[3] )
{
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void linalg::crossp (double a[3], double b[3], double c[3])
{
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

double linalg::angle (double x[3], double y[3] )
{
	return acos( dotp(x,y) );
}

double linalg::det( double M[3][3] )
{
	double xbit = M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]);
	double ybit = M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]);
	double zbit = M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
	return xbit - ybit + zbit;
}

void linalg::zero_mat( double M[3][3] ) 
{
	M[0][0] = M[0][1] = M[0][2] = 0.0;
	M[1][0] = M[1][1] = M[1][2] = 0.0;
	M[2][0] = M[2][1] = M[2][2] = 0.0;
}

void linalg::ident_mat( double M[3][3] ) 
{
	M[0][0] = M[0][1] = M[0][2] = 0.0;
	M[1][0] = M[1][1] = M[1][2] = 0.0;
	M[2][0] = M[2][1] = M[2][2] = 0.0;
}
/***************************************************************
 *                   EXIT WITH ERROR MESSAGE                   *
 ***************************************************************/
void die_now( const char * message, int ret_val )
{
	cerr << message << endl;
	exit(ret_val);
}

void print_err( const string & function_name, 
		const string & message)
{
	print_err( function_name.c_str(), message.c_str() );
}

void print_err( const char * function_name, 
		const char * message)
{
	const int wrap_width = 60;
	cerr << "\nA Message from " << function_name << endl;
	const char * i = &message[0];
	int cwritten = 0;
	while ( *i != '\0' )
	{
		if ( *i == '\t' || *i == '\n' || *i == ' ' )
		{
			i++;
			while ( *i == ' ' || *i == '\t' || *i == '\n' ) i++;
			if (cwritten > wrap_width)
			{
				cerr << '\n';
				cwritten = 0;
			}
			else
			{
				cerr << " ";
				cwritten++;
			}
		}
		else
		{
			cerr << *i++;
			cwritten++;
		}
	}
	cerr << endl;
}
