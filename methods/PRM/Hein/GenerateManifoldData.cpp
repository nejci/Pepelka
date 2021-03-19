#include "mex.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

/* uniform [0,1] random number generator
   developed by Pierre Lecuyer based on a clever
   and tested combination of two linear congruential
   sequences 
*/

/*
s1 and s2 are the seeds (nonnegative integers)
*/

double uni()
{
  static long s1 = 55555;
  static long s2 = 99999;
  static double factor = 1.0/2147483563.0;
  register long k,z;
  k= s1 /53668;
  s1 =40014*(s1%53668)-k*12211;
  if (s1 < 0) s1 += 2147483563;
  k=s2/52774;
  s2=40692*(s2%52774)-k*3791;
  if (s2 < 0) s2 += 2147483399;

  /*
  z = abs(s1 ^ s2);
  */
  z= (s1 - 2147483563) + s2;
  if (z < 1) z += 2147483562;

  return(((double)(z))*factor);
}


void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  //this part generates toy training data
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // choice of dataset
  // 0: sinusoid over the circle with high curvature
  // 1: a dim-1 dimensional sphere (uniformly sampled)
  // 2: a 3d-affine space in 5d
  // 3: a strange 4-dim figure in 5d (but very concentrated so it is ess 3d)
  // 4: a 4-dim manifold in 8 dimensions
  // 5: a 2d helix in 3d
  // 6: a 6-dim manifold in 36 dimensions
  // 7: "swiss roll"
  // 8: 12-dim manifold in 72 dimensions
  // 9: A 20-dim affine space in 100-dimensions
  //10: k-hypercube, uniformly sampled
  //11: Moebius band
  //12: Multivariate Gaussian
  //13: one-dimensional curve in dim dimensions

  // Test number of parameters.
  if (nrhs != 3 || nlhs != 1) {
    mexWarnMsgTxt("Usage: X = generate_Manifold_Data(Dataset,Dimension,Number_of_Points)");
    return;
  }

  int data_set= (int) (mxGetScalar(prhs[0]));
  int dim     = (int) (mxGetScalar(prhs[1]));
  int num     = (int) (mxGetScalar(prhs[2]));

  if(data_set < 0 || data_set>13) { printf("Dataset has to be between 0 and 12 !"); return; }
  if(dim <= 0) { printf("Dimension has to be larger than 0 !"); return; }
  if(num <= 0) { printf("The number of points has to be larger than 0 !"); return; }

  // output variable
  double* data;


  long i,j,k;

  double* radii;
  double* phi;
  double* theta;

  double norm;

  long* labeled;

  long nr_par;
  double** vectors;
  double** para;
  double* para1; double* para2; double* para3;
  double** para4;

  switch(data_set)
  {
    case 0:  //generates a (slightly distorted) sinusoid over the circle in R^3
      dim=3; // correct dimension is 1

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      // generate random angles and radii for the two moons
      radii=new double[num];
      phi  =new double[num];

      for(i=0;i<num;i++)
      {
        radii[i]=1;//+0.05*uni();
        phi[i]  =uni()*2*M_PI;  // random numbers in [0,2pi]
      }

      // the circle
      for(i=0;i<num;i++)
      {
        data[i*dim+0]=radii[i]*cos(phi[i]);
        data[i*dim+1]=radii[i]*sin(phi[i]);
        data[i*dim+2]=0.1*sin(150*phi[i]);
        //data[i][2]=0.2*uni()-0.1;
      }
      printf("Generated %i data points of a sinusoid on the circle in R^3\n",num);
      printf("The correct dimension of this submanifold is 1\n");

      delete radii;
      delete phi;
      break;

     case 1:  // generates a k-sphere  (uniformly sampled)
      nr_par=dim;
      if(nr_par %2 ==1) nr_par+=1;

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      para=new double*[nr_par];
      for(i=0;i<nr_par;i++)
       para[i]=new double[num];

      for(j=0;j<nr_par;j++)
       for(i=0;i<num;i++)
       {
         para[j][i]=uni();
         while(para[j][i]==0)
          para[j][i]=uni();
       }


      for(i=0;i<num;i++)
       { for(j=0;j<dim;j++) data[i*dim+j]=0; }
      k=dim;
      for(i=0;i<num;i++)
      {
        for(j=0;j<dim-2;j+=2)
        {
          data[i*dim+j  ]=sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i*dim+j+1]=sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        if(dim % 2==0)
        {
          data[i*dim+j  ]=sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i*dim+j+1]=sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        else
          data[i*dim+j  ]=sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);

        // now normalize these Gaussian variables
        norm=0;
        for(j=0;j<dim;j++)
         norm+=data[i*dim+j]*data[i*dim+j];
        norm=sqrt(norm);
        for(j=0;j<dim;j++)
         data[i*dim+j]/=norm;
      }

      printf("Generated %i data points of a dim-1 - dimensional sphere in R^dim\n",num);
      printf("The correct dimension of this submanifold is %i\n",dim-1);

      for(i=0;i<dim-1;i++)
       delete para[i];
      delete para;
      break;

    case 2: // generatesa 3d affine space in 5 dimensions
      dim=5;

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      para1=new double[num];
      para2=new double[num];
      para3=new double[num];

      for(i=0;i<num;i++)
      {
        para1[i]=4*uni();
        para2[i]=4*uni();
        para3[i]=4*uni();
      }

      // the affine space
      for(i=0;i<num;i++)
      {
        data[i*dim+0]= 1.2*para1[i]- 0.5*para2[i]+3;
        data[i*dim+1]= 0.5*para1[i]+ 0.9*para3[i]-1;
        data[i*dim+2]=-0.5*para1[i]- 0.2*para2[i] +   para3[i];
        data[i*dim+3]= 0.4*para1[i]- 0.9*para2[i] -  0.1*para3[i];
        data[i*dim+4]= 1.1*para1[i]- 0.3*para3[i]+8;
      }

      printf("Generated %i data points of a 3-dimensional affine subspace in R^5\n",num);
      printf("The correct dimension of this submanifold is 3\n");

      delete para1,para2,para3;
      break;

    case 3: // generates a strange 4-dim figure in 6 dimensions
      dim=6;

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      para=new double*[4];
      for(i=0;i<4;i++)
        para[i]=new double[num];

      for(i=0;i<4;i++)
       for(j=0;j<num;j++)
        para[i][j]=uni();

      //generate the figure
      for(i=0;i<num;i++)
      {
        data[i*dim+0]= para[1][i]*para[1][i]*cos(2*M_PI*para[0][i]);
        data[i*dim+1]= para[2][i]*para[2][i]*sin(2*M_PI*para[0][i]);
        data[i*dim+2]= para[1][i]+para[2][i]+pow(para[1][i]-para[3][i],2);
        data[i*dim+3]= para[1][i]-2*para[2][i]+pow(para[0][i]-para[3][i],2);
        data[i*dim+4]=-para[1][i]-2*para[2][i]+pow(para[2][i]-para[3][i],2);
        data[i*dim+5]= para[0][i]*para[0][i]-para[1][i]*para[1][i]
                   +para[2][i]*para[2][i]-para[3][i]*para[3][i];
      }

      for(i=0;i<4;i++)
       delete para[i];
      delete para;

      printf("Generated %i data points of a 4-dimensional strange submanifold in R^6\n",num);
      printf("The correct dimension of this submanifold is 4\n");

      break;

     case 4: // generates a 4-dim manifold in 8 dimensions
       dim=8;

       // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

       para=new double*[4];
       for(i=0;i<4;i++)
        para[i]=new double[num];

       for(i=0;i<4;i++)
        for(j=0;j<num;j++)
         para[i][j]=uni();

       //generate the figure
       for(i=0;i<num;i++)
       {
         data[i*dim+0]= para[1][i]*cos(2*M_PI*para[0][i]);
         data[i*dim+1]= para[1][i]*sin(2*M_PI*para[0][i]);
         data[i*dim+2]= para[2][i]*cos(2*M_PI*para[1][i]);
         data[i*dim+3]= para[2][i]*sin(2*M_PI*para[1][i]);
         data[i*dim+4]= para[3][i]*cos(2*M_PI*para[2][i]);
         data[i*dim+5]= para[3][i]*sin(2*M_PI*para[2][i]);
         data[i*dim+6]= para[0][i]*cos(2*M_PI*para[3][i]);
         data[i*dim+7]= para[0][i]*sin(2*M_PI*para[3][i]);
       }

       printf("Generated %i data points of a 4-dimensional  submanifold in R^8\n",num);
       printf("The correct dimension of this submanifold is 4\n");

       for(i=0;i<4;i++)
        delete para[i];
       delete para;
       break;

     case 5: // generates a 2d helix in 3d
      dim=3;

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      radii=new double[num];
      phi  =new double[num];

      for(i=0;i<num;i++)
      {
        radii[i]=2*M_PI*uni()-M_PI;
        phi[i]  =(uni())*2*M_PI;
      }

      // the helix
      for(i=0;i<num;i++)
      {
        data[i*dim+0]=radii[i]*sin(phi[i]);
        data[i*dim+1]=radii[i]*cos(phi[i]);
        data[i*dim+2]=phi[i];
      }

      printf("Generated %i data points of a 2-dimensional helix in R^3\n",num);
      printf("The correct dimension of this submanifold is 2\n");

      delete radii;
      delete phi;
      break;

    case 6: // generates a 6-dim manifold in 36- dimensions
       dim=36;

       // allocate memory for the output
       plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
       data = mxGetPr(plhs[0]);

       para=new double*[6];
       for(i=0;i<6;i++)
        para[i]=new double[num];

       for(i=0;i<6;i++)
        for(j=0;j<num;j++)
         para[i][j]=uni();

       //generate the figure
       for(i=0;i<num;i++)
       {
         for(j=0;j<10;j+=2)
         {
           data[i*dim+j  ]= para[j/2+1][i]*cos(2*M_PI*para[j/2][i]);
           data[i*dim+j+1]= para[j/2+1][i]*sin(2*M_PI*para[j/2][i]);
         }
         data[i*dim+10]= para[0][i]*cos(2*M_PI*para[5][i]);
         data[i*dim+11]= para[0][i]*sin(2*M_PI*para[5][i]);
       }
       for(i=0;i<num;i++)
       {
         for(j=0;j<12;j++)
         {
           data[i*dim+j+12]=data[i*dim+j];
           data[i*dim+j+24]=data[i*dim+j];
         }
       }

       printf("Generated %i data points of a 6-dimensional strange submanifold in R^36\n",num);
       printf("The correct dimension of this submanifold is 6\n");

       for(i=0;i<6;i++)
        delete para[i];
       delete para;
       break;

    case 7: //swiss roll
      dim=3;
      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      radii=new double[num];
      phi  =new double[num];

      for(i=0;i<num;i++)
      {
        radii[i]=2*M_PI*uni()-M_PI;
        phi[i]  =(uni())*2*M_PI;
      }

      // the swiss roll
      for(i=0;i<num;i++)
      {
        data[i*dim+0]=(phi[i])*sin(2.5*phi[i]);
        data[i*dim+1]=radii[i];
        data[i*dim+2]=(phi[i])*cos(2.5*phi[i]);

      }

      printf("Generated %i data points of the swiss roll in R^3\n",num);
      printf("The correct dimension of this submanifold is 2\n");

      delete radii;
      delete phi;
      break;



    case 8: // generates a 12-dim  in 72- dimensions
       dim=72;
       // allocate memory for the output
       plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
       data = mxGetPr(plhs[0]);

       para4=new double*[12];
       for(i=0;i<12;i++)
        para4[i]=new double[num];

       for(i=0;i<12;i++)
        for(j=0;j<num;j++)
         para4[i][j]=uni();

       //generate the figure
       for(i=0;i<num;i++)
       {
         for(j=0;j<22;j+=2)
         {
           data[i*dim+j  ]= para4[j/2+1][i]*cos(2*M_PI*para4[j/2][i]);
           data[i*dim+j+1]= para4[j/2+1][i]*sin(2*M_PI*para4[j/2][i]);
         }
         data[i*dim+22]= para4[0][i]*cos(2*M_PI*para4[11][i]);
         data[i*dim+23]= para4[0][i]*sin(2*M_PI*para4[11][i]);
       }
       for(i=0;i<num;i++)
       {
         for(j=0;j<24;j++)
         {
           data[i*dim+j+24]=data[i*dim+j];
           data[i*dim+j+48]=data[i*dim+j];
         }
       }

       printf("Generated %i data points of a 12-dimensional strange submanifold in R^72\n",num);
       printf("The correct dimension of this submanifold is 12\n");

       for(i=0;i<12;i++)
        delete para4[i];
       delete para4;
       break;

     case 9: // generates a 20-dim affine space  in 100- dimensions
       dim=20;

       // allocate memory for the output
       plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
       data = mxGetPr(plhs[0]);

       para=new double*[20];
       for(i=0;i<20;i++)
        para[i]=new double[num];

       for(i=0;i<20;i++)
        for(j=0;j<num;j++)
         para[i][j]=5.0*uni()-2.5;

       vectors=new double*[20];
       for(i=0;i<20;i++)
        vectors[i]=new double[dim];

       for(i=0;i<20;i++)
        for(j=0;j<dim;j++)
        {
          if(i==j)
           vectors[i][j]=1;
          else
           vectors[i][j]=0;//0.25*(uni()-0.125);
        }
         //vectors[i][j]=10*(uni()-0.5);

       //generate the figure
       for(i=0;i<num;i++)
       {
         for(j=0;j<dim;j++)
         {
           data[i*dim+j]=0;
           for(k=0;k<20;k++)
            data[i*dim+j]+=para[k][i]*vectors[k][j];
         }
       }

       printf("Generated %i data points of a 20-dimensional affine subspace in R^20\n",num);
       printf("The correct dimension of this submanifold is 20\n");

       for(i=0;i<20;i++)
        delete para[i];
       delete para;
       break;



    case 10: // generates a k-hypercube
      nr_par=dim-1;

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      para=new double*[nr_par];
      for(i=0;i<nr_par;i++)
       para[i]=new double[num];

      for(j=0;j<nr_par;j++)
       for(i=0;i<num;i++)
       {
         para[j][i]=uni();
         while(para[j][i]==0)
          para[j][i]=uni();
       }

      for(i=0;i<num;i++)
      {
        for(j=0;j<dim-1;j++)
        {
          data[i*dim+j]=para[j][i];
        }
        data[i*dim+dim-1]=0;
      }

      printf("Generated %i data points of a dim-1 dimensional hypercube R^dim\n",num);
      printf("The correct dimension of this submanifold is %i\n",dim-1);

      for(i=0;i<dim-1;i++)
       delete para[i];
      delete para;
      break;

    case 11: //generates a moebius-band in 3d
      dim=3;

      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      para = new double*[2];
      for(i=0;i<2;i++)
       para[i]=new double[num];

      for(i=0;i<num;i++)
      {
        para[0][i]=(uni())*2*M_PI;
        para[1][i]=2.0*(uni())-1.0;
      }

      // the moebius band
      for(i=0;i<num;i++)
      {
        data[i*dim+0]=(1+0.5*para[1][i]*cos(0.5*10.0*para[0][i]))*cos(para[0][i]);
        data[i*dim+1]=(1+0.5*para[1][i]*cos(0.5*10.0*para[0][i]))*sin(para[0][i]);
        data[i*dim+2]=0.5*para[1][i]*sin(0.5*10.0*para[0][i]);
      }

      printf("Generated %i data points of a 10 times tisted moebius strip in R^3\n",num);
      printf("The correct dimension of this submanifold is 2\n");

      for(i=0;i<2;i++)
       delete para[i];
      delete para;
      break;

    case 12: // generate a multivariate Gaussian
      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      para2=new double[dim];
      para2[0]=1.0;
      for(i=1;i<dim;i++)
       para2[i]=1.0;

      nr_par=dim;
      if(nr_par %2 ==1) nr_par+=1;

      para=new double*[nr_par];
      for(i=0;i<nr_par;i++)
       para[i]=new double[num];

      for(j=0;j<nr_par;j++)
       for(i=0;i<num;i++)
       {
         para[j][i]=uni();
         while(para[j][i]==0)
          para[j][i]=uni();
       }

      for(i=0;i<num;i++)
      {
        for(j=0;j<dim-2;j+=2)
        {
          data[i*dim+j  ]=para2[j]  *sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i*dim+j+1]=para2[j+1]*sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        if(dim % 2==0)
        {
          data[i*dim+j  ]=para2[j]  *sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);
          data[i*dim+j+1]=para2[j+1]*sqrt(-2*log(para[j][i]))*sin(2*M_PI*para[j+1][i]);
        }
        else
          data[i*dim+j  ]=para2[j]*sqrt(-2*log(para[j][i]))*cos(2*M_PI*para[j+1][i]);

      }

      printf("Generated %i data points of a dim dimensional Gaussian distribution with unit variance in R^dim\n",num);
      printf("The correct dimension of this submanifold is %i\n",dim);

      for(i=0;i<nr_par;i++)
       delete para[i];
      delete para;
      delete para2;
      break;

    case 13:  //generates a one-dimensional curve in dim dimensions
      // allocate memory for the output
      plhs[0] = mxCreateDoubleMatrix(dim, num, mxREAL);
      data = mxGetPr(plhs[0]);

      // parameter of the curve
      phi  =new double[num];
      for(i=0;i<num;i++)
       phi[i]  =uni()*2.0*M_PI;  // random numbers in [0,2pi]

      // the curve
      for(i=0;i<num;i++)
      {
        for(j=0;j<dim;j++)
        {
          data[i*dim+j]=phi[i]/(2.0*M_PI);
          for(k=0;k<j;k++)
           data[i*dim+j]+=sin((k+1)*phi[i]);
          data[i*dim+j]/=(1.0*(j+1));
        }
      }

      printf("Generated %i data points of a curve in R^dim\n",num);
      printf("The correct dimension of this submanifold is 1\n");
      delete phi;
      break;

  }
}

