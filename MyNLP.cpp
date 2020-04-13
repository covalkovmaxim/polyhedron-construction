// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.hpp"

#include <cassert>

typedef Number (*func)(int, const Number*);
using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP()
{}

MyNLP::~MyNLP()
{}
Number df(Number f(int,const Number*),int n,const Number*x,int i);
Number gf0(int n, const Number*x);
Number gf1(int n, const Number*x);
Number gf2(int n, const Number*x);
Number gf3(int n, const Number*x);
Number gf4(int n, const Number*x);

Number f(int n,const Number*x);
bool MyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in MyNLP.hpp has 2 variables, x1, & x2,
  n = 16;

  // one equality constraint,
  m = 5;

  // 2 nonzeros in the jacobian (one for x1, and one for x2),
  nnz_jac_g =80;

  // and 2 nonzeros in the hessian of the lagrangian
  // (one in the hessian of the objective for x2,
  //  and one in the hessian of the constraints for x1)
  nnz_h_lag = 40;

  // We use the standard fortran index style for row/col entries
  index_style = FORTRAN_STYLE;

  return true;
}

bool MyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 16);
  assert(m == 5);

  // x1 has a lower bound of -1 and an upper bound of 1
  //         x_l[0] = -1.0;
  //         x_u[0] = 1.0;

  // x2 has no upper or lower bound, so we set them to
  // a large negative and a large positive number.
  // The value that is interpretted as -/+infinity can be
  // set in the options, but it defaults to -/+1e19
  //        x_l[1] = -1.0e19;
  //        x_u[1] = +1.0e19;
  for(int i=0;i<16;i++)
  {
      x_l[i] = -1.0e19;
      x_u[i] = +1.0e19;

  }

  // we have one equality constraint, so we set the bounds on this constraint
  // to be equal (and zero).
  //        g_l[0] = g_u[0] = 0.0;
  g_l[0] = g_u[0] = 1.0;
  for(int i=1;i<5;i++)
  {
      g_l[i]=g_u[i]=0.0;
  }

  return true;
}

bool MyNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // we initialize x in bounds, in the upper right quadrant
  //    x[0] = 0.5;
  //    x[1] = 1.5;
  x[0]=1.0; x[1]=0.0; x[2]=1.0;
  x[3]=0.0; x[4]=0.0; x[5]=1.0;
  x[6]=0.0; x[7]=1.0; x[8]=1.0;
  x[9]=1.0; x[10]=1.0; x[11]=1.0;

  x[13]=x[14]=x[15]=0.0;
  x[14]=1.0;

  return true;
}

bool MyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
 /* Number x2 = x[1];
  Number x1=x[0];
  */

  //obj_value = -(x2 - 2.0) * (x1 - 2.0);
  obj_value=f(16,x);

  return true;
}

bool MyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{

    for(int i=0;i<16;i++)
    {
        grad_f[i]=df(f,16,x,i);
    }


  return true;
}

bool MyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{

    g[0]=gf0(16,x);
    g[1]=gf1(16,x);
    g[2]=gf2(16,x);
    g[3]=gf3(16,x);
    g[4]=gf4(16,x);
  return true;
}

bool MyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  func arr_g[5];
  arr_g[0]=gf0;
  arr_g[1]=gf1;
  arr_g[2]=gf2;
  arr_g[3]=gf3;
  arr_g[4]=gf4;
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

    // element at 1,1: grad_{x1} g_{1}(x)
    //      iRow[0] = 1;
    //      jCol[0] = 1;

    // element at 1,2: grad_{x2} g_{1}(x)
    //      iRow[1] = 1;
    //      jCol[1] = 2;
      for(int i=0;i<5;i++)
      {
          for(int j=0;j<16;j++)
          {
              iRow[i*16+j]=i+1;
              jCol[i*16+j]=j+1;
          }
      }

  }
  else {

      for(int i=0;i<5;i++)
          for(int j=0;j<16;j++)
          {

              values[i*16+j]=df(arr_g[i],16,x,j);
          }

  }

  return true;
}

bool MyNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  return true;
}

void MyNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{
    for(int i=0;i<4;i++)
    {
        printf("x=%f y=%f z=%f\n",(double) x[3*i+0],(double) x[3*i+1],(double) x[3*i+2]);
    }
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
}


Number f(int n, const Number*x)
{
    FILE*fp=fopen("input.txt","r");
    Number *X=new Number[4];
    double *input_x=new double[4];
    for(int i=0;i<4;i++)
    {
        X[i]=x[3*i];
    }
    Number *Y=new Number[4];
    double *input_y=new double[4];
    for(int i=0;i<4;i++)
    {
        Y[i]=x[3*i+1];
    }
    Number *Z=new Number[4];
    double *input_z=new double[4];
    for(int i=0;i<4;i++)
    {
        Z[i]=x[3*i+2];
    }
    Number A=x[12],B=x[13],C=x[14],D=x[15];
    for(int i=0;i<4;i++)
    {
        fscanf(fp,"%lf %lf %lf",&input_x[i],&input_y[i],&input_z[i]);
    }
    //obj_value = -(x2 - 2.0) * (x1 - 2.0);
    Number res=(X[0]-input_x[0])*(X[0]-input_x[0])+(Y[0]-input_y[0])*(Y[0]-input_y[0])+(Z[0]-input_z[0])*(Z[0]-input_z[0])
            +(X[1]-input_x[1])*(X[1]-input_x[1])+(Y[1]-input_y[1])*(Y[1]-input_y[1])+(Z[1]-input_z[1])*(Z[1]-input_z[1])
            +(X[2]-input_x[2])*(X[2]-input_x[2])+(Y[2]-input_y[2])*(Y[2]-input_y[2])+(Z[2]-input_z[2])*(Z[2]-input_z[2])
            +(X[3]-input_x[3])*(X[3]-input_x[3])+(Y[3]-input_y[3])*(Y[3]-input_y[3])+(Z[3]-input_z[3])*(Z[3]-input_z[3]);
    delete[]X;
    delete[]Y;
    delete[]Z;

    delete[]input_x;
    delete[]input_y;
    delete[]input_z;

    fclose(fp);
    return res;
}
Number df(Number f(int n,const Number*),int n,const Number*x,int i)
{
    Number*x1=new Number[n];
    Number*x2=new Number[n];
    for(int j=0;j<n;j++)
    {
        if(j!=i)
        {
            x1[j]=x2[j]=x[j];
        }
        else
        {
            //printf("\ninput:%lf \n",x[i]);
            x1[i]=x[i]-0.000001;
            x2[i]=x[i]+0.000001;
        }
    }
    //printf("\nup:%lf down:%lf\n",x2[i],x1[i]);
    Number res=(f(n,x2)-f(n,x1))/0.000002;

    delete[]x1;
    delete[]x2;
    return res;
}
Number gf0(int n, const Number*x)
{
    return x[14]*x[14]+x[12]*x[12]+x[13]*x[13];
}
Number gf1(int n, const Number*x)
{
    return x[12]*x[0]+x[13]*x[1]+x[14]*x[2]+x[15];
}
Number gf2(int n, const Number*x)
{
    return x[12]*x[3]+x[13]*x[4]+x[14]*x[5]+x[15];
}
Number gf3(int n, const Number*x)
{
    return x[12]*x[6]+x[13]*x[7]+x[14]*x[8]+x[15];
}
Number gf4(int n, const Number*x)
{
    return x[12]*x[9]+x[13]*x[10]+x[14]*x[11]+x[15];
}
