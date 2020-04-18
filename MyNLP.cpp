// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.hpp"
#include "polyhedron.hpp"
#include <cassert>


using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP()
{}

MyNLP::~MyNLP()
{}
Number df(std::function<Number(int,const Number*)> f,int n,const Number*x,int i);
std::vector<std::function<Number(int,const Number*)>> arr_g;
std::vector<std::vector<int>> support_index;
std::vector<int> part_support_index;
std::function<Number(int,const Number*)> my_functional;
polyhedron my_pol("start_model.txt");
std::vector<point> points_for_edges[2];
std::vector<double> coeffs;
int total_nonzero_jac=0;
bool MyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{

  for(auto edg=std::begin(my_pol.edges_list);edg!=std::end(my_pol.edges_list);++edg)
  {
      points_for_edges[0].push_back(my_pol.points_list[(*edg).coord[0]]);
      points_for_edges[1].push_back(my_pol.points_list[(*edg).coord[1]]);
  }
  for(int i=0;i<my_pol.edges_list.size();i++)
  {
      coeffs.push_back(0.1);
  }
  FILE*fp=fopen("corected_edges.txt","rw");
  FILE*fp1=fopen("draw_cor_edges.txt","w");
  int siz,num1,num2,tec_num;
  double xx,yy,zz;
  fscanf(fp,"%d",&siz);
  for(int i=0;i<siz;i++)
  {
      fscanf(fp,"%d %d",&num1,&num2);
      tec_num=my_pol.get_edge_num_by_two_facets(my_pol.facets_list[num1],my_pol.facets_list[num2]);
      fscanf(fp,"%lf %lf %lf",&xx,&yy,&zz);
      fprintf(fp1,"%f %f %f\n",xx,yy,zz);
      points_for_edges[0][tec_num]=point(xx,yy,zz);
      fscanf(fp,"%lf %lf %lf",&xx,&yy,&zz);
      fprintf(fp1,"%f %f %f   \n   \n\n",xx,yy,zz);
      points_for_edges[1][tec_num]=point(xx,yy,zz);
      coeffs[tec_num]=1.;
  }

  fclose(fp);
  fclose(fp1);
  std::vector<std::function<Number(int,const Number*)>> part_functional_vector;
  for(int i=0;i<points_for_edges[0].size();i++)
  {
      point tec_point1=points_for_edges[0][i];
      point tec_point2=points_for_edges[1][i];
      int tec_point_num1=my_pol.edges_list[i].coord[0];
      int tec_point_num2=my_pol.edges_list[i].coord[1];
      double my_coeff=coeffs[i];
      std::function<Number(int,const Number*)> tec_functional=
      [my_coeff,tec_point1,tec_point2,tec_point_num1,tec_point_num2](int n, const Number* x)
      {
          Number res=0.;
          Number numerator1_minus,numerator2_minus,numenator1_plus,numenator2_plus,denominator;
          //printf("%f %f %f %f %f %f\n",tec_point1.x,tec_point1.y,tec_point1.z,tec_point2.x,tec_point2.y,tec_point2.z);
          numerator1_minus=(tec_point2.x-tec_point1.x)*(x[tec_point_num1*3]-tec_point1.x)+
                     (tec_point2.y-tec_point1.y)*(x[tec_point_num1*3+1]-tec_point1.y)+
                     (tec_point2.z-tec_point1.z)*(x[tec_point_num1*3+2]-tec_point1.z);

          numenator1_plus=(x[tec_point_num1*3]-tec_point1.x)*(x[tec_point_num1*3]-tec_point1.x)+
                          (x[tec_point_num1*3+1]-tec_point1.y)*(x[tec_point_num1*3+1]-tec_point1.y)+
                          (x[tec_point_num1*3+2]-tec_point1.z)*(x[tec_point_num1*3+2]-tec_point1.z);

          numerator2_minus=(tec_point2.x-tec_point1.x)*(x[tec_point_num2*3]-tec_point1.x)+
                     (tec_point2.y-tec_point1.y)*(x[tec_point_num2*3+1]-tec_point1.y)+
                     (tec_point2.z-tec_point1.z)*(x[tec_point_num2*3+2]-tec_point1.z);

          numenator2_plus=(x[tec_point_num2*3]-tec_point1.x)*(x[tec_point_num2*3]-tec_point1.x)+
                          (x[tec_point_num2*3+1]-tec_point1.y)*(x[tec_point_num2*3+1]-tec_point1.y)+
                          (x[tec_point_num2*3+2]-tec_point1.z)*(x[tec_point_num2*3+2]-tec_point1.z);

          denominator=(tec_point2.x-tec_point1.x)*(tec_point2.x-tec_point1.x)+
                      (tec_point2.y-tec_point1.y)*(tec_point2.y-tec_point1.y)+
                      (tec_point2.z-tec_point1.z)*(tec_point2.z-tec_point1.z);

          res=(numenator1_plus+numenator2_plus)-(numerator1_minus*numerator1_minus+numerator2_minus*numerator2_minus)/denominator;


          return my_coeff*res;



      };
      part_functional_vector.push_back(tec_functional);

  }
  my_functional=[part_functional_vector](int n,const Number* x)
  {
    Number res=0.;
    for(auto funct=std::begin(part_functional_vector);funct!=std::end(part_functional_vector);++funct)
    {
        res=res+(*funct)(n,x);
    }
    return res;
  };
  int planes_index_start=3*(int)my_pol.points_list.size();
  for(int i=0;i<my_pol.facets_list.size();i++)
  {
      std::set<int> knowing_points;
      for(auto edg=std::begin(my_pol.facets_list[i].edges);edg!=std::end(my_pol.facets_list[i].edges);++edg)
      {
          int point_ind1,point_ind2;
          point_ind1=my_pol.edges_list[(*edg)].coord[0];
          point_ind2=my_pol.edges_list[(*edg)].coord[1];

          if(knowing_points.find(point_ind1)==std::end(knowing_points))
          {
              part_support_index={point_ind1*3,point_ind1*3+1,point_ind1*3+2,
                                  planes_index_start+i*4,planes_index_start+i*4+1,
                                  planes_index_start+i*4+2,planes_index_start+i*4+3};
              arr_g.push_back(
                                [point_ind1,i,planes_index_start](int n, const Number* x)
                                {
                                    return x[planes_index_start+i*4]*x[point_ind1*3]+
                                           x[planes_index_start+i*4+1]*x[point_ind1*3+1]+
                                           x[planes_index_start+i*4+2]*x[point_ind1*3+2]+
                                           x[planes_index_start+i*4+3];
                                }
                             );
              knowing_points.insert(point_ind1);
              support_index.push_back(part_support_index);
          }
          if(knowing_points.find(point_ind2)==std::end(knowing_points))
          {
              part_support_index={point_ind2*3,point_ind2*3+1,point_ind2*3+2,
                                  planes_index_start+i*4,planes_index_start+i*4+1,
                                  planes_index_start+i*4+2,planes_index_start+i*4+3};
              arr_g.push_back(
                                [point_ind2,i,planes_index_start](int n, const Number* x)
                                {
                                    return x[planes_index_start+i*4]*x[point_ind2*3]+
                                           x[planes_index_start+i*4+1]*x[point_ind2*3+1]+
                                           x[planes_index_start+i*4+2]*x[point_ind2*3+2]+
                                           x[planes_index_start+i*4+3];
                                }
                             );
              knowing_points.insert(point_ind2);
              support_index.push_back(part_support_index);
          }
      }
      part_support_index={planes_index_start+i*4,planes_index_start+i*4+1,planes_index_start+i*4+2};
      arr_g.push_back(
                        [i,planes_index_start](int n, const Number* x)
                        {
                            return x[planes_index_start+i*4]*x[planes_index_start+i*4]+
                                   x[planes_index_start+i*4+1]*x[planes_index_start+i*4+1]+
                                   x[planes_index_start+i*4+2]*x[planes_index_start+i*4+2]-1.;

                        }
                     );
      support_index.push_back(part_support_index);
  }
  for(auto ind=std::begin(support_index);ind!=std::end(support_index);++ind)
  {
      total_nonzero_jac+=(int)ind->size();
  }
  printf("total=%d\n",total_nonzero_jac);
  // The problem described in MyNLP.hpp has 2 variables, x1, & x2,
  n = 3*(int)my_pol.points_list.size()+4*(int)my_pol.facets_list.size();

  // one equality constraint,
  m = (int)arr_g.size();

  // 2 nonzeros in the jacobian (one for x1, and one for x2),
  nnz_jac_g =total_nonzero_jac;

  // and 2 nonzeros in the hessian of the lagrangian
  // (one in the hessian of the objective for x2,
  //  and one in the hessian of the constraints for x1)
  //nnz_h_lag = 40;

  // We use the standard fortran index style for row/col entries
  index_style = FORTRAN_STYLE;

  my_pol.print();
  return true;
}

bool MyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 3*(int)my_pol.points_list.size()+4*(int)my_pol.facets_list.size());
  assert(m == (int)arr_g.size());

  // x1 has a lower bound of -1 and an upper bound of 1
  //         x_l[0] = -1.0;
  //         x_u[0] = 1.0;

  // x2 has no upper or lower bound, so we set them to
  // a large negative and a large positive number.
  // The value that is interpretted as -/+infinity can be
  // set in the options, but it defaults to -/+1e19
  //        x_l[1] = -1.0e19;
  //        x_u[1] = +1.0e19;
  for(int i=0;i<n;i++)
  {
      x_l[i] = -1.0e19;
      x_u[i] = +1.0e19;

  }

  // we have one equality constraint, so we set the bounds on this constraint
  // to be equal (and zero).
  //        g_l[0] = g_u[0] = 0.0;

  for(int i=0;i<m;i++)
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

  for(int i=0;i<my_pol.points_list.size();i++)
  {
      x[i*3+0]=my_pol.points_list[i].x;
      x[i*3+1]=my_pol.points_list[i].y;
      x[i*3+2]=my_pol.points_list[i].z;
  }
  int planes_index_start=3*(int)my_pol.points_list.size();
  for(int i=0;i<my_pol.facets_list.size();i++)
  {
      x[planes_index_start+i*4+0]=my_pol.facets_list[i].A;
      x[planes_index_start+i*4+1]=my_pol.facets_list[i].B;
      x[planes_index_start+i*4+2]=my_pol.facets_list[i].C;
      x[planes_index_start+i*4+3]=my_pol.facets_list[i].D;
  }
  return true;
}

bool MyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
 /* Number x2 = x[1];
  Number x1=x[0];
  */

  //obj_value = -(x2 - 2.0) * (x1 - 2.0);
  obj_value=my_functional(n,x);

  return true;
}

bool MyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{

    for(int i=0;i<n;i++)
    {
        grad_f[i]=df(my_functional,n,x,i);

    }


  return true;
}

bool MyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    for(int i=0;i<m;i++)
    {
        g[i]=arr_g[i](n,x);
    }

  return true;
}

bool MyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  int tec_num=0;
  if (values == NULL) {
    // return the structure of the jacobian of the constraints

    // element at 1,1: grad_{x1} g_{1}(x)
    //      iRow[0] = 1;
    //      jCol[0] = 1;

    // element at 1,2: grad_{x2} g_{1}(x)
    //      iRow[1] = 1;
    //      jCol[1] = 2;
      tec_num=0;
      for(int i=0;i<m;i++)
      {
          for(auto gg=std::begin(support_index[i]);gg!=std::end(support_index[i]);++gg)
          {
              iRow[tec_num]=i+1;
              jCol[tec_num]=(*gg)+1;
              tec_num++;
          }
      }

  }
  else {
      tec_num=0;
      for(int i=0;i<m;i++)
      {
          for(auto gg=std::begin(support_index[i]);gg!=std::end(support_index[i]);++gg)
          {
              values[tec_num]=df(arr_g[i],n,x,(*gg));
              tec_num++;
          }
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
    /*for(int i=0;i<4;i++)
    {
        printf("x=%f y=%f z=%f\n",(double) x[3*i+0],(double) x[3*i+1],(double) x[3*i+2]);
    }*/
    /*for(int i=0;i<my_pol.points_list.size();i++)
    {
        printf("%f %f %f\n",x[i*3+0],x[i*3+1],x[i*3+2]);
        my_pol.points_list[i]=point(x[i*3+0],x[i*3+1],x[i*3+2]);

    }*/
    printf("%f\n",my_functional(n,x));
    point mass_center(0.,0.,0.);
    for(auto p=std::begin(my_pol.points_list);p!=std::end(my_pol.points_list);++p)
    {
        mass_center.x+=(*p).x;
        mass_center.y+=(*p).y;
        mass_center.z+=(*p).z;
    }
    mass_center.x/=(double)(my_pol.points_list.size());
    mass_center.y/=(double)(my_pol.points_list.size());
    mass_center.z/=(double)(my_pol.points_list.size());
    int planes_index_start=3*(int)my_pol.points_list.size();
    std::vector<plane> planes;
    //printf("center: %f %f %f\n",mass_center.x,mass_center.y,mass_center.z);
    for(int i=0;i<my_pol.facets_list.size();i++)
    {
        printf("plane: %f %f %f %f\n",x[planes_index_start+i*4+0],x[planes_index_start+i*4+1],x[planes_index_start+i*4+2],x[planes_index_start+i*4+3]);
        if(x[planes_index_start+i*4+0]*mass_center.x+
           x[planes_index_start+i*4+1]*mass_center.y+
           x[planes_index_start+i*4+2]*mass_center.z+
           x[planes_index_start+i*4+3]>0)
        {

            planes.push_back(plane(-x[planes_index_start+i*4+0],
                                   -x[planes_index_start+i*4+1],
                                   -x[planes_index_start+i*4+2],
                                   -x[planes_index_start+i*4+3]));

        }
        else
        {
            planes.push_back(plane(x[planes_index_start+i*4+0],
                                   x[planes_index_start+i*4+1],
                                   x[planes_index_start+i*4+2],
                                   x[planes_index_start+i*4+3]));
        }
    }
    my_pol=construct_polyhedron_by_planes_list(&planes);
    printf("%d %d %d\n",my_pol.points_list.size(),my_pol.edges_list.size(),my_pol.facets_list.size());
    my_pol.print();
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
}


Number df(std::function<Number(int,const Number*)> f,int n,const Number*x,int i)
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
            x1[i]=x[i]-1.e-6;
            x2[i]=x[i]+1.e-6;
        }
    }
    //printf("\nup:%lf down:%lf\n",x2[i],x1[i]);
    Number res=(f(n,x2)-f(n,x1))/2.e-6;

    delete[]x1;
    delete[]x2;
    return res;
}

