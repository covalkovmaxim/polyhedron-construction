#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>

#include "polyhedron.hpp"
#include "gnuplot.hpp"



std::vector<target_edge> read_target_edges(const char* filename,const polyhedron& pol,std::vector<point>* points_for_edges,std::vector<double>& coeffs);
void jordan(double*A,double*b);
point get_point_by_three_planes_intersection(const plane& p1,const plane& p2,const plane& p3);
void calc_gradient(const std::vector<target_edge>& target_edges_list, const std::vector<double>& init_vars, double* gradient);
double scalar_product(int n, double* old_dirrection, double* new_dirrection);
void plot_hist(const polyhedron& my_pol,const std::vector<target_edge>& target_edges_list,const std::vector<double>& init_vars,const std::vector<double>& result_vars,const std::string& image_name);
double myvfunc(const std::vector<double> &vars, std::vector<double> &grad, void *data);
double dist_func_1(std::vector<double>& vars, const std::vector<target_edge>& target_edges_list);

int main(int argc, char** argv)
{
    if(3!=argc)
    {
        printf("Correct arguments initpoly_filename corred_filename\n");
        return -1;
    }

    polyhedron my_pol(argv[1]);
    std::vector<point> points_for_edges[2];
    std::vector<double> coeffs((int)my_pol.edges_list.size(),0.1);
    std::function<double (double )> fun_for_lambda;

    for(auto edg=std::begin(my_pol.edges_list);edg!=std::end(my_pol.edges_list);++edg)
    {
        points_for_edges[0].push_back(my_pol.points_list[edg->coord[0]]);
        points_for_edges[1].push_back(my_pol.points_list[edg->coord[1]]);
    }
    std::vector<target_edge> target_edges_list=read_target_edges(argv[2],my_pol,points_for_edges,coeffs);

    for(auto target : target_edges_list)
    {
        if(my_pol.get_edge_num_by_two_facets(my_pol.facets_list[target.facet_numbers[0]],my_pol.facets_list[target.facet_numbers[1]])!=target.edge_num)
        {
            printf("Edge is a intersection of more than 2 planes. It's wrong\n");
        }

        facet face_1=my_pol.facets_list[target.facet_numbers[0]];
        //printf("%d %f %f %f %f\n", target.facet_numbers[0], face_1.A, face_1.B, face_1.C, face_1.D);
        facet face_2=my_pol.facets_list[target.facet_numbers[1]];
        //printf("%d %f %f %f %f\n", target.facet_numbers[1], face_2.A, face_2.B, face_2.C, face_2.D);
        point p_1=get_point_by_three_planes_intersection(plane(face_1.A,face_1.B,face_1.C,face_1.D),plane(face_2.A,face_2.B,face_2.C,face_2.D),target.normal_plane_1);
        point p_2=get_point_by_three_planes_intersection(plane(face_1.A,face_1.B,face_1.C,face_1.D),plane(face_2.A,face_2.B,face_2.C,face_2.D),target.normal_plane_2);

        if(target.coeff<1.)
        {

            if(p_1!=my_pol.points_list[my_pol.edges_list[target.edge_num].coord[0]])
            {
                printf("%f first val %f %f %f second val %f %f %f\n",target.coeff,p_1.x,p_1.y,p_1.z,my_pol.points_list[my_pol.edges_list[target.edge_num].coord[0]].x,
                        my_pol.points_list[my_pol.edges_list[target.edge_num].coord[0]].y,my_pol.points_list[my_pol.edges_list[target.edge_num].coord[0]].z);
                //printf("%d %d\n", target.facet_numbers[0], target.facet_numbers[1]);
            }
        }

    }

    std::vector<double> init_vars, next_vars, init_vars_copy;
    for(auto face : my_pol.facets_list)
    {
        init_vars.push_back(face.A);
        init_vars.push_back(face.B);
        init_vars.push_back(face.C);
        init_vars.push_back(face.D);
    }
    init_vars_copy=init_vars;
    int var_num = init_vars.size();

    nlopt::opt opt(nlopt::LD_LBFGS, var_num);
    opt.set_min_objective(myvfunc, (void*)&target_edges_list);
    opt.set_xtol_rel(1e-4);
    double minf;
    opt.optimize(init_vars,minf);
    printf("min=%f\n", minf);
    next_vars=init_vars;

    std::vector<plane> planes;

    point mass_center(0.,0.,0.);
    for(auto point : my_pol.points_list)
    {
        mass_center.x+=point.x;
        mass_center.y+=point.y;
        mass_center.z+=point.z;
    }
    mass_center.x/=(double)(my_pol.points_list.size());
    mass_center.y/=(double)(my_pol.points_list.size());
    mass_center.z/=(double)(my_pol.points_list.size());

    for(int i=0;i<var_num/4;i++)
    {
        if(next_vars[i*4]*mass_center.x+
           next_vars[i*4+1]*mass_center.y+
           next_vars[i*4+2]*mass_center.z+
           next_vars[i*4+3]>0.)
        {
            planes.push_back(plane(-next_vars[i*4],-next_vars[i*4+1],-next_vars[i*4+2],-next_vars[i*4+3]));
        }
        else
        {
            planes.push_back(plane(next_vars[i*4],next_vars[i*4+1],next_vars[i*4+2],next_vars[i*4+3]));
        }
    }

    for(auto plane : planes)
    {
        printf("plane: %f %f %f %f\n", plane.A, plane.B, plane.C, plane.D);
    }
    plot_hist(my_pol,target_edges_list,init_vars_copy,next_vars,"hists.png");
    my_pol=construct_polyhedron_by_planes_list(&planes);
    my_pol.print();
    return 0;
}

std::vector<target_edge> read_target_edges(const char* filename,const polyhedron& pol,std::vector<point>* points_for_edges,std::vector<double>& coeffs)
{
    FILE*fp=fopen(filename, "rw");
    FILE*fp1=fopen("draw_cor_edges.txt","w");

    int siz,num1,num2,tec_num;
    double xx,yy,zz;
    std::set<int> changed_facets_set;

    fscanf(fp,"%d",&siz);
    for(int i=0;i<siz;i++)
    {
        fscanf(fp,"%d %d",&num1,&num2);
        changed_facets_set.insert(num1);
        changed_facets_set.insert(num2);
        tec_num=pol.get_edge_num_by_two_facets(pol.facets_list[num1],pol.facets_list[num2]);
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

    std::set<int> changed_edges_set;

    for(auto facet_num=changed_facets_set.begin();facet_num!=changed_facets_set.end();++facet_num)
    {
        facet face=pol.facets_list[(*facet_num)];
        for(auto edge_num=face.edges.begin();edge_num!=face.edges.end();++edge_num)
        {
            changed_edges_set.insert(*edge_num);
        }
    }

    std::vector<int> changed_edges_vector(changed_edges_set.begin(),changed_edges_set.end());
    int len1=(int)changed_edges_vector.size();
    std::vector<target_edge> target_edges_list;

    for(int i=0;i<len1;i++)
    {
        std::vector<int> two_facets=pol.get_facet_nums_by_edge_num(changed_edges_vector[i]);
        int num=changed_edges_vector[i];
        target_edges_list.push_back(target_edge(num,coeffs[num],points_for_edges[0][num],points_for_edges[1][num],two_facets));
        if(2!=target_edges_list.back().facet_numbers.size())
        {
            printf("Very bad\n");
        }
    }

    return target_edges_list;

}

void jordan(double*A,double*b)
{
    double cur=-1.;
    int index;
    for(int t=0;t<3;t++)
    {
        cur=-1.;
        for(int i=t;i<3;i++)
        {
            if(fabs(A[i*3+t])>cur)
            {
                cur=fabs(A[i*3+t]);
                index=i;
            }
        }
        for(int j=0;j<3;j++)
        {
            cur=A[t*3+j];
            A[t*3+j]=A[index*3+j];
            A[index*3+j]=cur;
        }
        cur=b[t];
        b[t]=b[index];
        b[index]=cur;

        cur=1./A[t*3+t];
        for(int j=0;j<3;j++)
        {
            A[t*3+j]*=cur;
        }
        b[t]*=cur;

        for(int i=0;i<3;i++)
        {
            if(i==t)
            {
                continue;
            }
            for(int j=t+1;j<3;j++)
            {

                A[i*3+j]-=A[i*3+t]*A[t*3+j];
            }
            b[i]-=A[i*3+t]*b[t];
        }
    }
}

point get_point_by_three_planes_intersection(const plane& p1,const plane& p2,const plane& p3)
{
    double A[9], b[3];
    A[0]=p1.A; A[1]=p1.B; A[2]=p1.C; b[0]=-p1.D;
    A[3]=p2.A; A[4]=p2.B; A[5]=p2.C; b[1]=-p2.D;
    A[6]=p3.A; A[7]=p3.B; A[8]=p3.C; b[2]=-p3.D;
    //printf("\n%f %f %f\n%f %f %f\n%f %f %f\n\n",A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8]);
    jordan(A,b);
    return point(b[0],b[1],b[2]);
}

void calc_gradient(const std::vector<target_edge>& target_edges_list, const std::vector<double>& init_vars, double* gradient)
{
    int var_num = (int)init_vars.size();
    for(int i=0; i<var_num; i++)
    {
        gradient[i]=0.;

    }

    for(auto target : target_edges_list)
    {
        int num_1=target.facet_numbers[0];
        int num_2=target.facet_numbers[1];

        point p_1=get_point_by_three_planes_intersection(plane(init_vars[num_1*4],init_vars[num_1*4+1],init_vars[num_1*4+2],init_vars[num_1*4+3]),
                                                         plane(init_vars[num_2*4],init_vars[num_2*4+1],init_vars[num_2*4+2],init_vars[num_2*4+3]),
                                                         target.normal_plane_1);

        point p_2=get_point_by_three_planes_intersection(plane(init_vars[num_1*4],init_vars[num_1*4+1],init_vars[num_1*4+2],init_vars[num_1*4+3]),
                                                         plane(init_vars[num_2*4],init_vars[num_2*4+1],init_vars[num_2*4+2],init_vars[num_2*4+3]),
                                                         target.normal_plane_2);

        double loc_grad[3], loc_matr[9], A[9];

        loc_grad[0]=-p_1.x;
        loc_grad[1]=0.;
        loc_grad[2]=0.;

        loc_matr[0]=init_vars[num_1*4]; loc_matr[1]=init_vars[num_1*4+1]; loc_matr[2]=init_vars[num_1*4+2];
        loc_matr[3]=init_vars[num_2*4]; loc_matr[4]=init_vars[num_2*4+1]; loc_matr[5]=init_vars[num_2*4+2];
        loc_matr[6]=target.normal_plane_1.A; loc_matr[7]=target.normal_plane_1.B; loc_matr[8]=target.normal_plane_1.C;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                            2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                            2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;


        loc_grad[0]=-p_1.y;
        loc_grad[1]=0.;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+1]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;


        loc_grad[0]=-p_1.z;
        loc_grad[1]=0.;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+2]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;


        loc_grad[0]=-1.;
        loc_grad[1]=0.;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+3]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;


        loc_grad[0]=0.;
        loc_grad[1]=-p_1.x;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+0]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0.;
        loc_grad[1]=-p_1.y;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+1]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0.;
        loc_grad[1]=-p_1.z;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+2]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0.;
        loc_grad[1]=-1.;
        loc_grad[2]=0.;
        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+3]+=(2.*(p_1.x-target.target_point_1.x)*loc_grad[0]+
                              2.*(p_1.y-target.target_point_1.y)*loc_grad[1]+
                              2.*(p_1.z-target.target_point_1.z)*loc_grad[2])*target.coeff;

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        loc_matr[0]=init_vars[num_1*4]; loc_matr[1]=init_vars[num_1*4+1]; loc_matr[2]=init_vars[num_1*4+2];
        loc_matr[3]=init_vars[num_2*4]; loc_matr[4]=init_vars[num_2*4+1]; loc_matr[5]=init_vars[num_2*4+2];
        loc_matr[6]=target.normal_plane_2.A; loc_matr[7]=target.normal_plane_2.B; loc_matr[8]=target.normal_plane_2.C;

        loc_grad[0]=-p_2.x;
        loc_grad[1]=0.;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+0]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;


        loc_grad[0]=-p_2.y;
        loc_grad[1]=0.;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+1]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;


        loc_grad[0]=-p_2.z;
        loc_grad[1]=0.;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+2]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=-1.;
        loc_grad[1]=0.;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_1*4+3]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0;
        loc_grad[1]=-p_2.x;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+0]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0;
        loc_grad[1]=-p_2.y;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+1]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0;
        loc_grad[1]=-p_2.z;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+2]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;

        loc_grad[0]=0;
        loc_grad[1]=-1.;
        loc_grad[2]=0.;

        for(int i=0;i<9;i++)
        {
            A[i]=loc_matr[i];
        }
        jordan(A,loc_grad);

        gradient[num_2*4+3]+=(2.*(p_2.x-target.target_point_2.x)*loc_grad[0]+
                              2.*(p_2.y-target.target_point_2.y)*loc_grad[1]+
                              2.*(p_2.z-target.target_point_2.z)*loc_grad[2])*target.coeff;
    }

}

double scalar_product(int n, double* old_dirrection, double* new_dirrection)
{
    double res = 0.;
    for(int i=0; i<n; i++)
    {
        res += old_dirrection[i]*new_dirrection[i];
    }
    return  res;
}

void plot_hist(const polyhedron& my_pol,const std::vector<target_edge>& target_edges_list,const std::vector<double>& init_vars,const std::vector<double>& result_vars,const std::string& image_name)
{
    FILE*fp=fopen("hist_data.dat","w");
    int i=0;
    for(auto target : target_edges_list)
    {
        i++;
        int num_1=target.facet_numbers[0];
        int num_2=target.facet_numbers[1];
        double res_1,res_2;
        point p_1=get_point_by_three_planes_intersection(plane(init_vars[num_1*4],init_vars[num_1*4+1],init_vars[num_1*4+2],init_vars[num_1*4+3]),
                                                         plane(init_vars[num_2*4],init_vars[num_2*4+1],init_vars[num_2*4+2],init_vars[num_2*4+3]),
                                                         target.normal_plane_1);
        point p_2=get_point_by_three_planes_intersection(plane(init_vars[num_1*4],init_vars[num_1*4+1],init_vars[num_1*4+2],init_vars[num_1*4+3]),
                                                         plane(init_vars[num_2*4],init_vars[num_2*4+1],init_vars[num_2*4+2],init_vars[num_2*4+3]),
                                                         target.normal_plane_2);
        res_1=target.coeff*((target.target_point_1-p_1)*(target.target_point_1-p_1)+(target.target_point_2-p_2)*(target.target_point_2-p_2));

        p_1=get_point_by_three_planes_intersection(plane(result_vars[num_1*4],result_vars[num_1*4+1],result_vars[num_1*4+2],result_vars[num_1*4+3]),
                                                         plane(result_vars[num_2*4],result_vars[num_2*4+1],result_vars[num_2*4+2],result_vars[num_2*4+3]),
                                                         target.normal_plane_1);
        p_2=get_point_by_three_planes_intersection(plane(result_vars[num_1*4],result_vars[num_1*4+1],result_vars[num_1*4+2],result_vars[num_1*4+3]),
                                                         plane(result_vars[num_2*4],result_vars[num_2*4+1],result_vars[num_2*4+2],result_vars[num_2*4+3]),
                                                         target.normal_plane_2);
        res_2=target.coeff*((target.target_point_1-p_1)*(target.target_point_1-p_1)+(target.target_point_2-p_2)*(target.target_point_2-p_2));
        fprintf(fp,"%d %e %e\n",i,res_1,res_2);
    }
    fclose(fp);

    Gnuplot plot;

    plot("set terminal png size 1500, 800");
    plot("set output '"+image_name+"'");

    plot("set style data histograms");
    plot("set style histogram cluster");

    plot("plot 'hist_data.dat' u 2 fs solid 0.5 lt rgb 'red', '' u 3 fs solid 0.5 lt rgb 'green'");


}



double myvfunc(const std::vector<double> &vars, std::vector<double> &grad, void *data)
{
    std::vector<target_edge>* ptr=reinterpret_cast< std::vector<target_edge>*>(data);
    std::vector<target_edge> target_edges_list=*ptr;
    std::vector<double> norm_vars=vars;
    int var_num=vars.size();
    double res=dist_func_1(norm_vars,target_edges_list);

    double* gradient = new double[var_num];
    calc_gradient(target_edges_list,norm_vars,gradient);
    grad=std::vector<double>(gradient, gradient + var_num);
    delete [] gradient;
    return res;
}

double dist_func_1(std::vector<double>& vars, const std::vector<target_edge>& target_edges_list)
{
    double res=0.;
    int var_num=vars.size();
    for(int i=0;i<var_num/4;i++)
    {
        double norma = sqrt(pow(vars[i*4],2.)+pow(vars[i*4+1],2.)+pow(vars[i*4+2],2.));
        if(norma>1e-10)
        {
            vars[i*4]/=norma;
            vars[i*4+1]/=norma;
            vars[i*4+2]/=norma;
            vars[i*4+3]/=norma;
        }
    }
    for(auto target : target_edges_list)
    {
        int num_1=target.facet_numbers[0];
        int num_2=target.facet_numbers[1];
        point p_1=get_point_by_three_planes_intersection(plane(vars[num_1*4],vars[num_1*4+1],vars[num_1*4+2],vars[num_1*4+3]),
                                                         plane(vars[num_2*4],vars[num_2*4+1],vars[num_2*4+2],vars[num_2*4+3]),
                                                         target.normal_plane_1);
        point p_2=get_point_by_three_planes_intersection(plane(vars[num_1*4],vars[num_1*4+1],vars[num_1*4+2],vars[num_1*4+3]),
                                                         plane(vars[num_2*4],vars[num_2*4+1],vars[num_2*4+2],vars[num_2*4+3]),
                                                         target.normal_plane_2);
        res+=target.coeff*((target.target_point_1-p_1)*(target.target_point_1-p_1)+(target.target_point_2-p_2)*(target.target_point_2-p_2));
    }
    return res;
}
