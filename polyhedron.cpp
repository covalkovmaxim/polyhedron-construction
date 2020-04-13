#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <vector>

class point
{
public:
    double x,y,z;
    point(double x1,double y1,double z1)
    {
        x=x1;
        y=y1;
        z=z1;
    }
};
class edge
{
public:
    int coord[2];
    int del_flg;
    edge(int ind1,int ind2)
    {
        coord[0]=ind1;
        coord[1]=ind2;
        del_flg=0;
    }

};
class facet
{
public:
    double A, B, C, D;
    std::vector<int> edges;
    facet(double a,double b,double c, double d,int ind1, int ind2,int ind3,int ind4)
    {
        A=a; B=b;C=c;D=d;
        edges.push_back(ind1);
        edges.push_back(ind2);
        edges.push_back(ind3);
        edges.push_back(ind4);
    }
    facet(double a, double b, double c, double d, std::vector<int> vec)
    {
        A=a; B=b;C=c;D=d;
        edges=vec;
    }
};
class polyhedron
{
public:
    std::vector<point> points_list;
    std::vector<edge> edges_list;
    std::vector<facet> facets_list;
    polyhedron(polyhedron& init)
    {
        points_list=init.points_list;
        edges_list=init.edges_list;
        facets_list=init.facets_list;

    }
    polyhedron(double x)
    {
        points_list.push_back(point(-x,-x,-x));
        points_list.push_back(point(x,-x,-x));
        points_list.push_back(point(x,x,-x));
        points_list.push_back(point(-x,x,-x));
        points_list.push_back(point(-x,-x,x));
        points_list.push_back(point(x,-x,x));
        points_list.push_back(point(x,x,x));
        points_list.push_back(point(-x,x,x));

        edges_list.push_back(edge(0,1));//0
        edges_list.push_back(edge(0,3));//1
        edges_list.push_back(edge(0,4));//2
        edges_list.push_back(edge(1,2));//3
        edges_list.push_back(edge(1,5));//4
        edges_list.push_back(edge(2,3));//5
        edges_list.push_back(edge(2,6));//6
        edges_list.push_back(edge(3,7));//7
        edges_list.push_back(edge(4,5));//8
        edges_list.push_back(edge(4,7));//9
        edges_list.push_back(edge(5,6));//10
        edges_list.push_back(edge(6,7));//11

        facets_list.push_back(facet(0.,0.,1.,x,0,3,5,1));
        facets_list.push_back(facet(0.,1.,0.,x,0,4,8,2));
        facets_list.push_back(facet(1.,0.,1.,x,1,7,9,2));
        facets_list.push_back(facet(0.,1.,0.,-x,5,6,11,7));
        facets_list.push_back(facet(1.,0.,0.,-x,3,6,10,4));
        facets_list.push_back(facet(0.,0.,1.,-x,8,10,11,9));
    }

    void print()
    {
        FILE*fp=fopen("out.txt","w");

        int tec;

        for (int i=0;i<facets_list.size();i++)
        {
            facet tec_facet=facets_list[i];

            if(edges_list[tec_facet.edges[0]].coord[1]==edges_list[tec_facet.edges[1]].coord[0] || edges_list[tec_facet.edges[0]].coord[1]==edges_list[tec_facet.edges[1]].coord[1])
            {
               fprintf(fp,"%f %f %f\n",points_list[edges_list[tec_facet.edges[0]].coord[0]].x,points_list[edges_list[tec_facet.edges[0]].coord[0]].y,points_list[edges_list[tec_facet.edges[0]].coord[0]].z);
               tec=edges_list[tec_facet.edges[0]].coord[1];
            }
            else
            {
                fprintf(fp,"%f %f %f\n",points_list[edges_list[tec_facet.edges[0]].coord[1]].x,points_list[edges_list[tec_facet.edges[0]].coord[1]].y,points_list[edges_list[tec_facet.edges[0]].coord[1]].z);
                tec=edges_list[tec_facet.edges[0]].coord[0];
            }
            for(int j=1;j<tec_facet.edges.size();j++)
            {
                fprintf(fp,"%f %f %f\n",points_list[tec].x,points_list[tec].y,points_list[tec].z);
                if(tec==edges_list[tec_facet.edges[j]].coord[0])
                {
                    tec=edges_list[tec_facet.edges[j]].coord[1];
                }
                else
                {
                    tec=edges_list[tec_facet.edges[j]].coord[0];
                }
            }
            fprintf(fp,"%f %f %f\n  \n\n",points_list[tec].x,points_list[tec].y,points_list[tec].z);
        }
        fclose(fp);
    }
};
class plane
{
public:
    double A,B,C,D;

    plane(double a,double b,double c,double d)
    {
        A=a; B=b; C=c;D=d;
    }

};
void cross(polyhedron pol, plane space);
plane get_plane_by_three_points(double x1,double y1,double z1,double x2,double y2,double z2,double x3, double y3,double z3);
int main()
{
    double x=1.;
    polyhedron cub(x);
    cross(cub,plane(0,0,1,0));
    //cub.print();
    return 0;
}
void cross(polyhedron pol, plane space)
{
    double t,dx,dy,dz,x0,y0,z0,x,y,z;

    std::vector<std::vector<int>> points_for_new_facet;
    for (int i=0;i<pol.facets_list.size();i++)
    {
        std::vector<int> one_edge;
        facet tec_facet=pol.facets_list[i];
        for(int j=0;j<tec_facet.edges.size();j++)
        {
            edge tec_edge=pol.edges_list[tec_facet.edges[j]];

            x0=pol.points_list[tec_edge.coord[0]].x;
            dx=pol.points_list[tec_edge.coord[1]].x-pol.points_list[tec_edge.coord[0]].x;

            y0=pol.points_list[tec_edge.coord[0]].y;
            dy=pol.points_list[tec_edge.coord[1]].y-pol.points_list[tec_edge.coord[0]].y;

            z0=pol.points_list[tec_edge.coord[0]].z;
            dz=pol.points_list[tec_edge.coord[1]].z-pol.points_list[tec_edge.coord[0]].z;
            if(fabs(space.A*dx+space.B*dy+space.C*dz)<1.e-12)
            {
                if(space.A*x0+space.B*y0+space.C*z0+space.D>0)
                {
                    pol.edges_list[tec_facet.edges[j]].del_flg=1;
                }
                continue;
            }
            t=-(space.A*x0+space.B*y0+space.C*z0+space.D)/(space.A*dx+space.B*dy+space.C*dz);
            if(t<0||t>1)
            {
                if(space.A*x0+space.B*y0+space.C*z0+space.D>0)
                {
                    pol.edges_list[tec_facet.edges[j]].del_flg=1;
                }
                continue;
            }
            if(fabs(t)<1.e-12||fabs(t-1)<1.e-12)
            {
                if(fabs(t)<1.e-12)
                {
                    if(0==one_edge.size()||(1==one_edge.size()&&one_edge[0]!=tec_edge.coord[0]))
                    {
                        one_edge.push_back(tec_edge.coord[0]);
                    }
                    if(space.A*pol.points_list[tec_edge.coord[1]].x+
                       space.B*pol.points_list[tec_edge.coord[1]].y+
                       space.C*pol.points_list[tec_edge.coord[1]].z+
                       space.D>0)
                    {
                        pol.edges_list[tec_facet.edges[j]].del_flg=1;
                    }
                }
                else
                {
                    if(0==one_edge.size()||(1==one_edge.size()&&one_edge[0]!=tec_edge.coord[1]))
                    {
                        one_edge.push_back(tec_edge.coord[1]);
                    }
                    if(space.A*pol.points_list[tec_edge.coord[0]].x+
                       space.B*pol.points_list[tec_edge.coord[0]].y+
                       space.C*pol.points_list[tec_edge.coord[0]].z+
                       space.D>0)
                    {
                        pol.edges_list[tec_facet.edges[j]].del_flg=1;
                    }
                }
            }
            else
            {
                x=x0+dx*t;
                y=y0+dy*t;
                z=z0+dz*t;
                //printf("%f %f %f\n",x,y,z);
                pol.points_list.push_back(point(x,y,z));
                one_edge.push_back(pol.points_list.size()-1);
                if(pol.points_list[tec_edge.coord[0]].x*space.A+
                   pol.points_list[tec_edge.coord[0]].y*space.B+
                   pol.points_list[tec_edge.coord[0]].z*space.C+
                   space.D>0)
                {
                    pol.edges_list[tec_facet.edges[j]].coord[0]=pol.points_list.size()-1;
                }
                else
                {
                    pol.edges_list[tec_facet.edges[j]].coord[1]=pol.points_list.size()-1;
                }

            }

        }
        //printf("size=%d\n",one_edge.size());
        if(2==one_edge.size())
        {
            printf("%f %f %f, %f %f %f\n",  pol.points_list[one_edge[0]].x,
                                            pol.points_list[one_edge[0]].y,
                                            pol.points_list[one_edge[0]].z,
                                            pol.points_list[one_edge[1]].x,
                                            pol.points_list[one_edge[1]].y,
                                            pol.points_list[one_edge[1]].z);
            pol.edges_list.push_back(edge(one_edge[0],one_edge[1]));
            one_edge.push_back(pol.edges_list.size()-1);
            points_for_new_facet.push_back(one_edge);
            pol.facets_list[i].edges.push_back(one_edge[2]);
        }


        auto vec=std::begin(pol.facets_list[i].edges);
        while( vec!=std::end(pol.facets_list[i].edges))
        {
           if(pol.edges_list[*vec].del_flg==1)
           {
               vec=pol.facets_list[i].edges.erase(vec);
           }
           else
           {
               ++vec;

           }
        }

    }
    auto vec=std::begin(pol.facets_list);
    while(vec!=std::end(pol.facets_list))
    {
        if(0==(*vec).edges.size())
        {
            vec=pol.facets_list.erase(vec);
        }
        else
        {
            ++vec;
        }
    }
    pol.print();
}
