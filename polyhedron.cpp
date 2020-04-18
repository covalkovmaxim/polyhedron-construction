#include "polyhedron.hpp"

polyhedron::polyhedron(char *name)
{
    FILE*fp=fopen(name,"rw");
    int n,m,t1,t2;
    double x,y,z;
    fscanf(fp,"%d",&n);
    for(int i=0;i<n;i++)
    {
        fscanf(fp,"%lf%lf%lf",&x,&y,&z);
        points_list.push_back(point(x,y,z));
    }
    fscanf(fp,"%d",&n);
    for(int i=0;i<n;i++)
    {
        fscanf(fp,"%d%d",&t1,&t2);
        edges_list.push_back(edge(t1,t2));
    }
    fscanf(fp,"%d",&n);
    for(int i=0;i<n;i++)
    {
        std::vector<int> edges_num;
        fscanf(fp,"%d",&m);
        for(int j=0;j<m;j++)
        {
            fscanf(fp,"%d",&t1);
            edges_num.push_back(t1);
        }
        plane this_plane=this->get_plane_by_two_edges(edges_list[edges_num[0]],edges_list[edges_num[1]]);
        facets_list.push_back(facet(this_plane.A,this_plane.B,this_plane.C,this_plane.D,edges_num));
    }

    fclose(fp);
}
polyhedron::polyhedron(double x)
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
void polyhedron::print()
{
    FILE*fp=fopen("out.txt","w");

    int tec;

    for (int i=0;i<facets_list.size();i++)
    {
        facet tec_facet=facets_list[i];

        for(int j=0;j<tec_facet.edges.size();j++)
        {
            edge tec_edge=edges_list[tec_facet.edges[j]];
            fprintf(fp,"%f %f %f\n %f %f %f \n  \n\n",
                    points_list[tec_edge.coord[0]].x,points_list[tec_edge.coord[0]].y,points_list[tec_edge.coord[0]].z,
                    points_list[tec_edge.coord[1]].x,points_list[tec_edge.coord[1]].y,points_list[tec_edge.coord[1]].z);
        }

    }
    fclose(fp);
}
polyhedron cross(polyhedron pol, plane space)
{
    double t,dx,dy,dz,x0,y0,z0,x,y,z;
    int key1=0, key2=0;
    std::vector<int> points_for_new_facet;
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
                if(space.A*x0+space.B*y0+space.C*z0+space.D>0.)
                {
                    pol.edges_list[tec_facet.edges[j]].del_flg=1;
                }
                continue;
            }
            t=-(space.A*x0+space.B*y0+space.C*z0+space.D)/(space.A*dx+space.B*dy+space.C*dz);
            //printf("t=%e\n",t);

            if(fabs(t)<1.e-12||fabs(t-1.)<1.e-12)
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
            if(t<0||t>1)
            {
                if(space.A*x0+space.B*y0+space.C*z0+space.D>0.)
                {
                    pol.edges_list[tec_facet.edges[j]].del_flg=1;
                }
                continue;
            }
            if(fabs(t)>1.e-12&&fabs(t-1.)>1.e-12)
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
            key2=-1;
            for(int j=0;j<pol.edges_list.size();j++)
            {
                if((pol.edges_list[j].coord[0]==one_edge[0]&&pol.edges_list[j].coord[1]==one_edge[1])||
                   (pol.edges_list[j].coord[1]==one_edge[0]&&pol.edges_list[j].coord[0]==one_edge[1]))
                {
                    key2=j;
                }
            }
            if(-1==key2)
            {
                pol.edges_list.push_back(edge(one_edge[0],one_edge[1]));
                one_edge.push_back(pol.edges_list.size()-1);
                points_for_new_facet.push_back(pol.edges_list.size()-1);
                pol.facets_list[i].edges.push_back(one_edge[2]);
                key1=1;
            }
            else
            {
                one_edge.push_back(key2);
                points_for_new_facet.push_back(key2);
                pol.facets_list[i].edges.push_back(one_edge[2]);
            }

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
    if(1>=key1)
    {
        pol.facets_list.push_back(facet(space.A,space.B,space.C,space.D,points_for_new_facet));
    }
    auto vec=std::begin(pol.facets_list);
    while(vec!=std::end(pol.facets_list))
    {
        if(2>=(*vec).edges.size())
        {
            vec=pol.facets_list.erase(vec);
        }
        else
        {
            ++vec;
        }
    }

    std::vector<edge> new_edges_list;
    std::vector<int> new_edges_num;
    for(auto i=std::begin(pol.edges_list);i!=std::end(pol.edges_list);++i)
    {
        if (0==(*i).del_flg)
        {
            new_edges_num.push_back(new_edges_list.size());
            new_edges_list.push_back(*i);

        }
        else
        {
            new_edges_num.push_back(-1);
        }
    }
    pol.edges_list=new_edges_list;
    for(int i=0;i<pol.facets_list.size();i++)
    {
        for(int j=0;j<pol.facets_list[i].edges.size();j++)
        {
            pol.facets_list[i].edges[j]=new_edges_num[pol.facets_list[i].edges[j]];
        }
    }

    std::set<int> using_points;
    std::vector<point> new_points_list;
    std::vector<int> new_points_num;

    for(auto edg=std::begin(pol.edges_list);edg!=std::end(pol.edges_list);++edg)
    {
        using_points.insert((*edg).coord[0]);
        using_points.insert((*edg).coord[1]);
    }
    for(int i=0;i<pol.points_list.size();i++)
    {
        if(using_points.find(i)!=std::end(using_points))
        {
            new_points_num.push_back(new_points_list.size());
            new_points_list.push_back(pol.points_list[i]);
        }
        else
        {
            new_points_num.push_back(-1);
        }
    }
    pol.points_list=new_points_list;
    for(int i=0;i<pol.edges_list.size();i++)
    {
        pol.edges_list[i].coord[0]=new_points_num[pol.edges_list[i].coord[0]];
        pol.edges_list[i].coord[1]=new_points_num[pol.edges_list[i].coord[1]];
    }
    return pol;
    //pol.print();
}
plane get_plane_by_three_points(double x1,double y1,double z1,double x2,double y2,double z2,double x3, double y3,double z3)
{
    double A,B,C,D,norm;
    A=(y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
    B=(z2-z1)*(x3-x1)-(x2-x1)*(z3-z1);
    C=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
    D=-(x1*A+y1*B+z1*C);
    norm=sqrt(A*A+B*B+C*C);
    if(0==norm)
    {
        norm=1.;
    }
    A=A/norm;
    B=B/norm;
    C=C/norm;
    D=D/norm;
    plane res(A,B,C,D);
    return res;
}
plane polyhedron::get_plane_by_two_edges(edge ed1, edge ed2)
{
    plane res(points_list[ed1.coord[0]],points_list[ed1.coord[1]],points_list[ed2.coord[0]]);
    if(fabs(res.A)+fabs(res.B)+fabs(res.C)<1.e-9)
    {
        res=plane(points_list[ed1.coord[0]],points_list[ed1.coord[1]],points_list[ed2.coord[1]]);
    }
    return res;
}
int polyhedron::get_edge_num_by_two_facets(facet f1,facet f2)
{
    for(auto i=std::begin(f1.edges);i!=std::end(f1.edges);++i)
    {
        for(auto j=std::begin(f2.edges);j!=std::end(f2.edges);++j)
        {
            if((*i)==(*j))
            {
                return (*i);
            }
        }
    }
    return -1;
}
polyhedron construct_polyhedron_by_planes_list(std::vector<plane>* planes)
{
    double x=100.;
    polyhedron my_pol(x);
    int exit_key=1;
    std::vector<point> old_points_list;
    while(1==exit_key)
    {
        my_pol=polyhedron(x);
        old_points_list=my_pol.points_list;
        for(auto space=std::begin((*planes));space<std::end((*planes));++space)
        {
            //printf("%e %e %e %e\n",(*space).A,(*space).B,(*space).C,(*space).D);
            my_pol=cross(my_pol,(*space));
            printf("%d %d %d\n",my_pol.points_list.size(),my_pol.edges_list.size(),my_pol.facets_list.size());
            //break;
        }
        exit_key=0;
        for(auto old_point=std::begin(old_points_list);old_point!=std::end(old_points_list);++old_point)
        {

            for(auto new_point=std::begin(my_pol.points_list);new_point!=std::end(my_pol.points_list);++new_point)
            {
                if(fabs((*old_point).x-(*new_point).x)<1.e-12 &&
                   fabs((*old_point).y-(*new_point).y)<1.e-12 &&
                   fabs((*old_point).z-(*new_point).z)<1.e-12)
                {
                    exit_key=1;
                    break;
                }
            }
            if(1==exit_key)
            {
                break;
            }
        }
        x*=2.;
        //break;
    }
    return my_pol;
}
/*int main()
{
    double x=100.;

    polyhedron cub(x);
    plane pl(0.,0.,-1.,0.);
    printf("%d %d %d\n",(int)cub.facets_list.size(),(int)cub.edges_list.size(),(int)cub.points_list.size());
    cub=cross(cub,pl);
    printf("%d %d %d\n",(int)cub.facets_list.size(),(int)cub.edges_list.size(),(int)cub.points_list.size());
    //cub=cross(cub,plane(1,1,1,0));
    //printf("%d %d %d\n",(int)cub.facets_list.size(),(int)cub.edges_list.size(),(int)cub.points_list.size());
    //cub=cross(cub,plane(-1,1,1,0));
    //printf("%d %d %d\n",(int)cub.facets_list.size(),(int)cub.edges_list.size(),(int)cub.points_list.size());
    cub.print();

    return 0;
}*/
