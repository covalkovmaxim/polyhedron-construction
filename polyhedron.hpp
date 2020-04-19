#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <set>
#include <functional>
class plane;

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
    point operator+(const point& right_point) const;
    point operator-(const point& right_point) const;
    double operator*(const point& right_point) const;
    point operator*(const double& right_number) const;
    point operator/(const double& right_number) const;
    bool operator==(const point& right_point) const;

};
point operator*(const double& left_number, const point& right_point);
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
    bool operator==(const facet& right_facet) const;
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
    polyhedron(double x);
    polyhedron(char*name);
    plane get_plane_by_two_edges(edge ed1,edge ed2);
    int get_edge_num_by_two_facets(facet f1,facet f2);

    void print();

};
plane get_plane_by_three_points(double x1,double y1,double z1,double x2,double y2,double z2,double x3, double y3,double z3);

class plane
{
public:
    double A,B,C,D;

    plane(double a,double b,double c,double d)
    {
        A=a; B=b; C=c;D=d;
    }
    plane(point p1,point p2,point p3)
    {
        plane tec=get_plane_by_three_points(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,p3.x,p3.y,p3.z);
        A=tec.A;B=tec.B;C=tec.C;D=tec.D;
    }

};
polyhedron cross(polyhedron&& pol, plane&& space);
polyhedron construct_polyhedron_by_planes_list(std::vector<plane>* planes);
