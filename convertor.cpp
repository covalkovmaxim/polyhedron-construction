#include "polyhedron.hpp"


int main()
{
    FILE*in=fopen("big_initpoly.txt","r");
    FILE*out=fopen("new_big_initpoly.txt","w");
    int *index;
    std::pair<int,int>*edges;
    int n1,n2,n3,sides_num,tec_num,p1,p2,oth1,oth2;
    double x,y,z,other;
    fscanf(in,"%d%d%d",&n1,&n2,&n3);
    printf("%d %d %d\n",n1,n2,n3);
    fprintf(out,"%d\n",n1);
    for(int i=0;i<n1;i++)
    {
        fscanf(in,"%d%lf%lf%lf",&p1,&x,&y,&z);
        fprintf(out,"%f %f %f\n",x,y,z);
        printf("%f %f %f\n",x,y,z);
    }
    index=new int[n2];
    edges=new std::pair<int,int>[n2];
    std::vector<std::pair<int,int>> coord;
    int cor_num=0,key=0;
    for(int i=0;i<n2;i++)
    {
        fscanf(in,"%d%d%d%d",&p1,&p2,&oth1,&oth1);
        edges[i].first=p1;
        edges[i].second=p2;
        key=0;
        for(int j=0;j<i-1;j++)
        {
            if(edges[j].first==edges[i].second && edges[j].second==edges[i].first)
            {
                key=1;
                index[i]=index[j];
                break;
            }
        }
        if(0==key)
        {
            index[i]=cor_num;
            coord.push_back(edges[i]);
            cor_num++;
        }
    }
    fprintf(out,"%d\n",cor_num);
    for(int i=0;i<cor_num;i++)
    {
        fprintf(out,"%d %d\n",coord[i].first,coord[i].second);
    }
    fprintf(out,"%d\n",n3);
    for(int i=0;i<n3;i++)
    {
        fscanf(in,"%d%d%lf%lf%lf%lf",&oth1,&sides_num,&other,&other,&other,&other);
        fprintf(out,"%d ",sides_num);
        for(int j=0;j<sides_num;j++)
        {
            fscanf(in,"%d",&tec_num);
            fprintf(out,"%d ",index[tec_num]);
        }
        fprintf(out,"\n");
    }
    fclose(in);
    fclose(out);
    delete[]index;
    delete[]edges;
    return 0;
}
