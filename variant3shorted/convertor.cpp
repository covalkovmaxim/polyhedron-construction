#include "polyhedron.hpp"


int main()
{
    FILE*in=fopen("big_initpoly_8.txt","r");
    FILE*out=fopen("new_big_initpoly_8.txt","w");
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
        //printf("%f %f %f\n",x,y,z);
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
        for(int j=0;j<i;j++)
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
        int*tec_arr=new int [sides_num];

        for(int j=0;j<sides_num;j++)
        {
            fscanf(in,"%d",&tec_num);
            tec_arr[j]=tec_num;
            //fprintf(out,"%d ",index[tec_num]);
        }
        int tec_elem,next_elem;
        for(int j=0;j<sides_num;j++)
        {
            tec_elem=tec_arr[j];
            next_elem=tec_arr[(j+1)%sides_num];
            for(int k=0;k<n2;k++)
            {
                if(tec_elem==edges[k].first && next_elem==edges[k].second)
                {
                    fprintf(out,"%d ",index[k]);
                    //printf("%d ",index[k]);
                    break;
                }
            }
        }
        delete []tec_arr;
        fprintf(out,"\n");
    }
    printf("edges_num=%d\n",cor_num);
    polyhedron my_pol("new_big_initpoly_2.txt");
    printf("%d %d %d\n",my_pol.points_list.size(),my_pol.edges_list.size(),my_pol.facets_list.size());
    my_pol.print();
    fclose(in);
    fclose(out);
    delete[]index;
    delete[]edges;
    return 0;
}
