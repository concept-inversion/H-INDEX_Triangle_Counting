#include <iostream>
#include "graph.h"
#include "wtime.h"
#include <queue>
#include <set>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;


int binary_search(int start,int end,int value, int *arr)
{
    //printf("low:%d,high:%d,value:%f\n",start,end,value); 
    int low=start;
    int high=end;
    int index=start;
    while (low<=high)
    {
	index=((low+high)/2);
        if (value<arr[index])
		{
            //set high to index-1
            high= index-1;
	    //printf("high:%d\n",high);
        }
        else if (value>arr[index])
        {
            // set low to index+1
            low = index+1;
            //printf("low:%d\n",low);

	}
        else
        {
            break;
        } 
	}
	printf("Vlaue: %d,Found: %d\n",value,arr[index]);
    return index;
}


int start(int argc, char** argv)
{
	int total_process= atoi(argv[3]);
	int rank= atoi(argv[2]);
	string json_file 	= argv[1];
	 graph *graph_d 
	 	= new graph	(json_file); 
	//int N_GPUS=argv[1];
	int deviceCount=1;
	int vertex_count=	graph_d-> vert_count;
	int edge_count= graph_d-> edge_count;
	int edge_list_count= graph_d-> edge_list_count;
	int edges= edge_list_count/2;
	int ratio=2*(edges/total_process);
	int E_START=rank*ratio;
	//int E_END=E_START+ratio;
	//if(rank==(total_process-1)){E_END=edge_list_count;}
	int *counter=(int *)malloc(sizeof(int)*edges);
	//printf("Rank: %d, Devicecount: %d,  Start: %d, End: %d, Selected: %d\n",rank,deviceCount,E_START,E_END,(rank%deviceCount));
	double time_start=wtime();
	
	int *prefix=(int *)malloc(sizeof(int)*edges);
	int temp;
	for(int i=0;i<edge_list_count;i+=2)
	 {
	 	int vertex = graph_d->edge_list[i];
	 	int vertex1 = graph_d->edge_list[i+1];
	 	int N1_start=graph_d->beg_pos[vertex];
	 	int N1_end= graph_d->beg_pos[vertex+1];
	 	int L1= N1_end-N1_start;
	 	int N2_start= graph_d->beg_pos[vertex1];
	 	int N2_end= graph_d->beg_pos[vertex1+1];	
	 	int L2= N2_end-N2_start;
	 	int sum=L1+L2;
		temp =  sum +prefix[(i/2)-1];
		prefix[i/2]= temp;
	 	//printf("vertexA: %d, D1: %d, vertexB: %d, D2: %d, Degree: %d, prefix: %d\n",vertex,L1,vertex1,L2,L1+L2,temp);
	 }
	 int total_degree= temp;
	 int SIZE = (total_degree/total_process)*(rank+1);
	 int E_END= binary_search(0,edges,SIZE,prefix);
	 printf("total degree: %d,total edges: %d, E_END: %d, size: %d\n",temp,edges,E_END,SIZE);

}    


