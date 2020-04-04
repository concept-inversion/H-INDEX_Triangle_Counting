#include <iostream>
#include "graph.h"
#include "wtime.h"
#include <queue>
#include <set>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include "herror.h"
#include <math.h>
#include "TC.cuh"
#include <assert.h>

int BUCKET_SIZE=1000;
using namespace std;

__device__ 
void d_display(int *a, int column,int row,int start)
{
		printf("\n");
		for(int i=0; i<row; i++)
		{
			for(int j=0;j<column;j++)
			{
				printf("%d\t",a[i*column+j+start]);
			}
			printf("\n");

		}		  
}

__device__
void kogge_sum(int *A,int len, int start)
{
    /* We require enough threads for this method */
    int step=log2f(len);
    //printf("Len: %d, Steps: %d, start: %d\n",len,step,start);
    int pos,i;
    for(i=0;i<step;i++)
    {
        pos=powf(2,i);
        int j=pos+threadIdx.x+start;
        while(j<(len+start))
        {
            int temp=A[j-pos];
            __syncthreads;
            A[j]+=temp;
            //printf("Write:%d , Read:%d , Written: %d\n",j,j-pos,A[j]);
            j+=blockDim.x;
        }
        //if(threadIdx.x==0){printf("\n\n");}
        __syncthreads();
	}
}

__device__
int linear_search(int neighbor,int *partition1, int *bin_count, int bin, int BIN_OFFSET, int BIN_START,int BUCKETS)
{
	int len = bin_count[bin+BIN_OFFSET];
	//printf("\nPartStart: %d\n",BIN_START);
	int i = bin + BIN_START;
	int step=0;
	while(step<len)
	{
		int test=partition1[i];
		//printf("Neighbor: %d, Test: %d\n",neighbor,test);
		if(test==neighbor)
		{
			return 1;
		}
		else
		{
			i+=BUCKETS;
		}
		step+=1;
	}
	return 0;
}

__device__
int merge(int *A, int *B, int ai, int bi, int l1_e, int l2_e,int steps)
{
	/*Reminder: As the partition is coalesced, accessing next element in each partition would require: next element --> prev + Warpsize */
	int WARPSIZE=64;
	int count = 0;
	int steps_count=0;
	while((ai<=l1_e) && (bi<=l2_e))
	{
		steps_count+=1;
		//printf("\nAI: %d, value: %d \t",ai, A[ai]);
		//printf("BI: %d, value: %d \t",bi, B[bi]);
		if(A[ai]>B[bi])
		{
			bi+=WARPSIZE;	
		}
		else if(A[ai]<B[bi])
		{
			ai+=WARPSIZE;
		}
		else
		{
			count+=1;
			ai+=WARPSIZE;
			bi+=WARPSIZE;
		}
		//printf("\n");
		__syncthreads();
	}
	//printf("Thread: %d, count: %d \n",threadIdx.x,count);
	return count;
}

int binary_search(long long start,long long end,int value, long long *arr)
{
    //printf("low:%d,high:%d,value:%f\n",start,end,value); 
    long long low=start;
    long long high=end;
    long long index=start;
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
	//printf("Vaue: %d,Found: %d\n",value,arr[index]);
    return index;
}

__device__
int max_count(int *bin_count,int start,int end,int len)
{
	int max_count =bin_count[start];
	int min_count=bin_count[start];
	int zero_count=0;
	for (int i=start;i<end;i++)
	{
		if(bin_count[i]>max_count)
		{
			max_count=bin_count[i];			
		}
		if(bin_count[i]<min_count)
		{
			min_count=bin_count[i];
		}
		if(bin_count[i]==0)
		{
			zero_count+=1;
		}
	}
	// printf("%d,%d,%d\n",zero_count,max_count,len);
	return max_count;
}

void graph_reordering(graph *graph_temp)
{

}

__global__ void
warp_hash_count(vertex_t* adj_list, index_t* beg_pos, vertex_t* edge_list, int edge_count, int vertex_count,int edge_list_count, int *partition,unsigned long long *GLOBAL_COUNT,long long E_START, long long E_END, int device, int BUCKETS, int G_BUCKET_SIZE, int T_Group)
{
	// Uncomment the lines below and change partition to Gpartition for using shared version
	int *part;
	int S_BUCKET_SIZE=320;
	int tid=threadIdx.x+blockIdx.x*blockDim.x;
	int WARPSIZE=T_Group;
	int __shared__ bin_count[256*4];
	//int __shared__ partition[160*4];
	int PER_BLOCK_WARP=blockDim.x/WARPSIZE;
	int G_WARPID= tid/WARPSIZE;
	int WARPID = threadIdx.x/WARPSIZE;
	int __shared__ G_counter;
	G_counter=0;
	int P_counter=0;
	int BINsize = BUCKETS*G_BUCKET_SIZE;
	//int BINsize = BUCKETS*5;
	int BIN_START = G_WARPID*BINsize;
	//int BIN_START = WARPID*BINsize;
	long long i=G_WARPID*2;
	long long RANGE= E_END-E_START;
	int BIN_OFFSET= WARPID*BUCKETS;
	//for(int i=0;i<edge_list_count; i+=2)
	//TODO: Static assignment to dynamic assignment of edges 
	// unsigned long long TT=0,HT=0,IT=0;
	// unsigned long long __shared__ G_TT,G_HT,G_IT;
	// G_TT=0,G_HT=0,G_IT=0;
	while(i<( RANGE))
	{
		//if(threadIdx.x%32==0){printf("Warp:%d, G_WArp: %d,i: %d \n",WARPID,G_WARPID,i);}
		//if (device==1){printf("Device: %d, i: %d\n",device,i);}
		/* TODO: Divide edge list to multiple blocks*/
		// unsigned long long start_time=clock64();
		int destination = edge_list[i];
		int source = edge_list[i+1];
		int N1_start=beg_pos[destination];
		int N1_end= beg_pos[destination+1];	
		int L1= N1_end-N1_start;
		int N2_start=beg_pos[source];
		int N2_end= beg_pos[source+1];	
		int L2= N2_end-N2_start;
		
		// if ((L1==0))
		// {
		// 	//printf("continue %d\n",i);
		// 	continue;
		// }
		// // N2 is for hashing and N1 is lookup
		if(L1>L2)
		{
			int temp= N1_start;
			N1_start= N2_start;
			N2_start=temp;
			temp=N1_end;
			N1_end=N2_end;
			N2_end=temp;
			temp=L2;
			L2=L1;
			L1=temp;
		}
		
		// unsigned long long hash_start=clock64();
		int id=threadIdx.x%WARPSIZE+BIN_OFFSET;
		int end = BIN_OFFSET+BUCKETS;
		//if(threadIdx.x%32==0){printf("End: %d\n",end);}
		// We can remove this line
		
		__syncwarp();
		while(id<(end))
		{
			bin_count[id]=0;
			//printf("BIN: %d\n",id);
			id+=WARPSIZE;
		}
		int start=threadIdx.x%WARPSIZE + N2_start;
		// BIN_OFFSET is for count of number of element of each bin for all 4 warps
		
		__syncwarp();
		// Hash one list 
		while(start<N2_end)
		{
			int temp= adj_list[start];
			int bin=temp%BUCKETS;
			int index=atomicAdd(&bin_count[bin+BIN_OFFSET],1);
			partition[index*BUCKETS+ bin + BIN_START]=temp;
			//{printf("thread: %d,warp:%d, write: %d bin %d, index %d  at: %d\n",threadIdx.x,WARPID,temp,bin,index,(index*WARPSIZE+bin+BIN_START));}	
			start+=WARPSIZE;
		}
		__syncwarp();
		// unsigned long long hash_time=clock64()-hash_start;
		//int max_len_collision= max_count(bin_count,BIN_OFFSET,BIN_OFFSET+BUCKETS,L2);
		
		// unsigned long long intersection_start=clock64();
		start=threadIdx.x%WARPSIZE + N1_start;
		int count;
		//if(threadIdx.x==32){printf("start: %d, BIN_OFFSET: %d\n",start,BIN_OFFSET);}
		//P_counter=0;
		while(start<N1_end)
		{
			count=0;
			int neighbor=adj_list[start];
			int bin=neighbor%BUCKETS;
			count=linear_search(neighbor,partition,bin_count,bin,BIN_OFFSET,BIN_START,BUCKETS);
			P_counter+=count;
			start+=WARPSIZE;
			//printf("Tid: %d, Search:%d\n",threadIdx.x,neighbor);
		}
		//atomicAdd(&GLOBAL_COUNT[0],P_counter);
		
		__syncwarp();
		// unsigned long long intersection_time=clock64()-intersection_start;
		// if(threadIdx.x%32==0){printf("I: %d, Start:%d, End:%d, Count:%d\n",i,vertex,vertex1,G_counter);}
		i+=gridDim.x*PER_BLOCK_WARP*2;
		// unsigned long long total_time=clock64()-start_time;
		// if(threadIdx.x%32==0){
		// 	// printf("%d %d %d\n",total_time, hash_time, intersection_time);
		// 	TT+=total_time;
		// 	HT+=hash_time;
		// 	IT+=intersection_time;
		// }
	}
	atomicAdd(&G_counter,P_counter);
	// atomicAdd(&G_HT,HT);
	// atomicAdd(&G_TT,TT);
	// atomicAdd(&G_IT,IT);
	__syncthreads();
	if(threadIdx.x==0)
	{
		// printf("%d\n",G_TT);
		atomicAdd(&GLOBAL_COUNT[0],G_counter);
		// atomicAdd(&GLOBAL_COUNT[1],G_TT);
		// atomicAdd(&GLOBAL_COUNT[2],G_HT);
		// atomicAdd(&GLOBAL_COUNT[3],G_IT);
	}
	
	//if(threadIdx.x==0){printf("Device: %d, Count:%d\n",device,GLOBAL_COUNT[0]);}

}

__global__ void
CTA_hash_count(vertex_t* adj_list, index_t* beg_pos, vertex_t* edge_list, int edge_count, int vertex_count,int edge_list_count, int *partition,unsigned long long *GLOBAL_COUNT,int E_START, int E_END, int device, int BUCKETS, int BUCKET_SIZE,int T_Group)
{
	int tid=threadIdx.x+blockIdx.x*blockDim.x;
	int WARPSIZE=128;
	int __shared__ bin_count[512];
	int G_WARPID= blockIdx.x;
	int WARPID = blockIdx.x;
	int __shared__ G_counter;
	G_counter=0;
	int P_counter=0;
	int BINsize = BUCKETS*BUCKET_SIZE;
	int i=G_WARPID*2;
	int RANGE= E_END-E_START;
	int BIN_START = G_WARPID*BINsize;
	int divid=vertex_count/BUCKETS;
	int max_len_collision;
	//for(int i=0;i<edge_list_count; i+=2)
	//TODO: Static assignment to dynamic assignment of edges 
	
	while(i<( RANGE))
	{

		/* TODO: Divide edge list to multiple blocks*/
		int destination = edge_list[i];
		int source = edge_list[i+1];
		int N1_start=beg_pos[destination];
		int N1_end= beg_pos[destination+1];	
		int L1= N1_end-N1_start;
		int N2_start=beg_pos[source];
		int N2_end= beg_pos[source+1];	
		int L2= N2_end-N2_start;

		// N2 is for hashing and N1 is lookup
		if(L1>L2)
		{
			int temp= N1_start;
			N1_start= N2_start;
			N2_start=temp;
			temp=N1_end;
			N1_end=N2_end;
			N2_end=temp;
			temp=L2;
			L2=L1;
			L1=temp;
		}
		

		int id=threadIdx.x;
		int end = BUCKETS;

		while(id<(end))
		{
			bin_count[id]=0;
			id+=blockDim.x;
		}
		__syncthreads();
		int start=threadIdx.x + N2_start;
		
		// Hash one list 
		while(start<N2_end)
		{
			int temp= adj_list[start];
			int bin=temp%BUCKETS;
			int index=atomicAdd(&bin_count[bin],1);
			partition[index*BUCKETS+ bin + BIN_START]=temp;
			//{printf("thread: %d,warp:%d, write: %d bin %d, index %d  at: %d\n",threadIdx.x,WARPID,temp,bin,index,(index*WARPSIZE+bin+BIN_START));}	
			start+=blockDim.x;
		}
		__syncthreads();
		if (threadIdx.x==0)
		{
			max_len_collision= max_count(bin_count,0,0+BUCKETS,L2);
			// printf("max_len_collision: %d\n",max_len_collision );
		}
		__syncthreads();
		start=threadIdx.x + N1_start;
		int count;
		//if(threadIdx.x==32){printf("start: %d, BIN_OFFSET: %d\n",start,BIN_OFFSET);}
		//P_counter=0;
		while(start<N1_end)
		{
			count=0;
			int neighbor=adj_list[start];
			int bin=neighbor%BUCKETS;
			count=linear_search(neighbor,partition,bin_count,bin,0,BIN_START,BUCKETS);
			P_counter+=count;
			start+=blockDim.x;
			//printf("Tid: %d, Search:%d\n",threadIdx.x,neighbor);
		}
		//atomicAdd(&GLOBAL_COUNT[0],P_counter);
		
		//if(threadIdx.x%32==0){printf("I: %d, Start:%d, End:%d, Count:%d\n",i,vertex,vertex1,G_counter);}
		i+=gridDim.x*2;
		
	}
	atomicAdd(&G_counter,P_counter);
	__syncthreads();
	if(threadIdx.x==0){atomicAdd(&GLOBAL_COUNT[0],max_len_collision);}
	
	//if(threadIdx.x==0){printf("Device: %d, Count:%d\n",device,GLOBAL_COUNT[0]);}
}

struct arguments Triangle_count(int rank, char name[100], struct arguments args, int total_process,int n_threads , int n_blocks, int BUCKETS, int select_thread_group, int select_partition)
{

	// printf("---------------Here----------------");
	int T_Group= 32;
	int PER_BLOCK_WARP= n_threads/T_Group;
	int total=n_blocks*PER_BLOCK_WARP*BUCKETS*BUCKET_SIZE;
    unsigned long long *counter=(unsigned long long *)malloc(sizeof(unsigned long long)*10);
	string json_file 	= name;
	 graph *graph_d 
		 = new graph	(json_file);
	 graph *graph_b 
	     = new graph	(json_file);

	//printf("Graph Adj Read: %d",graph_d->adj_list[10]); 
	//int N_GPUS=argv[1];
	int deviceCount;
	HRR(cudaGetDeviceCount(&deviceCount));
	// cout<<deviceCount<<endl;
	//fprintf(stderr,"----------------Device count: %d\n",deviceCount);
	//cudaSetDevice();
	HRR(cudaSetDevice(rank%deviceCount));
	//cudaDeviceProp devProp;
	//HRR(cudaGetDeviceProperties(&devProp, rank));
	index_t vertex_count=	graph_d-> vert_count;
	index_t edge_count= graph_d-> edge_count;
	index_t edge_list_count= graph_d-> edge_list_count;
	index_t edges= edge_list_count>>1;
	/* Preprocessing Step to calculate the ratio */
	long long *prefix=(long long *)malloc(sizeof(long long)*edges);
	long long temp;
	for(int i=0;i<edge_list_count;i+=2)
	 {
	 	int destination = graph_d->edge_list[i];
	 	int source = graph_d->edge_list[i+1];
	 	int N1_start=graph_d->beg_pos[destination];
	 	int N1_end= graph_d->beg_pos[destination+1];
	 	int L1= N1_end-N1_start;
	 	int N2_start= graph_d->beg_pos[source];
	 	int N2_end= graph_d->beg_pos[source+1];	
	 	int L2= N2_end-N2_start;
	 	int sum=L1+L2;
		if(i==0)
		{
			temp=0;
		}
		else
		{
			temp =  sum +prefix[(i>>1)-1];
		}
		prefix[i>>1]= temp;
	 	//printf("vertexA: %d, D1: %d, vertexB: %d, D2: %d, Degree: %d, prefix: %d\n",vertex,L1,vertex1,L2,L1+L2,temp);
	 }
	//  cout<<"edge_list_count OK"<<endl;
	int total_degree= temp;
	 //printf("total degree: %d,total edges: %d, E_END: %d,E_start:%d, size: %d, rank: %d\n",temp,edges,E_END,E_START,SIZE,rank);
	int SIZE,ratio;
	long long E_END,E_START;
	//-------------------------------------------//
	if(select_partition==1)
	{
		SIZE = (total_degree/total_process);
		E_END= binary_search(0,edge_count,SIZE*(rank+1),prefix);
		E_START= binary_search(0,edge_count,SIZE*rank,prefix);
		E_END=E_END<<1;
		E_START=E_START<<1;
	}
	//--------------------------------------------//
	
	else
	{
		ratio=2*(edges/total_process);
		E_START=rank*ratio;
		E_END=E_START+ratio;
		SIZE= prefix[E_END/2]-prefix[E_START/2];
	}
	//--------------------------------------------//
	// E_START=0;E_END=edge_list_count;
	// cout<<"partition OK!"<<endl;
	// cout<<E_START<<' '<<E_END<<endl;
	assert(E_END>E_START);
	//fprintf(stderr,"Rank: %d, Devicecount: %d,  Start: %d, End: %d, Selected: %d\n",rank,deviceCount,E_START,E_END,(rank%deviceCount));
	if(rank==(total_process-1)){E_END=edge_list_count;}
	int *hash,* BIN_MEM;
	unsigned long long *GLOBAL_COUNT;
	index_t *d_beg_pos;
	vertex_t *d_adj_list,*d_edge_list;
	float memory_req = (sizeof(int)*total + sizeof(index_t)*(vertex_count+1)+ sizeof(vertex_t)*(edge_count)+sizeof(vertex_t)*(E_END-E_START+1))/(1024*1024);
	HRR(cudaMalloc((void **) &GLOBAL_COUNT,sizeof(unsigned long long)*10));
	HRR(cudaMalloc((void **) &BIN_MEM,sizeof(int)*total));
	HRR(cudaMalloc((void **) &d_beg_pos,sizeof(index_t)*(vertex_count+1)));
	HRR(cudaMalloc((void **) &d_adj_list,sizeof(vertex_t)*(edge_count)));
	HRR(cudaMalloc((void **) &d_edge_list,sizeof(vertex_t)*(E_END-E_START+1))); // Swap edge list count with Eend - estart; --> gives error; may add some more
	
	
	HRR(cudaMemcpy(d_edge_list,graph_d->edge_list+E_START,sizeof(vertex_t)*(E_END-E_START+1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_beg_pos,graph_d->beg_pos,sizeof(index_t)*(vertex_count+1), cudaMemcpyHostToDevice));
	HRR(cudaMemcpy(d_adj_list,graph_d->adj_list,sizeof(vertex_t)*edge_count, cudaMemcpyHostToDevice));

	double t1=wtime();
	double cmp_time;
	if(select_thread_group==1)
	{
		double time_start=wtime();
		CTA_hash_count<<<n_blocks,n_threads>>>(d_adj_list, d_beg_pos, d_edge_list, edge_count, vertex_count,edge_list_count, BIN_MEM,GLOBAL_COUNT,E_START,E_END,rank,BUCKETS,BUCKET_SIZE, T_Group);
		HRR(cudaDeviceSynchronize());
    	cmp_time = wtime()-time_start;	
	}
	else
	{
		double time_start=wtime();
		warp_hash_count<<<n_blocks,n_threads>>>(d_adj_list, d_beg_pos, d_edge_list, edge_count, vertex_count,edge_list_count, BIN_MEM,GLOBAL_COUNT,E_START,E_END,rank,BUCKETS,BUCKET_SIZE, T_Group);
		HRR(cudaDeviceSynchronize());
    	cmp_time = wtime()-time_start;
	}
	HRR(cudaMemcpy(counter,GLOBAL_COUNT,sizeof(unsigned long long)*10, cudaMemcpyDeviceToHost));
	//printf("Edges: %d,Start: %d, End: %d, Rank: %d,ratio:%d, Triangle: %d, Time: %f\n",E_END-E_START,E_START,E_END,rank,SIZE,counter[0],cmp_time);
	args.time=cmp_time;
	args.count=counter[0];
	args.edge_count=edges;
	args.degree= SIZE;
	args.vertices= vertex_count-1;
	return args;
}    
