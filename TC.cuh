
struct arguments{
    int edge_count;
    int count;
    double time;
    int degree;
    int vertices;
};

struct arguments Triangle_count(int rank, char input[100], struct arguments args,int total_process, int threads, int blocks, int buckets, int SELECT_thread_group, int SELECT_partition);