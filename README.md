## Compile the code:
    make

## Dataset:
    For generating the binary of datasets, we use the converter from:
        https://github.com/huyang1988/TC/blob/master/README.md
    Some example datasets are provided in the data folder of this repository. For using different repository, provide the path to the dataset.
    i.e. for p2p08 dataset, "data/p2p08/input"

## Multi-GPU:
    We use jsrun to run multi-gpu version of the code on Summit.
    https://www.olcf.ornl.gov/for-users/system-user-guides/summit/summit-user-guide/#running-jobs


## Input arguments: 
1. Folder name containing the binary file of dataset. 
2. Total number of process 
3. Number of threads 
4. Number of Blocks 
5. Number of Buckets for hashing (Limit: 256) 
6. Block-based (1) or Warp-based(0) 
7. Degree based Balance (1) or Index based partition(0)

## Output Format: 
1. Name of the dataset 
2. Edge count 
3. Triangles count 
4. Max Time 
5. Min Time 
6. Total number of process 
7. Number of threads 
8. Number of block 
9. Number of buckets used 
10. TEPS rate based upon max time
