exe=trianglecounting.bin
N=1
cucc= "$(shell which nvcc)"
cc= "$(shell which mpicxx)"
commflags=-lcudart -L"$(shell dirname $(cucc))"/../lib64
cuflags= --compiler-options -v -Xcudafe  -\# --resource-usage 

objs	= $(patsubst %.cu,%.o,$(wildcard *.cu)) \
	$(patsubst %.cpp,%.o,$(wildcard *.cpp))

deps	= $(wildcard ./*.cuh) \
	$(wildcard ./*.hpp) \
	$(wildcard ./*.h) \


%.o:%.cu $(deps)
	$(cucc) -c $(cuflags) $< -o $@

%.o:%.cpp $(deps)
	$(cc) -c  $< -o $@

$(exe):$(objs)
	$(cc) $(objs) $(commflags) -o $(exe)


test:$(exe)
	@#mpirun -n $(N) cuda-memcheck $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/p2p-Gnutella08_adj/ $(N) 32 1 32 0 0
	mpirun -n $(N) $(exe) data/p2p08/input $(N) 32 1 32 0 0 
	@#mpirun -n $(N) $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/email-EuAll_adj/ $(N) 32 1 32 0 0


clean:
	rm -rf *.o ${exe}
