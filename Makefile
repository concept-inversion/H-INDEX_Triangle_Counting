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
	mpirun -n $(N) cuda-memcheck $(exe) data/p2p/input $(N) 16 1 16 0 0 
#	mpirun -n $(N) $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/email-EuAll_adj/ $(N)


clean:
	rm -rf *.o ${exe}
