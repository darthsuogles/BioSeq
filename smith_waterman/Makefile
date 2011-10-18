
# use CUDA to compile the 
gpu_src=$(shell ls *.cu)
gpu_obj=$(gpu_src:.cu=.cu.o)

cpu_src=$(shell ls *.cpp)
cpu_obj=$(cpu_src:.cpp=.o)

all: test

$(gpu_obj): %.cu.o: %.cu
	@echo compiling gpu source code
	nvcc -o $@ -c $< $(CUDA_MAKEFILE_I)

$(cpu_obj): %.o: %.cpp
	@echo compiling cpu source code
	g++ -c $< $(CUDA_MAKEFILE_I)

test: $(gpu_obj) $(cpu_obj)
	@echo compiling the data file
	g++ -fPIC -o $@ $^ $(CUDA_MAKEFILE_L) -lcuda -lcudart -lshrutil_x86_64 -lcutil_x86_64

clean:
	rm test *.o