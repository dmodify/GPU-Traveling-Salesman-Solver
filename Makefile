NVCC=/usr/local/cuda-6.5/bin/nvcc
NVCCFLAGS= --use_fast_math -O3 -Xptxas -dlcm=cg  --relocatable-device-code=true
CUFILES= two_opt.cu util.cu
NVCCARCH= -arch=sm_30
all:
	$(NVCC) $(NVCCFLAGS) $(CUFILES) -o two_opt $(NVCCARCH) -lm
clean:
	rm -f two_opt *~ *.core out
