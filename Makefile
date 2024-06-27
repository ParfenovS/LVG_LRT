LVG_LRT.exe :
	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -llapacke -llapack -lblas -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -llapacke -llapack -lblas -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	clang++ -std=c++17 -Wall -O3 -mavx2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	clang++ -std=c++17 -Wall -O3 -mavx2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
# AMD AOCC, AOCL
#
#	clang++ -std=c++17 -Wall -O3 -march=native -mavx2 -flto -zopt -fremap-arrays -mllvm -reduce-array-computations=3 -fnt-store -DUSE_LAPACKE -lflame -lblis -laoclutils -lalm -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	clang++ -std=c++17 -Wall -O3 -march=native -mavx2 -flto -zopt -fremap-arrays -mllvm -reduce-array-computations=3 -fnt-store -DUSE_LAPACKE -lflame -lblis -laoclutils -lalm -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -lflame -lblis -laoclutils -lalm -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -lflame -lblis -laoclutils -lalm -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
# INTEL MKL, compilers, see https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html
#
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -DMKL_ILP64 -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -DMKL_ILP64 -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -DMKL_ILP64 -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -DUSE_LAPACKE -DMKL_ILP64 -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icpc -std=c++17 -Wall -O3 -xCORE-AVX2 -DUSE_LAPACKE -DMKL_ILP64 -I"${MKLROOT}/include" -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	icpc -std=c++17 -Wall -O3 -xCORE-AVX2 -DUSE_LAPACKE -DMKL_ILP64 -I"${MKLROOT}/include" -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icpx -std=c++17 -Wall -O3 -xCORE-AVX2 -flto -DUSE_LAPACKE -DMKL_ILP64 -I"${MKLROOT}/include" -L${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	icpx -std=c++17 -Wall -O3 -xCORE-AVX2 -flto -DUSE_LAPACKE -DMKL_ILP64 -I"${MKLROOT}/include" -L${MKLROOT}/lib -lmkl_rt -lpthread -lm -ldl -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icpx -std=c++17 -Wall -O3 -xCORE-AVX2 -flto -DUSE_LAPACKE -DMKL_ILP64 -qmkl-ilp64=sequential -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	icpx -std=c++17 -Wall -O3 -xCORE-AVX2 -flto -DUSE_LAPACKE -DMKL_ILP64 -qmkl-ilp64=sequential -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
