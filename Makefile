LVG_LRT.exe :
	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
	c++ -std=c++17 -Wall -O3 -mtune=native -mavx2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icc -std=c++17 -Wall -O3 -xCORE-AVX2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	icc -std=c++17 -Wall -O3 -xCORE-AVX2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icpx -std=c++17 -Wall -O3 -xCORE-AVX2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	icpx -std=c++17 -Wall -O3 -xCORE-AVX2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	clang++ -std=c++17 -Wall -O3 -mavx2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	clang++ -std=c++17 -Wall -O3 -mavx2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
# AMD AOCC
#	clang++ -std=c++17 -Wall -O3 -march=native -mavx2 -flto -zopt -fremap-arrays -mllvm -reduce-array-computations=3 -fnt-store -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp src/cubature/hcubature.c
#	clang++ -std=c++17 -Wall -O3 -march=native -mavx2 -flto -zopt -fremap-arrays -mllvm -reduce-array-computations=3 -fnt-store -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
