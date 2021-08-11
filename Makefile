LVG_LRT.exe :
	c++ -std=c++17 -O3 -mtune=native -mavx2 -Wall -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp
	c++ -std=c++17 -O3 -mtune=native -mavx2 -Wall -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icc -std=c++17 -O3 -Wall -xCORE-AVX2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp
#	icc -std=c++17 -O3 -Wall -xCORE-AVX2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
#
#	icpx -std=c++17 -O3 -Wall -xCORE-AVX2 -o LVG_LRT.exe src/LVG_LRT.cpp src/stdafx.cpp
#	icpx -std=c++17 -O3 -Wall -xCORE-AVX2 -o pupmpit.exe src/pupmpit.cpp src/stdafx.cpp
