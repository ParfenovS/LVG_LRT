# LVG_LRT
The code for radiative transfer calculations in LVG (large velocity gradient) approximation and analysis of the maser pumping cycles.
The details of the radiative transfer modelling can be found in [Sobolev & Deguchi (1994)](https://ui.adsabs.harvard.edu/abs/1994A%26A...291..569S/abstract), [Sobolev et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26A...324..211S/abstract), and [Cragg et al. (2002)](https://ui.adsabs.harvard.edu/abs/2002MNRAS.331..521C/abstract). The procedure of the maser pumping analysis is described by [Sobolev (1986)](https://ui.adsabs.harvard.edu/abs/1986AZh....63..674S/abstract) and [Sobolev & Deguchi (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...433..719S/abstract).

The code can be compiled with a compiler supporting C++17 under Linux (see Makefile) and Windows (LVG_LRT.sln is a Visual Studio solution file).

To run the LVG code type in a terminal:
```
./LVG_LRT.exe
```
The code reads input model parameters from files in Parameters directory. Please check also src/hiddenParameters.h file containing additional parameters for calculations. The output will be stored in files in Results directory.

To run the code for the maser pumping analysis type in a terminal:

```
./pupmpit.exe 1 3
```
where the first and second command arguments denote the indexes of, respectively, starting and final energy levels between which one wants to search for pumping cycles (these are usually the lower and upper levels of a maser transition). The level indexes begin from 1 and correspond to level indexes given in the input file with spectroscopic data in [LAMDA](https://home.strw.leidenuniv.nl/~moldata/) format. The code will read the input parameters from Parameters directory and the output files produced by LVG_LRT.exe. The input parameters given in this repository are in the range considered by [van der Walt (2014)](https://ui.adsabs.harvard.edu/abs/2014A%26A...562A..68V/abstract) (see their Figure 6) so one can compare the results obtained with two different independent codes. The output will be into a terminal and will contain the information on found cycles. For each cycle in the output, there will be the number of levels in the cycle, the number of times the code found the cycle (the code search for cycles with Monte-Carlo method and can encounter the same cycle many times), the relative efficiency of the cycle, the list of level indexes within the cycle.

The repository also contains examples of python scripts which can be used to calculate grids of models for a range of input parameters (see Scripts directory). There is also the python script Scripts/plot_pops_vs_time.py to plot the dependence of level populations on time from the binary file that is produced by LVG_LRT.exe when it is configured to integrate kinetic equations for level populations (by default LVG_LRT.exe solves statistical equilibrium equations).

An example of results produced with the code can be found in [Chen et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...877...90C/abstract).

The source code from [LINPACK](https://people.sc.fsu.edu/~jburkardt/cpp_src/linpack/linpack.html) and  [https://github.com/ttk592/spline/](https://github.com/ttk592/spline/) is used.
