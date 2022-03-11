#pragma once
#include "math.h"
#include <float.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include "physconsts.h"

using namespace std;

inline double planck_function(const double & T, const double & nu)
{
    constexpr double inv_squared_c = 2. / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    const double exp_factor = (PLANK_CONSTANT * nu) / (BOLTZMANN_CONSTANT * T);
    if (fabs(exp_factor) > 1.e-4) return PLANK_CONSTANT * (nu*nu*nu*inv_squared_c) / ( exp(exp_factor) - 1.0 );
    else return (nu*nu*inv_squared_c) * BOLTZMANN_CONSTANT * T;
}

template <typename T>
inline T oneMinusExp(const T & tau) // = 1 - exp(-tau)
{
    if (tau < 1.e-5) return tau - tau * tau * 0.5;
    else return 1 - exp(-tau);
}

template <typename T, typename Tfin>
inline T readline(Tfin & fin)       // auxiliary function to read lines with empty spaces after symbols and before \n
{
    T a;
    string str;
    istringstream sfin;
    getline(fin, str);
    sfin.str(str);
    sfin >> a;
    sfin.clear();
    return a;
}

inline string trim(string & str)    // remove left and right spaces in a string
{
    str.erase(0, str.find_first_not_of(' '));       // prefixing spaces
    str.erase(str.find_last_not_of(' ') + 1);       // surfixing spaces
    return str;
}


// see https://gist.github.com/HViktorTsoi/58eabb4f7c5a303ced400bcfa816f6f5
template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] > array[right];
              });

    return indices;
}
