#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

// Simple placeholder code for BLAS calls - replace with calls to a real BLAS library

#include <vector>

namespace BLAS {

// dot products ==============================================================

inline double dot(const std::vector<double> &x, const std::vector<double> &y) { 
    //return cblas_ddot((int)x.size(), &x[0], 1, &y[0], 1); 

    double sum = 0;
    for(size_t i = 0; i < x.size(); i++) {
        sum += x[i] * y[i];
    }
    return sum;
}

// inf-norm (maximum absolute value: index of max returned) ==================

inline int indexAbsMax(const std::vector<double> &x) { 
    //return cblas_idamax((int)x.size(), &x[0], 1); 

    int maxind = 0;
    double maxvalue = 0;
    for(int i = 0; i < (int)x.size(); i++) {
        if (fabs(x[i]) > maxvalue) {
            maxvalue = fabs(x[i]);
            maxind = i;
        }
    }
    return maxind;
}

// inf-norm (maximum absolute value) =========================================
// technically not part of BLAS, but useful

inline double absMax(const std::vector<double> &x) { 
    return std::fabs(x[indexAbsMax(x)]); 
}

// saxpy (y=alpha*x+y) =======================================================

inline void addScaled(double alpha, const std::vector<double> &x, std::vector<double> &y) { 
    //cblas_daxpy((int)x.size(), alpha, &x[0], 1, &y[0], 1); 

    for (size_t i = 0; i < x.size(); i++) {
        y[i] += alpha * x[i];
    }
}

}
#endif
