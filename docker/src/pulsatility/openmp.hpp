//
//  openmp.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 02/06/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef openmp_hpp
#define openmp_hpp

#include <iostream>
#include <cstddef>
#include <omp.h>


namespace mci {
/**
 * Simple wrapper class to encapsulate OpenMP.
 */
    class OpenMP
    {
    public:
        /**
         * Returns the number of available threads.
         *
         */
        static inline std::size_t getMaxThreads() {
            return omp_get_max_threads();
        }
        
        /**
         * Returns the thread-identifier or 0 if OpenMP is disabled.
         */
        static std::size_t getThreadNum() {
            return omp_get_thread_num();
        }
    };
}
#endif /* openmp_hpp */
