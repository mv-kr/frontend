//
//  stochastic.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 22/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef stochastic_hpp
#define stochastic_hpp

#include "openmp.hpp"
//#include <omp.h>

#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <Eigen/Dense>


namespace mci {
    
    using rng_type = boost::random::mt19937;
    
    class stochastic {
        
        boost::random::normal_distribution<double> normal_dist;
        //Uniform distribution object
        boost::random::uniform_01<> U_dist;
        //bernulli object
        boost::random::bernoulli_distribution<> bernoulli_dist;
 
    public:
        
        
        // vector of thread-private rng generators
        static std::vector<rng_type> generators;

        /** sample a normally distributed random number**/
        const  double randn() ;
        /** sample a uniformly random number in [0, 1]**/
        const   double rand();
        /** draw n samples from a discete distribution with probabilities given in W**/
        const Eigen::MatrixXd rand_sample(const Eigen::VectorXd& W, const int n) ;
        /** draw one sample from a discete distribution with probabilities given in W**/
        const int rand_sample(const Eigen::VectorXd& W) ;
        /** draw one sample from a discete distribution with probabilities given in W**/
        const int rand_sample(const Eigen::VectorXd& W,const Eigen::VectorXd& W2,const Eigen::VectorXd& W3) ;
        
        const Eigen::VectorXd rand_v(size_t N);
        
        const rng_type get_generator();
        
        

        /**sample form multivariate normal**/
        const Eigen::VectorXd rand_mvn(const Eigen::VectorXd& m, const Eigen::MatrixXd& normTransform);
        
        /**sample form multivariate normal**/
        const Eigen::VectorXd rand_mixed(const Eigen::VectorXd& m, const Eigen::MatrixXd& normTransform);
        
    };
    
}

#endif /* stochastic_hpp */
