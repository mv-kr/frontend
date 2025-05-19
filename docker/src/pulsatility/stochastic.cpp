//
//  stochastic.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 22/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#include "stochastic.hpp"



using namespace mci;

// define static member
std::vector<mci::rng_type> mci::stochastic::generators(0);

/** sample a normally distributed random number**/
const  double stochastic::randn() {
    //return normal_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    return normal_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    
}

/** sample a uniformly random number in [0, 1]**/
const   double stochastic::rand() {
    return U_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
}

/** sample from a discete distribution with probabilities given in W**/
const Eigen::MatrixXd stochastic::rand_sample(const Eigen::VectorXd& W, const int n) {
    
    Eigen::MatrixXd ind(n, 1);
    std::vector<double> WW (W.data(), W.data()  + W.rows()* W.cols() );
    boost::random::discrete_distribution<> discrete_dist( WW.begin(), WW.end() );
    
    for (size_t i = 0; i<n; i++) {
        ind(i, 0) = discrete_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    }
    
    return ind;
    
}


/** sample from a discete distribution with probabilities given in weigths **/
const int stochastic::rand_sample(const Eigen::VectorXd& weigths) {
    Eigen::VectorXd W = weigths;
    
    bool all_nan = true;
    for (size_t i = 0;i<W.rows();i++ ) {
        all_nan = all_nan && std::isnan(W(i));
        if (std::isnan(W(i)) ) {
            W(i) = 0;
        }
        
        if (std::isinf(W(i))) {
            //std::cout<< std::endl << W.transpose() << std::endl;
            throw;
        }
    }
    
    std::vector<double> WW (W.data(), W.data()  + W.rows()* W.cols() );
    boost::random::discrete_distribution<> discrete_dist( WW.begin(), WW.end() );
    
    if (all_nan) {
        
        //std::cout<< all_nan << " " << weigths.transpose() << std::endl;
        return 0;
    }
    else {
        return discrete_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    }
    
}

/** sample from a discete distribution with probabilities given in weigths **/
const int stochastic::rand_sample(const Eigen::VectorXd& weigths, const Eigen::VectorXd& W2, const Eigen::VectorXd& W3) {
    Eigen::VectorXd W = weigths;
    
    bool all_nan = true;
    for (size_t i = 0;i<W.rows();i++ ) {
        all_nan = all_nan && std::isnan(W(i));
        if (std::isnan(W(i)) ) {
            W(i) = 0;
        }
        
        if (std::isinf(W(i))) {
            //std::cout<< std::endl << W2.transpose()<< std::endl;
            //std::cout<< std::endl << W3.transpose()<< std::endl;
            //std::cout<< std::endl << W3.maxCoeff() << std::endl;
            throw;
        }
    }
    
    std::vector<double> WW (W.data(), W.data()  + W.rows()* W.cols() );
    boost::random::discrete_distribution<> discrete_dist( WW.begin(), WW.end() );
    
    if (all_nan) {
        
        //std::cout<< all_nan << " " << weigths.transpose() << std::endl;
        return 0;
    }
    else {
        return discrete_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    }
    
}

const Eigen::VectorXd stochastic::rand_mvn(const Eigen::VectorXd& mean, const Eigen::MatrixXd& normTransform) {
    
    size_t size = normTransform.rows(); // Dimension (rows)
    size_t nn=1;             // How many samples (columns) to draw
    Eigen::VectorXd m_tmp(size,nn);
    
    
    for (size_t i=0; i<size; i++) {
        m_tmp(i) = normal_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    }
    
    
    // get sample
    Eigen::VectorXd sample = (normTransform * m_tmp) + mean;
    
    //return
    return sample;
}

const Eigen::VectorXd stochastic::rand_v(size_t N) {
    Eigen::VectorXd V(N);
    for (size_t i=0; i<N; i++) {
        //V(i) =  0.5+bernoulli_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
        V(i) =  rand();
    }
    
    return V;
}

const Eigen::VectorXd stochastic::rand_mixed(const Eigen::VectorXd& mean, const Eigen::MatrixXd& normTransform) {
    
    size_t size = normTransform.rows(); // Dimension (rows)
    size_t nn=1;             // How many samples (columns) to draw
    Eigen::VectorXd m_tmp(size,nn);
    Eigen::VectorXd mean2(size);
    
    double new_K_on = 20;//mean(size-1);
    double new_K_off = 4;//mean(size-1);
    if (this->rand() < -1.) {
        if (this->rand() < .5)
            new_K_on = new_K_on - 1 ;
        else
            new_K_on = new_K_on + 1;
    }
    
    for (size_t i=0; i<size; i++) {
        m_tmp(i) = normal_dist(mci::stochastic::generators[OpenMP::getThreadNum()]);
    }
    
    
    // make sure you perturb K_on/k_on
    mean2 = mean;
    mean2(0) = new_K_on/pow(10,mean(0));
    mean2(1) = new_K_off/pow(10,mean(1));
    // get sample
    Eigen::VectorXd sample = (normTransform * m_tmp) + mean2;
    if (sample(0) < 1.)
        sample(0) = 1;
    if (sample(1) < 1.)
        sample(1) = 1;
    
    sample(0)      =  - log10( sample(0) / new_K_on )  ;
    sample(1)      =  - log10( sample(1) / new_K_off )  ;
    //sample(size-1) = new_K;
    
    //std::cout << mean << std::endl;
    //std::cout << sample << std::endl;
    //return
    return sample;
}
