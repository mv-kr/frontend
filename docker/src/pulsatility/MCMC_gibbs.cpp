//
//  MCMC_gibbs.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 24/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#include "MCMC_gibbs.hpp"
#include <fstream>
#include <iostream>

using namespace mci::algorithms;
using namespace mci::models;


/** const static members definitions**/
const std::string mci::algorithms::MCMC_gibbs::PARAM_iterations = "iterations";
const std::string mci::algorithms::MCMC_gibbs::PARAM_theta_dim = "theta_dim";
const std::string mci::algorithms::MCMC_gibbs::PARAM_report_every = "report_every";
const std::string mci::algorithms::MCMC_gibbs::NAME = "MCMC_gibbs";


MCMC_gibbs::MCMC_gibbs(const cMCMC_gibbs_parameter_t& p, mci::models::basic_model* m, cSMC_AS* csmc, SMC* smc,
                       const Eigen::MatrixXd& d, const Eigen::MatrixXd& c,
                       const stochastic& r,
                       const Eigen::VectorXd& i_g,
                       const Eigen::VectorXd& l_limit,
                       const Eigen::VectorXd& u_limit,
                       const Eigen::VectorXd& p_blocks) {
    
    /*set up internal state of the algorithm*/
    parameters = p;
    model_ptr = m;
    data = d;
    rng = r;
    cSMC_AS_ptr = csmc;
    SMC_ptr = smc;
    cov = c;
    i_guess = i_g;
    lower_limit = l_limit;
    upper_limit = u_limit;
    parameter_blocks = p_blocks;
    
    parameters[MCMC_gibbs::PARAM_theta_dim] = i_guess.size();
    
    Eigen::LLT<Eigen::MatrixXd> cholSolver(cov);
    
    // We can only use the cholesky decomposition if
    // the covariance matrix is symmetric, pos-definite.
    // But a covariance matrix might be pos-semi-definite.
    // In that case, we'll go to an EigenSolver
    if (cholSolver.info()==Eigen::Success) {
        // Use cholesky solver
        normTransform = cholSolver.matrixL();
    } else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(cov);
        normTransform = eigenSolver.eigenvectors()
        * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    
    
    
}

/** generate intial vector of parameters () **/
const void MCMC_gibbs::generate_init_params()  {
    
    for (size_t i = 0 ;i<this->lower_limit.size(); i++) {
        this->i_guess(i) = this->lower_limit(i) + (this->upper_limit(i)-this->lower_limit(i))*rng.rand();
    }
}

/** run the algorithm and store results in res **/
const void MCMC_gibbs::run(chain_doubles_t& ch) {
    
    unsigned int iterations = parameters[MCMC_gibbs::PARAM_iterations];
    unsigned int report_every = parameters[MCMC_gibbs::PARAM_report_every];
    size_t theta_dim =  parameters[MCMC_gibbs::PARAM_theta_dim];
    size_t model_state_dim = model_ptr->get_dim();
    size_t model_timepoints = round( (data(data.rows()-1,data.cols()-1) - data(data.rows()-1,0))/ model_ptr->get_dt() ) + 1;
    size_t start_block_i = 0;
    size_t size_block_i  = 0;
    size_t block_no = parameter_blocks.size()-1;
    /** variables hodling chain results **/
    Eigen::MatrixXd param_chain(theta_dim, iterations);
    // vector holding target distribution values
    Eigen::VectorXd likelihoods(iterations);
    // vector holding accept/reject events
    //Eigen::Matrix<size_t, Eigen::Dynamic, 1> accepts(iterations, 1);
    Eigen::MatrixXd accepts(iterations, block_no);
    
    /** aux variables **/
    Eigen::VectorXd theta_current(theta_dim), theta_proposal(theta_dim);
    Eigen::MatrixXd theta_tmp(1, theta_dim);
    Eigen::MatrixXd X_sample(model_state_dim, model_timepoints);
    Eigen::MatrixXd X_sample_1(model_state_dim, model_timepoints);
    
    double zz = 0.0;
    double L = 0.0;
    double log_prior_ratio = 0.0;
    double log_proposal_density_ratio = 0.0;
    double L1 = 0.0;
    double L2 = 0.0;
    double ll = 0.0;
    double ll_c = 0.0;
    double acceptance_ratio;
    double time;
    bool cond;
    size_t mm =0;
    
    std::vector<Eigen::MatrixXd> trajectories_res_tmp(iterations/report_every, X_sample) ;
    
    /*model_ptr->update_i_parameters( this->i_guess );
     trajectory_t Xfinal = model_ptr->run(model_ptr->generate_ic(-1), -1,480);
     std::cout<< Xfinal <<std::endl;
     */
    /*** Initialise ***/
    
    std::cout<< "Initialising using vanilla SMC" << std::endl;
    
   // generate_init_params();
    model_ptr->update_i_parameters( this->i_guess );
    std::cout<<this->i_guess<<std::endl;
    SMC_ptr->run(X_sample);
    model_ptr->log_joint_density(X_sample, data, &ll);
    
    param_chain.col(0) = this->i_guess;
    accepts(0)         = 0;
    likelihoods(0)     = ll;
    
    
    std::cout<< "Starting MCMC Gibbs ... @ theta (" << param_chain.col(0).transpose() << ")" << std::endl;
    double  mm_v  = ll;
    
    // MCMC main body
    time = (omp_get_wtime() );
    for (unsigned int i=1 ; i<iterations; i++) {
        
        // - current parameter vector
        theta_current = param_chain.col(i-1);
        X_sample_1 = X_sample;
        // - evaluate log-likelihood ratio L
        model_ptr->update_i_parameters( theta_current.transpose() );
        L1 = model_ptr->log_joint_density(X_sample_1, data, &ll_c);
        
        accepts(i) = 0;
        for (size_t pb_i=0 ; pb_i<parameter_blocks.size()-1; pb_i++) {
            
            /* [0, 4, 8] *** last element is the size of the entire parameter dataset
             start_block_i = 4
             size_bloc_i   = 7-4*/
            start_block_i = parameter_blocks(pb_i);
            size_block_i   = parameter_blocks(pb_i+1) - parameter_blocks(pb_i);
            // - draw a new sample vector for the current block of parameters using the proposal distribution
            theta_proposal = theta_current;
            theta_proposal.middleRows(start_block_i,size_block_i) = rng.rand_mvn(
                                                                                 theta_current.middleRows(start_block_i, size_block_i),
                                                                                 normTransform.block(start_block_i, start_block_i, size_block_i, size_block_i)
                                                                                 );
            
            //if proposal within bounds
            cond = (theta_proposal.array() < this->lower_limit.array()).any() || (theta_proposal.array() > this->upper_limit.array()).any();
            
            if ( ! cond  ) {
                
                // - evaluate log-likelihood ratio L
                //model_ptr->update_i_parameters( theta_current.transpose() );
                //L1 = model_ptr->log_joint_density(X_sample_1, data, &ll_c);
                
                model_ptr->update_i_parameters(theta_proposal.transpose());
                L2 = model_ptr->log_joint_density(X_sample_1, data, &ll);
                //std::cout<<L2<<std::endl;
                L = - L1 + L2;
                
                /* calculate the acceptance ratio */
                log_prior_ratio = model_ptr->log_prior_density(theta_proposal)
                - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
                //curently only symmetric proposal densities are supported
                log_proposal_density_ratio = model_ptr->log_proposal_density_ratio( theta_proposal, theta_current); // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
                
                acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
                
                
                if (std::isnan(acceptance_ratio)) {
                    acceptance_ratio = 0;
                  //  throw std::runtime_error("Error evaluating acceptance ratio.");
                }
                
                // - draw a random number between 0 and 1
                // - accept or reject by comparing the random number to the
                // Metropolis-Hastings ratio (acceptance probability); if accept,
                // change the current value of theta to the proposed theta,
                // update the current value of the target and keep track of the
                // number of accepted proposals
                if (rng.rand() < std::min(1., acceptance_ratio)) {
                    
                    param_chain.col(i) = theta_proposal;
                    likelihoods(i)     = ll;
                    accepts(i,pb_i)         = 1;
                    
                    theta_current      = theta_proposal;
                    L1 = L2;
                    ll = ll_c;
                    
                } else {
                    param_chain.col(i) = theta_current;
                    likelihoods(i)     = ll_c;
                    accepts(i,pb_i)        = 0;
                    
                }
                
            } else {
                
                param_chain.col(i) = theta_current;
                likelihoods(i)     = ll_c;
                accepts(i,pb_i)        = 0;
                
            }
            
            
            theta_current = param_chain.col(i);
        }
        
        // - sample  trajectory given set of parameters and store in X_sample
        model_ptr->update_i_parameters(param_chain.col(i));
        mm_v = cSMC_AS_ptr->run(X_sample_1, X_sample);
        
        
        
        
        /* report progress every report_every iterations */
        if ((i+1)%report_every  == 0 ){
            
            time = (omp_get_wtime() - time);
            std::cout << "thread " << OpenMP::getThreadNum()
            << " --- iteration " << i+1
            << " log_lik " << likelihoods(i)
            << " acceptance prob. " << (accepts.block(i+1-report_every,0,report_every,block_no).colwise().sum()/(report_every))
            << " time " << time
            << "| theta (" << param_chain.col(i).transpose()
            << ") | " <<std::endl;
            
            trajectories_res_tmp[mm]= X_sample ;
            mm++;
            time = (omp_get_wtime() );
            
            //std::cout<< X_sample <<std::endl;
        }
        
        
        
    }//end iterations
    
    
    
    
    this->chain_res = param_chain;
    this->likelihood_res = likelihoods;
    this->trajectories_res = trajectories_res_tmp;
    
    
}//end run function

/** run the algorithm using parameters in params and store results in res **/
const void MCMC_gibbs::run(const cMCMC_gibbs_parameter_t& p, chain_doubles_t& ch) {}


const void MCMC_gibbs::print_chain_parameters() {
    std::cout << this->chain_res << std::endl;
}



const void MCMC_gibbs::print_chain_parameters(std::string str_file) {
    
    //remove( str_file.c_str() );
    
    std::ofstream file_tmp;
    file_tmp.open(str_file);
    file_tmp << this->chain_res << std::endl;
    file_tmp.close();
    
}


const void MCMC_gibbs::print_chain_likelihood(std::string str_file) {
    
    //remove( str_file.c_str() );
    
    std::ofstream file_tmp;
    file_tmp.open(str_file);
    file_tmp << this->likelihood_res << std::endl;
    file_tmp.close();
    
}


const void MCMC_gibbs::print_chain_trajectories(std::string str_file) {
    
    //remove( str_file.c_str() );
    
    std::ofstream file_tmp;
    file_tmp.open(str_file);
    for (size_t mm = 0; mm< this->trajectories_res.size(); mm++) {
        file_tmp << this->trajectories_res[mm] << std::endl;
    }
    file_tmp.close();
    
}

