//
//  MCMC_gibbs_hmc.cpp
//  mci
//
//  Created by Margaritis Voliotis on 20/09/2020.
//  Copyright © 2020 Margaritis Voliotis. All rights reserved.
//

#include "MCMC_gibbs_sMMALA.hpp"

#include <fstream>

using namespace mci::algorithms;
using namespace mci::models;


/** const static members definitions**/
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iterations = "iterations";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_warmup_pct = "warmup_pct";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_theta_dim = "theta_dim";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_hmc_theta_dim = "hmc_theta_dim";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_epsilon = "epsilon";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_steps = "steps";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iter_adapt = "iter_adapt";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_delta = "delta";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_report_every = "report_every";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_blocks_no = "blocks_no";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::PARAM_mmala_blocks_no = "mmala_blocks_no";
const std::string mci::algorithms::MCMC_gibbs_sMMALA::NAME = "MCMC_gibbs_sMMALA";


MCMC_gibbs_sMMALA::MCMC_gibbs_sMMALA(const MCMC_gibbs_hmc_parameter_t& p, mci::models::basic_model* m, cSMC_AS* csmc, SMC* smc,
                   const Eigen::MatrixXd& d, const Eigen::MatrixXd& c, const Eigen::MatrixXd& M,
                   const stochastic& r,
                   const Eigen::VectorXd& i_g,
                   const Eigen::VectorXd& l_limit,
                   const Eigen::VectorXd& u_limit,
                   const Eigen::VectorXd& bl_i) {
    
    parameters = p;
    model_ptr = m;
    data = d;
    rng = r;
    cSMC_AS_ptr = csmc;
    SMC_ptr = smc;
    cov = c;
    mass = M;
    i_guess = i_g;
    lower_limit = l_limit;
    upper_limit = u_limit;
    blocks_i = bl_i;
    
    parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim] = i_guess.size();
    
    /**int iterations = parameters["iterations"];
     int theta_dim =  parameters["theta_dim"];
     size_t model_timepoints = parameters["timepoints"];
     size_t model_state_dim = model_ptr->get_dim();
     Eigen::MatrixXd X_sample(model_state_dim, model_timepoints);**/
    
    //chain_res(theta_dim,iterations);
    //likelihood_res(iterations);
    //trajectories_res();
    
    
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
    
    Eigen::LLT<Eigen::MatrixXd> cholSolver2(mass);
    
    // We can only use the cholesky decomposition if
    // the covariance matrix is symmetric, pos-definite.
    // But a covariance matrix might be pos-semi-definite.
    // In that case, we'll go to an EigenSolver
    if (cholSolver2.info()==Eigen::Success) {
        // Use cholesky solver
        M_normTransform = cholSolver2.matrixL();
    } else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(mass);
        M_normTransform = eigenSolver.eigenvectors()
        * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    
    
    mass_inv = M.inverse();
    mass_det = M.determinant();
}

Eigen::MatrixXd MCMC_gibbs_sMMALA::leapfrog( Eigen::VectorXd& theta, Eigen::VectorXd &p, const Eigen::MatrixXd& X_sample, double epsilon, bool *cond) {
    
    int steps = parameters[MCMC_gibbs_sMMALA::PARAM_steps];
    //std::min(std::max(1., std::round(parameters[MCMC_gibbs_hmc::PARAM_steps]/epsilon)),10000.);//;////parameters[MCMC_gibbs_hmc::PARAM_steps];;//
    
   
    
    size_t number_of_i_parameters = parameters[MCMC_gibbs_sMMALA::PARAM_hmc_theta_dim];
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd grad_term = Eigen::VectorXd::Zero(theta_dim) ;
    Eigen::VectorXd theta_star = Eigen::VectorXd::Zero(theta_dim) ;
    Eigen::MatrixXd X_sample_1 = X_sample ;
    Eigen::MatrixXd X_sample_2 = X_sample ;
    //Eigen::MatrixXd mass_inv_sub = mass_inv.topLeftCorner(number_of_i_parameters,number_of_i_parameters)
    //Eigen::MatrixXd theta_seq(steps+1, theta_dim);
     Eigen::MatrixXd theta_seq(1, theta_dim);
    bool exit_cond = false;
    
    //theta_seq.row(0) = theta;
    
    //grad_term = get_grad_AIS(theta, X_sample_1, &ll_star, &joint_Lik,epsilon);
    model_ptr-> update_i_parameters( theta );
    grad_term.topRows(number_of_i_parameters) = model_ptr->grad_log_joint_density(X_sample, data);//

    //grad_term = get_grad_fast(theta, X_sample, &ll_star, &joint_Lik,epsilon);
    // std::cout<<" leapfrog grad1:" << k << " " <<grad_term.transpose()<<std::endl;
    grad_term.topRows(number_of_i_parameters) += model_ptr->grad_log_prior_density(theta);
    // std::cout<<" leapfrog grad2:" <<grad_term.transpose()<<std::endl;
    
    
    for (size_t k = 1; k<=steps; k++) {
        //std::cout<<" leapfrog step "<< k << "theta: "  <<theta.transpose()<<std::endl;
        
        p     = p  + epsilon/2.0*grad_term;
        
        
        theta = theta + epsilon*mass_inv*p;
     
        //check if propose theta out of bounds and break if indeed
        exit_cond = (theta.array() < this->lower_limit.array()).any() || (theta.array() > this->upper_limit.array()).any() || isnan(theta.array()).any();// || theta(4)<theta(3);
        if (exit_cond) {
            
      //      theta = theta_seq.row(k-1).transpose();
            
            break;
        }
        
        // update model parameters
        //model_ptr->update_i_parameters(theta);
        //cSMC_AS_ptr->run(X_sample_1, X_sample_2);
        
        /*** estimate grad. In principle one would like to get the grad log P(obs|theta) =
         1/(2e) [logP(obs|theta+e)-ogP(obs|theta-e)] **/
        //grad_term = get_grad_AIS(theta, X_sample_2, &ll_star, &joint_Lik,epsilon);
        //grad_term = get_grad_fast(theta, X_sample, &ll_star, &joint_Lik,epsilon);
        model_ptr->update_i_parameters( theta );
        //grad_term = get_grad_fast(theta, X_sample, &ll_star, &joint_Lik,epsilon);
        //std::cout<<" leapfrog grad1:" << k << " " <<grad_term.transpose()<<std::endl;
        
        grad_term.topRows(number_of_i_parameters) = model_ptr->grad_log_joint_density(X_sample, data);
        //std::cout<<" leapfrog grad2:" << k << " " <<grad_term.transpose()<<std::endl;
        grad_term.topRows(number_of_i_parameters) += model_ptr->grad_log_prior_density(theta);
        
        p     = p + epsilon/2.0*grad_term;
        
  
     //   theta_seq.row(k) = theta.transpose();
        //X_sample_1 = X_sample_2;
    }
    
    if (cond)
        *cond = exit_cond;
    
   // std::cout<<" ---> " << theta.transpose()<<std::endl;
    return theta_seq;
}



const void MCMC_gibbs_sMMALA::run(chain_doubles_t& ch) {
    
    
    
    // model parameters
    size_t theta_dim =  parameters[PARAM_theta_dim]; // number of inferred parameters
    //size_t hmc_theta_dim =  parameters[PARAM_hmc_theta_dim]; // number of parameters in the hmc step
    
    size_t model_state_dim = model_ptr->get_dim(); // number of state variables in the model
    size_t model_timepoints = round( (data(data.rows()-1,data.cols()-1) - data(data.rows()-1,0))/ model_ptr->get_dt() ) + 1;
    
    size_t parameter_blocks = parameters[PARAM_blocks_no];
    size_t parameter_mmala_blocks =  parameters[PARAM_mmala_blocks_no]; // number of parameters in the hmc step
    /** aux variables **/
    Eigen::VectorXd theta(theta_dim);
    Eigen::MatrixXd X_sample(model_state_dim, model_timepoints);
    Eigen::MatrixXd X_sample_1(model_state_dim, model_timepoints);
    double acceptance_prob = 0.0;
    double joint_llik = 0.0;
    double marginal_llik = 0.0;
    /** output parameters & variables **/
    unsigned int iterations = parameters[PARAM_iterations];
   // unsigned int warmup_iterations = round(iterations*parameters[PARAM_warmup_pct]);
    unsigned int report_every = parameters[PARAM_report_every];
    unsigned int div = 0;
    size_t start_block_i, size_block_i;
    
    size_t mm =0;
    std::vector<Eigen::MatrixXd> trajectories_res_tmp(iterations/report_every, X_sample) ;
    double time;
    
    
    
    //
    Eigen::VectorXd mu_theta = Eigen::VectorXd::Zero(theta_dim );
    Eigen::MatrixXd cov_theta = Eigen::MatrixXd::Zero(theta_dim, theta_dim);
    
    // mcmc parameters
    
    double epsilon = parameters[MCMC_gibbs_sMMALA::PARAM_epsilon];
    unsigned int iter_adapt = parameters[MCMC_gibbs_sMMALA::PARAM_iter_adapt];
    double eps_0;
    double mu;
    double H_0_bar = 0;
    double gamma_param = 0.05;
    double delta_param = parameters[MCMC_gibbs_sMMALA::PARAM_delta];
    double t_0 = 10;
    double kappa = 0.75;
    Eigen::MatrixXd epsilons = Eigen::MatrixXd::Zero(iter_adapt,parameter_blocks);
    Eigen::MatrixXd epsilons_bar = Eigen::MatrixXd::Zero(iter_adapt,parameter_blocks);
    Eigen::MatrixXd H_bar = Eigen::MatrixXd::Zero(iter_adapt,parameter_blocks);
    
    
    /*** Initialise ***/
    std::cout<< "Initialising latent variables using vanilla SMC and theta_0 = ( " << i_guess.transpose() << " )" << std::endl;
    
    model_ptr->update_i_parameters( this->i_guess );
    SMC_ptr->run(X_sample);
    
    /** initalise method parameters **/
    std::cout<< "Initialising sMMALA parameters ..." << std::endl;
    
    eps_0 =  1;//find_epsilon(this->i_guess, X_sample);
    mu = log(10.*eps_0);
    for (size_t bli = 0; bli<parameter_blocks; bli++) {
        epsilons(0,bli) = log(eps_0);
        epsilons_bar(0,bli) = log(eps_0);
        H_bar(0,bli) = H_0_bar;
    }
    epsilon = eps_0;
    
    //report parameters
    std::cout<< " --- epsilon = " << eps_0 << std::endl;
    std::cout<< " --- mu = " << mu << std::endl;
    std::cout<< " --- gamma = " << gamma_param << std::endl;
    std::cout<< " --- kappa = " << kappa << std::endl;
    std::cout<< " --- t_0 = " << t_0 << std::endl;
    std::cout<< " --- delta = " << delta_param << std::endl;
   // std::cout<< " --- steps = " << steps << std::endl;
    
    

    cSMC_AS_ptr->run(X_sample, X_sample_1);
    X_sample =X_sample_1;
    joint_llik = model_ptr->log_joint_density(X_sample, data, &marginal_llik);
    
    
    
    // ------------------------- warm-up --------------------
    /** variables hodling chain results **/
    Eigen::MatrixXd param_chain(theta_dim, iterations);
    // vector holding marginal log likelihood
    Eigen::VectorXd likelihoods(iterations);
    // vector holding acceptance probabilities
    Eigen::MatrixXd accepts = Eigen::MatrixXd::Zero(iterations, parameter_blocks);
    
    
    std::cout<< "Starting MCMC hmc warmup ... (" << iter_adapt << " iterations)" << std::endl;
    param_chain.col(0) = this->i_guess;
    likelihoods(0)     = marginal_llik;
    
    time = (omp_get_wtime() ); // get time
    for (unsigned int i=1 ; i<iter_adapt; i++) {
        // - get current parameter vector
        theta = param_chain.col(i-1);
        
        
        // - 1. sample trajectory (latent variables) given all parameters
        model_ptr->update_i_parameters(theta); // update model parameters
        cSMC_AS_ptr->run(X_sample, X_sample_1); // run cSMC to get a new trajecotry conditioned on the current
        X_sample = X_sample_1; // update trajectory
        
        
        
        
        // - 2 sample using simplified MMALA
        for (size_t pb_i=0 ; pb_i<parameter_mmala_blocks; pb_i++) {
        // size_t pb_i=0;
            start_block_i  = blocks_i(pb_i); //+ pb_i;
            size_block_i   = blocks_i(pb_i+1)-blocks_i(pb_i);
            
            epsilon = exp( epsilons(i-1,pb_i) );
            sMMALA_sampler_block(theta,
                           X_sample,
                           epsilon,
                           acceptance_prob,
                           joint_llik,
                           marginal_llik,
                           div,
                           start_block_i,size_block_i);
            
            param_chain.col(i) = theta;
            likelihoods(i)     = marginal_llik;
            accepts(i,pb_i)       = acceptance_prob;
            
            
            
            
            // - adapt epsilon parameter
            
            H_bar(i,pb_i) = (1.-1./(i+t_0))*H_bar(i-1) + 1./(i+t_0) * (delta_param - accepts(i,pb_i)  );
            epsilons(i,pb_i) = mu - sqrt(i)/gamma_param * H_bar(i,pb_i);
            epsilons_bar(i,pb_i) = pow(i,-kappa) * epsilons(i,pb_i) + (1- pow(i,-kappa)) * epsilons_bar(i-1,pb_i);
        
        }
        
        
        
        
        // - 3. sample using MH
        
        for (size_t pb_i=parameter_mmala_blocks ; pb_i<parameter_blocks; pb_i++) {
            
            start_block_i  = blocks_i(pb_i); //+ pb_i;
            size_block_i   = blocks_i(pb_i+1)-blocks_i(pb_i);
            epsilon  =  exp( epsilons(i-1,pb_i) );
           
            MH_sampler_single_parameter( theta,
                                        X_sample,
                                        epsilon,
                                        acceptance_prob,
                                        joint_llik,
                                        marginal_llik,
                                        start_block_i,size_block_i);
                
            param_chain.col(i) = theta;
            likelihoods(i)     = marginal_llik;
            accepts(i,pb_i)    = acceptance_prob;
         
            
            
            // - adapt epsilon parameters
            H_bar(i,pb_i) = (1.-1./(i+t_0))*H_bar(i-1,pb_i) + 1./(i+t_0) * (delta_param - accepts(i,pb_i)  );
            epsilons(i,pb_i) = mu - sqrt(i)/gamma_param * H_bar(i,pb_i);
            epsilons_bar(i,pb_i) = pow(i,-kappa) * epsilons(i,pb_i) + (1- pow(i,-kappa)) * epsilons_bar(i-1,pb_i);
        
            
        }
        
      
        // - report progress every report_every iterations //
        if ((i+1)%report_every  == 0 ){
            
            time = (omp_get_wtime() - time);
            std::cout << "thread " << OpenMP::getThreadNum() << " --- iteration (adapt) " << i+1 << "/" << (iter_adapt)
            << "(time " << time
            << "): log_lik " << marginal_llik
            << " div. " << (((double)div)/report_every)
            << " | theta (" << theta.transpose()
            << ") | epsilon: " << epsilons.row(i)
            << std::endl;
            
            div = 0;
            
            time = (omp_get_wtime() );
        }
        
        
        
    }//end adaptive iterations
    
    
    std::cout<< std::endl << "Starting sampling current theta_0 = (" << param_chain.col(iter_adapt-1).transpose() << ")"<< std::endl;
    
    
    // --------------------------- MCMC gibbs main body ------------------------
    
    param_chain.col(0) = param_chain.col(iter_adapt-1);
    likelihoods(0)     =  likelihoods(iter_adapt-1) ;
    
    Eigen::VectorXd epsilons_f = epsilons_bar.row(iter_adapt-1);
    
    time = (omp_get_wtime() ); // get time
    
    
    for (unsigned int i=1 ; i<iterations; i++) {
        // - get current parameter vector
        theta = param_chain.col(i-1);
        
        
        // - 1. sample trajectory (latent variables) given all parameters
        model_ptr->update_i_parameters(theta); // update model parameters
        cSMC_AS_ptr->run(X_sample, X_sample_1); // run cSMC to get a new trajecotry conditioned on the current
        X_sample = X_sample_1; // update trajectory
        
        
        
        
        // - 2 sample using simplified MMALA
        for (size_t pb_i=0 ; pb_i<parameter_mmala_blocks; pb_i++) {
        // size_t pb_i=0;
            start_block_i  = blocks_i(pb_i); //+ pb_i;
            size_block_i   = blocks_i(pb_i+1)-blocks_i(pb_i);
            
            epsilon = exp( epsilons_f(pb_i) );
            sMMALA_sampler_block(theta,
                           X_sample,
                           epsilon,
                           acceptance_prob,
                           joint_llik,
                           marginal_llik,
                           div,
                           start_block_i,size_block_i);
            
            param_chain.col(i) = theta;
            likelihoods(i)     = marginal_llik;
            accepts(i,pb_i)       = acceptance_prob;
            
            
            
            
            // - adapt epsilon parameter (hmc)
           /* if (i <= iter_adapt) {
                H_bar(i,pb_i) = (1.-1./(i+t_0))*H_bar(i-1) + 1./(i+t_0) * (delta_param - accepts(i,pb_i)  );
                
                epsilons(i,pb_i) = mu - sqrt(i)/gamma_param * H_bar(i,pb_i);
                epsilons_bar(i,pb_i) = pow(i,-kappa) * epsilons(i,pb_i) + (1- pow(i,-kappa)) * epsilons_bar(i-1,pb_i);
            } else {
                epsilons(i,pb_i) = epsilons_bar(iter_adapt,pb_i);
            }*/
        }
        
        
        
        
        // - 3. sample using MH
        
        for (size_t pb_i=parameter_mmala_blocks ; pb_i<parameter_blocks; pb_i++) {
            
            start_block_i  = blocks_i(pb_i); //+ pb_i;
            size_block_i   = blocks_i(pb_i+1)-blocks_i(pb_i);
            epsilon  =  exp( epsilons_f(pb_i) );
           
            MH_sampler_single_parameter( theta,
                                        X_sample,
                                        epsilon,
                                        acceptance_prob,
                                        joint_llik,
                                        marginal_llik,
                                        start_block_i,size_block_i);
                
            param_chain.col(i) = theta;
            likelihoods(i)     = marginal_llik;
            accepts(i,pb_i)    = acceptance_prob;
         
            
            
            // - adapt epsilon parameters
           /* if (i <= iter_adapt) {
                H_bar(i,pb_i) = (1.-1./(i+t_0))*H_bar(i-1,pb_i) + 1./(i+t_0) * (delta_param - accepts(i,pb_i)  );
                
                epsilons(i,pb_i) = mu - sqrt(i)/gamma_param * H_bar(i,pb_i);
                epsilons_bar(i,pb_i) = pow(i,-kappa) * epsilons(i,pb_i) + (1- pow(i,-kappa)) * epsilons_bar(i-1,pb_i);
            } else {
                epsilons(i,pb_i) = epsilons_bar(iter_adapt,pb_i);
            }*/
            
        }
        
      
        // - report progress every report_every iterations //
        if ((i+1)%report_every  == 0 ){
            //double acc_prob = (((double)accepts.middleRows(i+1-report_every,report_every).sum())/(report_every));
            time = (omp_get_wtime() - time);
            std::cout << "thread " << OpenMP::getThreadNum() << " --- iteration " << i+1 << "/" << (iterations)
            << "(time " << time
            << "): log_lik " << likelihoods(i)
            << " acceptance prob. " << (accepts.block(i+1-report_every,0,report_every,accepts.cols()).colwise().sum()/(report_every))
            << " div. " << (((double)div)/report_every)
            << " | theta (" << param_chain.col(i).transpose()
            << ") | epsilon: " << epsilons_f.transpose()
            << std::endl;
            
            div = 0;
            
            trajectories_res_tmp[mm]= X_sample ;
            mm++;
            time = (omp_get_wtime() );
        }
        
        
        
    }//end iterations
    
    
    
    
    this->chain_res = param_chain;
    this->likelihood_res = likelihoods;
    this->trajectories_res = trajectories_res_tmp;
    
    
}//end run function


const void MCMC_gibbs_sMMALA::MH_sampler_single_parameter(Eigen::VectorXd& theta,
                                const Eigen::MatrixXd& X_sample,
                                const double& epsilon,
                                double& acceptance_prob,
                                double& joint_llik,
                                double& marginal_llik,
                                const size_t& start_block_i,
                                const size_t& size_block_i) {
    
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd theta_current(theta_dim), theta_star(theta_dim);
    double joint_Lik_star,  ll_star, L, acceptance_ratio;
    bool  cond = false;
    
    theta_current = theta;
    // - draw a new sample vector for the current block of parameters using the proposal distribution
    theta_star = theta_current;
    
    theta_star.segment(start_block_i,size_block_i) = rng.rand_mvn(    theta_current.segment(start_block_i, size_block_i),
                                                               epsilon*Eigen::MatrixXd::Identity(size_block_i,size_block_i)
                                                                         );
    
    //if proposal within bounds
    //if (start_block_i <5)
        cond = (theta_star.segment(start_block_i,size_block_i).array() < this->lower_limit.segment(start_block_i,size_block_i).array()).any() ||        (theta_star.segment(start_block_i,size_block_i).array() > this->upper_limit.segment(start_block_i,size_block_i).array()).any();
    /*else if (start_block_i == 5) {
        
        double tmp = theta_star(5) ;
        double lb = this->lower_limit(5);
        double ub = this->upper_limit(5);
        
        cond = ( tmp < lb || tmp > ub) ;
        
        
    } else if (start_block_i == 7) {
        
        double tmp = theta_star(7) ;
        double lb = this->lower_limit(7);
        double ub = this->upper_limit(7);
        
        cond = ( tmp < lb || tmp > ub) ;
        
    } else if (start_block_i == 6) {
        
      //  double tmp = theta_star(5) + theta_star(6);//log10(1 + theta_star(6) );
      //  double lb = this->lower_limit(5);
      //  double ub = this->upper_limit(5);
        
       // double tmp2 = theta_star(6) ;
       // double lb2 = this->lower_limit(6);
       // double ub2 = this->upper_limit(6);
        
        cond = false;//( tmp < lb || tmp > ub);// || ( tmp2 < lb2 || tmp2 > ub2);
        
    } else if (start_block_i == 8) {
      //  double tmp = theta_star(7) +  theta_star(8) ;//log10(1 + theta_star(8) );
      //  double lb = this->lower_limit(7);
      //  double ub = this->upper_limit(7);
        
      //  double tmp2 = theta_star(8) ;
      //  double lb2 = this->lower_limit(8);
      //  double ub2 = this->upper_limit(8);
        
        //std::cout<< theta_current.transpose() <<std::endl;
        //std::cout<< theta_star.transpose()    <<std::endl;
        
        cond = false;//( tmp < lb || tmp > ub);// || ( tmp2 < lb2 || tmp2 > ub2);
    }*/
        
    
    acceptance_prob = 0.0;
    
    if ( ! cond  ) {
        
        model_ptr->update_i_parameters(theta_star);
        joint_Lik_star = model_ptr->log_joint_density(X_sample, data, &ll_star);
        
        //L = - L1 + L2;
        L = - joint_llik + joint_Lik_star;
        
        /* calculate the acceptance ratio */
        double log_prior_ratio = model_ptr->log_prior_density(theta_star)
        - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
        //curently only symmetric proposal densities are supported
        double log_proposal_density_ratio = model_ptr-> log_proposal_density_ratio( theta_star, theta_current); // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
        
        acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
        
        
        if (std::isnan(acceptance_ratio)) {
            acceptance_prob = 0.0;
            //  throw std::runtime_error("Error evaluating acceptance ratio.");
        } else {
            acceptance_prob = std::min(1., acceptance_ratio);
            // - draw a random number between 0 and 1
            // - accept or reject by comparing the random number to the
            // Metropolis-Hastings ratio (acceptance probability); if accept,
            // change the current value of theta to the proposed theta,
            // update the current value of the target and keep track of the
            // number of accepted proposals
            if (rng.rand() < acceptance_prob) {
                
                theta             = theta_star;
                marginal_llik     = ll_star;
                joint_llik        = joint_Lik_star;
                
               
            }
        }
    } // end out of bounds if
}



/** simplified MMALA implementation **/
const void MCMC_gibbs_sMMALA::sMMALA_sampler(Eigen::VectorXd& theta,
                                       const Eigen::MatrixXd& X_sample,
                                       const double& epsilon,
                                       double& acceptance_prob,
                                       double& joint_llik,
                                       double& marginal_llik,
                                       unsigned int& div,
                                       unsigned int iter) {
    
    size_t number_of_i_parameters = parameters[PARAM_hmc_theta_dim];
    size_t theta_dim =  parameters[PARAM_theta_dim];
    Eigen::VectorXd theta_current(theta_dim) , theta_star(theta_dim) ;
    
    Eigen::VectorXd mu(number_of_i_parameters);
    Eigen::VectorXd grad_term(number_of_i_parameters);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters);
    Eigen::MatrixXd G_inv(number_of_i_parameters,number_of_i_parameters), sigma(number_of_i_parameters,number_of_i_parameters);
    
    double H_star_to_current = 0.0, H_current_to_star = 0.0, joint_Lik_star = 0.0 , ll_star = 0.0 , L = 0.0, acceptance_ratio= 0.0;
    bool  cond = false;
    
    // - initialise current parameter vector and proposal
    theta_current = theta;
    theta_star = theta_current;
    
    // - set up model parameters
    model_ptr->update_i_parameters(theta_current);
    // - gather info to calculate proposal
    // - 1. gradient of log target density
    grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_llik, &marginal_llik, &G);//
    grad_term += (model_ptr->grad_log_prior_density(theta_current));
    model_ptr->add_prior_information(theta_current,&G);
    
    G_inv =  G.inverse();
    
    mu = theta_current.topRows(number_of_i_parameters) + pow(epsilon,2)*G_inv*grad_term/2.;
    
    Eigen::LLT<Eigen::MatrixXd> cholSolver( pow(epsilon,2)*G_inv );
    if (cholSolver.info()==Eigen::Success) {
        // Use cholesky solver
        sigma = cholSolver.matrixL();
    } else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver( pow(epsilon,2)*G_inv);
        sigma = eigenSolver.eigenvectors()
        * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    
    // - sample proposal
    theta_star.topRows(number_of_i_parameters)  = rng.rand_mvn(mu , sigma);
    
    //if proposal within bounds
    cond = (theta_star.array() < this->lower_limit.array()).any() || (theta_star.array() > this->upper_limit.array()).any();
    acceptance_prob = 0.0;
    
    if ( ! cond  ) {

        H_current_to_star= -1./(2.*pow(epsilon,2))*(theta_star.topRows(number_of_i_parameters)-mu).transpose()*G*(theta_star.topRows(number_of_i_parameters)-mu) - 1./2.*log(G_inv.determinant());
        
       /* std::cout<< "--------------- "<<std::endl;
        std::cout<< sigma  << std::endl << "--" << H_current_to_star <<std::endl;
         std::cout<< G_inv  << std::endl << "--"<< epsilon <<std::endl;
        std::cout << mu.transpose() <<std::endl;
        std::cout<< grad_term.transpose() <<std::endl;
         std::cout<<  theta_star.transpose() <<std::endl;
        std::cout<<  theta_current.transpose() <<std::endl;*/
        
        // gather info to form proposal
        model_ptr->update_i_parameters(theta_star);
        grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_Lik_star, &ll_star, &G);//
        grad_term += (model_ptr->grad_log_prior_density(theta_star));
        model_ptr->add_prior_information(theta_star,&G);
        
        G_inv   = G.inverse();
        mu = theta_star.topRows(number_of_i_parameters) + pow(epsilon,2)*G_inv*grad_term/2.;
        H_star_to_current= -1./(2.*pow(epsilon,2))*(theta_current.topRows(number_of_i_parameters)-mu).transpose()*G*(theta_current.topRows(number_of_i_parameters)-mu)
                            + 1./2.*log(G.determinant());
        
        L = - joint_llik + joint_Lik_star;
        /* calculate the acceptance ratio */
        double log_prior_ratio = model_ptr->log_prior_density(theta_star)
        - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
        //curently only symmetric proposal densities are supported
        double log_proposal_density_ratio = H_star_to_current-H_current_to_star ; // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
        
        acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
        
        
        if (std::isnan(acceptance_ratio)) {
            acceptance_prob = 0.0;
            div++;
            //  throw std::runtime_error("Error evaluating acceptance ratio.");
        } else {
            acceptance_prob = std::min(1., acceptance_ratio);
            // - draw a random number between 0 and 1
            // - accept or reject by comparing the random number to the
            // Metropolis-Hastings ratio (acceptance probability); if accept,
            // change the current value of theta to the proposed theta,
            // update the current value of the target and keep track of the
            // number of accepted proposals
            if (rng.rand() < acceptance_prob ) {
                
                theta             = theta_star;
                marginal_llik     = ll_star;
                joint_llik        = joint_Lik_star;
                
                
            }
        }
    } // end out of bounds if
    
}


/** MMALA implementation**/
const void MCMC_gibbs_sMMALA::MMALA_sampler(Eigen::VectorXd& theta,
                                          const Eigen::MatrixXd& X_sample,
                                          const double& epsilon,
                                          double& acceptance_prob,
                                          double& joint_llik,
                                          double& marginal_llik,
                                          unsigned int& div,
                                          unsigned int iter) {
    size_t number_of_i_parameters = parameters[MCMC_gibbs_sMMALA::PARAM_hmc_theta_dim];
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd theta_current(theta_dim) , theta_star(theta_dim) ;
    
    Eigen::VectorXd mu(number_of_i_parameters);
    Eigen::VectorXd grad_term(number_of_i_parameters);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters);
    std::vector<Eigen::MatrixXd> G_k = {Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters), Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters), Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters)};
    Eigen::MatrixXd G_inv(number_of_i_parameters,number_of_i_parameters), sigma(number_of_i_parameters,number_of_i_parameters);
    double  H_star_to_current = 0.0, H_current_to_star = 0.0, joint_Lik_star = 0.0 , ll_star = 0.0 , L = 0.0, acceptance_ratio= 0.0;
    bool  cond = false;
    
    // - current parameter vector
    theta_current = theta;
    theta_star = theta_current;
    // - set up model parameters
    model_ptr->update_i_parameters(theta_current);
    // - gather info to calculate proposal
    // - 1. grad
    
    
    grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_llik, &marginal_llik, &G, &G_k);//
    grad_term += (model_ptr->grad_log_prior_density(theta_current));
    model_ptr->add_prior_information(theta_current,&G, &G_k);
    G_inv =  G.inverse();
    //G_inv = (G_inv+G_inv.transpose())/2 + 1e-3*Eigen::MatrixXd::Identity(3,3);
    mu = theta_current.topRows(number_of_i_parameters) + pow(epsilon,2)*G_inv*grad_term/2. ;
    for (size_t k=0;k<number_of_i_parameters;k++) {
        mu -= pow(epsilon,2) * (G_inv*G_k[k])*G_inv.col(k);
        mu += pow(epsilon,2)/2. * G_inv.col(k)*(G_inv*G_k[k]).trace();
    }
    
    
    Eigen::LLT<Eigen::MatrixXd> cholSolver( pow(epsilon,2)*G_inv );
    if (cholSolver.info()==Eigen::Success) {
        // Use cholesky solver
        sigma = cholSolver.matrixL();
    } else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver( pow(epsilon,2)*G_inv);
        sigma = eigenSolver.eigenvectors()
        * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    
    
    
    // - sample proposal
    theta_star.topRows(number_of_i_parameters)  = rng.rand_mvn(mu , sigma);
    
    
    
    //if proposal within bounds
    cond = (theta_star.array() < this->lower_limit.array()).any() || (theta_star.array() > this->upper_limit.array()).any();
    acceptance_prob = 0.0;
    
    if ( ! cond  ) {
        
        H_current_to_star= -1./(2.*pow(epsilon,2))*(theta_star.topRows(number_of_i_parameters)-mu).transpose()*G*(theta_star.topRows(number_of_i_parameters)-mu) - 1./2.*log(G_inv.determinant());
        
        /* std::cout<< "--------------- "<<std::endl;
         std::cout<< sigma  << std::endl << "--" << H_current_to_star <<std::endl;
         std::cout<< G_inv  << std::endl << "--"<< epsilon <<std::endl;
         std::cout << mu.transpose() <<std::endl;
         std::cout<< grad_term.transpose() <<std::endl;
         std::cout<<  theta_star.transpose() <<std::endl;
         std::cout<<  theta_current.transpose() <<std::endl;*/
        
        // gather info to form proposal
        model_ptr->update_i_parameters(theta_star);
        grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_Lik_star, &ll_star, &G, &G_k);//
        grad_term += (model_ptr->grad_log_prior_density(theta_star));
        model_ptr->add_prior_information(theta_star,&G, &G_k);
        
        G_inv   = G.inverse();
        mu = theta_star.topRows(number_of_i_parameters) + pow(epsilon,2)*G_inv*grad_term/2.;
        for (size_t k=0;k<number_of_i_parameters;k++) {
            mu -= pow(epsilon,2) * (G_inv*G_k[k])*G_inv.col(k);
            mu += pow(epsilon,2)/2. * G_inv.col(k)*(G_inv*G_k[k]).trace();
        }
        
        H_star_to_current= -1./(2.*pow(epsilon,2))*(theta_current.topRows(number_of_i_parameters)-mu).transpose()*G*(theta_current.topRows(number_of_i_parameters)-mu)
        + 1./2.*log(G.determinant());
        
        L = - joint_llik + joint_Lik_star;
        /* calculate the acceptance ratio */
        double log_prior_ratio = model_ptr->log_prior_density(theta_star)
        - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
        //curently only symmetric proposal densities are supported
        double log_proposal_density_ratio = H_star_to_current-H_current_to_star ; // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
        
        acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
        
        
        if (std::isnan(acceptance_ratio)) {
            acceptance_prob = 0.0;
            div++;
            //  throw std::runtime_error("Error evaluating acceptance ratio.");
        } else {
            acceptance_prob = std::min(1., acceptance_ratio);
            // - draw a random number between 0 and 1
            // - accept or reject by comparing the random number to the
            // Metropolis-Hastings ratio (acceptance probability); if accept,
            // change the current value of theta to the proposed theta,
            // update the current value of the target and keep track of the
            // number of accepted proposals
            if (rng.rand() < acceptance_prob ) {
                
                theta             = theta_star;
                marginal_llik     = ll_star;
                joint_llik        = joint_Lik_star;
                
                
            }
        }
    } // end out of bounds if
    
}

/** function HMC implementation**/
const void MCMC_gibbs_sMMALA::hmc_sampler(Eigen::VectorXd& theta,
                                       const Eigen::MatrixXd& X_sample,
                                       const double& epsilon,
                                       double& acceptance_prob,
                                       double& joint_llik,
                                       double& marginal_llik,
                                       unsigned int& div) {
    
    size_t number_of_i_parameters = parameters[MCMC_gibbs_sMMALA::PARAM_hmc_theta_dim];;
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd p_current(theta_dim), p_star(theta_dim), theta_current(theta_dim), theta_star(theta_dim);
    Eigen::VectorXd p_zero = Eigen::VectorXd::Zero(theta_dim);
    double H_current, H_star, joint_Lik_c, joint_Lik_star, ll_c, ll_star, L, acceptance_ratio;
    bool  cond = false;
    
    theta_current = theta;
    // - sample momentum vector
    p_current = p_zero;
    p_current.topRows(number_of_i_parameters)     = rng.rand_mvn(p_zero.topRows(number_of_i_parameters),
                                                 M_normTransform.block(0, 0, number_of_i_parameters, number_of_i_parameters));
    
    
    // - calculate contribution of the momentum vector to the hamiltonian
    H_current = 1./2.*p_current.transpose()*mass_inv*p_current;
    
    /** - Evaluate log prob density of observed and latent variables given parameters (joint_Lik)
     Evaluate likelihood of observed variables given parameters and latent variables (ll_c) **/
    model_ptr->update_i_parameters(theta_current);
    joint_Lik_c = model_ptr->log_joint_density(X_sample, data, &ll_c);
    //joint_Lik_c = joint_llik;
    //ll_c = marginal_llik;
    
    //make proposal using leapfrog (X_sample or X_sample_1)
    p_star = p_current;
    theta_star = theta_current;
    leapfrog(theta_star, p_star, X_sample, epsilon, &cond);
    
    if (! cond) {
        //  - update model with theta_star and calculate the various terms
        model_ptr->update_i_parameters(theta_star);
        joint_Lik_star = model_ptr->log_joint_density(X_sample, data, &ll_star);
        H_star = 1./2.*p_star.transpose()*mass_inv*p_star;
    
    
        // - calculate the acceptance ratio
        L = - joint_Lik_c + joint_Lik_star;
        acceptance_ratio = exp( L +
                                       model_ptr->log_prior_density(theta_star) - model_ptr->log_prior_density(theta_current)
                                       + H_current - H_star);
    
    /* std::cout<< "----------------- " << i <<" ----------------- " <<std::endl;
     std::cout<<"theta_c:" <<theta_current.transpose() <<std::endl;
     std::cout<<"theta_s:" <<theta_star.transpose() <<std::endl;
     std::cout<<"a_ratio Hc-Hs:" <<H_current<< " - " << H_star << " = " << H_current - H_star <<std::endl;
     std::cout<<"a_ratio -L_c+L_s:" << - joint_Lik_c <<" + " << joint_Lik_star<< " = "<<L<< " - "  <<std::endl;
     std::cout<<"a_ratio prior:" <<model_ptr->log_prior_density(theta_star)<< " - " << model_ptr->log_prior_density(theta_current) << " = " << model_ptr->log_prior_density(theta_star) - model_ptr->log_prior_density(theta_current) <<std::endl;
     std::cout<<"a_ratio cond:" << !cond << std::endl;
     std::cout<<"a_ratio f:" << acceptance_ratio << std::endl;
     std::cout<<"ll " << ll <<std::endl;*/
    
    
    
    
        if (std::isnan(acceptance_ratio)) { // divergence
            acceptance_prob = 0;
            div += 1;
            std::cout<< " D -- " << theta_star.transpose() <<std::endl;
            marginal_llik  = ll_c;
            joint_llik = joint_Lik_c;
            
        } else { //calculate the acceptance prob
            acceptance_prob = std::min(1., acceptance_ratio);
            
            
            if (rng.rand() < acceptance_prob) {
                //model_ptr->log_joint_density(X_sample, data, &ll);
                //std::cout<<"--- | "<<theta_star.transpose()<<std::endl;
                theta = theta_star;
                marginal_llik  = ll_star;
                joint_llik = joint_Lik_star;
                
                
            } else {
                marginal_llik  = ll_c;
                joint_llik = joint_Lik_c;
            }
        }
    } else {
        if (isnan(theta_star.array()).any())
            div += 1;
        //std::cout<< " OoB -- " << theta_star.transpose() <<std::endl;
        acceptance_prob = 0.0;
        marginal_llik  = ll_c;
        joint_llik = joint_Lik_c;
    }

}


/** run the algorithm using parameters in params and store results in res **/
const void MCMC_gibbs_sMMALA::run(const MCMC_gibbs_hmc_parameter_t& p, chain_doubles_t& ch) {}


const void MCMC_gibbs_sMMALA::print_chain_parameters() {
    std::cout << this->chain_res << std::endl;
}

const void MCMC_gibbs_sMMALA::AIS_sampler(Eigen::VectorXd& theta,
                                 Eigen::MatrixXd& X_sample,
                                 const double& epsilon,
                                 double& acceptance_prob,
                                 double& joint_llik,
                                 double& marginal_llik,
                                 const size_t& start_block_i,
                                 const size_t& size_block_i) {
   
   
    
    Eigen::MatrixXd& X_sample_1 = X_sample;
    Eigen::MatrixXd& X_sample_2 = X_sample;
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd theta_current(theta_dim), theta_star(theta_dim), theta_tmp(theta_dim);
    double L1, L2;
    double ll_star =0.0;
    double L, acceptance_ratio;
    bool  cond = false;
    int   K = 5;
    double zz ;
   
    
    theta_current = theta;
    // - draw a new sample vector for the current block of parameters using the proposal distribution
    theta_star = theta_current;
    
    theta_star.middleRows(start_block_i,size_block_i) = rng.rand_mvn(
                                                          theta_current.middleRows(start_block_i, size_block_i),
                                                          epsilon*normTransform.block(start_block_i, start_block_i, size_block_i, size_block_i)
                                                          );
    
    
    //if proposal within bounds
    cond = (theta_star.array() < this->lower_limit.array()).any() || (theta_star.array() > this->upper_limit.array()).any();
    acceptance_prob = 0.0;
    
    
    L = 0.0 ;
    if ( ! cond  ) {
        
        
        // - evaluate first pair in log-likelihood ratio L
        model_ptr->update_i_parameters( theta_current );
        cSMC_AS_ptr->run(X_sample, X_sample_1);
        //X_sample = X_sample_1;
        L1 = model_ptr->log_joint_density(X_sample_1, data);
        
        
        // update parameter vector
        zz                    = 1./(K + 1.);
        theta_tmp = theta_current*(1-zz) + theta_star*zz;
        model_ptr->update_i_parameters(theta_tmp);
        L2 = model_ptr->log_joint_density(X_sample_1, data);
        
        L = - L1 + L2;

        for (size_t k = 1; k<K+1; k++) {
            // - calculate theta_k
            //zz                    = (k)/(K + 1.);
            //theta_tmp = theta_current*(1-zz) + theta_star*zz;
            
            
            // - sample u_k
            //model_ptr->update_i_parameters(theta_tmp);
            cSMC_AS_ptr->run(X_sample_1, X_sample_2);
            X_sample_1 = X_sample_2;
            
            
            // - calculate kth term in L ratio
            L1 = model_ptr->log_joint_density(X_sample_1, data);
            
            zz                    = (k+1.)/(K + 1.);
            theta_tmp = theta_current*(1.-zz) + theta_star*zz;
            model_ptr->update_i_parameters(theta_tmp);
            L2 = model_ptr->log_joint_density(X_sample_1, data, &ll_star);
            
            L = L - L1 + L2;
            
            //make current trajectory the ref trajectory for the next k
            
            
            
        }
        
        /* calculate the acceptance ratio */
        double log_prior_ratio = model_ptr->log_prior_density(theta_star)
        - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
        //curently only symmetric proposal densities are supported
        double log_proposal_density_ratio = model_ptr->log_proposal_density_ratio( theta_star, theta_current); // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
        
        acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
        
        if (std::isnan(acceptance_ratio)) {
            acceptance_prob = 0.0;
            //  throw std::runtime_error("Error evaluating acceptance ratio.");
        } else {
            acceptance_prob = std::min(1., acceptance_ratio);
            // - draw a random number between 0 and 1
            // - accept or reject by comparing the random number to the
            // Metropolis-Hastings ratio (acceptance probability); if accept,
            // change the current value of theta to the proposed theta,
            // update the current value of the target and keep track of the
            // number of accepted proposals
            if (rng.rand() < acceptance_prob) {
                
                theta             = theta_star;
                marginal_llik     = ll_star;
                joint_llik        = L2;
                X_sample = X_sample_1;
                
                
                //std::cout << "Accepted " << marginal_llik << " " << L << std::endl;
                
            }
        }
    }
    
    
    
    
    
    
}



double MCMC_gibbs_sMMALA::find_epsilon(const Eigen::VectorXd& theta, const Eigen::MatrixXd& X_sample) {
    double epsilon = 1;
    double ratio = 0.0;
    double a = 0.0;
    bool cond;
    double ll =0.0;
    double ll_c = 0;
    double L = 0.0;
    double L1 = 0.0;
    double L2 = 0.0;
    double H_current = 0.0;
    double H_star = 0.0;
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    size_t steps         =  parameters[MCMC_gibbs_sMMALA::PARAM_steps];
    Eigen::MatrixXd X_sample_1 = X_sample;
    Eigen::VectorXd p_current(theta_dim);
    Eigen::VectorXd p_star(theta_dim);
    Eigen::VectorXd theta_star(theta_dim);
    Eigen::MatrixXd theta_seq(steps+1, theta_dim);
    
    
    
    p_current = rng.rand_mvn(Eigen::VectorXd::Zero(theta_dim), M_normTransform);
    p_current(3) = 0.;
    p_current(4) = 0.;
    
    p_star = p_current;
    theta_star = theta;
    
    // update parameters
    model_ptr->update_i_parameters(theta_star);
    L1 = model_ptr->log_joint_density(X_sample, data, &ll_c);
    
    //theta_seq = leapfrog(theta_star, p_star, X_sample, epsilon, &cond);
    leapfrog(theta_star, p_star, X_sample, epsilon, &cond);
    
    if (!cond) {
        // update parameters
        model_ptr->update_i_parameters(theta_star);
        L2 = model_ptr->log_joint_density(X_sample, data, &ll);
        
        L = - L1 + L2;
        H_star = 1./2.*p_star.transpose()*mass_inv*p_star;
        H_current = 1./2.*p_current.transpose()*mass_inv*p_current;
        ratio = exp( L
                    + model_ptr->log_prior_density(theta_star) - model_ptr->log_prior_density(theta)
                    + H_current - H_star);
        
        if (std::isnan(ratio)) {
            ratio = 0;
        }
        
    } else {
        ratio = 0;
    }
        
    
    
    a = (ratio>0.5) ? 1 : -1;
    
    
    while (pow(ratio,a) > pow(2,-a)) {
        p_star = p_current;
        theta_star = theta;
        epsilon = pow(2,a)*epsilon;
        
        //theta_seq = leapfrog(theta_star, p_star, X_sample, epsilon, &cond);
        leapfrog(theta_star, p_star, X_sample, epsilon, &cond);
        
        // if theta within bounds calculate acceptance ratio
        if (!cond) {
            // update parameters
            model_ptr->update_i_parameters(theta_star);
            L2 = model_ptr->log_joint_density(X_sample, data, &ll);
            
            L = -L1 + L2;
            H_star = 1./2.*p_star.transpose()*mass_inv*p_star;
            H_current = 1./2.*p_current.transpose()*mass_inv*p_current;
            ratio = exp( L
                        + model_ptr->log_prior_density(theta_star) - model_ptr->log_prior_density(theta)
                        + H_current - H_star);
            
            if (std::isnan(ratio)) {
                ratio = 0.;
            }
            
            
        } else {
            ratio = 0.;
        }
      
        
        
    }
    
    return epsilon;
}



const void MCMC_gibbs_sMMALA::print_chain_parameters(std::string str_file) {
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
     
    std::ofstream file_tmp;
    file_tmp.open(str_file);
    file_tmp << this->chain_res.format(CSVFormat) << std::endl;
    file_tmp.close();
    
}


const void MCMC_gibbs_sMMALA::print_chain_likelihood(std::string str_file) {
    
    std::ofstream file_tmp;
    file_tmp.open(str_file);
    file_tmp << this->likelihood_res << std::endl;
    file_tmp.close();
    
}


const void MCMC_gibbs_sMMALA::print_chain_trajectories(std::string str_file) {
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    
    
    std::ofstream file_tmp;
    file_tmp.open(str_file);
    for (size_t mm = 0; mm< this->trajectories_res.size(); mm++) {
        file_tmp << this->trajectories_res[mm].format(CSVFormat) << std::endl;
    }
    file_tmp.close();
    
}


/** simplified MMALA implementation **/
const void MCMC_gibbs_sMMALA::sMMALA_sampler_block(Eigen::VectorXd& theta,
                                          const Eigen::MatrixXd& X_sample,
                                          const double& epsilon,
                                          double& acceptance_prob,
                                          double& joint_llik,
                                          double& marginal_llik,
                                          unsigned int& div,
                                          const size_t& start_index,
                                          const size_t& block_size) {
    
    size_t number_of_i_parameters = parameters[PARAM_hmc_theta_dim];
    size_t theta_dim =  parameters[PARAM_theta_dim];
    Eigen::VectorXd theta_current(theta_dim) , theta_star(theta_dim) ;
    Eigen::MatrixXd G_t = Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters);
    
    
    Eigen::VectorXd mu(block_size);
    Eigen::VectorXd grad_term(block_size);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(block_size,block_size);
    Eigen::MatrixXd G_inv(block_size,block_size), sigma(block_size,block_size);
    
    // aux variables
    double H_star_to_current = 0.0, H_current_to_star = 0.0, joint_Lik_star = 0.0 , ll_star = 0.0 , L = 0.0, acceptance_ratio= 0.0;
    bool  cond = false;
    
    // - initialise current parameter vector and proposal
    theta_current = theta;
    theta_star = theta_current;
    
    // - set up model parameters
    model_ptr->update_i_parameters(theta_current);
    // - gather info to calculate proposal
    // - 1. gradient of log target density
    grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_llik, &marginal_llik, &G_t).segment(start_index,block_size);//
    grad_term += (model_ptr->grad_log_prior_density(theta_current)).segment(start_index,block_size);
    model_ptr->add_prior_information(theta_current,&G_t);
    G = G_t.block(start_index,start_index,block_size,block_size);
    
    G_inv =  G.inverse();
    
    mu = theta_current.segment(start_index,block_size) + pow(epsilon,2)*G_inv*grad_term/2.;
    
    Eigen::LLT<Eigen::MatrixXd> cholSolver( pow(epsilon,2)*G_inv );
    if (cholSolver.info()==Eigen::Success) {
        // Use cholesky solver
        sigma = cholSolver.matrixL();
    } else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver( pow(epsilon,2)*G_inv);
        sigma = eigenSolver.eigenvectors()
        * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    
    // - sample proposal
    theta_star.segment(start_index,block_size)  = rng.rand_mvn(mu , sigma);
    
    //check if proposal within bounds
    cond = (theta_star.array() < this->lower_limit.array()).any() || (theta_star.array() > this->upper_limit.array()).any();
    
    acceptance_prob = 0.0;
    
    if ( ! cond  ) {
        
        H_current_to_star= -1./(2.*pow(epsilon,2))*(theta_star.segment(start_index,block_size)-mu).transpose()*G*(theta_star.segment(start_index,block_size)-mu) - 1./2.*log(G_inv.determinant());
        
        /* std::cout<< "--------------- "<<std::endl;
         std::cout<< sigma  << std::endl << "--" << H_current_to_star <<std::endl;
         std::cout<< G_inv  << std::endl << "--"<< epsilon <<std::endl;
         std::cout << mu.transpose() <<std::endl;
         std::cout<< grad_term.transpose() <<std::endl;
         std::cout<<  theta_star.transpose() <<std::endl;
         std::cout<<  theta_current.transpose() <<std::endl;*/
        
        // gather info to form proposal
        model_ptr->update_i_parameters(theta_star);
        grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_Lik_star, &ll_star, &G_t).segment(start_index,block_size);//
        grad_term += model_ptr->grad_log_prior_density(theta_star).segment(start_index,block_size);
        model_ptr->add_prior_information(theta_star,&G_t);
        G = G_t.block(start_index,start_index,block_size,block_size);
        
        G_inv   = G.inverse();
        mu = theta_star.segment(start_index,block_size) + pow(epsilon,2)*G_inv*grad_term/2.;
        H_star_to_current= -1./(2.*pow(epsilon,2))*(theta_current.segment(start_index,block_size)-mu).transpose()*G*(theta_current.segment(start_index,block_size)-mu)
        + 1./2.*log(G.determinant());
        
        L = - joint_llik + joint_Lik_star;
        /* calculate the acceptance ratio */
        double log_prior_ratio = model_ptr->log_prior_density(theta_star)
        - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
        //curently only symmetric proposal densities are supported
        double log_proposal_density_ratio = H_star_to_current-H_current_to_star ; // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
        
        acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
        
        
        if (std::isnan(acceptance_ratio)) {
            acceptance_prob = 0.0;
            div++;
            //  throw std::runtime_error("Error evaluating acceptance ratio.");
        } else {
            acceptance_prob = std::min(1., acceptance_ratio);
            // - draw a random number between 0 and 1
            // - accept or reject by comparing the random number to the
            // Metropolis-Hastings ratio (acceptance probability); if accept,
            // change the current value of theta to the proposed theta,
            // update the current value of the target and keep track of the
            // number of accepted proposals
            if (rng.rand() < acceptance_prob ) {
                
                theta             = theta_star;
                marginal_llik     = ll_star;
                joint_llik        = joint_Lik_star;
                
                
            }
        }
    } // end out of bounds if
    
}



/** MMALA implementation**/
const void MCMC_gibbs_sMMALA::MMALA_sampler_block(Eigen::VectorXd& theta,
                                         const Eigen::MatrixXd& X_sample,
                                         const double& epsilon,
                                         double& acceptance_prob,
                                         double& joint_llik,
                                         double& marginal_llik,
                                         unsigned int& div,
                                         const size_t& start_index,
                                         const size_t& block_size) {
    
    size_t number_of_i_parameters = parameters[MCMC_gibbs_sMMALA::PARAM_hmc_theta_dim];
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd theta_current(theta_dim) , theta_star(theta_dim) ;
    Eigen::MatrixXd G_t = Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters);
    std::vector<Eigen::MatrixXd> G_k = {Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters), Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters), Eigen::MatrixXd::Zero(number_of_i_parameters,number_of_i_parameters)};
    
    Eigen::VectorXd mu(block_size);
    Eigen::VectorXd grad_term(block_size);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(block_size,block_size);
    Eigen::MatrixXd G_inv(block_size,block_size), sigma(block_size,block_size);
    Eigen::MatrixXd tmp_G_k = Eigen::MatrixXd::Zero(block_size,block_size);
    
    
    
    
    
    
    double  H_star_to_current = 0.0, H_current_to_star = 0.0, joint_Lik_star = 0.0 , ll_star = 0.0 , L = 0.0, acceptance_ratio= 0.0;
    bool  cond = false;
    
    // - current parameter vector
    theta_current = theta;
    theta_star = theta_current;
    // - set up model parameters
    model_ptr->update_i_parameters(theta_current);
    // - gather info to calculate proposal
    // - 1. grad
    
    
    grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_llik, &marginal_llik, &G_t, &G_k).segment(start_index,block_size);//
    grad_term += (model_ptr->grad_log_prior_density(theta_current)).segment(start_index,block_size);
    model_ptr->add_prior_information(theta_current,&G_t, &G_k);
    G = G_t.block(start_index,start_index,block_size,block_size);
    
    G_inv =  G.inverse();
    //G_inv = (G_inv+G_inv.transpose())/2 + 1e-3*Eigen::MatrixXd::Identity(3,3);
    mu = theta_current.segment(start_index,block_size) + pow(epsilon,2)*G_inv*grad_term/2. ;
    for (size_t k=0;k<block_size;k++) {
        tmp_G_k = G_k[start_index+k].block(start_index,start_index,block_size,block_size);
        mu -= pow(epsilon,2) * (G_inv*tmp_G_k)*G_inv.col(k);
        mu += pow(epsilon,2)/2. * G_inv.col(k)*(G_inv*tmp_G_k).trace();
    }
    
    
    Eigen::LLT<Eigen::MatrixXd> cholSolver( pow(epsilon,2)*G_inv );
    if (cholSolver.info()==Eigen::Success) {
        // Use cholesky solver
        sigma = cholSolver.matrixL();
    } else {
        // Use eigen solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver( pow(epsilon,2)*G_inv);
        sigma = eigenSolver.eigenvectors()
        * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }
    
    
    
    // - sample proposal
    theta_star.segment(start_index,block_size)   = rng.rand_mvn(mu , sigma);
    
    
    
    //if proposal within bounds
    cond = (theta_star.array() < this->lower_limit.array()).any() || (theta_star.array() > this->upper_limit.array()).any();
    acceptance_prob = 0.0;
    
    if ( ! cond  ) {
        
        H_current_to_star= -1./(2.*pow(epsilon,2))*(theta_star.segment(start_index,block_size)-mu).transpose()*G*(theta_star.segment(start_index,block_size)-mu) - 1./2.*log(G_inv.determinant());
        
        /* std::cout<< "--------------- "<<std::endl;
         std::cout<< sigma  << std::endl << "--" << H_current_to_star <<std::endl;
         std::cout<< G_inv  << std::endl << "--"<< epsilon <<std::endl;
         std::cout << mu.transpose() <<std::endl;
         std::cout<< grad_term.transpose() <<std::endl;
         std::cout<<  theta_star.transpose() <<std::endl;
         std::cout<<  theta_current.transpose() <<std::endl;*/
        
        // gather info to form proposal
        model_ptr->update_i_parameters(theta_star);
        grad_term = model_ptr->grad_log_joint_density(X_sample, data, &joint_Lik_star, &ll_star, &G_t, &G_k).segment(start_index,block_size);//
        grad_term += model_ptr->grad_log_prior_density(theta_star).segment(start_index,block_size);
        model_ptr->add_prior_information(theta_star,&G_t, &G_k);
        G = G_t.block(start_index,start_index,block_size,block_size);
        
        G_inv   = G.inverse();
        mu = theta_star.segment(start_index,block_size)  + pow(epsilon,2)*G_inv*grad_term/2.;
        
        for (size_t k=0;k<block_size;k++) {
            tmp_G_k = G_k[start_index+k].block(start_index,start_index,block_size,block_size);
            mu -= pow(epsilon,2) * (G_inv*tmp_G_k)*G_inv.col(k);
            mu += pow(epsilon,2)/2. * G_inv.col(k)*(G_inv*tmp_G_k).trace();
        }
        
        H_star_to_current= -1./(2.*pow(epsilon,2))*(theta_current.segment(start_index,block_size)-mu).transpose()*G*(theta_current.segment(start_index,block_size)-mu)
        + 1./2.*log(G.determinant());
        
        L = - joint_llik + joint_Lik_star;
        /* calculate the acceptance ratio */
        double log_prior_ratio = model_ptr->log_prior_density(theta_star)
        - model_ptr->log_prior_density(theta_current); // log( pi(theta_proposal) / pi(theta_current) )
        //curently only symmetric proposal densities are supported
        double log_proposal_density_ratio = H_star_to_current-H_current_to_star ; // p(theta_proposal -> theta_current)/p(theta_current -> theta_proposal)
        
        acceptance_ratio = exp(L + log_prior_ratio + log_proposal_density_ratio);
        
        
        if (std::isnan(acceptance_ratio)) {
            acceptance_prob = 0.0;
            div++;
            //  throw std::runtime_error("Error evaluating acceptance ratio.");
        } else {
            acceptance_prob = std::min(1., acceptance_ratio);
            // - draw a random number between 0 and 1
            // - accept or reject by comparing the random number to the
            // Metropolis-Hastings ratio (acceptance probability); if accept,
            // change the current value of theta to the proposed theta,
            // update the current value of the target and keep track of the
            // number of accepted proposals
            if (rng.rand() < acceptance_prob ) {
                
                theta             = theta_star;
                marginal_llik     = ll_star;
                joint_llik        = joint_Lik_star;
                
                
            }
        }
    } // end out of bounds if
    
}


/** HMC implementation always accepting proposal if **/
const void MCMC_gibbs_sMMALA::hmc_sampler_warmup(Eigen::VectorXd& theta,
                                       const Eigen::MatrixXd& X_sample,
                                       const double& epsilon,
                                       double& acceptance_prob,
                                       double& joint_llik,
                                       double& marginal_llik,
                                       unsigned int& div) {
    
    size_t number_of_i_parameters = parameters[MCMC_gibbs_sMMALA::PARAM_hmc_theta_dim];;
    size_t theta_dim =  parameters[MCMC_gibbs_sMMALA::PARAM_theta_dim];
    Eigen::VectorXd p_current(theta_dim), p_star(theta_dim), theta_current(theta_dim), theta_star(theta_dim);
    Eigen::VectorXd p_zero = Eigen::VectorXd::Zero(theta_dim);
    double H_current, H_star, joint_Lik_c, joint_Lik_star, ll_c, ll_star, L, acceptance_ratio;
    bool  cond = false;
    
    theta_current = theta;
    // - sample momentum vector
    p_current = p_zero;
    p_current.topRows(number_of_i_parameters)     = rng.rand_mvn(p_zero.topRows(number_of_i_parameters),
                                                                 M_normTransform.block(0, 0, number_of_i_parameters, number_of_i_parameters));
    
    
    // - calculate contribution of the momentum vector to the hamiltonian
    H_current = 1./2.*p_current.transpose()*mass_inv*p_current;
    
    /** - Evaluate log prob density of observed and latent variables given parameters (joint_Lik)
     Evaluate likelihood of observed variables given parameters and latent variables (ll_c) **/
    model_ptr->update_i_parameters(theta_current);
    joint_Lik_c = model_ptr->log_joint_density(X_sample, data, &ll_c);
    //joint_Lik_c = joint_llik;
    //ll_c = marginal_llik;
    
    //make proposal using leapfrog (X_sample or X_sample_1)
    p_star = p_current;
    theta_star = theta_current;
    leapfrog(theta_star, p_star, X_sample, epsilon, &cond);
    
    if (! cond) {
        //  - update model with theta_star and calculate the various terms
        model_ptr->update_i_parameters(theta_star);
        joint_Lik_star = model_ptr->log_joint_density(X_sample, data, &ll_star);
        H_star = 1./2.*p_star.transpose()*mass_inv*p_star;
        
        
        // - calculate the acceptance ratio
        L = - joint_Lik_c + joint_Lik_star;
        acceptance_ratio = ( L +
                               model_ptr->log_prior_density(theta_star) - model_ptr->log_prior_density(theta_current)
                           );
        
        
        
        if (std::isnan(acceptance_ratio)) { // divergence
            acceptance_prob = 0;
            div += 1;
            std::cout<< " D -- " << theta_star.transpose() <<std::endl;
            marginal_llik  = ll_c;
            joint_llik = joint_Lik_c;
            
        } else { //calculate the acceptance prob
            //acceptance_prob = std::min(1., acceptance_ratio);
            
            
            if (acceptance_ratio > 0) {
                
                theta = theta_star;
                marginal_llik  = ll_star;
                joint_llik = joint_Lik_star;
                
                
            } else {
                marginal_llik  = ll_c;
                joint_llik = joint_Lik_c;
            }
        }
    } else {
        if (isnan(theta_star.array()).any())
            div += 1;
        //std::cout<< " OoB -- " << theta_star.transpose() <<std::endl;
        acceptance_prob = 0.0;
        marginal_llik  = ll_c;
        joint_llik = joint_Lik_c;
    }
    
}

