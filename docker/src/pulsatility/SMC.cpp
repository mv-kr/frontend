//
//  SMC.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 24/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//


#include "SMC.hpp"

using namespace mci;
using namespace mci::algorithms;


const std::string mci::algorithms::SMC::PARAM_nparticles1 = "nparticles";
const std::string mci::algorithms::SMC::NAME1 = "SMC";

SMC::SMC(mci::models::basic_model* m, const Eigen::MatrixXd& d, const stochastic& r) {
    
    model_obj = m;
    data = d;
    rng = r;
    
}

/** run the algorithm and store results in res_traj **/
const void SMC::run(trajectory_t& res_traj) {
    
    // number of timepoints in dataset
    size_t obs_timepoints = data.cols();
    // start time
    double t_0   = data(data.rows()-1,  0);
    // end time
    double t_end = data(data.rows()-1,  data.cols()-1);
    // number of timepoints considered while running the model;
    size_t model_timepoints = round( (t_end - t_0)/ (model_obj->get_dt() ) ) + 1;
    // number of particles
    unsigned int nparticles = parameters[SMC::PARAM_nparticles1];
    size_t model_timesteps_per_observation_period = round((data(data.rows()-1, 1)-data(data.rows()-1, 0))/model_obj->get_dt());
    /*** Initialisation of the algorithm  ***/
    //std::cout <<" --- SMC init --- " << std::endl;

    std::vector<Eigen::MatrixXd> particles(nparticles) ; // particles
    Eigen::MatrixXd A(nparticles, obs_timepoints) ; // anscestral indices
    Eigen::MatrixXd weights(nparticles, obs_timepoints) ; //weights
    Eigen::VectorXd WW; // tmp vector for weigths in each iteration
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> model_times_ind(obs_timepoints,1); //keeping track of the timepoints
    double next_time, current_time;
    size_t model_t_index;
    Eigen::VectorXd obs_value;
    
    
    
    /***  First timepoint (t = t_0) ***/
    current_time  = data(data.rows()-1, 0);
    obs_value    = data.col(0).head(data.rows()-1);
    
    model_times_ind(0) = 0;
    
    // initialise particles
    for (unsigned int pi=0;pi<nparticles; pi++) {
        particles[pi].resize(model_obj->get_dim(), model_timepoints);
        particles[pi].col(0) = model_obj->generate_ic(current_time) ;
        weights(pi, 0) = model_obj->log_obs_density(particles[pi].col(0), obs_value);
        A(pi, 0) = pi;
    }
    
    
    /***  loop over timepoints ***/
    for (size_t ti = 1; ti<obs_timepoints; ti++ ) {
        //std::cout<<" --- SMC at timepoint " << ti << std::endl;
        
        //get observed value
        obs_value = data.col(ti).head(data.rows()-1);
        
        //set interval
        next_time = data(data.rows()-1, ti);
        
        //
        model_t_index = model_times_ind(ti-1);
        
        // get weights
        WW = weights.col(ti-1);
        WW.array() -= WW.maxCoeff();
        // for each particle
        Eigen::Matrix<size_t,Eigen::Dynamic, 1> tmp(nparticles,1);
        for (unsigned int pi = 0; pi< nparticles; pi++) {
            
            //sample anscestral particle X_ansc
            unsigned int ansc_i = rng.rand_sample(WW.array().exp()) ;
            mci::models::state_t X_ansc    = particles[ansc_i].col(model_t_index);
            // Propagate the particle from current observation time to the next one
            trajectory_t Xfinal = model_obj->run(X_ansc, current_time, next_time,
                                                 data.middleCols<2>(ti-1));
            
            // save trajectory
            particles[pi].middleCols(model_t_index+1, Xfinal.cols()) = Xfinal;
            
            // Weight the particle with the likelihood of the observed data point
            weights(pi, ti) = model_obj->log_obs_density(Xfinal.rightCols(1), obs_value);
            // update A matrix
            A(pi, ti) = ansc_i;
            tmp(pi) = Xfinal.cols();
        }//End for() loop over particles
        
        model_times_ind(ti)  = model_t_index + tmp(0) ;
        current_time = next_time;
        
    }
    
    
    //Compute and return the marginal log-likelihood
    //sum of the log of the mean of the weights at each observation time
    //marg_ll = (( ( (weights.array().exp() ).colwise().sum() ).array().log() ).rowwise().sum() )(0);
    
    
    
    // sample trajectory according to particle weights
    WW = weights.col(obs_timepoints-1);
    WW.array() -= WW.maxCoeff();
    
    res_traj.resize(model_obj->get_dim(), model_timepoints);
    unsigned int pi = rng.rand_sample(WW.array().exp()) ;
    int ti = (int)obs_timepoints-1;
    
    //std::cout<< particles[0].col(model_times_ind(ti)) << std::endl;
    res_traj.middleCols(model_times_ind(ti-1)+1,model_times_ind(ti)-model_times_ind(ti-1)) =
           particles[pi].middleCols(model_times_ind(ti-1)+1,model_times_ind(ti)-model_times_ind(ti-1));
    
    
    for (ti = (int)obs_timepoints-2; ti>=0; ti-- ) {
        pi = A(pi, ti+1);
        if (ti>0)
            res_traj.middleCols(model_times_ind(ti-1)+1,model_times_ind(ti)-model_times_ind(ti-1)) = particles[pi].middleCols(model_times_ind(ti-1)+1,model_times_ind(ti)-model_times_ind(ti-1));
        else
            res_traj.col(0) = particles[pi].col(0);
    }
    
    
    
    
}

/** run the algorithm using parameters in params and store results in res **/
const void SMC::run(const SMC_parameter_t& p, trajectory_t& res_traj) {
    
    this->set_parameters(p);
    this->run(res_traj);
    
}


