//
//  cSMC_AS.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 23/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

// Updates TODO list
// 1.In the resampling step: simulate model in [t_current t_end] from rather than [t_0 t_end]


#include "cSMC_AS.hpp"


using namespace mci;
using namespace mci::algorithms;


const std::string mci::algorithms::cSMC_AS::PARAM_nparticles = "nparticles";
const std::string mci::algorithms::cSMC_AS::PARAM_DT = "DT";
const std::string mci::algorithms::cSMC_AS::PARAM_AS_probability = "AS_probability";
const std::string mci::algorithms::cSMC_AS::NAME = "cSMC_AS";


cSMC_AS::cSMC_AS(mci::models::basic_model* m, const Eigen::MatrixXd& d, const stochastic& r) {
    
    model_obj = m;
    data = d;
    rng = r;
    
}

/**  **/
const void cSMC_AS::run(trajectory_t& res_traj) {}
/**  **/
const void cSMC_AS::run(const cSMC_AS_parameter_t& p, trajectory_t& res_traj) {}



/** run the algorithm and store results in res_traj **/
const double cSMC_AS::run(const trajectory_t& ref_traj, trajectory_t& res_traj) {
    
    // number of timepoints in dataset
    size_t obs_timepoints = data.cols();
    // number of timepoints considered while running the model;
    size_t model_timepoints = ref_traj.cols();
    // number of particles
    unsigned int nparticles = parameters[cSMC_AS::PARAM_nparticles];
    // AS probability
    double as_prob = parameters[cSMC_AS::PARAM_AS_probability];
    
    
    /*** Initialisation of the algorithm  ***/
    //std::cout <<" --- SMC_AC init --- " << std::endl;
    
    std::vector<Eigen::MatrixXd> particles(nparticles) ; // particles
    std::vector<Eigen::MatrixXd> particles_tmp(nparticles) ; //
    std::vector<Eigen::MatrixXd> putative_states(nparticles) ; // AS
    //Eigen::MatrixXd A(nparticles, obs_timepoints) ; // anscestral indices
    Eigen::MatrixXd weights(nparticles, obs_timepoints) ; //weights
    Eigen::MatrixXd weights_tmp(nparticles, obs_timepoints) ; //weights
    Eigen::VectorXd WW(nparticles); // tmp vector for weigths in each iteration
    Eigen::VectorXd WW_as(nparticles); // vector for weigths in each AS step
    Eigen::VectorXd WW_as2(nparticles); // vector for weigths in each AS step
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> model_times_ind(obs_timepoints,1); //keeping track of the timepoints
    double next_time, current_time;
    size_t model_t_index;
    Eigen::VectorXd obs_values;
    //size_t model_timesteps_per_observation_period = round((data(data.rows()-1, 1)-data(data.rows()-1, 0))/model_obj->get_dt());
    
    weights.setZero(nparticles, obs_timepoints);
    weights_tmp.setZero(nparticles, obs_timepoints);
    
    /***  First timepoint (t = t_0; column index 0) ***/
    current_time  = data(data.rows()-1, 0);
    obs_values    = data.col(0).head(data.rows()-1) ;
    model_times_ind(0) = 0;
    
    // set the first particle to be the ref_trajerence particle
    particles[0].resize(model_obj->get_dim(), model_timepoints);
    
    trajectory_t ref_traj_update = model_obj->run(ref_traj, ref_traj.col(0), data);
    particles[0] = ref_traj_update;
    weights(0, 0) = model_obj->log_obs_density(particles[0].col(0), obs_values);
    
    // initialise remaining particles
    for (size_t pi=1;pi<nparticles; pi++) {
        particles[pi].resize(model_obj->get_dim(), model_timepoints);
        particles[pi].col(0) = model_obj->generate_ic(current_time) ;
        weights(pi, 0) = model_obj->log_obs_density(particles[pi].col(0), obs_values);
        //A(0, pi) = pi;
    }
    
  
    /***  loop over timepoints ***/
    for (size_t ti = 1; ti<obs_timepoints; ti++ ) {
        //std::cout<<" --- SMC-AC at timepoint " << ti << std::endl;
        
        //get observed value
        obs_values = data.col(ti).head(data.rows()-1) ;
        
        //set interval
        next_time = data(data.rows()-1, ti);
        
        //
        model_t_index = model_times_ind(ti-1);
        
        // get weights
        weights_tmp.setZero(nparticles, obs_timepoints);
        WW = weights.col(ti-1);
        WW.array() -= WW.maxCoeff();
        
        // for each particle
        Eigen::Matrix<size_t,Eigen::Dynamic, 1> tmp(nparticles,1);
        particles_tmp = particles;
        for (unsigned int pi = 1; pi< nparticles; pi++) {
            
            //sample anscestral particle X_ansc
            unsigned int ansc_i = rng.rand_sample(WW.array().exp()) ;
            mci::models::state_t X_ansc    = particles[ansc_i].col(model_t_index);
            // Propagate the particle from current observation time to the next one
            trajectory_t Xfinal = model_obj->run(X_ansc, current_time, next_time,
                                                 data.middleCols<2>(ti-1));
            
            
            
            // save trajectory[t_0 next_time] in tmp data_structure;
            particles_tmp[pi].leftCols(model_t_index+1) = particles[ansc_i].leftCols(model_t_index+1);
            //particles_tmp[pi].middleCols(model_t_index+1, Xfinal.cols()) = Xfinal;
            particles_tmp[pi].middleCols(model_t_index+1, Xfinal.cols()) = Xfinal;
            
            // Weight the particle with the likelihood of the observed data point
            weights_tmp.row(pi) = weights.row(ansc_i);
            weights_tmp(pi, ti) = model_obj->log_obs_density(Xfinal.rightCols<1>(), obs_values);
            
            // update A matrix
            //A(pi, ti) = ansc_i;
            tmp(pi) = Xfinal.cols();
        }//End for() loop
        
        
        /** Ancestral sampling for ref_trajectory particle **/
        if (rng.rand() <= as_prob) {
            
            WW_as = weights.col(ti-1);
            size_t end_t = obs_timepoints;
            
            
            for (unsigned int pi = 0; pi< nparticles; pi++) {
                // combine current particle trajectories [t_0 current_time]
                // with future reference trajectory [next_time t_end]
                trajectory_t X_ref = particles[0];
                X_ref.leftCols(model_t_index+1) = particles[pi].leftCols(model_t_index+1);
                
                // Run using reference trajecotry
                trajectory_t Xfinal = model_obj->run(X_ref, X_ref.col(0), data);
                                
                //store the
                //putative_states[pi] = Xfinal.leftCols(model_t_index+1+tmp(1));
                putative_states[pi] = Xfinal;
                
                // calculate weight
                double g_p;
                Eigen::VectorXd obs_values_2;
                double obs_time_2;
                size_t model_t_index_2;
                for (size_t tii = ti; tii<end_t; tii++   ) {
                
                    //set interval
                    obs_time_2 = data(data.rows()-1, tii);
                    obs_values_2 = data.col(tii).head(data.rows()-1) ;

                    model_t_index_2 = round(obs_time_2/model_obj->get_dt());
                    //size_t model_t_index_3 = tii*model_timesteps_per_observation_period; /////FIX THIS
                    
                    //calculate g
                    g_p = model_obj->log_obs_density(Xfinal.col(model_t_index_2), obs_values_2);
                    WW_as(pi) = WW_as(pi) +  g_p ;
                    
                    if (std::isnan(WW_as(pi)) ) {
                        //std::cout<< std::endl<<Xfinal.col(model_t_index_2).transpose() << std::endl;
                        //std::cout<<" "<< obs_values_2 << " " << g_p;
                    }
                    
                }
                
                WW_as(pi) = WW_as(pi) + model_obj->log_transition_density(Xfinal.col(model_t_index+1), Xfinal.col(model_t_index), model_t_index);
                
                if (std::isnan(WW_as(pi)) ) {
                    //std::cout<< std::endl<<model_obj->log_transition_density(Xfinal.col(model_t_index+1), Xfinal.col(model_t_index), model_t_index);
                }
                
            }
            

            // sample accoring to AS weight
            WW_as2 =WW_as;

            
            WW_as.array() -= WW_as.maxCoeff();
            unsigned int ancestral_ind = rng.rand_sample(WW_as.array().exp(), WW_as , WW_as2) ;
            //store
            //particles_tmp[0].leftCols(model_t_index+1+tmp(1)) = putative_states[ancestral_ind];
            particles_tmp[0] = putative_states[ancestral_ind];
            // Weight the particle with the likelihood of the observed data point
            
            double obs_time_2;
            Eigen::VectorXd obs_values_2 ;
            size_t model_t_index_2;
            for (size_t tii = 0; tii<ti; tii++   ) {
                obs_time_2 = data(data.rows()-1, tii);
                obs_values_2 = data.col(tii).head(data.rows()-1) ;
                model_t_index_2 = round(obs_time_2/model_obj->get_dt());
                //size_t model_t_index_2 = tii*model_timesteps_per_observation_period; /////FIX THIS
                
                weights_tmp(0, tii) = model_obj->log_obs_density(putative_states[ancestral_ind].col(model_t_index_2), obs_values_2);
            }
            
            /*if( abs(weights_tmp(0,ti-1) - weights(ancestral_ind, ti-1)) > 0.1 ) {
                std::cout << ti << " " << ancestral_ind << std::endl;
                
                trajectory_t Xfinal2 = model_obj->run(particles[ancestral_ind].leftCols(model_t_index+1), particles[ancestral_ind].col(0), data);
                std::cout <<particles[ancestral_ind].leftCols(model_t_index+1).row(6) <<std::endl;
                std::cout <<Xfinal2.row(6) <<std::endl;
                std::cout <<Xfinal2.row(14) <<std::endl;
            }*/
            
            weights_tmp.row(0) = weights.row(ancestral_ind);
            //weights_tmp(0, ti) = model_obj->log_obs_density(putative_states[ancestral_ind].rightCols<1>(), obs_values);
            weights_tmp(0, ti) = model_obj->log_obs_density(putative_states[ancestral_ind].col(model_t_index+tmp(1)), obs_values);
            
        } else {
            //
            trajectory_t X_ref = particles[0].leftCols(model_t_index+1+tmp(1));
            // Run with ref trajectory
            trajectory_t X_ref_updated = model_obj->run(X_ref, X_ref.col(0), data);
            
            //store
            particles_tmp[0].middleCols(model_t_index+1, tmp(1)) = X_ref_updated.middleCols(model_t_index+1, tmp(1));
            // Weight the particle with the likelihood of the observed data point
            weights_tmp.row(0) = weights.row(0);
            weights_tmp(0, ti) = model_obj->log_obs_density(X_ref_updated.col(model_t_index+tmp(1)), obs_values);
            
        }
        
        weights = weights_tmp;
        particles = particles_tmp;
        
        model_times_ind(ti)  = model_t_index + tmp(1) ;
        current_time = next_time;

         
    }
        
        
    //Compute and return the marginal log-likelihood
    //sum of the log of the mean of the weights at each observation time
    //double logc = (( ( (weights.array().exp() ).colwise().sum() ).array().log() ).rowwise().sum() )(0);
    
    // sample trajectory
    WW = weights.col(obs_timepoints-1);
    WW.array() -= WW.maxCoeff();
    unsigned int pi = rng.rand_sample(WW.array().exp()) ;
    //res_traj.resize(model_obj->get_dim(), model_timepoints);
    res_traj = particles[pi];
    
    //model quality control
    //std::cout<< "SMC " << model_obj->log_joint_density(res_traj, data);
    //std::cout<< " --" << weights.rowwise().sum()(pi) << std::endl;
    
    return weights.rowwise().sum()(pi);
}

/** run the algorithm using parameters in params and store results in res **/
const void cSMC_AS::run(const cSMC_AS_parameter_t& p, const trajectory_t& ref_traj, trajectory_t& res_traj) {
    
    //this->set_parameters(p);
    //this->run2(ref_traj, res_traj);
    
}


