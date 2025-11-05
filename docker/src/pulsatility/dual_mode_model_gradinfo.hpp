//
//  dual_mode_model_gradinfo.hpp
//  mci
//
//  Created by Margaritis Voliotis on 06/10/2020.
//  Copyright © 2020 Margaritis Voliotis. All rights reserved.
//

#ifndef dual_mode_model_gradinfo_hpp
#define dual_mode_model_gradinfo_hpp

#include <iostream>

#include "basic_model.hpp"
#include "stochastic.hpp"


#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/random/mersenne_twister.hpp>



using boost::math::pdf; // Probability mass function.
using boost::math::cdf; // Cumulative density function.

namespace mci {
    namespace models {
        
        
        class dual_mode_model_gradinfo  : public basic_model{
            
            
        private:
            
            /** index of observed variable **/
            Eigen::Vector2d obs_variable_index;
            
            
            /** names of inferred parameters**/
            std::vector<std::string> iparam_names;
            
            /** object for sampling random numbers (thread safe)**/
            mci::stochastic rng;
            
            /** cache**/
            std::map< size_t, double> cache_haz_off;
            std::map< size_t, double> cache_haz_on;
            
        public:
            
            const static std::string PARAM_dt;
            const static std::string PARAM_B_0;
            const static std::string PARAM_B_1;
            const static std::string PARAM_mu_tau_on;
            const static std::string PARAM_mu_tau_off;
            const static std::string PARAM_std_tau_on;
            const static std::string PARAM_std_tau_off;
            const static std::string PARAM_f;
            const static std::string PARAM_alpha_1;
            const static std::string PARAM_alpha_2;
            const static std::string PARAM_tau_1;
            const static std::string PARAM_tau_2;
            const static std::string PARAM_LH_0;
            const static std::string PARAM_sigma_obs;
            const static std::string PARAM_observed_variable_index_0;
            const static std::string PARAM_dim;
            const static std::string PARAM_grad_dim;
            const static std::string NAME;
            const static std::string PARAM_prior_mean_B_1 ;
            const static std::string PARAM_prior_std_B_1 ;
            const static std::string PARAM_prior_mean_alpha_1 ;
            const static std::string PARAM_prior_std_alpha_1 ;
            
            /** Create a model object **/
            dual_mode_model_gradinfo(const model_params_t& p, const mci::stochastic& r, const std::vector<std::string>& pi) ;
            
            /** Get model parameters **/
            const model_params_t get_parameters();
            
            /** Get named model parameter  **/
            const double get_parameter(const std::string& str);
            
            /** Set all model parameters **/
            void set_parameters(const model_params_t& p);
            
            /** update the inferred model parameters**/
            void update_i_parameters(const Eigen::VectorXd& pp);
            
            /** return logarithm of parameter prior density \pi(theta) **/
            const double log_prior_density(const Eigen::VectorXd& theta);
            
            /** return gradient of log prior density **/
            const Eigen::VectorXd grad_log_prior_density(const Eigen::VectorXd&);
            
            /** return logarithm of proposal density p(theta_proposed | theta_current) **/
            const double log_proposal_density(const Eigen::VectorXd& theta_current, const Eigen::VectorXd& theta_proposed);
            
            /** return values of inferred model parameters**/
            const Eigen::VectorXd get_i_parameters() ;
            
            /** generate initial S **/
            const state_t generate_ic(const double t_0);
            
            /** return logarithm of conditional observation density g(X_t,Y_t) **/
            const double log_obs_density(const state_t& X, const Eigen::VectorXd& Y);
            
            /** calculate and return the FIM **/
            Eigen::MatrixXd FIM( const state_t & X);
            
            /** return grad of the logarithm of the conditional observation density grad_theta g(X_t,Y_t) wiht re**/
            const Eigen::VectorXd grad_log_obs_density(const state_t& X, const Eigen::VectorXd& Y);
            
            /** return logartihm of transtion density f(X_t2 | X_t1)**/
            const double log_transition_density(const state_t& X_t1, const state_t& X_t2, const size_t& t);
            
            /** return the logarithm of the joint density P(X,Y)**/
            const double log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* ll=0);
            
            /****/
            const Eigen::VectorXd grad_log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* jll=0, double* mll=0, Eigen::MatrixXd* FIM=0, std::vector<Eigen::MatrixXd>* G_k=0);
            
            /**add prior inofmration to FIM and its derivatives**/
            const void add_prior_information(const Eigen::VectorXd& theta, Eigen::MatrixXd* G, std::vector<Eigen::MatrixXd>* G_k = 0);
            
            /** run model from t0 to t1 using X0 as the initial condition **/
            const trajectory_t run(const state_t& X0, const double& t_start, const double& t_end, const trajectory_t& data);
            
            /** run model using X_ref_traj as a reference trajectory**/
            const trajectory_t run(const trajectory_t& X_ref_traj, const state_t& X0, const trajectory_t& data) ;
            
            /*double event_propensity(const negative_binomial& dist, const double& x) ;
             double event_propensity( const double x, const double r, const double p) ;
             **/
            
            
            double propensity(const double x, const double m, const double s) ;
            double get_quantile(const double x, const double m, const double s) ;
            
            double propensity_on(const double x, const double m, const double s) ;
            double get_quantile_on(const double x, const double m, const double s) ;
            
            const size_t get_dim();
            const double get_dt() ;
            
        };
        
        
    }
}

#endif /* dual_mode_model_gradinfo_hpp */

