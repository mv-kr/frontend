//
//  basic_model.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 23/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef basic_model_hpp
#define basic_model_hpp

#include <iostream>
#include <map>
#include <Eigen/Dense>

namespace mci {
    
    using trajectory_t = Eigen::MatrixXd;
    
    namespace models {
    
        
        using state_t = Eigen::VectorXd ;
        using model_params_t = std::map<std::string, double>;
        
        class basic_model  {
         

            
        protected:
            
            /** model parameters **/
            model_params_t params;
            
            /** dimension (numder of state variables + time) **/
            size_t variables_num;
            
        public:
     
            virtual ~basic_model() {};
            
            /** update the inferred model parameters**/
            virtual void update_i_parameters(const Eigen::VectorXd& pp) = 0;
            
            /** calculate the log prior density: log \pi(theta) **/
            virtual const double log_prior_density(const Eigen::VectorXd& theta) = 0;
            
            /** calculate the log proposal density: log p(theta_proposed | theta_current) **/
            virtual const double log_proposal_density(const Eigen::VectorXd& theta_proposed, const Eigen::VectorXd& theta_current) = 0;
            
            /** calculate the log ratio of proposal densities: p(theta_proposed | theta_current)/ p(theta_current  | theta_proposed) **/
            const double log_proposal_density_ratio(const Eigen::VectorXd& theta_current, const Eigen::VectorXd& theta_proposed) ;
            
            
            /** return values of inferred model parameters**/
            virtual const Eigen::VectorXd get_i_parameters() = 0;
            
            /** Get the model parameters **/
            virtual const model_params_t get_parameters() = 0;
            
            /** Get the model parameters **/
            virtual const double get_parameter(const std::string& str) = 0;
            virtual void set_parameters(const model_params_t& p) = 0;
            /** generate initial S **/
            virtual const state_t generate_ic(const double t_0) = 0;
            
            /** return logarithm of conditional observation density g(X_t,Y_t) **/
            virtual const double log_obs_density(const state_t& X, const Eigen::VectorXd& Y) = 0;
            
            /** return logartihm of transtion density f(X_t2 | X_t1)**/
            virtual  const double log_transition_density(const state_t& X_t1, const state_t& X_t2, const size_t& t) = 0;
            
            /** return the logarithm of the joint density P(X,Y)**/
            virtual const double log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* ll=0) = 0;
            
            /** run model from t0 to t1 using X0 as the initial condition **/
            virtual const trajectory_t run(const state_t& X0, const double& t_start, const double& t_end, const trajectory_t& data) = 0;
            
            /** run model using X_ref_traj as a reference trajectory**/
            virtual const trajectory_t run(const trajectory_t& X_ref_traj, const state_t& X0, const trajectory_t& data) = 0;
            
            
            virtual const size_t get_dim() = 0;
            virtual const double get_dt() = 0;
            
            
            /** return gradient of log prior density **/
            virtual const Eigen::VectorXd grad_log_prior_density(const Eigen::VectorXd&) = 0;
            
            /** return grad of the logarithm of the conditional observation density g(X_t,Y_t) **/
            virtual const Eigen::VectorXd grad_log_obs_density(const state_t& X, const Eigen::VectorXd& Y) = 0;
            
            /** return the grad of the logarithm of the joint density P(X,Y)**/
            virtual const Eigen::VectorXd grad_log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* jll=0, double* mll=0, Eigen::MatrixXd* FIM=0, std::vector<Eigen::MatrixXd>* G_k=0) = 0;
            
            /**add prior inofmration to FIM and its derivatives**/
            virtual const void add_prior_information(const Eigen::VectorXd& theta, Eigen::MatrixXd* G, std::vector<Eigen::MatrixXd>* G_k=0) = 0;
            
            /** return the grad of the logarithm of the joint density P(X,Y)**/
          //  virtual const Eigen::VectorXd grad_log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* jll=0, double* mll=0, Eigen::MatrixXd* FIM=0, std::vector<Eigen::MatrixXd>* G_k=0) = 0;
            
        };
    }
}

#endif /* basic_model_hpp */
