//
//  basic_model_with_grad_info.hpp
//  mci
//
//  Created by Margaritis Voliotis on 06/10/2020.
//  Copyright © 2020 Margaritis Voliotis. All rights reserved.
//

#ifndef basic_model_with_grad_info_hpp
#define basic_model_with_grad_info_hpp

#include <iostream>
#include <map>
#include <Eigen/Dense>
#include "basic_model.hpp"

namespace mci {
    
     //using trajectory_t = Eigen::MatrixXd;
    
    namespace models {
        
       // using state_t = Eigen::VectorXd ;
       // using model_params_t = std::map<std::string, double>;
        
        class basic_model_with_grad_info  : public basic_model {
            
            
        public:
            
            /** return gradient of log prior density **/
            virtual const Eigen::VectorXd grad_log_prior_density(const Eigen::VectorXd&) = 0;
            
            /** return grad of the logarithm of the conditional observation density g(X_t,Y_t) **/
            virtual const Eigen::VectorXd grad_log_obs_density(const state_t& X, const Eigen::VectorXd& Y) = 0;
            
            /** return the grad of the logarithm of the joint density P(X,Y)**/
            virtual const Eigen::VectorXd grad_log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* ll=0) = 0;
            
            
        };
        
    }
}
#endif /* basic_model_with_grad_info_hpp */
