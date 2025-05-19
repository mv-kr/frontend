//
//  SMC.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 24/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef SMC_hpp
#define SMC_hpp

#include <iostream>


#include "basic_algorthm.hpp"
#include "basic_model.hpp"
#include "stochastic.hpp"



#include <cmath>
#include <map>
#include <Eigen/Dense>


namespace mci {
    namespace algorithms {
        
        using SMC_parameter_t = std::map<std::string, double>;
        
        
        
        class SMC: public basic_algorithm<SMC_parameter_t, mci::trajectory_t> {
            
    
            
            
        protected:
            
            mci::models::basic_model* model_obj;
            Eigen::MatrixXd data;
            stochastic rng;
            
        public:
            
            const static std::string PARAM_nparticles1;
            const static std::string NAME1;
            
            /** constructor **/
            SMC(mci::models::basic_model* m, const Eigen::MatrixXd& d, const stochastic& r);
            
            /** run the algorithm and store results in res **/
            const void run(trajectory_t& res_traj) ;
            
            /** run the algorithm using parameters in params and store results in res **/
            const void run(const SMC_parameter_t& p, trajectory_t& res_traj) ;
            
            
            
        };
        
        
    }
}



#endif /* SMC_hpp */
