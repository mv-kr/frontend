//
//  cSMC_AS.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 23/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef cSMC_AS_hpp
#define cSMC_AS_hpp

#include "basic_algorthm.hpp"
#include "basic_model.hpp"
#include "stochastic.hpp"



#include <cmath>
#include <map>
#include <Eigen/Dense>




namespace mci {
    namespace algorithms {
        
        using cSMC_AS_parameter_t = std::map<std::string, double>;
        

        
        class cSMC_AS : public basic_algorithm<cSMC_AS_parameter_t, mci::trajectory_t> {
        
        protected:
            
            //cSMC_AS_parameter_t parameters;
            mci::models::basic_model* model_obj;
            Eigen::MatrixXd data;
            stochastic rng;
            
        public:
            
            const static std::string PARAM_nparticles;
            const static std::string PARAM_DT;
            const static std::string PARAM_AS_probability;
            const static std::string NAME;
            
            cSMC_AS(mci::models::basic_model* m, const Eigen::MatrixXd& d, const stochastic& r) ;
            
            
            const void run(mci::trajectory_t& res_traj) ;
            
            const void run(const cSMC_AS_parameter_t& p, mci::trajectory_t& res_traj) ;
                     
            /** run the algorithm and store results in res **/
            const double run(const mci::trajectory_t& ref_traj, mci::trajectory_t& res_traj) ;
            
            /** run the algorithm using parameters in params and store results in res **/
            const void run(const cSMC_AS_parameter_t& p, const mci::trajectory_t& ref_traj, mci::trajectory_t& res_traj) ;
            
            /** run the algorithm and store results in res **/
            const void run_2(const mci::trajectory_t& ref_traj, mci::trajectory_t& res_traj) ;
            
            
        };
        
        
    }
}


#endif /* cSMC_AS_hpp */
