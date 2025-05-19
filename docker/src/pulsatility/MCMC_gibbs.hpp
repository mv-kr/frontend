//
//  MCMC_gibbs.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 24/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef MCMC_gibbs_hpp
#define MCMC_gibbs_hpp


#include "basic_algorthm.hpp"
#include "cSMC_AS.hpp"
#include "SMC.hpp"
#include "basic_model.hpp"
#include "stochastic.hpp"

#include <iostream>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include <vector>
#include <omp.h>

namespace mci {
    namespace algorithms {
        
        using cMCMC_gibbs_parameter_t = std::map<std::string, double>;
        using chain_doubles_t         = std::vector<double> ;
        
        class MCMC_gibbs: public basic_algorithm<cMCMC_gibbs_parameter_t, chain_doubles_t> {
            
        protected:
            
            mci::models::basic_model* model_ptr;
            Eigen::MatrixXd data;
            Eigen::MatrixXd cov;
            Eigen::MatrixXd normTransform;
            cSMC_AS* cSMC_AS_ptr;
            SMC* SMC_ptr;
            stochastic rng;
            
            Eigen::VectorXd i_guess;
            Eigen::VectorXd lower_limit;
            Eigen::VectorXd upper_limit;
            Eigen::VectorXd parameter_blocks;
            
            Eigen::MatrixXd chain_res;
            Eigen::VectorXd likelihood_res;
            std::vector<Eigen::MatrixXd> trajectories_res;
            
            /** generate intial vector of parameters () **/
            const void generate_init_params() ;
            
        public:
            
            /** const static members declarations**/
            const static std::string PARAM_iterations;
            const static std::string PARAM_theta_dim;
            const static std::string PARAM_K;
            const static std::string PARAM_report_every;
            const static std::string NAME;
            
            /** constructor **/
            MCMC_gibbs(const cMCMC_gibbs_parameter_t& p,
                       mci::models::basic_model* m, cSMC_AS* csmc, SMC* smc,
                       const Eigen::MatrixXd& d, const Eigen::MatrixXd& c,
                       const stochastic& r,
                       const Eigen::VectorXd& i_g,
                       const Eigen::VectorXd& l_limit,
                       const Eigen::VectorXd& u_limit,
                       const Eigen::VectorXd& p_blocks  );
            
            /** run the algorithm and store results in res **/
            const void run(chain_doubles_t& res_traj) ;
            
            
             
            /** run the algorithm using parameters in params and store results in res **/
            const void run(const cMCMC_gibbs_parameter_t& p, chain_doubles_t& res_traj) ;
            
            const void print_chain_parameters();
            const void print_chain_parameters(std::string str_file);
            const void print_chain_likelihood(std::string str_file);
            const void print_chain_trajectories(std::string str_file);
        };
        
        
    }
}




#endif /* MCMC_gibbs_hpp */

