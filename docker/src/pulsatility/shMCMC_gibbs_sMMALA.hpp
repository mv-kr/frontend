//
//  shMCMC_gibbs_sMMALA.hpp
//  mci
//
//  Created by Margaritis Voliotis on 20/09/2020.
//  Copyright © 2020 Margaritis Voliotis. All rights reserved.
//

#ifndef shMCMC_gibbs_sMMALA_hpp
#define shMCMC_gibbs_sMMALA_hpp



#include "basic_algorthm.hpp"
#include "cSMC_AS.hpp"
#include "SMC.hpp"
#include "basic_model.hpp"
#include "basic_model_with_grad_info.hpp"
#include "stochastic.hpp"

#include <iostream>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include <vector>
#include <omp.h>

namespace mci {
    namespace algorithms {
        
        using MCMC_gibbs_hmc_parameter_t = std::map<std::string, double>;
        using  chain_doubles_t = std::vector<double> ;
        
        class shMCMC_gibbs_sMMALA: public basic_algorithm<MCMC_gibbs_hmc_parameter_t, chain_doubles_t> {
            
        protected:
            
            mci::models::basic_model* model_ptr;
            Eigen::MatrixXd data;
            Eigen::MatrixXd cov;
            Eigen::MatrixXd normTransform;
            Eigen::MatrixXd mass;
            Eigen::MatrixXd M_normTransform;
            Eigen::MatrixXd mass_inv;
            double mass_det;
            
            cSMC_AS* cSMC_AS_ptr;
            SMC* SMC_ptr;
            stochastic rng;
            Eigen::VectorXd blocks_i;
            
            Eigen::VectorXd i_guess;
            Eigen::VectorXd lower_limit;
            Eigen::VectorXd upper_limit;
            
            Eigen::MatrixXd chain_res;
            Eigen::VectorXd likelihood_res;
            std::vector<Eigen::MatrixXd> trajectories_res;
            
        public:
            
            /** const static members declarations**/
            const static std::string PARAM_iterations;
            const static std::string PARAM_warmup_pct;
            const static std::string PARAM_theta_dim;
            const static std::string PARAM_hmc_theta_dim;
            const static std::string PARAM_epsilon;
            const static std::string PARAM_M;
            const static std::string PARAM_iter_adapt;
            const static std::string PARAM_delta;
            const static std::string PARAM_steps;
            const static std::string PARAM_blocks_no;
            const static std::string PARAM_mmala_blocks_no;
            const static std::string PARAM_report_every;
            const static std::string NAME;
            const static std::string PARAM_adapt_runs;
            
            /** constructor **/
            shMCMC_gibbs_sMMALA(const MCMC_gibbs_hmc_parameter_t& p,
                     mci::models::basic_model* m, cSMC_AS* csmc, SMC* smc,
                     const Eigen::MatrixXd& d, const Eigen::MatrixXd& c, const Eigen::MatrixXd& M,
                     const stochastic& r,
                     const Eigen::VectorXd& i_g,
                     const Eigen::VectorXd& l_limit,
                     const Eigen::VectorXd& u_limit,
                     const Eigen::VectorXd& bl_i);
            
            /** run the algorithm and store results in res **/
            const void run(chain_doubles_t& res_traj) ;
            
            /** run the algorithm using parameters in params and store results in res **/
            const void run(const MCMC_gibbs_hmc_parameter_t& p, chain_doubles_t& res_traj) ;
            
            Eigen::MatrixXd leapfrog( Eigen::VectorXd&, Eigen::VectorXd&, const Eigen::MatrixXd&, double, bool*);
        
            double find_epsilon(const Eigen::VectorXd& , const Eigen::MatrixXd&);
            const void hmc_sampler(Eigen::VectorXd& theta ,
                                   const Eigen::MatrixXd& X_sample,
                                   const double& epsilon,
                                   double& acceptance_prob,
                                   double& joint_llik,
                                   double& marginal_llik,
                                   unsigned int& div);
            
            
            const double MH_sampler_single_parameter(Eigen::VectorXd& theta,
                                                   const Eigen::MatrixXd& X_sample,
                                                   const double& epsilon,
                                                   double& acceptance_prob,
                                                   double& joint_llik,
                                                   double& marginal_llik,
                                                   const size_t& start_block_i,
                                                   const size_t& size_block_i);
            
            const void sMMALA_sampler(Eigen::VectorXd& theta,
                                                      const Eigen::MatrixXd& X_sample,
                                                      const double& epsilon,
                                                      double& acceptance_prob,
                                                      double& joint_llik,
                                                      double& marginal_llik,
                                                      unsigned int& div,
                                                       unsigned int iter  );
            
            const void MMALA_sampler(Eigen::VectorXd& theta,
                                      const Eigen::MatrixXd& X_sample,
                                      const double& epsilon,
                                      double& acceptance_prob,
                                      double& joint_llik,
                                      double& marginal_llik,
                                      unsigned int& div,
                                      unsigned int iter  );
            
            const void AIS_sampler( Eigen::VectorXd& theta,
                                      Eigen::MatrixXd& X_sample,
                                     const double& epsilon,
                                     double& acceptance_prob,
                                     double& joint_llik,
                                     double& marginal_llik,
                                     const size_t& start_block_i,
                                     const size_t& size_block_i) ;
            
            
            const void sMMALA_sampler_block(Eigen::VectorXd& theta,
                                                            const Eigen::MatrixXd& X_sample,
                                                            const double& epsilon,
                                                            double& acceptance_prob,
                                                            double& joint_llik,
                                                            double& marginal_llik,
                                                            unsigned int& div,
                                                            const size_t& start_index,
                                                            const size_t& block_size);
            
            
            const void MMALA_sampler_block(Eigen::VectorXd& theta,
                                                           const Eigen::MatrixXd& X_sample,
                                                           const double& epsilon,
                                                           double& acceptance_prob,
                                                           double& joint_llik,
                                                           double& marginal_llik,
                                                           unsigned int& div,
                                                           const size_t& start_index,
                                                           const size_t& block_size);
            
            const void hmc_sampler_warmup(Eigen::VectorXd& theta ,
                                   const Eigen::MatrixXd& X_sample,
                                   const double& epsilon,
                                   double& acceptance_prob,
                                   double& joint_llik,
                                   double& marginal_llik,
                                   unsigned int& div);
            
            const void print_chain_parameters();
            const void print_chain_parameters(std::string str_file);
            const void print_chain_likelihood(std::string str_file);
            const void print_chain_trajectories(std::string str_file);
            
            
        };
        
        
    }
}


#endif /* shMCMC_gibbs_sMMALA */
