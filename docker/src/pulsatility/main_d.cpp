//
//  main_d.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 25/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#include "MCMC_gibbs.hpp"
#include "MCMC_gibbs_sMMALA.hpp"
#include "cSMC_AS.hpp"
#include "SMC.hpp"


#include "stochastic.hpp"
#include "LH_secretion_model.hpp"
#include "LH_secretion_across_studies_model.hpp"
#include "discrete_time_model.hpp"
#include "parameter_parser.hpp"

#include <iostream>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <filesystem>

using namespace mci;




int main(int argc, char *argv[]) {
    
    /* parse the arguments */
    if (argc < 3) {
        std::string name = argv[0] ;
        std::cout << "Usage:" +name + " data_file parameter_file [output_folder]" <<std::endl;
        return -1;
    }
    
    //parser object
    mci::utils::parameter_parser p;
    
    try {
        //parse parameters
        p.read_parameters(argv[2]);
        //parse data
        p.read_data(argv[1]);
    } catch (std::invalid_argument& e) {
        std::cout << e.what() <<std::endl;
        return -1;
    }
    
    //load parsed parameters
    std::map<std::string, mci::utils::parameter_set> res = p.resD; //contains scalar-valued parameters
    std::map<std::string, Eigen::VectorXd> resV = p.resV; //contains vector-vlaued parameters
    std::map<std::string, Eigen::MatrixXd> resM = p.resM; //contains matrix-valued parameters
    std::map<std::string, std::string> resS_s = p.resS_s; //contains strings
    std::map<std::string, std::vector<std::string>> resV_s = p.resV_s; // contains vector of strings
    
    //load parsed data
    Eigen::MatrixXd my_data_all = p.data;
    
    // get outputfolder
    std::string output_folder = "";
    if(argc >= 4) {
        output_folder = argv[3] ;
        std::filesystem::create_directory(output_folder);
        output_folder = output_folder +  "/";
    } else {
        output_folder = resS_s["save_location_str"];
    }
    
    //std::cout<<output_folder <<std::endl;
    
    //get number of threads
    int num_threads_used = res["general"]["threads_no"];
    if (num_threads_used < 1 && num_threads_used != -1) {
        num_threads_used = 1;
        std::cout << "Setting parameter threads_no == "<< num_threads_used << std::endl;
    }
    if (num_threads_used > OpenMP::getMaxThreads() || num_threads_used == -1) {
        num_threads_used = (int) (OpenMP::getMaxThreads());
        std::cout << "Setting parameter threads_no == "<< num_threads_used << std::endl;
    }
    
    
    //set up rng object for each thread
    mci::stochastic rng;
    mci::rng_type master_rng;
    unsigned int seed = res["general"]["seed"];
    if (seed == 0 )
        seed = (unsigned int)std::time(0);
    master_rng.seed(seed);
    boost::random::uniform_int_distribution<> dist;
    for (size_t i=0; i < num_threads_used; i++) {
        mci::rng_type gnr; //create rng
        mci::stochastic::generators.push_back(gnr); //add it to the list
        mci::stochastic::generators[i].seed(dist(master_rng)); //initialise seed using master_rng
        
    }
    
    size_t chain_repeats = res["general"]["chain_repeats"];
    if (chain_repeats < 1) {
        chain_repeats = 1;
        std::cout << "Setting parameter chain_repeats == "<< chain_repeats << std::endl;
    }
    
    size_t rows_per_dataset = res["general"]["rows_per_dataset_entry"];
    if (rows_per_dataset==0) rows_per_dataset=2;
    size_t start_dataset_entry = res["general"]["start_dataset_entry"];
    size_t end_dataset_entry = res["general"]["end_dataset_entry"];
    
    
    size_t dataset_no = (my_data_all.rows())/rows_per_dataset;
    if (start_dataset_entry > 0 ) {
        start_dataset_entry--;
    }
    if (end_dataset_entry > 0 ) {
        end_dataset_entry--;
        if( end_dataset_entry < start_dataset_entry || end_dataset_entry >= dataset_no ) {
            end_dataset_entry = dataset_no-1;
        }
    } else {
        end_dataset_entry = dataset_no-1;
    }
    
    // compile list of jobs as pairs: (dataset id, chain id)
    int number_of_jobs = (int) (chain_repeats*(end_dataset_entry - start_dataset_entry + 1));
    Eigen::MatrixXd job_order_matrix(number_of_jobs,2);
    int k =0;
    for (size_t i =start_dataset_entry; i <= end_dataset_entry ; i++) {
        for (size_t ci= 0; ci<chain_repeats; ci++) {
            
            job_order_matrix(k,0)= i;
            job_order_matrix(k,1)= ci;
            k++;
        }
    }

    // start processing all job in parallel
    #pragma omp parallel for if (num_threads_used>1) num_threads(num_threads_used) schedule(static)
    for (size_t k =0; k < job_order_matrix.rows() ; k++) {
        //std::cout<< job_order_matrix(k,0) << "-" << job_order_matrix(k,1) << std::endl;
        size_t ci = job_order_matrix(k,1); //chain index
        size_t i = job_order_matrix(k,0); // dataset index
        size_t current_dataset_i = rows_per_dataset*i;
    
        //get start time, final time and number of timepoints
        size_t final_index = my_data_all.cols()-1;
        while (my_data_all(current_dataset_i,final_index) <= 0 ) {
            final_index--;
        }
        double t_0    = my_data_all(current_dataset_i,0);
        double t_end  = my_data_all(current_dataset_i,final_index);
        double exp_dt = my_data_all(current_dataset_i,1) - t_0;
        size_t timepoints_no = final_index+1;
    
        //load current dataset into a MatrixXd object;
        Eigen::MatrixXd data(rows_per_dataset, timepoints_no);
        for (size_t jj = 0; jj < rows_per_dataset; jj++) {
            
            if (jj%2 == 0)
                data.row(jj) = my_data_all.block(current_dataset_i+rows_per_dataset-1-jj, 0, 1, timepoints_no);
            else
                data.row(jj) = my_data_all.block(current_dataset_i+rows_per_dataset-1-jj, 0, 1, timepoints_no) - Eigen::MatrixXd::Constant(1,timepoints_no,t_0);
        }
        
        t_end = t_end-t_0;
        t_0 = 0;
        
        
        
        
        //get CV for first observed variable
        double cv;
        if (i <= resV["experiment_CV"].size()-1)
            cv = resV["experiment_CV"](i);
        else
            cv = resV["experiment_CV"](resV["experiment_CV"].size()-1);
        
        // the model and MCMC algorithm that will be used
        std::string model_str = resS_s["model_code_str"];
        std::string mcmc_algo_str = resS_s["alg_code_str"];
        
        /** model parameters **/
        std::map<std::string, double> model_params;
        model_params = res[model_str];
        
        /** SMC parameters **/
        std::map<std::string, double> SMC_params;
        SMC_params = res[mci::algorithms::SMC::NAME1];
    
        /** cSMC_AS parameters**/
        std::map<std::string, double> cSMC_AS_params;
        cSMC_AS_params = res[mci::algorithms::cSMC_AS::NAME];
    
        /** MCMC parameters**/
        std::map<std::string, double> MCMC_params;
        MCMC_params = res[mcmc_algo_str];
        
    
        /** names of model parameter to be inferred **/
        std::vector<std::string> i_param_names = resV_s["i_param_names"];
        /** initial values for inferred parameters **/
        Eigen::VectorXd initial_values = resV["initial_values"];
        /** ranges for inferred parameters **/
        Eigen::VectorXd l_limit = resV["min_values"];
        Eigen::VectorXd u_limit = resV["max_values"];
        /** covariance matrix **/
        Eigen::MatrixXd cov = resM["cov"];
    
        /** A few parameter checks **/
        
        //make sure the hasn't been an error with  parameter cv
        if (cv <= 0 )
            throw std::invalid_argument("The value of parameter cv is out of range.");

        //make sure user defined value for parameter dt is positive too small or too large
        double dt = model_params["dt"];
        if (dt <= 0 || dt >= exp_dt ) {
            if (data(1,1) > 10)
                model_params["dt"] = 1;
            else
                model_params["dt"] = data(1,1)/10.0;
            
            std::cout << "Reset parameter dt value: " << dt << " --> " <<  model_params["dt"] <<std::endl;
            
            dt = model_params["dt"];
        }
        
        
        
        //make sure user has provided the right number of values for lower and upper parameter limits
        size_t i_param_number = i_param_names.size();
        if ( initial_values.size() != i_param_number)
            throw std::invalid_argument("The number of initial values provided is not consistent with the number of parameters.");
        if (l_limit.size() != i_param_number )
            throw std::invalid_argument("The number of max values is not consistent with the number of parameters.");
        if (u_limit.size() != i_param_number )
            throw std::invalid_argument("The number of min values is not consistent with the number of parameters.");
        //  if  (cov.cols()  != i_param_number  ||
        //      cov.rows() != i_param_number) {
        //      throw std::invalid_argument("The dimension of the cov matrix is not consistent with the number of parameters.");
        //  }
        
        if ( (initial_values.array() <= l_limit.array()).any() ||
                (initial_values.array() >= u_limit.array()).any() ) {
            throw std::invalid_argument("Initial parameter values are outside the range");
        }
            
        MCMC_params["timepoints"] = round( (t_end - t_0)/ dt ) + 1; //CAN I FIX THIS
        
        /*  ---------------------------------------------------------------------------------------- */
        std::cout<< "Thread " << omp_get_thread_num() << " started. "
        << "Analysing dataset "  << i+1 << "/" << dataset_no << " (chain_id = "<< ci << ")"<<std::endl;
        
        std::cout<< " ----- Data summary ----- " << std::endl;
        std::cout<< " Dataset       : " << i+1 << std::endl;
        std::cout<< " Starting time : " << t_0 << std::endl;
        std::cout<< " Final time    : " << t_end << std::endl;
        std::cout<< " time points   : " << timepoints_no << std::endl;
        std::cout<< " 1st sampling interval   : " << exp_dt << std::endl;
        std::cout<< " mean sampling interval   : " << (t_end-t_0)/(timepoints_no-1) << std::endl;
        std::cout<< " CV            : " << cv*100 << "%" << std::endl;
        std::cout<< std::endl;
        std::cout<< " ----- Data -------------- " << std::endl;
        std::cout<< data.transpose() << std::endl;
        
        std::cout<< " ----- Model summary -----" << std::endl;
        std::cout<< " time resolution (dt): " << dt << std::endl;
        std::cout<< " timepoints          : " << round( (t_end - t_0)/ dt ) + 1 << std::endl;
        std::cout<< " model code          : " << model_str << std::endl;
        std::cout<< " MCMC alg code       : " << mcmc_algo_str << std::endl;
        std::cout<< std::endl;
        
        std::cout<< " ----- Initial param (log) values and ranges -----" << std::endl;
        size_t max_str_length=0;
        for (size_t pi = 0 ; pi < i_param_names.size();pi++) {
            if (max_str_length<i_param_names[pi].size())
                max_str_length =i_param_names[pi].size();
        }
        for (size_t pi = 0 ; pi < i_param_names.size();pi++) {
            std::cout << " " << i_param_names[pi] << std::string(max_str_length -i_param_names[pi].size(), ' ' ) <<": " << initial_values(pi)
            << "(" <<l_limit(pi) << "," << u_limit(pi) <<")"
            << std::endl;
        }
        std::cout<< std::endl;
        
        //   std::cout<< " ----- Covariance matrix -----" << std::endl;
        //   std::cout<<cov<<std::endl;
        //   std::cout<< " -----------------------------" << std::endl;
        /*  ---------------------------------------------------------------------------------------- */
        
    
    
        // rng utility object
        mci::stochastic rng;
        
        //set up model object
        mci::models::basic_model *model_obj;
    
        if  ( model_str.compare(mci::models::LH_secretion_model::NAME) == 0  ) { //LH_secretion_model
            
            model_params[mci::models::LH_secretion_model::PARAM_LH_0]         =  data(0,0);
            
            // LH index
            model_params[mci::models::LH_secretion_model::PARAM_observed_variable_index_0]  =  6;
            if (model_params[mci::models::LH_secretion_model::PARAM_tau_1]  == 0)
                model_params[mci::models::LH_secretion_model::PARAM_tau_1]  =  (t_end-t_0)/(timepoints_no-1);
            // LH clearance rate prior distribution parameters (default values)
            if (model_params[mci::models::LH_secretion_model::PARAM_prior_std_alpha_1] <=0 ) {
                model_params[mci::models::LH_secretion_model::PARAM_prior_mean_alpha_1] = 80.0;
                model_params[mci::models::LH_secretion_model::PARAM_prior_std_alpha_1] = 9.3;
            }
            
            model_params[mci::models::LH_secretion_model::PARAM_sigma_obs]    =  cv;
            
            if (model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_alpha] <= 0)
                model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_alpha] = 1;
            if (model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_beta] <= 0)
                model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_beta] = 1;
            // LH secretion prior distribution parameters (default values)
            if (model_params[mci::models::LH_secretion_model::PARAM_prior_std_B_1] <= 0) {
                model_params[mci::models::LH_secretion_model::PARAM_prior_std_B_1] = 5;
                model_params[mci::models::LH_secretion_model::PARAM_prior_mean_B_1] = 0;
            }
            
            model_obj = new mci::models::LH_secretion_model(model_params, rng, i_param_names ) ;
            
        }  else if (model_str.compare(mci::models::discrete_time_model::NAME) == 0) { //discrete_time_model
            
            model_params[mci::models::discrete_time_model::PARAM_LH_0]         =  data(0,0);
            
            // LH index
            model_params[mci::models::discrete_time_model::PARAM_observed_variable_index_0]  =  6;
            
            // LH clearance rate prior distribution parameters (default values)
            if (model_params[mci::models::discrete_time_model::PARAM_prior_std_alpha_1] <=0 ) {
                model_params[mci::models::discrete_time_model::PARAM_prior_mean_alpha_1] = 80.0;
                model_params[mci::models::discrete_time_model::PARAM_prior_std_alpha_1] = 9.3;
            }
            
            model_params[mci::models::discrete_time_model::PARAM_sigma_obs]    =  cv;
            
            if (model_params[mci::models::discrete_time_model::PARAM_f_Beta_prior_alpha] <= 0)
                model_params[mci::models::discrete_time_model::PARAM_f_Beta_prior_alpha] = 1;
            if (model_params[mci::models::discrete_time_model::PARAM_f_Beta_prior_beta] <= 0)
                model_params[mci::models::discrete_time_model::PARAM_f_Beta_prior_beta] = 1;
            // LH secretion prior distribution parameters (default values)
            if (model_params[mci::models::discrete_time_model::PARAM_prior_std_B_1] <= 0) {
                model_params[mci::models::discrete_time_model::PARAM_prior_std_B_1] = 5;
                model_params[mci::models::discrete_time_model::PARAM_prior_mean_B_1] = 0;
            }
            
            model_obj = new mci::models::discrete_time_model(model_params, rng, i_param_names ) ;
            
            
        }else if ( model_str.compare(mci::models::LH_secretion_across_studies_model::NAME) == 0 ) { //LH_secretion_across_studies_model
            
            model_params[mci::models::LH_secretion_across_studies_model::PARAM_LH_0_b]         =  data(0,0);
            model_params[mci::models::LH_secretion_across_studies_model::PARAM_LH_0_a]         =  data(2,0);
            model_params[mci::models::LH_secretion_across_studies_model::PARAM_sigma_obs_a]    =  cv;
            model_params[mci::models::LH_secretion_across_studies_model::PARAM_sigma_obs_b]    =  cv;
            model_params[mci::models::LH_secretion_across_studies_model::PARAM_observed_variable_index_0]  =  6;
            if (model_params[mci::models::LH_secretion_across_studies_model::PARAM_tau_1]  == 0)
                model_params[mci::models::LH_secretion_across_studies_model::PARAM_tau_1]  =  (t_end-t_0)/(timepoints_no-1);
            
            // LH clearance rate prior distribution parameters (default values)
            if (model_params[mci::models::LH_secretion_model::PARAM_prior_std_alpha_1] <=0 ) {
                model_params[mci::models::LH_secretion_model::PARAM_prior_mean_alpha_1] = 80.0;
                model_params[mci::models::LH_secretion_model::PARAM_prior_std_alpha_1] = 9.3;
            }
            if (model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_alpha] <= 0)
                model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_alpha] = 1;
            if (model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_beta] <= 0)
                model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_beta] = 1;
            // LH secretion prior distribution parameters (default values)
            if (model_params[mci::models::LH_secretion_model::PARAM_prior_std_B_1] <= 0) {
                model_params[mci::models::LH_secretion_model::PARAM_prior_std_B_1] = 5;
                model_params[mci::models::LH_secretion_model::PARAM_prior_mean_B_1] = 0;
            }
            
            model_obj = new mci::models::LH_secretion_across_studies_model(model_params, rng, i_param_names ) ;
            
            
        } else { //initialise the model_obj using the LH_secretion_model model
            
            model_params[mci::models::LH_secretion_model::PARAM_LH_0]         =  data(0,0);
            model_params[mci::models::LH_secretion_model::PARAM_sigma_obs]    =  cv;
            // LH index
            model_params[mci::models::LH_secretion_model::PARAM_observed_variable_index_0]  =  6;
            if (model_params[mci::models::LH_secretion_model::PARAM_tau_1]  == 0)
                model_params[mci::models::LH_secretion_model::PARAM_tau_1]  =  (t_end-t_0)/(timepoints_no-1);
            
            // LH clearance rate prior distribution parameters (default values)
            model_params[mci::models::LH_secretion_model::PARAM_prior_mean_alpha_1] = 80.0;
            model_params[mci::models::LH_secretion_model::PARAM_prior_std_alpha_1] = 9.3;
            // pulsatility parameter prior distribution parameters (default values)
            model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_alpha] = 1;
            model_params[mci::models::LH_secretion_model::PARAM_f_Beta_prior_beta] = 1;
            
            model_params[mci::models::LH_secretion_model::PARAM_prior_mean_B_1] = 0;
            model_params[mci::models::LH_secretion_model::PARAM_prior_std_B_1] = 5;
            
            
            
            model_obj = new mci::models::LH_secretion_model(model_params, rng, i_param_names ) ;
            
        }

        // set up SMC objs
        mci::algorithms::SMC SMC_obj(model_obj, data.topRows(rows_per_dataset), rng);
        SMC_obj.set_parameters(SMC_params);
        
        // set up cSMC obj
        mci::algorithms::cSMC_AS cSMC_obj(model_obj, data.topRows(rows_per_dataset), rng);
        cSMC_obj.set_parameters(cSMC_AS_params);
        
        std::string str_file_prefix = output_folder  + resS_s["identifier_str"] + "_" + std::to_string(i+1) + "_chain_" + std::to_string(ci) ;
        
        
        if ( mcmc_algo_str.compare(mci::algorithms::MCMC_gibbs_sMMALA::NAME) == 0 ) {
        
            Eigen::VectorXd blocks_i = resV["block_start_indices"];
            
            if (MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iterations] == 0)
                MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iterations] = 30000;
            if (MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_epsilon] == 0)
                MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_epsilon] = 0.01;
            if (MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iter_adapt] == 0)
                MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iter_adapt] = round(MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_iterations]/10);
            if (MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_delta] == 0)
                MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_delta] = 0.234;
                
            if (MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_report_every] == 0)
                MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_report_every] = 100;
            if (MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_report_every] == 0)
                MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_report_every] = 20;
            
            
            
            MCMC_params[mci::algorithms::MCMC_gibbs_sMMALA::PARAM_blocks_no] = blocks_i.size()-1;
            
            
                
                
            //set up MCMC_gibbs_AIC obj
            mci::algorithms::MCMC_gibbs_sMMALA alg(MCMC_params, model_obj, &cSMC_obj, &SMC_obj ,
                                                data,
                                                cov,
                                                cov,
                                                rng, initial_values, l_limit, u_limit, blocks_i);
            
            
            //run
            std::vector<double> tmp;
            //alg.warmup_run(tmp);
            alg.run(tmp);
            
            //alg.print_chain_parameters();
            std::cout <<std::endl << "Thread " << OpenMP::getThreadNum() << ": writing output files  " <<std::endl ;
            std::cout << "Thread " << OpenMP::getThreadNum() << ": ... " << str_file_prefix + "_parameters.csv" <<std::endl;
            alg.print_chain_parameters(str_file_prefix+"_parameters.csv") ;
            std::cout << "Thread " << OpenMP::getThreadNum() << ": ... " << str_file_prefix + "_lik.csv" <<std::endl;
            alg.print_chain_likelihood(str_file_prefix+"_lik.csv") ;
            std::cout << "Thread " << OpenMP::getThreadNum() << ": ... " << str_file_prefix + "_trajectories.csv" <<std::endl;
            alg.print_chain_trajectories(str_file_prefix + "_trajectories.csv") ;
            
            
        } else { //use MCMC_gibbs algorithm
        
            Eigen::VectorXd p_blocks = resV["parameter_blocks"];
            //set up MCMC obj
            mci::algorithms::MCMC_gibbs alg(MCMC_params, model_obj, &cSMC_obj, &SMC_obj , data, cov, rng, initial_values, l_limit, u_limit, p_blocks);
            
            //run
            std::vector<double> tmp;
            alg.run(tmp);
            
            //alg.print_chain_parameters();
            alg.print_chain_parameters(str_file_prefix+"_parameters.txt") ;
            alg.print_chain_likelihood(str_file_prefix+"_lik.txt") ;
            alg.print_chain_trajectories(str_file_prefix + "_trajectories.txt") ;
            
        
        }
        
    }
    return 0;
    
}

