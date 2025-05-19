//
//  composite_model.cpp
//  mci
//
//  Created by Margaritis Voliotis on 19/10/2020.
//  Copyright © 2020 Margaritis Voliotis. All rights reserved.
//

#include "LH_secretion_across_studies_model.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/assert.hpp>




using namespace mci;
using namespace mci::models;
using namespace boost::assign; // bring 'operator+=()' into scope

/** static members definitions**/
const std::string mci::models::LH_secretion_across_studies_model::PARAM_dt = "dt";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_B_1 = "k";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_mu_tau_on =  "mu_tau_on";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_mu_tau_off = "mu_tau_off";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_f = "f";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_alpha_1 = "d";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_tau_1 = "tau_1";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_LH_0_a = "LH_0_a";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_LH_0_b = "LH_0_b";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_sigma_obs_a = "sigma_obs_a";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_sigma_obs_b = "sigma_obs_b";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_diff_B_1 = "diff_B_1" ;
const std::string mci::models::LH_secretion_across_studies_model::PARAM_diff_mu_tau_on = "diff_mu_tau_on"  ;
const std::string mci::models::LH_secretion_across_studies_model::PARAM_diff_mu_tau_off = "diff_mu_tau_off"  ;
const std::string mci::models::LH_secretion_across_studies_model::PARAM_diff_f = "diff_f" ;

const std::string mci::models::LH_secretion_across_studies_model::PARAM_observed_variable_index_0 = "observed_variable_index_0";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_dim = "dim";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_dim_a = "dim_a";
const std::string mci::models::LH_secretion_across_studies_model::PARAM_dim_b = "dim_b";
const std::string mci::models::LH_secretion_across_studies_model::NAME = "LH_secretion_across_studies_model";


/** constructor **/
LH_secretion_across_studies_model::LH_secretion_across_studies_model(const model_params_t& p, const mci::stochastic& r, const std::vector<std::string>& pi) {
    params = p;

    
    obs_variable_index(0) = params[LH_secretion_across_studies_model::PARAM_observed_variable_index_0];
    iparam_names = pi;
    
    /** construct first model**/
    std::vector<std::string> iparam_names_a;
    iparam_names_a.push_back( iparam_names[0]);
    iparam_names_a.push_back( iparam_names[1]);
    iparam_names_a.push_back( iparam_names[2]);
    iparam_names_a.push_back( iparam_names[5]);
    iparam_names_a.push_back( iparam_names[7]);
  
    // -- construct first model
    model_params_t params_a = params;
    params_a[mci::models::LH_secretion_model::PARAM_LH_0] = params[PARAM_LH_0_a];
    params_a[mci::models::LH_secretion_model::PARAM_sigma_obs] = params[PARAM_sigma_obs_a];
    //params_a[mci::models::dual_mode_model_gradinfo::PARAM_dim] = params[PARAM_dim_a];
    model_obj_a = new mci::models::LH_secretion_model(params_a, r, iparam_names_a);
    //model_obj_a(params_a, r, iparam_names_a);
    
    
    // -- construct second model
    std::vector<std::string> iparam_names_b = pi;
    model_params_t params_b = params;
    params_b[mci::models::LH_secretion_model_interv::PARAM_LH_0] = params[PARAM_LH_0_b];
    params_b[mci::models::LH_secretion_model_interv::PARAM_sigma_obs] = params[PARAM_sigma_obs_b];
    //params_b[mci::models::LH_secretion_model_interv::PARAM_dim] = params[PARAM_dim_b];
    model_obj_b = new mci::models::LH_secretion_model_interv(params_b, r, iparam_names_b);
    
    
    variables_num = model_obj_b->get_dim() + model_obj_a->get_dim();
    
    
}

/** Get all model parameters **/
const model_params_t LH_secretion_across_studies_model::get_parameters() {
    return params;
}

/** Get parameter with identifier str **/
const double LH_secretion_across_studies_model::get_parameter(const std::string& str) {
    return params[str];
}

const size_t LH_secretion_across_studies_model::get_dim() {
    return this->variables_num;
}

const double LH_secretion_across_studies_model::get_dt() {
    return this->params[LH_secretion_across_studies_model::PARAM_dt];
}


/** Set the model parameters **/
void LH_secretion_across_studies_model::set_parameters(const model_params_t& p) {
    params = p;
}


/** update the inferred model parameters**/
void LH_secretion_across_studies_model::update_i_parameters(const Eigen::VectorXd& pp) {
    
    size_t k = 0;
    for (auto str : iparam_names) {
        params[str] = pp(k);
        k++;
        
        
    }
    
    Eigen::VectorXd p_model_a(5), p_model_b(9) ;
    Eigen::VectorXd p_d(4) ;
    break_i_parameters(pp, p_model_a, p_d);
    
    // update parameters in model 1
    model_obj_a->update_i_parameters(p_model_a);
    
    // update parameters in model 2
    model_obj_b->update_i_parameters(pp);
    
}

/** return values of inferred model parameters**/
const Eigen::VectorXd LH_secretion_across_studies_model::get_i_parameters() {
    
    
    Eigen::VectorXd i_params(iparam_names.size());
    size_t k = 0;
    for (auto str : iparam_names) {
        i_params(k) = params[str];
        k++;
    }
    return i_params;
    
}


/** return logarithm of proposal density p(theta_proposed | theta_current) **/
const double LH_secretion_across_studies_model::log_proposal_density(const Eigen::VectorXd& theta_current, const Eigen::VectorXd& theta_proposed) {
    return 0.0;
}

/** return logarithm of parameter prior density \pi(theta) **/
const double LH_secretion_across_studies_model::log_prior_density(const Eigen::VectorXd& theta) {
    
    double val = 0.0;
    Eigen::VectorXd theta_model(5);
    Eigen::VectorXd theta_d(4) ;
    break_i_parameters(theta, theta_model, theta_d);

    // get contribution for model 1 parameters
    val += model_obj_a->log_prior_density(theta_model);
    
    // get contribution for model 2 parameters
    val += model_obj_b->log_prior_density(theta_d);
    
    return val;
}

/** return logarithm of parameter prior density \pi(theta) **/
const Eigen::VectorXd LH_secretion_across_studies_model::grad_log_prior_density(const Eigen::VectorXd& theta) {
    
    size_t N = 5;
    Eigen::VectorXd prior_grad = Eigen::VectorXd::Zero(N);
    
    Eigen::VectorXd theta_model(5), theta_d(4) ;
    break_i_parameters(theta, theta_model, theta_d);
    
    // get values for model 1 parameters
    prior_grad.topRows<3>() = model_obj_a->grad_log_prior_density(theta_model);
    
    // get values for model 2 parameters
    prior_grad.bottomRows<2>() = model_obj_b->grad_log_prior_density(theta_d);
    
    return prior_grad;
    
}

/**add prior inofmration to FIM and its derivatives**/
const void LH_secretion_across_studies_model::add_prior_information(const Eigen::VectorXd& theta, Eigen::MatrixXd* G, std::vector<Eigen::MatrixXd>* G_k) {
    
    
    Eigen::MatrixXd G1 = G->topLeftCorner<3,3>();
    Eigen::MatrixXd G2 = G->bottomRightCorner<2,2>();
    
    Eigen::VectorXd theta_model(5);
    Eigen::VectorXd theta_d(4) ;
    break_i_parameters(theta, theta_model, theta_d);
    
    // get values for model 1 parameters
    model_obj_a->add_prior_information(theta_model, &G1);
    
    // get values for model 1 parameters
    model_obj_b->add_prior_information(theta_d, &G2);
    
    G->topLeftCorner<3,3>() = G1;
    G->bottomRightCorner<2,2>() = G2;
    
    
}

/** return logarithm of conditional observation density g(X_t,Y_t) **/
const double LH_secretion_across_studies_model::log_obs_density( const state_t & X, const Eigen::VectorXd& Y) {
    
    state_t X_a, X_b;
    Eigen::VectorXd Y_a, Y_b;
    double val = 0.0;
    
    break_up_state(X,X_a,X_b);
    break_up_data(Y, Y_a, Y_b);
    
    val += model_obj_a->log_obs_density(X_a,Y_a);
    val += model_obj_b->log_obs_density(X_b,Y_b);
    return val;
}

/** return gradient of the logarithm of conditional observation density g(X_t,Y_t) with respect to the parameters**/
const Eigen::VectorXd LH_secretion_across_studies_model::grad_log_obs_density( const state_t & X, const Eigen::VectorXd& Y) {
    //composite models shouldn't worry about this function
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(9);
    
    return grad;
}

/** return the grad logarithm of the joint density P(X,Y)**/
const Eigen::VectorXd LH_secretion_across_studies_model::grad_log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data,double* jll, double* mll, Eigen::MatrixXd* FIM_ret,std::vector<Eigen::MatrixXd>* G_k) {
    
    
    size_t N = 5;
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(N);
    //Eigen::MatrixXd FIM = Eigen::MatrixXd::Zero(N,N);
    double mlog_lik_a=0.0, mlog_lik_b=0.0,  jlog_lik_a=0.0, jlog_lik_b=0.0 ;
    
    trajectory_t data_a, data_b;
    trajectory_t X_a, X_b;
    Eigen::MatrixXd FIM_res_a = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd FIM_res_b = Eigen::MatrixXd::Zero(2,2);
    
    
    break_up_data(data, data_a, data_b);
    break_up_trajectory(X_ref_traj, X_a, X_b);
    
    grad.topRows<3>() = model_obj_a->grad_log_joint_density(X_a, data_a,  &jlog_lik_a, &mlog_lik_a, &FIM_res_a);
    grad.bottomRows<2>() = model_obj_b->grad_log_joint_density(X_b, data_b,  &jlog_lik_b, &mlog_lik_b, &FIM_res_b);
    
    
   
    if (jll) {
        *jll = jlog_lik_a+jlog_lik_b;
    }
    if (mll) {
        *mll = mlog_lik_a+mlog_lik_b;
    }
    if (FIM_ret) {
        
        *FIM_ret =  Eigen::MatrixXd::Zero(5,5);
        //std::cout << *FIM_ret << std::endl;
        //std::cout << "--------" << std::endl;
        FIM_ret->topLeftCorner<3,3>() = FIM_res_a;
        //std::cout << *FIM_ret << std::endl;
        //std::cout << "--------" << std::endl;
        FIM_ret->bottomRightCorner<2,2>() = FIM_res_b;
        
    }
    
    return grad;
    
}


/** return logartihm of transtion density f(X_t2 | X_t1)**/
const double LH_secretion_across_studies_model::log_transition_density(const state_t& X_t2, const state_t& X_t1, const size_t& t) {
    double val = 0.0;
    state_t X_a_t2, X_a_t1, X_b_t2, X_b_t1;
    
    break_up_state(X_t2, X_a_t2, X_b_t2);
    break_up_state(X_t1, X_a_t1, X_b_t1);
    
    val += model_obj_a->log_transition_density(X_a_t2,X_a_t1, t);
    val += model_obj_b->log_transition_density(X_b_t2,X_b_t1, t);
    
    return val;
    
}



/** return the logarithm of the joint density P(X,Y)**/
const double LH_secretion_across_studies_model::log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* ll) {
    
    double log_lik_a=0.0, log_lik_b = 0.0;
    double ret_val=0.0;
    trajectory_t data_a, data_b;
    trajectory_t X_a, X_b;
    
    break_up_data(data, data_a, data_b);
    break_up_trajectory(X_ref_traj, X_a, X_b);
    
    ret_val += model_obj_a->log_joint_density(X_a, data_a,  &log_lik_a);
    ret_val += model_obj_b->log_joint_density(X_b, data_b,  &log_lik_b);
    
    
    if (ll) {
        *ll = log_lik_a+log_lik_b;
    }
    
    return ret_val;
}

/** generate initial state **/
const state_t LH_secretion_across_studies_model::generate_ic(const double t_0) {
    
    
    // initialise state vector
    state_t ic(this->variables_num);
    
    //ic.topRows(params[composite_model::PARAM_dim_a]) = model_obj_a->generate_ic(t_0);
    //ic.bottomRows(params[composite_model::PARAM_dim_b]) = model_obj_b->generate_ic(t_0);
    ic.topRows<23>() = model_obj_a->generate_ic(t_0);
    ic.bottomRows<23>() = model_obj_b->generate_ic(t_0);
    return ic;
}


/** run model from t0 to t1 using X_0 as the initial condition **/
const trajectory_t LH_secretion_across_studies_model::run(const state_t& X0, const double& t_start, const double& t_end, const trajectory_t& data) {
    
    trajectory_t data_a, data_b;
    state_t X0_a, X0_b;
    
    break_up_data(data, data_a, data_b);
    break_up_state(X0, X0_a, X0_b);
    
    trajectory_t X_a = model_obj_a->run(X0_a, t_start, t_end, data_a);
    trajectory_t X_b = model_obj_b->run(X0_b, t_start, t_end, data_b);
    
    return join_trajectories(X_a,X_b); // return
}

/** run model from t0 to t1 using X_ref_traj as a reference trajectory**/
const trajectory_t LH_secretion_across_studies_model::run(const trajectory_t& X_ref_traj, const state_t& X0, const trajectory_t& data) {
    
    trajectory_t data_a, data_b;
    trajectory_t X_ref_a, X_ref_b;
    state_t X0_a, X0_b;
    
    break_up_data(data, data_a, data_b);
    break_up_state(X0, X0_a, X0_b);
    break_up_trajectory(X_ref_traj, X_ref_a, X_ref_b);
    
    trajectory_t X_a = model_obj_a->run(X_ref_a, X0_a, data_a);
    trajectory_t X_b = model_obj_b->run(X_ref_b, X0_b, data_b);
    
    return join_trajectories(X_a,X_b); // return
    
}


inline void LH_secretion_across_studies_model::break_up_data(const trajectory_t& data,  trajectory_t& data_a,  trajectory_t& data_b) {
 
    //data for model B are in the top rows of data
    data_b = data.topRows<2>();
    
    //data for model A are in the bottom rows
    data_a = data.bottomRows<2>();

}

inline void LH_secretion_across_studies_model::break_up_data(const Eigen::VectorXd& Y,  Eigen::VectorXd& Y_a,  Eigen::VectorXd& Y_b) {
    
    Y_b =Y.topRows<1>();
    Y_a =Y.bottomRows<1>();
}



inline void LH_secretion_across_studies_model::break_up_trajectory(const trajectory_t& X,  trajectory_t& X_a,  trajectory_t& X_b) {
   
    X_a = X.topRows<23>();//(params[composite_model::PARAM_dim_a]);
    X_b = X.bottomRows<23>();//(params[composite_model::PARAM_dim_b]);
}

inline trajectory_t LH_secretion_across_studies_model::join_trajectories(const trajectory_t& X_a,  const trajectory_t& X_b) {
    
    trajectory_t X(46, X_a.cols());
    X.topRows<23>() = X_a;
    X.bottomRows<23>() = X_b;
    
    return X;
}

inline void LH_secretion_across_studies_model::break_up_state(const state_t& X,  state_t& X_a,  state_t& X_b) {
    
    //std::cout << X << std::endl;
    
    X_a = X.topRows<23>();//(params[composite_model::PARAM_dim_a]);
    X_b = X.bottomRows<23>();//(params[composite_model::PARAM_dim_b]);
}


inline void LH_secretion_across_studies_model::break_i_parameters(const Eigen::VectorXd& theta, Eigen::VectorXd& theta_model, Eigen::VectorXd&  theta_d) {
    
    theta_model(0) = theta(0); //B_1
    theta_model(1) = theta(1); //alpha
    theta_model(2) = theta(2); //f
    theta_model(3) = theta(5); //mu_tau_on
    theta_model(4) = theta(7); //mu_tau_off
    
    theta_d(0) = theta(3); //diff_B_1
    theta_d(1) = theta(4); //diff_f
    theta_d(2) = theta(6); // diff_mu_tau_on
    theta_d(3) = theta(8); // diff_mu_tau_off
   
}
