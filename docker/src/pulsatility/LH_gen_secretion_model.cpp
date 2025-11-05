//
//  LH_gen_secretion_model.cpp
//  pulsatility
//
//  Created by Margaritis Voliotis on 30/10/2024.
//

#include "LH_gen_secretion_model.hpp"


#include <boost/math/special_functions/gamma.hpp>


using namespace mci;
using namespace mci::models;

/** static members definitions**/
const std::string mci::models::LH_gen_secretion_model::PARAM_dt = "dt";
const std::string mci::models::LH_gen_secretion_model::PARAM_B_1 = "k";
const std::string mci::models::LH_gen_secretion_model::PARAM_mu_tau_on =  "mu_tau_on";
const std::string mci::models::LH_gen_secretion_model::PARAM_mu_tau_off = "mu_tau_off";
const std::string mci::models::LH_gen_secretion_model::PARAM_std_tau_on = "std_tau_on";
const std::string mci::models::LH_gen_secretion_model::PARAM_std_tau_off ="std_tau_off";

const std::string mci::models::LH_gen_secretion_model::PARAM_f = "f";
const std::string mci::models::LH_gen_secretion_model::PARAM_alpha_1 = "d";
const std::string mci::models::LH_gen_secretion_model::PARAM_tau_1 = "tau_1";
const std::string mci::models::LH_gen_secretion_model::PARAM_LH_0 = "LH_0";
const std::string mci::models::LH_gen_secretion_model::PARAM_sigma_obs = "sigma_obs";
const std::string mci::models::LH_gen_secretion_model::PARAM_observed_variable_index_0 = "observed_variable_index_0";
//const std::string mci::models::dual_mode_model_gradinfo::PARAM_dim = "dim";
//const std::string mci::models::dual_mode_model_gradinfo::PARAM_grad_dim = "grad_dim";
const std::string mci::models::LH_gen_secretion_model::NAME = "LH_gen_secretion_model";

const std::string mci::models::LH_gen_secretion_model::PARAM_prior_mean_B_1 = "k_normal_prior_mean";
const std::string mci::models::LH_gen_secretion_model::PARAM_prior_std_B_1 = "k_normal_prior_std";
const std::string mci::models::LH_gen_secretion_model::PARAM_prior_mean_alpha_1 = "d_normal_prior_mean";
const std::string mci::models::LH_gen_secretion_model::PARAM_prior_std_alpha_1 = "d_normal_prior_std";
const std::string mci::models::LH_gen_secretion_model::PARAM_f_Beta_prior_alpha = "f_Beta_prior_alpha";
const std::string mci::models::LH_gen_secretion_model::PARAM_f_Beta_prior_beta = "f_Beta_prior_beta";
/** constructor **/
LH_gen_secretion_model::LH_gen_secretion_model(const model_params_t& p,
                                       const mci::stochastic& r,
                                       const std::vector<std::string>& pi) {
    params = p;
    rng = r;
    variables_num = 23; //params[dual_mode_model_gradinfo::PARAM_dim];
    grad_num = 3;
    obs_variable_index(0) = params[LH_gen_secretion_model::PARAM_observed_variable_index_0];
    
    std::cout<< "s " << obs_variable_index(0) <<std::endl;
    
    iparam_names = pi;
    iparam_id = {0,1,2,3,4};
    size_t k = 0;
    
    
    /*for (size_t t = 0;t<=4;t++) {
        std::cout <<  iparam_names[t] << " " ;
    }
    
    std::cout<<std::endl;*/
    
    for (auto str : iparam_names) {
        
        if ( str.compare(PARAM_B_1) == 0 )
            iparam_id[0] = k;
        else if ( str.compare(PARAM_alpha_1) == 0 )
            iparam_id[1] = k;
        else if ( str.compare(PARAM_f) == 0 )
            iparam_id[2] = k;
        else if ( str.compare(PARAM_mu_tau_on) == 0 )
            iparam_id[3] = k;
        else if ( str.compare(PARAM_mu_tau_off) == 0 )
            iparam_id[4] = k;
        
        k++;
        
        
    }
   
    //std::cout << iparam_id[0] << " " << iparam_id[1] <<" " << " " << iparam_id[2];
    
}

/** Get all model parameters **/
const model_params_t LH_gen_secretion_model::get_parameters() {
    return params;
}

/** Get parameter with identifier str **/
const double LH_gen_secretion_model::get_parameter(const std::string& str) {
    return params[str];
}

const size_t LH_gen_secretion_model::get_dim() {
    return this->variables_num;
}

const double LH_gen_secretion_model::get_dt() {
    return this->params[LH_gen_secretion_model::PARAM_dt];
}


/** Set the model parameters **/
void LH_gen_secretion_model::set_parameters(const model_params_t& p) {
    params = p;
}


/** return logarithm of parameter prior density \pi(theta) **/
const double LH_gen_secretion_model::log_prior_density(const Eigen::VectorXd& theta) {
    
    double val = 0.0;
    // uniform priors suppported currently
    //for (auto str : iparam_names) {
    //}
    /*val = (theta(0) + log10(1+theta(7)) + theta(1) + log10(1+theta(8)) )*log(10)
     - 1./(2.*pow(0.25,2))*pow( theta(2) - 1,2)
     - 1./(2.*pow(0.3,2))*pow( theta(3) - 1.5,2)
     - 1./(2.*pow(0.3,2))*pow( theta(4) + log10(1+theta(11)) - 1,2)
     - 1./(2.*pow(0.3,2))*pow( theta(5) + log10(1+theta(12))- 0.9322,2)//1./(2.*pow(??,2))*pow(pow(10,theta(5)) - 0.1169,2) //basal secretion
     - 1./(2.*pow(9.3,2))*pow( log(2)/pow(10,theta(6)) - 80,2); // clearance */
    
    
    val = + ( theta(iparam_id[3]) + theta(iparam_id[4]) )*log(10)
    //- theta(3) - 2*log(1 + exp(-theta(3)) )
    //- theta(4) - 2*log(1 + exp(-theta(4)) )
    // - 1./(2.*pow(0.25,2))*pow( theta(2) - 1.,2)
    // - 1./(2.*pow(0.3,2))*pow( theta(3) - 1.5,2)
    //- 1./(2.*pow(0.3,2))*pow( theta(0) - 2.6,2)
    - 1./(2.*pow(params[PARAM_prior_std_B_1],2))*pow( theta(iparam_id[0]) - params[PARAM_prior_mean_B_1],2)
    //- 1./(2.*pow(0.3,2))*pow( theta(5) - 0.9322,2)//1./(2.*pow(??,2))*pow(pow(10,theta(5)) - 0.1169,2) //basal secretion
    //- 1./(2.*pow(9.3,2))*pow( log(2)/pow(10,theta(1)) - 80.0,2); // clearance
    - 1./(2.*pow(params[PARAM_prior_std_alpha_1],2))*pow( log(2.)/pow(10.,theta(iparam_id[1])) - params[PARAM_prior_mean_alpha_1],2) - theta(iparam_id[1])*log(10.)
    //- 1./(2.*pow(params[PARAM_prior_std_alpha_1],2))*pow( log(2.)/pow(10.,theta(iparam_id[1])) - 80.0,2)-theta(iparam_id[1])*log(10.)
    //  - 1./(2.*pow(0.3,2))*pow( theta(1) + 2.0623,2) // clearance
    -(params[PARAM_f_Beta_prior_alpha]-1.0)*log(1+exp(-theta(iparam_id[2])))
    -(params[PARAM_f_Beta_prior_beta]-1.0)*log(1+exp(+theta(iparam_id[2])))
    - theta(iparam_id[2]) - 2*log(1 + exp(-theta(iparam_id[2])) ) ;
    
    
    return val;
}

/** return grad of logarithm of parameter prior density \pi(theta) **/
const Eigen::VectorXd LH_gen_secretion_model::grad_log_prior_density(const Eigen::VectorXd& theta) {
    
        size_t N = grad_num;
        Eigen::VectorXd prior_grad = Eigen::VectorXd::Zero(N);
        

        prior_grad(0) = - 1./(pow(params[PARAM_prior_std_B_1],2))*( theta(iparam_id[0]) - params[PARAM_prior_mean_B_1]);
        prior_grad(1) = log(10.) * (-1 + 1./(pow(params[PARAM_prior_std_alpha_1],2))*log(2.)/pow(10,theta(iparam_id[1]))*( log(2.)/pow(10,theta(iparam_id[1])) - params[PARAM_prior_mean_alpha_1] ) ); // clearance
    
        //prior_grad(iparam_id[1]) = log(10.) * (-1 + 1./(pow(params[PARAM_prior_std_alpha_1],2))*log(2.)/pow(10,theta(iparam_id[1]))*( log(2.)/pow(10,theta(iparam_id[1])) - 80.0) ); // clearance
    
        prior_grad(2) = - 1 + 2./(1+exp(theta(iparam_id[2])))
                                   + (params[PARAM_f_Beta_prior_alpha]-1.0) /(1.0+exp(theta(iparam_id[2])))
                                   - (params[PARAM_f_Beta_prior_beta]-1.0) /(1.0+exp(-theta(iparam_id[2])));
    
   
        return prior_grad;
        
    }

/**add prior inofmration to FIM and its derivatives**/
const void LH_gen_secretion_model::add_prior_information(const Eigen::VectorXd& theta, Eigen::MatrixXd* G, std::vector<Eigen::MatrixXd>* G_k) {
    
    size_t N = grad_num;
    Eigen::VectorXd prior_grad = Eigen::VectorXd::Zero(N);
    
    
    (*G)(0,0) -=  - 1./(pow(params[PARAM_prior_std_B_1],2));
    (*G)(1,1) -=  log(2)*log(10.)*log(10.)/(pow(params[PARAM_prior_std_alpha_1],2)) * ( params[PARAM_prior_mean_alpha_1]*pow(10,-theta(iparam_id[1])) - log(4)*pow(10,-2*theta(iparam_id[1]) ) );
    //(*G)(iparam_id[1],iparam_id[1]) -=  log(2)*log(10.)*log(10.)/(pow(params[PARAM_prior_std_alpha_1],2)) * ( 80*pow(10,-theta(iparam_id[1])) - log(4)*pow(10,-2*theta(iparam_id[1]) ) );
    (*G)(2,2) -=  - 2. * exp(theta(iparam_id[2])) /pow(1+exp(theta(iparam_id[2])),2)
                                        - (params[PARAM_f_Beta_prior_alpha]-1.0) * exp(theta(iparam_id[2])) /pow(1+exp(theta(iparam_id[2])),2)
                                        - (params[PARAM_f_Beta_prior_beta]-1.0) * exp(-theta(iparam_id[2])) /pow(1+exp(-theta(iparam_id[2])),2);

    
}


/** return logarithm of proposal density p(theta_proposed | theta_current) **/
const double LH_gen_secretion_model::log_proposal_density(const Eigen::VectorXd& theta_current, const Eigen::VectorXd& theta_proposed) {
    return 0.0;
}

/** update the inferred model parameters**/
void LH_gen_secretion_model::update_i_parameters(const Eigen::VectorXd& pp) {
    
    size_t k = 0;
    for (auto str : iparam_names) {
        if ( str.compare(PARAM_f) == 0 ) {
            params[str] = 1./(1.+exp(-pp(k)));
      //  }else if ( str.compare(PARAM_mu_tau_off) == 0  ) {
      //      params[str] = 5. + (240.-5.)*1./(1.+exp(-pp(k)));//pow(10, 0.6990 + (2.3802-0.6990)*1./(1.+exp(-pp(k)))   )//
      //  } else if (str.compare(PARAM_mu_tau_on) == 0 ) {
      //      params[str] = 1. + (240.-1.)*1./(1.+exp(-pp(k)));//pow(10, 0.6990 + (2.3802-0.6990)*1./(1.+exp(-pp(k)))   )//
        } else
            params[str] = pow(10,pp(k));
        k++;
        
        
    }
    
  //  std::cout << params[PARAM_mu_tau_off] << " "<< std::endl;
    
    cache_haz_off.clear();
    cache_haz_on.clear();
}

/** return values of inferred model parameters**/
const Eigen::VectorXd LH_gen_secretion_model::get_i_parameters() {
    
    
    Eigen::VectorXd i_params(iparam_names.size());
    size_t k = 0;
    for (auto str : iparam_names) {
        if ( str.compare(PARAM_f) == 0 )
            i_params(k) = log(params[str]/(1-params[str]));
        else
            i_params(k) = log10(params[str]);
        
        k++;
    }
    return i_params;
    
}





/** return logarithm of conditional observation density g(X_t,Y_t) **/
const double LH_gen_secretion_model::log_obs_density( const state_t & X, const Eigen::VectorXd& Y) {
    
    double sigma_obs = params[LH_gen_secretion_model::PARAM_sigma_obs];
    
    
    double s_LH  = sigma_obs*abs(X(obs_variable_index(0),0));
    
    double val  = - log(2*M_PI)/2
                  - log(s_LH)
                  - 1./(2.*pow(s_LH,2)) * pow( Y(0) - X(obs_variable_index(0),0),2);
    
    return val;
}

/** return gradient of the logarithm of conditional observation density g(X_t,Y_t) with respect to the parameters**/
Eigen::MatrixXd LH_gen_secretion_model::FIM( const state_t & X) {
    
    size_t N = grad_num;
    Eigen::MatrixXd FIM = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(N);
    double alpha_1 = params[LH_gen_secretion_model::PARAM_alpha_1];
    double B_1 = params[LH_gen_secretion_model::PARAM_B_1];
    double f = params[LH_gen_secretion_model::PARAM_f];
    double sigma_obs = params[LH_gen_secretion_model::PARAM_sigma_obs];
    
    grad(0) =  X(18,0) * log(10) * B_1;
    grad(1) =  X(19,0) * log(10) * alpha_1;
    grad(2) =  X(21,0) * f * (1-f);
    
    //double s_LH  = pow(0.15,2);
    //FIM = ( grad * grad.transpose() ) / s_LH;

    // double s_LH  =  sigma_obs*abs(X(obs_variable_index(0),0));
     FIM = ( grad * grad.transpose() ) / pow(X(obs_variable_index(0),0),2) * (2 + 1./pow(sigma_obs,2));
    
    return FIM;
}

/** return gradient of the logarithm of conditional observation density g(X_t,Y_t) with respect to the parameters**/
const Eigen::VectorXd LH_gen_secretion_model::grad_log_obs_density( const state_t & X, const Eigen::VectorXd& Y) {
    
    size_t N = grad_num;
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(N);
    double alpha_1 = params[LH_gen_secretion_model::PARAM_alpha_1];
    double B_1 = params[LH_gen_secretion_model::PARAM_B_1];
    double f = params[LH_gen_secretion_model::PARAM_f];
    double sigma_obs = params[LH_gen_secretion_model::PARAM_sigma_obs];
    
    
    /*double s_LH  = pow(0.15,2);//
    grad(0) = (Y(0) - X(obs_variable_index(0),0))/s_LH * X(18,0) * log(10) * B_1;
    grad(1) = (Y(0) - X(obs_variable_index(0),0))/s_LH * X(19,0) * log(10) * alpha_1;
    grad(2) = (Y(0) - X(obs_variable_index(0),0))/s_LH * X(21,0) * f * (1-f);
    */
    
    double s_LH  = sigma_obs*abs(X(obs_variable_index(0),0)) ;
    grad(0) =  - 1./s_LH*sigma_obs*X(18,0) * log(10) * B_1
            + 1./(pow(s_LH,3))*pow(Y(0) - X(obs_variable_index(0),0),2) * sigma_obs * X(18,0) * log(10) * B_1
            + 1./(pow(s_LH,2))*(Y(0) - X(obs_variable_index(0),0)) * X(18,0) * log(10) * B_1;
    grad(1) = - 1./s_LH*sigma_obs*X(19,0) * log(10) * log(10) * alpha_1
            + 1./(pow(s_LH,3))*pow(Y(0) - X(obs_variable_index(0),0),2) * sigma_obs * X(19,0) * log(10) * alpha_1
            + 1./(pow(s_LH,2))*(Y(0) - X(obs_variable_index(0),0)) * X(19,0) * log(10) * alpha_1;
    grad(2) = - 1./s_LH*sigma_obs*X(21,0) * f * (1-f)
            + 1./(pow(s_LH,3))*pow(Y(0) - X(obs_variable_index(0),0),2) * sigma_obs * X(21,0) * f * (1-f)
            + 1./(pow(s_LH,2))*(Y(0) - X(obs_variable_index(0),0)) * X(21,0) * f * (1-f);
    
  
    
    return grad;
}

/** return logartihm of transtion density f(X_t2 | X_t1)**/
const double LH_gen_secretion_model::log_transition_density(const state_t& X_t2, const state_t& X_t1, const size_t& t) {
    
    double val = 0.0  - 0.5*pow( X_t2(15),2) - log(2*M_PI)/2;
    return val;
    
}


/** return the logarithm of the joint density P(X,Y)**/
const double LH_gen_secretion_model::log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data, double* ll) {
    
    double dt = params[LH_gen_secretion_model::PARAM_dt];
    
    trajectory_t Xfinal = run( X_ref_traj, X_ref_traj.col(0), data); //run using reference trajecotry
    
    size_t timepoints_data = data.cols();
    size_t timepoints_model = X_ref_traj.cols();
    
    double log_lik = 0;
    double f_t = 0;
    double ret_val=0;
    
    //first timepoint
    Eigen::VectorXd obs_v = data.col(0).head(data.rows()-1);
    double obs_t = data(data.rows()-1, 0) ;
    log_lik = log_obs_density(Xfinal.col(0), obs_v);
    
    //loop over timepoints
    
    // get time and value of next observation
    size_t k = 1;
    obs_v = data.col(k).head(data.rows()-1);
    obs_t = data(data.rows()-1, k) ;
    for ( size_t tii = 1;  tii<timepoints_model ; tii++ ) {
        
        double model_t = X_ref_traj(X_ref_traj.rows()-1, tii);
        
        if ( abs(obs_t - model_t) <dt/2.) {
            
            // - calculate cond obs density
            log_lik = log_lik + log_obs_density(Xfinal.col(tii), obs_v);
            k = k + 1;
            // get time and value of next observation
            if (k<timepoints_data) {
                obs_v = data.col(k).head(data.rows()-1);
                obs_t = data(data.rows()-1, k) ;
            }
        }
        
        f_t = f_t +
        log_transition_density(Xfinal.col(tii), Xfinal.col(tii-1), tii-1);
        
        
    }
    
    ret_val = log_lik +  f_t;
    if (ll) {
        *ll = log_lik;
    }
    
    return ret_val;
    
}



                                                                       
                                                                       
/** return the grad logarithm of the joint density P(X,Y)**/
const Eigen::VectorXd LH_gen_secretion_model::grad_log_joint_density(const trajectory_t& X_ref_traj, const trajectory_t& data,double* jll, double* mll, Eigen::MatrixXd *FIM_ret, std::vector<Eigen::MatrixXd>* G_k) {
    
    size_t N = grad_num;
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd FIM_res = Eigen::MatrixXd::Zero(N,N);
    double dt = params[LH_gen_secretion_model::PARAM_dt];
    
    trajectory_t Xfinal = run( X_ref_traj, X_ref_traj.col(0), data); //run using reference trajecotry
    
    size_t timepoints_data = data.cols();
    size_t timepoints_model = X_ref_traj.cols();
    
    double log_lik = 0;
    double f_t = 0;
    double joint_log_lik=0;
    
    //first timepoint
    Eigen::VectorXd obs_v = data.col(0).head(data.rows()-1);
    double obs_t = data(data.rows()-1, 0) ;
    log_lik = log_obs_density(Xfinal.col(0), obs_v);
    grad = grad_log_obs_density(Xfinal.col(0), obs_v);
    FIM_res = FIM(Xfinal.col(0));
    //loop over timepoints
     //std::cout<<"0 --- " <<grad.transpose() << std::endl;
    // get time and value of next observation
    size_t k = 1;
    obs_v = data.col(k).head(data.rows()-1);
    obs_t = data(data.rows()-1, k) ;
    for ( size_t tii = 1;  tii<timepoints_model ; tii++ ) {
        
        double model_t = X_ref_traj(X_ref_traj.rows()-1, tii);
        
        if ( abs(obs_t - model_t) <dt/2.) {
            
            // - get cond obs density
            log_lik = log_lik + log_obs_density(Xfinal.col(tii), obs_v);
            // - get the grad info
            grad += grad_log_obs_density(Xfinal.col(tii), obs_v);
            
            FIM_res += FIM(Xfinal.col(tii));
            //std::cout<<tii<< "   "<<grad.transpose() << std::endl;
            k = k + 1;
            // get time and value of next observation
            if (k<timepoints_data) {
                obs_v = data.col(k).head(data.rows()-1);
                obs_t = data(data.rows()-1, k) ;
            }
        }
        
        f_t = f_t +
        log_transition_density(Xfinal.col(tii), Xfinal.col(tii-1), tii-1);
        
        
    }
    
    joint_log_lik = log_lik +  f_t;
    if (jll) {
        *mll = log_lik;
        *jll = joint_log_lik;
    }
    if (FIM_ret) {
        *FIM_ret = FIM_res;
        
    }
    
    return grad;
    
}




/** generate initial state **/
const state_t LH_gen_secretion_model::generate_ic(const double t_0) {
    
    double dt           = params[LH_gen_secretion_model::PARAM_dt];
    double mu_tau_on    = params[LH_gen_secretion_model::PARAM_mu_tau_on] ;
    double mu_tau_off   = params[LH_gen_secretion_model::PARAM_mu_tau_off] ;
    double std_tau_on   = params[LH_gen_secretion_model::PARAM_std_tau_on] ;
    double std_tau_off  = params[LH_gen_secretion_model::PARAM_std_tau_off] ;
    
    //std::cout<< (1+params[dual_mode_model_gradinfo::PARAM_diff_std_tau_on]) << " " << (1+params[dual_mode_model_gradinfo::PARAM_diff_std_tau_on])<< std::endl;
    double tau_1 = params[LH_gen_secretion_model::PARAM_tau_1];
    double B_1 = params[LH_gen_secretion_model::PARAM_B_1] ;
    double alpha_1 = params[LH_gen_secretion_model::PARAM_alpha_1];
    double LH_0         = params[LH_gen_secretion_model::PARAM_LH_0];
    double f = params[LH_gen_secretion_model::PARAM_f];//1./(exp(params[dual_mode_model_gradinfo::PARAM_f])+1);
    
    double prob_off_event = 0.0;
    double prob_on_event = 0.0;
    double tmp = 0.;
    double u_sample = 0.0;
    double n_sample = 0.0;
    double n_sample_BS = 0.0;
    double u_sample_progress = 0.0;
    // initialise state vector
    state_t ic(this->variables_num);
    
    u_sample = rng.rand();
    n_sample = rng.rand();
    n_sample_BS = rng.randn();
    u_sample_progress = rng.rand();
    
    //initialize generator randomly
    u_sample < mu_tau_on/(mu_tau_on + mu_tau_off) ? ic(4) =  1 : ic(4) =  0 ;
    
    if (ic(4) == 1) { //if generator is on
        
        tmp = floor(get_quantile_on(u_sample_progress, mu_tau_on, std_tau_on)/dt )*dt;
        //tmp = floor(get_quantile_on(u_sample_progress, mu_tau_on, std_tau_on))*dt;
        
        prob_off_event  = dt * propensity_on(tmp, mu_tau_on, std_tau_on);
        prob_on_event = 0;
        
        ic(2) = tmp;
        ic(3) = tmp;
        
        
    }
    else { //if generator is off
        
        tmp = floor(get_quantile(u_sample_progress, mu_tau_off, std_tau_off)/dt )*dt;
        //tmp = floor(get_quantile(u_sample_progress, mu_tau_off, std_tau_off))*dt;
        
        prob_on_event = dt * propensity(tmp, mu_tau_off, std_tau_off);
        prob_off_event = 0;
        
        ic(2)  =  tmp;
        ic(3)  =  tmp;
        
    }
    
    ic(14) = n_sample_BS;//1/(1+exp(-1*ic(15)));
    double aaa = (ic(4)*f + 1/(1+exp(-1*ic(14)))*(1-f))/tau_1;//1;
    
    ic(5) = aaa*exp(-tau_1*tmp) +   (ic(4)*f + 1/(1+exp(-1*ic(14)))*(1-f) )  * (1 - exp(-tau_1*tmp))/tau_1;
    ic(6) = 5*LH_0*n_sample*exp(-alpha_1*tmp) + B_1*(ic(4)*f + 1/(1+exp(-1*ic(14)))*(1-f) ) /tau_1
    + exp(-alpha_1*tmp)*B_1/(tau_1/alpha_1-1)*(aaa-(ic(4)*f + 1/(1+exp(-1*ic(14)))*(1-f) ) /alpha_1)
    + exp(-tau_1*tmp)*B_1/(tau_1/alpha_1-1)*((ic(4)*f + 1/(1+exp(-1*ic(14)))*(1-f) )/tau_1-aaa);
    
    /*** grad infomration ***/
    //dLH/dB_1
    ic(18) = ((1 - f)/(exp(-ic(14)) + 1) + ic(4)*f)/tau_1
    +(exp(-alpha_1*dt)*(aaa + ((f - 1)/(exp(-ic(14)) + 1) - ic(4)*f)/alpha_1))/(tau_1/alpha_1 - 1)
    -(exp(-dt*tau_1)*(aaa + ((f - 1)/(exp(-ic(14)) + 1) - ic(4)*f)/tau_1))/(tau_1/alpha_1 - 1);
    //dLH/dalpha_1
    ic(19) =  - dt*5*LH_0*n_sample*exp(-alpha_1*dt)
    + (B_1*exp(-alpha_1*dt)*((f - 1)/(exp(-ic(14)) + 1) - ic(4)*f))/(alpha_1*(alpha_1 - tau_1)) + (B_1*tau_1*exp(-alpha_1*dt)*(aaa + ((f - 1.0)/(exp(-ic(14)) + 1) - ic(4)*f)/alpha_1))/pow(alpha_1 - tau_1,2) + (B_1*alpha_1*dt*exp(-alpha_1*dt)*(aaa + ((f - 1.0)/(exp(-ic(14)) + 1) - ic(4)*f)/alpha_1))/(alpha_1 - tau_1)
    - (B_1*tau_1*exp(-dt*tau_1)*(aaa + ((f - 1)/(exp(-ic(14)) + 1) - ic(4)*f)/tau_1))/pow(alpha_1 - tau_1,2);
    //dPG/df
    ic(20) =  - ((ic(4) - 1/(exp(-ic(14)) + 1))*(exp(-dt*tau_1) - 1))/tau_1;
    //dLH/df
    ic(21) = (B_1*(ic(4) - 1/(exp(-ic(14)) + 1)))/tau_1
    -(B_1*exp(-alpha_1*dt)*((ic(4) - 1/(exp(-ic(14)) + 1))/alpha_1 - ic(20)))/(tau_1/alpha_1 - 1)
    + (B_1*exp(-dt*tau_1)*((ic(4) - 1/(exp(-ic(14)) + 1))/tau_1 - ic(20)))/(tau_1/alpha_1 - 1);
    
    
    ic(7) = rng.rand(); // sample random number
    ic(8) = rng.rand(); // sample random number
    
    ic(7) < prob_on_event  ? ic(0)=1 : ic(0)=0 ; //on event indicator
    ic(8) < prob_off_event ? ic(1)=1 : ic(1)=0 ; //off event indicator
    
    ic(9)  = prob_off_event;
    ic(10) = prob_on_event;
    
    ic(11) = u_sample; //rng for generator initial state
    ic(12) = n_sample; // rng for LH initial state
    ic(13) = u_sample_progress; //rng for  time since previous event
    ic(15) = n_sample_BS;//(n_sample_BS > 0.5) ? 10 : -10; // rng for BS initial state
    ic(ic.size()-1) = t_0; //initial time
    
    return ic;
}

/** run model from t0 to t1 using X_0 as the initial condition **/
const trajectory_t LH_gen_secretion_model::run(const state_t& X0, const double& t_start, const double& t_end, const trajectory_t& data) {
    
    double dt = params[LH_gen_secretion_model::PARAM_dt];
    double mu_tau_on    = params[LH_gen_secretion_model::PARAM_mu_tau_on] ;
    double mu_tau_off   = params[LH_gen_secretion_model::PARAM_mu_tau_off] ;
    double std_tau_on   = params[LH_gen_secretion_model::PARAM_std_tau_on] ;
    double std_tau_off  = params[LH_gen_secretion_model::PARAM_std_tau_off] ;
    double tau_1 = params[LH_gen_secretion_model::PARAM_tau_1];
    
    
    double B_1 = params[LH_gen_secretion_model::PARAM_B_1] ;
    double alpha_1 = params[LH_gen_secretion_model::PARAM_alpha_1];
    double f = params[LH_gen_secretion_model::PARAM_f];
    
    double prob_on_event = 0;
    double prob_off_event = 0;
    
    // length of trajecotry
    int timepoints = round( (t_end - t_start)/dt ) + 1;
    
    //initialise trajecotry
    size_t dim = this->variables_num;
    trajectory_t Xc(dim, timepoints);
    
    Xc.col(0) = X0;
    
    
    //loop
    for (size_t k = 0; k< timepoints-1 ; k++) {
        
        Xc(dim-1, k+1) = Xc(dim-1, k) + dt; //update time
        
        Xc(11, k+1) = Xc(11, k);
        Xc(12, k+1) = Xc(12, k);
        Xc(7, k+1) = rng.rand();
        Xc(8, k+1) = rng.rand();
        Xc(15, k+1) = rng.randn();
        
        if (Xc(4, k) == 1) { // if generator 'on'
            
            if ( Xc(1, k) == 1 ) { //if turn-off event occured
                Xc(4, k+1) = 0;  //change state of generator to off
                Xc(2, k+1) = 0;  //store time of event
                Xc(3, k+1) = 0;  //store time of event
                Xc(12, k+1)  = 1;
                
                
                //calculate prob 'on' event occurs
                prob_off_event  = 0;
                prob_on_event   = dt * propensity(Xc(3,k+1), mu_tau_off, std_tau_off);
                
                Xc(9,k+1)  = prob_off_event;
                Xc(10,k+1) = prob_on_event;
                
            }
            else {
                Xc(4, k+1) = Xc(4, k);
                Xc(2, k+1) = Xc(2, k)+dt;
                Xc(3, k+1) = Xc(3, k)+dt;
                
                
                //calculate prob 'off' event occurs
                prob_off_event  = dt * propensity_on(Xc(2,k+1), mu_tau_on, std_tau_on);
                prob_on_event   = 0;
                
                Xc(9,k+1)  = prob_off_event;
                Xc(10,k+1) = prob_on_event;
                
            }
            
        } else { // if generator 'off'
            
            if (Xc(0, k) == 1 ) {//if turn-on event occurent
                Xc(4, k+1) = 1; //change state of generator to on
                Xc(2, k+1) = 0;//store time of event
                Xc(3, k+1) = 0;//store time of event
                Xc(12, k+1)  = 1;
                
                
                //calculate prob 'off' event occurs
                prob_off_event  = dt * propensity_on(Xc(2,k+1), mu_tau_on, std_tau_on);
                prob_on_event   = 0;
                
                Xc(9,k+1) =prob_off_event;
                Xc(10,k+1) =prob_on_event;
                
            }
            else {
                Xc(4, k+1) = Xc(4, k);
                Xc(2, k+1) = Xc(2, k) + dt;
                Xc(3, k+1) = Xc(3, k) + dt;
                
                //calculate prob 'on' event occurs
                prob_off_event  = 0;
                prob_on_event   = dt * propensity(Xc(3,k+1), mu_tau_off, std_tau_off);
                
                Xc(9,k+1) =prob_off_event;
                Xc(10,k+1) =prob_on_event;
                
            }
            
        }
        
        Xc(7, k+1) < prob_on_event   ? Xc(0, k+1)=1 : Xc(0 ,k+1)=0 ;
        Xc(8, k+1) < prob_off_event  ? Xc(1, k+1)=1 : Xc(1, k+1)=0 ;
        
        Xc(5, k+1)      = Xc(5, k)*exp(-tau_1*dt)   +   (Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1-f)  ) * (1 - exp(-tau_1*dt))/tau_1;
        
        
        Xc(6, k+1) = Xc(6, k)*exp(-alpha_1*dt) + (B_1)*(Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1.0-f)  )/tau_1
        + exp(-alpha_1*dt)*(B_1)/(tau_1/alpha_1-1)*(Xc(5, k)-(Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1-f)  )/alpha_1)
        + exp(-tau_1*dt)*(B_1)/(tau_1/alpha_1-1)*((Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1-f)  )/tau_1-Xc(5, k));
        
        
        
        Xc(14, k+1) = Xc(14, k)  + dt/2 * (-  Xc(14, k) ) +   sqrt(dt) * Xc(15, k);
       
        /*** grad infomration ***/
        //dLH/dB_1
        Xc(18, k+1) = Xc(18, k)*exp(-alpha_1*dt) +  ((1.0 - f)/(exp(-Xc(14,k)) + 1) + Xc(4,k)*f)/tau_1
                     +(exp(-alpha_1*dt)*(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14,k)) + 1) - Xc(4,k)*f)/alpha_1))/(tau_1/alpha_1 - 1)
                     -(exp(-dt*tau_1)  *(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14,k)) + 1) - Xc(4,k)*f)/tau_1))/(tau_1/alpha_1 - 1);
        //dLH/dalpha_1
        Xc(19, k+1) = Xc(19, k)*exp(-alpha_1*dt) - dt*Xc(6, k)*exp(-alpha_1*dt)
        + (B_1*exp(-alpha_1*dt)*((f - 1)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f))/(alpha_1*(alpha_1 - tau_1))
        + (B_1*tau_1*exp(-alpha_1*dt)*(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/alpha_1))/pow(alpha_1 - tau_1,2)
        + (B_1*alpha_1*dt*exp(-alpha_1*dt)*(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/alpha_1))/(alpha_1 - tau_1)
        - (B_1*tau_1*exp(-dt*tau_1)*(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/tau_1))/pow(alpha_1 - tau_1,2);
        
    
        //dPG/df
        Xc(20, k+1) = Xc(20, k)*exp(-tau_1*dt) - ((Xc(4, k) - 1/(exp(-Xc(14, k)) + 1))*(exp(-dt*tau_1) - 1))/tau_1;
        //dLH/df
        Xc(21, k+1) = Xc(21, k)*exp(-alpha_1*dt) + (B_1*(Xc(4, k) - 1/(exp(-Xc(14, k)) + 1)))/tau_1
           -(B_1*exp(-alpha_1*dt)*((Xc(4, k) - 1/(exp(-Xc(14, k)) + 1))/alpha_1 - Xc(20, k)))/(tau_1/alpha_1 - 1)
           + (B_1*exp(-dt*tau_1)*((Xc(4, k) - 1/(exp(-Xc(14, k)) + 1))/tau_1 - Xc(20, k)))/(tau_1/alpha_1 - 1);
    }
    
    ///Xc.erase(Xc.begin()); // remove first entty
    
    return Xc.middleCols(1, timepoints-1); // return
}

/** run model from t0 to t1 using X_ref_traj as a reference trajectory**/
const trajectory_t LH_gen_secretion_model::run(const trajectory_t& X_ref_traj, const state_t& X0, const trajectory_t& data) {
    
    size_t timepoints = X_ref_traj.cols();
    
    double dt = params[LH_gen_secretion_model::PARAM_dt];
    double mu_tau_on    = params[LH_gen_secretion_model::PARAM_mu_tau_on] ;
    double mu_tau_off   = params[LH_gen_secretion_model::PARAM_mu_tau_off] ;
    double std_tau_on   = params[LH_gen_secretion_model::PARAM_std_tau_on] ;
    double std_tau_off  = params[LH_gen_secretion_model::PARAM_std_tau_off] ;
    double tau_1 = params[LH_gen_secretion_model::PARAM_tau_1];
    
    
    double B_1 = params[LH_gen_secretion_model::PARAM_B_1] ;
    double alpha_1 = params[LH_gen_secretion_model::PARAM_alpha_1];
    
    double f = params[LH_gen_secretion_model::PARAM_f];
    double LH_0         = params[LH_gen_secretion_model::PARAM_LH_0];
    
    
    double tmp = 0;
    
    double prob_on_event = 0;
    double prob_off_event = 0;
    
    //initialise trajecotry
    size_t dim = this->variables_num;
    trajectory_t Xc(dim, timepoints);
    Xc.col(0)  = X_ref_traj.col(0);
    
    Xc(11,0) < mu_tau_on/(mu_tau_on + mu_tau_off) ? Xc(4,0) =  1 : Xc(4,0) =  0;
    
    
    if (Xc(4, 0) == 1) { // if generator 'on'
        
        tmp = floor(get_quantile_on(Xc(13,0), mu_tau_on, std_tau_on)/dt )*dt;
        //tmp = floor(get_quantile_on(Xc(13,0), mu_tau_on, std_tau_on))*dt;
        
        prob_off_event  = dt * propensity_on(tmp, mu_tau_on, std_tau_on);
        prob_on_event = 0;
        
        Xc(2,0) = tmp;
        Xc(3,0) = tmp;
        
        Xc(9, 0) = prob_off_event;
        Xc(10,0) = prob_on_event;
        
        
        
    } else { // if generator 'off'
        
        
        tmp = floor(get_quantile(Xc(13,0), mu_tau_off, std_tau_off)/dt )*dt;
        //tmp = floor(get_quantile(Xc(13,0), mu_tau_off, std_tau_off))*dt;
        
        prob_off_event = 0;
        prob_on_event  = dt * propensity(tmp, mu_tau_off, std_tau_off);
        
        Xc(2,0) = tmp;
        Xc(3,0) = tmp;
        
        Xc(9,0)  = prob_off_event;
        Xc(10,0) = prob_on_event;
  
        
    }
    
    Xc(14, 0) = Xc(15,0);//1/(1+exp(-1*Xc(15,0)));
    
    double aaa = (Xc(4, 0)*f +1/(1+exp(-1*Xc(14,0)))*(1-f) ) /tau_1;//1
    Xc(5,0) = aaa*exp(-tau_1*tmp) +   (Xc(4, 0)*f + 1/(1+exp(-1*Xc(14,0)))*(1-f) )* (1 - exp(-tau_1*tmp))/tau_1;
    
    Xc(6, 0) = 5*LH_0*Xc(12,0)*exp(-alpha_1*tmp) + (B_1)*(Xc(4, 0)*f + 1/(1+exp(-1*Xc(14,0)))*(1-f) )/tau_1
    + exp(-alpha_1*tmp)*(B_1)/(tau_1/alpha_1-1)*(aaa-(Xc(4, 0)*f + 1/(1+exp(-1*Xc(14,0)))*(1-f) )/alpha_1)
    + exp(-tau_1*tmp)*(B_1)/(tau_1/alpha_1-1)*((Xc(4, 0)*f + 1/(1+exp(-1*Xc(14,0)))*(1-f) )/tau_1-aaa);
    
    //dLH/dB_1
    Xc(18,0) = ((1 - f)/(exp(-Xc(14, 0)) + 1) + Xc(4, 0)*f)/tau_1
    +(exp(-alpha_1*dt)*(aaa + ((f - 1)/(exp(-Xc(14, 0)) + 1) - Xc(4, 0)*f)/alpha_1))/(tau_1/alpha_1 - 1)
    -(exp(-dt*tau_1)*(aaa + ((f - 1)/(exp(-Xc(14, 0)) + 1) - Xc(4, 0)*f)/tau_1))/(tau_1/alpha_1 - 1);
    //dLH/dalpha_1
    Xc(19,0) =  - dt*5*LH_0*Xc(12,0)*exp(-alpha_1*dt)
    + (B_1*exp(-alpha_1*dt)*((f - 1)/(exp(-Xc(14, 0)) + 1) - Xc(4, 0)*f))/(alpha_1*(alpha_1 - tau_1))
    + (B_1*tau_1*exp(-alpha_1*dt)*(aaa + ((f - 1)/(exp(-Xc(14, 0)) + 1) - Xc(4, 0)*f)/alpha_1))/pow(alpha_1 - tau_1,2)
    + (B_1*alpha_1*dt*exp(-alpha_1*dt)*(aaa + ((f - 1)/(exp(-Xc(14, 0)) + 1) - Xc(4, 0)*f)/alpha_1))/(alpha_1 - tau_1)
    - (B_1*tau_1*exp(-dt*tau_1)*(aaa + ((f - 1)/(exp(-Xc(14, 0)) + 1) - Xc(4, 0)*f)/tau_1))/pow(alpha_1 - tau_1,2);
    //dPG/df
    Xc(20,0) =  - ((Xc(4, 0) - 1/(exp(-Xc(14, 0)) + 1))*(exp(-dt*tau_1) - 1))/tau_1;
    //dLH/df
    Xc(21,0) =  (B_1*(Xc(4, 0) - 1/(exp(-Xc(14, 0)) + 1)))/tau_1
    -(B_1*exp(-alpha_1*dt)*((Xc(4, 0) - 1/(exp(-Xc(14, 0)) + 1))/alpha_1 - Xc(20, 0)))/(tau_1/alpha_1 - 1)
    + (B_1*exp(-dt*tau_1)*((Xc(4, 0) - 1/(exp(-Xc(14, 0)) + 1))/tau_1 - Xc(20, 0)))/(tau_1/alpha_1 - 1);
    
    
    Xc(7,0) < prob_on_event   ? Xc(0,0)=1 : Xc(0,0)=0 ;
    Xc(8,0) < prob_off_event  ? Xc(1,0)=1 : Xc(1,0)=0 ;
    
    
    //loop
    for (size_t k = 0; k< timepoints-1 ; k++) { //1:1:size(Xc,1)-1) {
       // std::cout<<dim-1<<" " << Xc(21, k) << " "<<Xc(22, k)<< std::endl;
        Xc(dim-1, k+1) = Xc(dim-1, k) + dt; //time
        
        Xc(11, k+1) = Xc(11, k) ;
        Xc(12, k+1) = Xc(12, k) ;
        Xc(7, k+1) = X_ref_traj(7, k+1); //get noise form reference traj
        Xc(8, k+1) = X_ref_traj(8, k+1);
        Xc(15, k+1) = X_ref_traj(15, k+1);
        
        if (Xc(4, k) == 1) { // if generator 'on'
            
            
            if ( Xc(1, k) == 1 ) { // if turn-off event occured (t_k, t_{k+1}]
                Xc(4, k+1) = 0; //change state of generator to off
                Xc(2, k+1) = 0; //store time of event -> for calculating probability of event in [t_{k+1}, t_{k+2}]
                Xc(3, k+1) = 0; //store time of event
                
                
                //calculate prob 'on' event occurs
                prob_off_event  = 0;
                prob_on_event   = dt * propensity(Xc(3,k+1), mu_tau_off, std_tau_off);
                
                Xc(9,k+1) = prob_off_event; //reset cdf cache
                Xc(10,k+1) = prob_on_event; //reset cdf cache
                Xc(12, k+1) = 1;
                
            }
            else {
                Xc(4, k+1) = Xc(4, k);
                Xc(2, k+1) = Xc(2, k) + dt;
                Xc(3, k+1) = Xc(3, k) + dt;
                
                
                //calculate prob 'off' event occurs
                prob_off_event  = dt * propensity_on(Xc(2,k+1), mu_tau_on, std_tau_on);
                prob_on_event   = 0;
                
                
                Xc(9, k+1) = prob_off_event;
                Xc(10, k+1) = prob_on_event;
            }
            
        } else { // if generator 'off'
            
            if (Xc(0, k) == 1 ) {//if turn-on event occurent
                Xc(4, k+1) = 1;  //change state of generator to on
                Xc(2, k+1) = 0; //store time of event
                Xc(3, k+1) = 0; //store time of event
                
                
                //calculate prob 'off' event occurs
                prob_off_event  = dt * propensity_on(Xc(2,k+1), mu_tau_on, std_tau_on);
                prob_on_event   = 0;
                
                Xc(9,k+1) = prob_off_event; //reset cdf cache
                Xc(10,k+1) = prob_on_event; //reset cdf cache
                Xc(12, k+1) = 1;
                
            }
            else {
                Xc(4, k+1) = Xc(4, k);
                Xc(2, k+1) = Xc(2, k) + dt;
                Xc(3, k+1) = Xc(3, k) + dt;
                
                
                //calculate prob 'on' event occurs
                prob_off_event  = 0;
                prob_on_event   = dt * propensity(Xc(3,k+1), mu_tau_off, std_tau_off);
                
                Xc(9, k+1) = prob_off_event;
                Xc(10, k+1) = prob_on_event;
                
            }
            
        }
        
        
        Xc(7, k+1) < prob_on_event   ? Xc(0, k+1)=1 : Xc(0 ,k+1)=0 ;
        Xc(8, k+1) < prob_off_event  ? Xc(1, k+1)=1 : Xc(1, k+1)=0 ;
        
        
        Xc(5, k+1) = Xc(5,k)*exp(-tau_1*dt)   +   (Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1.0-f)) * (1 - exp(-tau_1*dt))/tau_1;
        
        Xc(6, k+1) = Xc(6,k)*exp(-alpha_1*dt) + (B_1)*(Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1-f))/tau_1
        + exp(-alpha_1*dt)*(B_1)/(tau_1/alpha_1-1)*(Xc(5, k)-(Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1-f))/alpha_1)
        + exp(-tau_1*dt)*(B_1)/(tau_1/alpha_1-1)*((Xc(4,k)*f + 1/(1+exp(-1*Xc(14, k)))*(1-f))/tau_1-Xc(5, k));
        
        Xc(14, k+1) = Xc(14, k)  + dt/2 * (-  Xc(14, k) ) +   sqrt(dt) * Xc(15, k);
    
        
        /*** grad infomration ***/
        //dLH/dB_1
        Xc(18, k+1) = Xc(18, k)*exp(-alpha_1*dt) + ((1.0 - f)/(exp(-Xc(14,k)) + 1) + Xc(4,k)*f)/tau_1
        +(exp(-alpha_1*dt)*(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/alpha_1))/(tau_1/alpha_1 - 1)
        -(exp(-dt*tau_1)*(Xc(5, k) + ((f - 1.0)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/tau_1))/(tau_1/alpha_1 - 1);
        
        
        //dLH/dalpha_1
        Xc(19, k+1) = Xc(19, k)*exp(-alpha_1*dt) - dt*Xc(6, k)*exp(-alpha_1*dt)
        + (B_1*exp(-alpha_1*dt)*((f - 1)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f))/(alpha_1*(alpha_1 - tau_1))
        + (B_1*tau_1*exp(-alpha_1*dt)*(Xc(5, k) + ((f - 1)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/alpha_1))/pow(alpha_1 - tau_1,2)
        + (B_1*alpha_1*dt*exp(-alpha_1*dt)*(Xc(5, k) + ((f - 1)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/alpha_1))/(alpha_1 - tau_1)
        - (B_1*tau_1*exp(-dt*tau_1)*(Xc(5, k) + ((f - 1)/(exp(-Xc(14, k)) + 1) - Xc(4, k)*f)/tau_1))/pow(alpha_1 - tau_1,2);
        
        //dPG/df
        Xc(20, k+1) = Xc(20, k)*exp(-tau_1*dt) - ((Xc(4, k) - 1/(exp(-Xc(14, k)) + 1))*(exp(-dt*tau_1) - 1))/tau_1;
        //dLH/df
        Xc(21, k+1) = Xc(21, k)*exp(-alpha_1*dt) + (B_1*(Xc(4, k) - 1/(exp(-Xc(14, k)) + 1)))/tau_1
        -(B_1*exp(-alpha_1*dt)*((Xc(4, k) - 1/(exp(-Xc(14, k)) + 1))/alpha_1 - Xc(20, k)))/(tau_1/alpha_1 - 1)
        + (B_1*exp(-dt*tau_1)*((Xc(4, k) - 1/(exp(-Xc(14, k)) + 1))/tau_1 - Xc(20, k)))/(tau_1/alpha_1 - 1);
        
    }
    
    return Xc; // return
}



double LH_gen_secretion_model::propensity(const double x, const double m, const double s) {
    
    double dt = params[LH_gen_secretion_model::PARAM_dt];
    size_t t = x/dt;
    double haz =0.0;
    if ( cache_haz_off.find(t) == cache_haz_off.end() ) {
        double p1 = 1/pow(s,2) ;
        double p2 = m/p1 ;
        boost::math::gamma_distribution<> my_ig(p1, p2);
        
        //double p1 = 1./m ;
        //boost::math::exponential my_ig(p1);
        
   //     try {
            haz = hazard(my_ig, x);
            cache_haz_off[t] = haz;
   //     } catch (std::overflow_error ee) {
   //         if (t>1e-3)
   //             throw ee;
            //std::cout<<x<<" "<<t<<" ---- error caught ----- "<<std::endl;
   //     }
        
        
    } else {
        
        haz = cache_haz_off.at(t);
        
        
    }
    
    return  haz;
    
    
}

double LH_gen_secretion_model::get_quantile(const double x, const double m, const double s) {
    //double dt = params[LH_gen_secretion_model::PARAM_dt];
    
    
    
    double p1 = 1/pow(s,2) ;
    double p2 = m/p1 ;
    boost::math::gamma_distribution<> my_ig(p1, p2);
    
    //double p1 = 1./m ;
    //boost::math::exponential my_ig(p1);
    
    
    double prop = 0.;
    //if (x<0.9)
    prop = quantile(my_ig, x);
    
    
    
    
    if (std::isnan(prop))
    {
        
        
        double dt = params[LH_gen_secretion_model::PARAM_dt];
        //if (m>10) {
        prop = 0.;
        double cdfv = cdf(my_ig, 0);
        while (cdfv<x) {
            prop += 0.01*dt;
            cdfv = cdf(my_ig, prop);
        }
        
        
    }
    return prop;
    
}


double LH_gen_secretion_model::propensity_on(const double x, const double m, const double s) {
    
    double dt = params[LH_gen_secretion_model::PARAM_dt];
    size_t t = x/dt;
    double haz =0.0;
    if ( cache_haz_on.find(t) == cache_haz_on.end() ) {
        double p1 = 1/pow(s,2) ;
        double p2 = m/p1 ;
        boost::math::gamma_distribution<> my_ig(p1, p2);
        
        //double p1 = 1./m ;
        //boost::math::exponential my_ig(p1);
        
        
        
   //     try {
            haz = hazard(my_ig, x);
            cache_haz_on[t] = haz;
  //      } catch (std::overflow_error ee) {
  //          if (t>1e-3)
   //             throw ee;
            //std::cout<<x<<" "<<t<<" ---- error caught ----- "<<std::endl;
   //     }
        
        
    } else {
        
        haz = cache_haz_on.at(t);
        
    }
    
    return  haz;
    
    
}

double LH_gen_secretion_model::get_quantile_on(const double x, const double m, const double s) {
    
    
    double p1 = 1/pow(s,2) ;
    double p2 = m/p1 ;
    boost::math::gamma_distribution<> my_ig(p1, p2);
    
    //double p1 = 1./m ;
    //boost::math::exponential my_ig(p1);
    
    
    double prop = 0.;
    //if (x<0.9)
    prop = quantile(my_ig, x);
    
    if (std::isnan(prop))
    {
        
        
        double dt = params[LH_gen_secretion_model::PARAM_dt];
        //if (m>10) {
        prop = 0.;
        double cdfv = cdf(my_ig, 0);
        while (cdfv<x) {
            prop += 0.01*dt;
            cdfv = cdf(my_ig, prop);
        }
        
        
    }
    return prop;
    
}
