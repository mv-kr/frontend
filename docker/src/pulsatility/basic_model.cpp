//
//  basic_model.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 23/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#include "basic_model.hpp"



/** calculate the log ratio of proposal densities: p(theta_proposed | theta_current)/ p(theta_current  | theta_proposed) **/
const double mci::models::basic_model::log_proposal_density_ratio(const Eigen::VectorXd& theta_proposed, const Eigen::VectorXd& theta_current) {
    return log_proposal_density( theta_proposed, theta_current) -
            log_proposal_density(  theta_current, theta_proposed);
}


