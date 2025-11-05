//
//  parameter_parser.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 29/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef parameter_parser_hpp
#define parameter_parser_hpp

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

namespace mci {
    
    namespace utils {
        
        // A parameter_set is a map containing a string description of the parameter (key) and its value
        typedef std::map<std::string, double> parameter_set;
        
        class parameter_parser {
            
        
            
        public:
            
            std::map<std::string, mci::utils::parameter_set> resD;
            std::map<std::string, Eigen::VectorXd> resV;
            std::map<std::string, Eigen::MatrixXd> resM;
            std::map<std::string, std::vector<std::string>> resV_s;
            std::map<std::string, std::string> resS_s;
            Eigen::MatrixXd data;
            std::map<std::string, Eigen::MatrixXd> data_map;

            parameter_parser();
            
            /** read parameters from a file ***/
            void read_parameters(const std::string& str_file);
            
            /** read parameters from a file ***/
            void read_data(const std::string& str_file, const size_t rows_per_entry);

            std::vector<double> parse_line_to_vector( std::stringstream& line_stream, char delimiter);
            
            
        };
        
    }
    
    
}

#endif /* parameter_parser_hpp */
