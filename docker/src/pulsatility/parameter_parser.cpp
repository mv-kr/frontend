//
//  parameter_parser.cpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 29/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#include "parameter_parser.hpp"
#include <fstream>
#include <sstream> 
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <exception>
using namespace mci::utils;


parameter_parser::parameter_parser() {}

/** helper function that reads line into a vector */
std::vector<double> mci::utils::parameter_parser::parse_line_to_vector( std::stringstream& line_stream, char delimiter) {
    std::vector<double> result;
    std::string item;
    while (std::getline(line_stream, item, delimiter)) {

        result.push_back(std::stod(item));
            
        
    }
    return result;
}

/** read data from str_file **/
void parameter_parser::read_data(const std::string& str_file, const size_t rows_per_entry) {
    
    // aux variables
    std::string line, val_str;
    std::stringstream row_stream;
    
    
    //open file
    std::ifstream in_file;
    in_file.open(str_file);

    if(in_file) {
        //read data line by line and store in a map of eigen::matrixXd objects
        std::vector<std::vector<double>> m(0);
        size_t max_dim = 0;

        std::string key;
        

        while ( std::getline(in_file, line ) ) {
            
            bool dataset_id_present = true;
            
            row_stream.clear();
            row_stream.str(line);
            
            std::string raw_first;
            std::getline(row_stream, raw_first, ',' ) ; // get dataset_id

            // Try to convert to number
            try {
                std::stod(raw_first);
                // If successful, this is numeric data without ID
                key = "dataset_" + std::to_string(data_map.size() + 1);
                // Reset stream to include first value in data
                row_stream.clear();
                row_stream.str(line);
                dataset_id_present = false;
            } catch (std::exception&) {
                // If conversion fails, treat as string ID
                key = raw_first;
            }    


            std::vector<double> v_tmp(0);
            try{
                //read data
                v_tmp = parse_line_to_vector(row_stream, ',' );
            } catch (std::exception& e) {
                std::cout<<row_stream.str()<<std::endl;
                throw std::invalid_argument("Error parsing data for dataset " + key + ": Numeric values should be separated by ',' .");
            }
            

            m.push_back(v_tmp);

            for (int r=1; r<rows_per_entry; r++) {
                // read next row of the same dataset
                std::getline(in_file, line );
                row_stream.clear();
                row_stream.str(line);
                //skip dataset_id
                if (dataset_id_present)
                    std::getline(row_stream, key, ',' );
                //read data
                std::vector<double> v_tmp2(0);
                try {
                    v_tmp2 = parse_line_to_vector(row_stream, ',' );    
                } catch (std::exception& e) {
                    std::cout<<row_stream.str()<<std::endl;
                    throw std::invalid_argument("Error parsing data for dataset " + key + ": Numeric values should be separated by ',' .");
                }
                
                //append to matrix
                m.push_back(v_tmp2);
            }

            //convert to eigen::matrixXd
            Eigen::MatrixXd eig_mat(m.size(),m[0].size());
            for (size_t i =0 ; i<m.size(); i++ ) {
                for (size_t j =0 ; j<m[i].size(); j++ )
                    eig_mat(i, j) = m[i][j];
            }
            

            data_map[key] = eig_mat;
            m.clear();
            
        }
    in_file.close();
 //std::cout<<eig_mat<<std::endl;
    
    } else {
        throw std::invalid_argument("Couldn't find data file: " + str_file );
    }
}






/** read parameters from a file ***/
void parameter_parser::read_parameters(const std::string& str_file) {
    
    //aux variables
    std::string line;
    std::string alg_name, param_name, type, val_str, val_row;
    std::stringstream line_stream;
    std::stringstream row_stream;
    
    //open file
    std::ifstream in_file;
    in_file.open(str_file);
    if (in_file) {
        //read lines
        while (std::getline(in_file, line)) {
            
            //conver line to stringstream object
            
            line_stream.clear();
            line_stream.str(line);
            
            //read in the following order : 1) name of alogrithm 2) name of parameter and 3) type
            line_stream >> alg_name >> param_name >> type;
            if (!line_stream.eof()) {
                
                //if type is scalar
                if (type.compare("Sd") == 0) {
                    //read value
                    line_stream >> val_str;
                    
                    try {
                        //store in the corresponding map
                        resD[alg_name][param_name] = std::stod(val_str);
                    } catch (std::exception& e) {
                        throw std::invalid_argument("Error parsing parameter " + param_name + ": Numeric value expected.");
                    }
                    
                    
                } else if (type.compare("Vd") == 0) { //if vector
                    
                    std::vector<double> v_tmp(0);
                    try {
                        //read values
                        while ( std::getline(line_stream, val_str, ',' )   )
                            v_tmp.push_back(std::stod(val_str));
                    } catch (std::exception& e) {
                        throw std::invalid_argument("Error parsing vector " + param_name + ": Numeric values should be separated by ',' expected.");
                    }
                        
                    //convert vector to Eigen::VectorXd
                    Eigen::VectorXd eig_v(v_tmp.size());
                    for (size_t i =0 ; i<v_tmp.size(); i++ )
                        eig_v(i) = v_tmp[i];
                    
                    //store in resV map
                    resV[param_name] = eig_v;
                        
                    
                    
                    
                    
                }  else if (type.compare("Ss") == 0) { //if  string
                    
                    //read parameter
                    line_stream >> val_str;
                    //store in resS_s map
                    resS_s[param_name] = val_str;
                    
                }
            
                else if (type.compare("Vs") == 0) { //if vector of strings
                    
                    
                    //read vector of values
                    std::vector<std::string> v_tmp(0);
                    while ( std::getline(line_stream, val_str, ',' )   ) {
                        boost::trim ( val_str );
                        v_tmp.push_back(val_str);
                    }
                    
                    
                    //store in resV_s map
                    resV_s[param_name] = v_tmp;
                    
                    
                    
                } else if (type.compare("Md") == 0) { //if matrix
                 
                    
                    
                    //read matrix data
                    std::vector<std::vector<double>> m(0);
                    try{
                        while (std::getline(line_stream, val_row, ';' ) ) {
                            
                            row_stream.clear();
                            row_stream.str(val_row);
                            
                            //read vector of values separated by ','
                            std::vector<double> v_tmp(0);
                            while (std::getline(row_stream, val_str, ',' ) ) {
                                v_tmp.push_back(std::stod(val_str));
                            }
                            
                            m.push_back(v_tmp);
                            
                        }
                    }  catch (std::exception& e) {
                        throw std::invalid_argument("Error parsing Matrix " + param_name + ": Values in each row should be separated by ',' and rows should be separated by ';'.");
                    }
                    
                    //convert to eigen::matrixXd
                    Eigen::MatrixXd eig_mat(m.size(),m[0].size());
                    
                    for (size_t i =0 ; i<m.size(); i++ ) {
                        if (m[i].size() != m[0].size())
                            throw std::invalid_argument("Error parsing Matrix " + param_name + ": Check all rows have the same number of elements.");
                        for (size_t j =0 ; j<m[i].size(); j++ )
                            eig_mat(i, j) = m[i][j];
                    }
                    
                        
                    
                    //store in resM map
                    resM[param_name] = eig_mat;
                        
                    
                }
            
            }
        }
        in_file.close();
    } else {
        throw std::invalid_argument("Couldn't find parameter file: " + str_file);
    }
}

