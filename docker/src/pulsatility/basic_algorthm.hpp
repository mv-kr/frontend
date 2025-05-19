//
//  basic_algorthm.hpp
//  pmcmc
//
//  Created by Margaritis Voliotis on 23/11/2018.
//  Copyright © 2018 Margaritis Voliotis. All rights reserved.
//

#ifndef basic_algorthm_hpp
#define basic_algorthm_hpp

#include <iostream>

namespace mci {
    namespace algorithms {
        
        template <typename P, typename R>
        class basic_algorithm {
        
        protected:
            /** parameters **/
            P parameters;
            /** defualt parameters **/
            static P default_parameters;
            
        public:
            
            virtual ~basic_algorithm() {};
            
            /** run the algorithm and store results in res **/
            virtual const void run(R& res) = 0;
            
            /** run the algorithm using parameters in params and store results in res **/
            virtual const void run(const P& p, R& res) = 0;
            
            /** get default parameters **/
            const P get_default_parameters() {
                return basic_algorithm::default_parameters;
            }
            
            /** get parameters **/
            const P get_parameters() {
                return parameters;
            }
            
            
            /** get default parameters **/
            void set_parameters(const P& p) {
                parameters = p;
            }
            
            
        };
    }
    
}

#endif /* basic_algorthm_hpp */
