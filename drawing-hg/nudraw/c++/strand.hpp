#ifndef STRAND_HPP
#define STRAND_HPP

#include <vector>

#include "base.hpp"

using namespace std;

namespace nudraw
{
    
    class strand
    {
        public:
        
        explicit strand();
            
        private:
         
        vector<base *> contents;
    };
    
}

#endif        //  #ifndef STRAND_HPP

