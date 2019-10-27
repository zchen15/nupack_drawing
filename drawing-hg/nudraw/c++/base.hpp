#ifndef BASE_HPP
#define BASE_HPP

#include "globals.hpp"
#include "point.hpp"

#define POINT_NUMERIC_TYPE float

namespace nudraw
{
    class base
    {
        public:
        
        // possible base types
        
        explicit base(constants::base_type initial_type, int initial_index);
        
        // getters
        const constants::base_type & get_type()                     const;
        const int & get_index()                                     const;
        const nudraw::point<POINT_NUMERIC_TYPE> & get_position()    const;
        
        // setters
        void set_base_type(const constants::base_type  & new_type);
        void set_index(const int & new_index);
        
        private:
            
        // data members
        constants::base_type type;
        int index;
        
        point<POINT_NUMERIC_TYPE> position;
        
        bool bonded_left;
        bool bonded_right;
    };
}

#endif        //  #ifndef BASE_HPP

