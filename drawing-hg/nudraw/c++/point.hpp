#ifndef POINT_HPP
#define POINT_HPP

namespace nudraw
{
    template <class numeric>
    class point
    {
        public:
        
        //explicit point(const numeric & x, const numeric & y, const numeric & z);
        explicit point();
        
        // getters
        const numeric & get_x() const;
        const numeric & get_y() const;
        const numeric & get_z() const;
        
        private:
        
        numeric x;
        numeric y;
        numeric z;
    };
}

#include "point_implementation.hpp"

#endif        //  #ifndef POINT_HPP
