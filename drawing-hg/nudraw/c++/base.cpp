#include "base.hpp"

nudraw::base::base(constants::base_type initial_type, int initial_index)
: type (initial_type), index (initial_index)
{
    // nothing to do for now
}

const nudraw::constants::base_type & nudraw::base::get_type() const
{
    return type;
}

const int & nudraw::base::get_index() const
{
    return index;
}

const nudraw::point<POINT_NUMERIC_TYPE> & nudraw::base::get_position() const
{
    return position;
}
