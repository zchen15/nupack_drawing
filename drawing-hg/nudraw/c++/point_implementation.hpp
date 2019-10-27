#ifndef POINT_IMPLEMENTATION_HPP
#define POINT_IMPLEMENTATION_HPP

/*nudraw::point::point<numeric>(  const numeric & initial_x,
                                const numeric & initial_y,
                                const numeric & initial_z)
: x(initial_x), y(initial_y), z(initial_z)
{
    // no op
};*/

template <class numeric>
nudraw::point<numeric>::point()
: x(0), y(0), z(0)
{
    // no op
};

template <class numeric>
const numeric & nudraw::point<numeric>::get_x() const
{
    return x;
}

template <class numeric>
const numeric & nudraw::point<numeric>::get_y() const
{
    return y;
}

template <class numeric>
const numeric & nudraw::point<numeric>::get_z() const
{
    return z;
}

#endif        //  #ifndef POINT_IMPLEMENTATION_HPP

