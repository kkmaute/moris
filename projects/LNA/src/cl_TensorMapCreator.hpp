#ifndef SRC_LINALG_CL_TENSORMAPCREATOR_HPP_
#define SRC_LINALG_CL_TENSORMAPCREATOR_HPP_

// MORIS library header files.
#include "typedefs.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src

namespace moris {
    /**
     * Helper to the Tensor Class
     *
     * Creates a specialized map for implementing symmetric and
     * unsymmetric Voigt notation for 2D and 3D tensors. For more,
     * see @ref TensorClass.
     *
     */
    template< int Dim, bool Sym >
    class TensorMapCreator;

    /// @cond
    class TensorMapCreator< 2, true >;
    class TensorMapCreator< 2, false>;
    class TensorMapCreator< 3, true >;
    class TensorMapCreator< 3, false>;
    /// @endcond

}

template< int Dim, bool Sym >
class moris::TensorMapCreator
{
};

/// @cond
class moris::TensorMapCreator< 2, true>
{
public:

    /**
     * This is the first map.
     */
    static
    moris::Mat< moris::uint >
    makeMap()
    {
        moris::Mat< moris::uint > temp( { { 0, 2}, {2, 1} } );
        return temp;
    }
};


class moris::TensorMapCreator< 2, false>
{
public:

    /**
     * This is the second map.
     */
    static
    moris::Mat< moris::uint >
    makeMap()
    {
        moris::Mat< moris::uint > temp( { { 0, 3}, {2, 1} } );
        return temp;
    }
};

class moris::TensorMapCreator< 3, true>
{
public:

    /**
     * This is the third map.
     */
    static
    moris::Mat< moris::uint >
    makeMap()
    {
        moris::Mat< moris::uint > temp( { { 0, 5, 4}, {5, 1, 3}, {4, 3, 2} } );
        return temp;
    }
};


class moris::TensorMapCreator< 3, false>
{
public:

    /**
     * This is the fourth map.
     */
    static
    moris::Mat< moris::uint >
    makeMap()
    {
        moris::Mat< moris::uint > temp( { { 0, 5, 4}, {8, 1, 3}, {7, 6, 2} } );
        return temp;
    }
};
/// @endcond

#endif /* SRC_LINALG_CL_TENSORMAPCREATOR_HPP_ */
