/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_TOL_Capacities.hpp
 *
 */

#ifndef SRC_TOOLS_FN_TOL_CAPACITIES_HPP_
#define SRC_TOOLS_FN_TOL_CAPACITIES_HPP_

#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
namespace moris
{
    /*!
     * @brief Calculates the internal data structure capacity of the moris::cell
     * the internal class must have a capacity function defined.
     */
    template<typename Class_With_Capacity_FN>
    inline
    size_t
    internal_capacity(moris::Cell<Class_With_Capacity_FN> & aCell)
    {
        size_t tInternalCapacity = 0;

        for(moris::uint i = 0; i < aCell.size(); i++)
        {
            tInternalCapacity = aCell(i).capacity();
        }

        // return the calculated internal memory usage,
        // add the capacity ofthe vector itself
        return tInternalCapacity + aCell.capacity();
    }

     /*!
     * @brief Calculates the internal data structure capacity of the moris::cell
     * the internal class must have a capacity function defined.
     */

    template<>
    inline
    size_t
    internal_capacity(moris::Cell<std::string> & aCell)
    {
        size_t tInternalCapacity = 0;

        for(moris::uint i = 0; i < aCell.size(); i++)
        {
            tInternalCapacity = aCell(i).length()*sizeof(char);
        }

        // return the calculated internal memory usage,
        // add the capacity ofthe vector itself
        return tInternalCapacity + aCell.capacity();
    }

    template<typename Class_With_Capacity_FN>
    inline
    size_t
    internal_capacity_nested(moris::Cell<moris::Cell<Class_With_Capacity_FN>> & aCell)
    {
        size_t tInternalCapacity = 0;

        for (moris::uint i = 0; i < aCell.size(); i++)
        {
            for (moris::uint j = 0; j < aCell(i).size(); j++)
            {
                tInternalCapacity = aCell(i)(j).capacity();
            }
        }

        // return the calculated internal memory usage,
        // add the capacity ofthe vector itself
        return tInternalCapacity + aCell.capacity();
    }

    //------------------------------------------------------------------------------

    /**
     * @brief compute the cpacity of moris::Cell<mori::Cell<DataTyp>>
     * where DataType is a native data type that does not have a capacity
     *
     * @tparam DataType
     * @param aCell
     * @return size_t
     */
    template< typename DataType >
    inline size_t
    internal_capacity_nested_datatype( moris::Cell< moris::Cell< DataType > > &aCell )
    {
        size_t tInternalCapacity = 0;

        for ( const auto &iInnerCells : aCell )
        {
            tInternalCapacity += iInnerCells.capacity() * ( sizeof( DataType ) + 1 );
        }

        // return the calculated internal memory usage,
        // add the capacity ofthe vector itself
        return tInternalCapacity + aCell.capacity();
    }

     //------------------------------------------------------------------------------
    template <typename Class_PTR_With_Capacity_FN>
    inline
    size_t
    internal_capacity_ptr(moris::Cell<Class_PTR_With_Capacity_FN> &aCell)
    {
        size_t tInternalCapacity = 0;

        for (moris::uint i = 0; i < aCell.size(); i++)
        {
            tInternalCapacity = aCell(i)->capacity();
        }

        // return the calculated internal memory usage,
        // add the capacity ofthe vector itself
        return tInternalCapacity + aCell.capacity();
    }

    template <typename Class_PTR_With_Capacity_FN>
    inline
    size_t
    internal_capacity_nested_ptr(moris::Cell<moris::Cell<Class_PTR_With_Capacity_FN>> &aCell)
    {
        size_t tInternalCapacity = 0;

        for (moris::uint i = 0; i < aCell.size(); i++)
        {
            for (moris::uint j = 0; j < aCell(i).size(); j++)
            {
                tInternalCapacity = aCell(i)(j)->capacity();
            }
        }

        // return the calculated internal memory usage,
        // add the capacity ofthe vector itself
        return tInternalCapacity + aCell.capacity();
    }
}

#endif
