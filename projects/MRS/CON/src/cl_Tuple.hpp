/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Tuple.hpp
 *
 */

#ifndef MORIS_CONTAINERS_CL_TUPLE_HPP_
#define MORIS_CONTAINERS_CL_TUPLE_HPP_

// C++ header files.
#include <tuple>

namespace moris
{
    /**
     * @brief Tuple wrapper class
     *
     * moris::Tuple is a template wrapper around std::tuple.
     *
     * Think of a tuple as a "vector" of "things". This is in contrast
     * with moris::Mat, for example, whose elements are always the same type.
     * When declaring a tuple, you must know the types of each argument. For
     * example, here is a tuple with an integer, a moris::real, and a string:
     *
     * @include CON/src/cl_tuple.inc
     */
    template< typename... Types >
    class Tuple
    {
    public:

        /**
         *
         * @brief default constructor
         *
         * tuple default constructor
         */
        Tuple() = default;

        /**
         * tuple initialization constructor
         *
         * @param[in] aArgs list of input arguments
         */
        explicit
        Tuple(
                Types const &... aArgs )
            : mTuple( aArgs... )
        {
        }

        /**
         * Wraps a standard tuple into a moris::Tuple
         *
         * @param[in] stdTuple Standard tuple to be wrapped
         */
        Tuple(
                std::tuple<Types...> stdTuple )
            : mTuple(stdTuple)
        {
        }

        /**
         * tuple copy constructor
         *
         * @param[in] aTuple Given tuple to be copied
         */
        template< class ... UTypes >
        Tuple(
                moris::Tuple< UTypes... > const & aTuple )
            : mTuple( aTuple.data() )
        {
        }

        /**
         * tuple default destructor
         */
        ~Tuple() = default;

        template<typename... UTypes>
        void
        operator=(
                moris::Tuple<UTypes...> aTuple )
        {
            mTuple = aTuple.data();
        }

        /**
         * tuple non-constant accessor
         *
         * @return const reference to underlying standard tuple
         */
        const std::tuple< Types... > &
        data()
        {
            return mTuple;
        }

        /**
         * tuple constant accessor
         *
         * @return const reference to underlying standard tuple
         */
        const std::tuple< Types... > &
        data() const
        {
            return mTuple;
        }

        /**
         * Gets a component of the tuple
         *
         * @note N is not given as a parameter, but as a type for the get
         * function.For example, to get the 2nd component, use myTuple.get<1>()
         * instead of myTuple.get(1)!
         */
        template< moris::size_t N >
        typename std::tuple_element<N, std::tuple< Types...>>::type&
        get()
        {
            return std::get< N >( mTuple );
        }

        /**
         *  Compares every element of the tuple lhs with
         *  the corresponding element of the tuple rhs.
         *  Note that both Tuples must be the same size.
         */
        template< class ... UTypes >
        bool
        operator==(
                moris::Tuple< UTypes... > const & aTuple )
        {
            return mTuple == aTuple.data();
        }

    private:

        std::tuple< Types... > mTuple;
    };
}

namespace moris
{
    /**
     * moris::tie is used to unpack a Tuple returned by some functions
     *
     * @param[in] aArgs List of arguments to be unpacked
     *
     * Using moris::tie requires knowledge of the order in which returns are
     * organized. For example, moris::qr accepts a matrix as an argument and
     * returns a Tuple of two matrices, Q and R, in that order. Thus,
     * moris::tie is used as follows: moris::tie(Q,R) = moris::qr(A). Q and R
     * must be initialized before using moris::tie.
     *
     * Additionally, moris::tie can be used to ignore a particular component
     * of the output Tuple. For example, moris::tie(Q,std::ignore) = moris::qr(A)
     * would ignore the R, which does not need to be initialized.
     */
    template< typename... Types>
    moris::Tuple<Types &... >
    tie(Types &... aArgs)
    {
        return moris::Tuple<Types &... >(aArgs...);
    }
}

namespace moris
{
    /**
     * Makes a moris::Tuple from a list of arguments.
     *
     * @param[in] aArgs List of arguments to be packed into a Tuple.
     */
    template< typename... Types>
    moris::Tuple<Types...>
    make_tuple(Types &... aArgs)
    {
        return moris::Tuple<Types...>(aArgs...);
    }
}

namespace moris
{
    /**
     * Renaming of std::ignore to moris::tilde.
     */
//    decltype(std::ignore)& tilde = std::ignore;
    decltype(std::ignore) tilde;
}

#endif /* CL_TUPLE_HPP_ */

