/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Vector.hpp
 *
 */
#ifndef MORIS_CONTAINERS_CL_VECTOR_HPP_
#define MORIS_CONTAINERS_CL_VECTOR_HPP_

// C++ header files.
#include <vector>
#include <algorithm>    // for unique
#include <iostream>
#include <iterator>

// MORIS library header files.
#include "moris_typedefs.hpp"    // COR/src
#include "assert.hpp"

namespace moris
{
    //------------------------------------------------------------------
    template< typename T >
    class Vector
    {
      private:
        /**
         * MORIS vector
         */
        std::vector< T > mVector;

#ifdef CHECK_MEMORY
        uint mNumResizeCalls     = 0;
        uint mNumReserveCalls    = 0;
        uint mNumImplicitResizes = 0;
        bool mPushBackMemWarn    = true;
#endif

      public:
        using value_type = T;

        /**
         * Vector constructor
         */
        Vector() = default;

        //------------------------------------------------------------------
        /**
         * Vector constructor
         *
         * @param[in] aVector A Vector
         */

        Vector(
                std::initializer_list< T > const & aVector )
                : mVector( aVector.begin(), aVector.end() )
        {
        }

        //------------------------------------------------------------------
        /**
         *
         * Vector copy constructor
         */
        Vector( const Vector< T >& aVector )
        {
            mVector = aVector.data();
        }

        //------------------------------------------------------------------

        /**
         * Vector constructor
         *
         * @param[in] aSize Size of the Vector to be initialized
         * @param[in] aValue Value of the elements of the initialized Vector
         */

        template< typename A >
        Vector(
                moris::uint const aSize,
                A const           aValue )
                : mVector( aSize, aValue )
        {
            MORIS_CHECK_MEMORY( sizeof( T ) * aSize < MORIS_MAX_CELL_CAPACITY,
                    "Vector: Maximum allowable capacity exceeded: %f MB.\n",
                    sizeof( T ) * aSize / 1e6 );
        }

        //------------------------------------------------------------------

        Vector(
                moris::uint const aSize )
                : mVector( aSize )
        {
            MORIS_CHECK_MEMORY( sizeof( T ) * aSize < MORIS_MAX_CELL_CAPACITY,
                    "Vector: Maximum allowable capacity exceeded: %f MB.\n",
                    sizeof( T ) * aSize / 1e6 );
        }

        //------------------------------------------------------------------

        /**
         * Vector destructor
         */
#ifdef CHECK_MEMORY
        ~Vector()
        {
            MORIS_CHECK_MEMORY(
                    sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? this->size() > MORIS_CELL_UTILIZATION_FRACTION_LIMIT * this->capacity() : true,
                    "Vector::~Vector:: At destruction memory used is less than 75 percent of capacity of large matrix: size %zu capacity %zu\n",
                    this->size(),
                    this->capacity() );
        }
#else
        ~Vector() = default;    // 'default' tells the compiler to automatically
                              // delete the underlying Vector
#endif

        //------------------------------------------------------------------

        /**
         * @brief Operator to replace the contents of the container
         *
         * @param[in] val Contents to be copied into the container
         */

        const Vector< T >&
        operator=( Vector< T > const & aVector )
        {
            mVector = aVector.data();

            return *this;
        }

        //------------------------------------------------------------------

        /**
         * @brief Operator to replace the contents of the container
         *
         * @param[in] val Replacement Vector
         *
         * @return Assignment.
         */

        template< typename A >
        const Vector< T >&
        operator=(
                const A& val )
        {
            mVector = val;

            return *this;
        }

        //------------------------------------------------------------------

        /**
         * Replaces the contents of the container
         *
         * @param[in] aCount number of entries
         * @param[in] aval Value to set as contents of container
         */

        void
        assign(
                moris::size_t const aCount,
                T const &           aval )
        {
            mVector.assign( aCount, aval );
        }

        //------------------------------------------------------------------

        /**
         * @brief Operator for accessing a specified Vector element.
         *
         * @param[in] i_index Position of an element in the Vector.The function
         *            at() automatically checks whether an element is within the
         *            bounds of valid elements in the container, throwing an
         *            out_of_range exception if it is not. This is in contrast
         *            with member operator[], that does not check against bounds
         *            and hence is faster.
         */

        auto
        operator()(
                moris::size_t const i_index )
                -> decltype(
#ifdef MORIS_HAVE_DEBUG
                        ( mVector.at( i_index ) )
#else
                        ( mVector[ i_index ] )
#endif
                )
        {
#ifdef MORIS_HAVE_DEBUG
            return ( mVector.at( i_index ) );
#else
            return ( mVector[ i_index ] );
#endif
        }

        //------------------------------------------------------------------

        /**
         * const version of above
         */

        auto
        operator()(
                moris::size_t const i_index ) const
                -> decltype(
#ifdef MORIS_HAVE_DEBUG
                        ( mVector.at( i_index ) )
#else
                        ( mVector[ i_index ] )
#endif
                )
        {
#ifdef MORIS_HAVE_DEBUG
            return ( mVector.at( i_index ) );
#else
            return ( mVector[ i_index ] );
#endif
        }

        //------------------------------------------------------------------

        /**
         * Direct access to the underlying array.
         *
         * @return Direct access to the underlying array.
         */

        std::vector< T > const &
        data() const
        {
            return mVector;
        }

        std::vector< T >&
        data()
        {
            return mVector;
        }

        //------------------------------------------------------------------

        /**
         * Direct access to the memory of the underlying array.
         *
         * @return Pointer to memory of the underlying array.
         */

        T const *
        memptr() const
        {
            return mVector.data();
        }

        T*
        memptr()
        {
            return mVector.data();
        }

        //------------------------------------------------------------------

        /**
         * Checks whether the Vector container is empty.
         *
         * @return bool value indicating whether the Vector container is empty.
         */

        bool
        empty() const
        {
            return mVector.empty();
        }

        //------------------------------------------------------------------

        /**
         * Returns the last element of the Vector.
         *
         * @return the last element of the Vector.
         */

        auto
        back()
                -> decltype( mVector.back() )
        {
            return mVector.back();
        }

        //------------------------------------------------------------------

        /**
         * Returns the number of elements in the Vector.
         *
         * @return Number of elements in the Vector.
         */

        auto
        size()
                -> decltype( mVector.size() )
        {
            return mVector.size();
        }

        //------------------------------------------------------------------

        /**
         * const version of above
         */

        auto
        size() const
                -> decltype( mVector.size() )
        {
            return mVector.size();
        }

        //------------------------------------------------------------------

        /**
         * Returns the maximum possible number of elements in the Vector.
         *
         * @return Maximum possible number of elements in the Vector.
         */

        auto
        max_size()
                -> decltype( mVector.max_size() )
        {
            return mVector.max_size();
        }

        //------------------------------------------------------------------

        /**
         * @brief Increase the capacity of the container to a value that's
         *        greater or equal to new_cap
         *
         * @param[in] new_cap The new capacity of the Vector
         */

        void
        reserve( moris::size_t const & new_cap )
        {
            MORIS_CHECK_MEMORY( sizeof( T ) * new_cap < MORIS_MAX_CELL_CAPACITY,
                    "Vector::reserve: Maximum allowable capacity exceeded: %f MB.\n",
                    sizeof( T ) * new_cap / 1e6 );

            MORIS_CHECK_MEMORY( mNumReserveCalls++ < MORIS_CELL_RESERVE_CALL_LIMIT,
                    "Vector::reserve: number of reserve calls exceeds limit.\n" );

            mVector.reserve( new_cap );
        }

        //------------------------------------------------------------------

        /**
         * @brief Returns the number elements that can be held in currently
         *        allocated storage.
         *
         * @return Number elements that can be held in currently allocated storage.
         */

        auto
        capacity() const
                -> decltype( mVector.capacity() )
        {
            return mVector.capacity();
        }

        //------------------------------------------------------------------

        /**
         * @brief Requests the removal of unused capacity
         */
        void
        shrink_to_fit()
        {
            // check that resize on large cells does not indicate overly conservative initial size allocation
            MORIS_CHECK_MEMORY(
                    sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? this->size() > MORIS_CELL_RESIZE_FRACTION_LIMIT * this->capacity() : true,
                    "Vector::shrink_to_fit: Shrink to less than 10 percent of capacity of large matrix: size %zu capacity %zu.\n",
                    this->size(),
                    this->capacity() );

            // check that number of resize + shrink_to_fit calls does not exceed limit
            MORIS_CHECK_MEMORY( mNumResizeCalls++ != MORIS_CELL_RESIZE_CALL_LIMIT,
                    "Vector::shrink_to_fit: number of resize + shrink_to_fit calls exceeds limit.\n" );

            mVector.shrink_to_fit();
        }

        //------------------------------------------------------------------

        /**
         * @brief Requests the removal of unused capacity of outer and inner cells
         */
        void
        shrink_to_fit_all()
        {
            // unless template specialization holds call shrink_to_fit on outer cell
            shrink_to_fit();
        }

        //------------------------------------------------------------------

        /**
         * Clears the contents of the Vector.
         */
        void
        clear()
        {
            mVector.clear();
        }

        //------------------------------------------------------------------

        /**
         * @brief Inserts elements at the specified location in the container
         *
         * @param[in] pos Iterator at which the content will be inserted.
         * @param[in] value Element value to insert
         */

        void
        insert(
                moris::size_t const & pos,
                T const &             value )
        {
            MORIS_CHECK_MEMORY( pos >= mVector.capacity() ? mNumResizeCalls++ != MORIS_CELL_RESIZE_CALL_LIMIT : true,
                    "Vector::insert: number of resize calls exceeds limit.\n" );

            mVector.insert( mVector.begin() + pos, value );

            MORIS_CHECK_MEMORY(
                    sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? this->size() > MORIS_CELL_UTILIZATION_FRACTION_LIMIT * this->capacity() : true,
                    "Vector::insert: After insert memory used is less than 75 percent of capacity of large matrix: size %zu capacity %zu\n",
                    this->size(),
                    this->capacity() );
        }

        //------------------------------------------------------------------

        /**
         * Appends a Vector, similar to how push_back appends a single value.
         *
         * @param[in] aVector Vector to be appended
         */

        Vector&
        append(
                Vector< T > const & aVector )
        {
            MORIS_CHECK_MEMORY( mVector.size() + aVector.size() > mVector.capacity() ? mNumResizeCalls++ != MORIS_CELL_RESIZE_CALL_LIMIT : true,
                    "Vector::append: number of resize calls exceeds limit.\n" );

            mVector.insert( mVector.end(), aVector.data().begin(), aVector.data().end() );

            MORIS_CHECK_MEMORY(
                    sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? this->size() > MORIS_CELL_UTILIZATION_FRACTION_LIMIT * this->capacity() : true,
                    "Vector::append: After append memory used is less than 75 percent of capacity of large matrix: size %zu capacity %zu\n",
                    this->size(),
                    this->capacity() );

            return *this;
        }

        //------------------------------------------------------------------

        /**
         * @brief Removes specified element from the container
         *
         * @param[in] pos Position (not index) of element to be remove
         */

        void
        erase(
                moris::size_t const & pos )
        {
            mVector.erase( mVector.begin() + pos );
        }

        //------------------------------------------------------------------

        /**
         * @brief Removes specified element from the container lying in the range
         *        [ pos1, pos2 ).
         *
         * @param[in] pos1 Position(not index) of first element to be removed
         * @param[in] pos2 Position(not index) of last element to be removed + 1
         *
         */

        void
        erase_range(
                moris::size_t const & pos1,
                moris::size_t const & pos2 )
        {
            mVector.erase( mVector.begin() + pos1, mVector.begin() + pos2 );
        }

        //------------------------------------------------------------------

        void
        remove( T const & value )
        {
            mVector.erase( std::remove( mVector.begin(), mVector.end(), value ), mVector.end() );
        }

        //------------------------------------------------------------------

        /**
         * @brief Appends the given element value to the end of the container
         *
         * @param[in] value The value of the element to append
         */

        void
        push_back(
                T const & value )
        {
#ifdef CHECK_MEMORY
            const moris::size_t tOldCapacity = mVector.capacity();

            mVector.push_back( value );

            MORIS_CHECK_MEMORY( mVector.capacity() != tOldCapacity && sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? mNumImplicitResizes++ != MORIS_CELL_IMPLICIT_RESIZE_CALL_LIMIT : true,
                    "Vector::push_back: number of implicit resize calls exceeds limit.\n" );

            if ( sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY && this->size() < MORIS_CELL_UTILIZATION_FRACTION_LIMIT * this->capacity() && mPushBackMemWarn )
            {
                MORIS_CHECK_MEMORY( false,
                        "Vector::push_back: After push_back memory used is less than 75 percent of capacity of large matrix: size %zu capacity %zu\n",
                        this->size(),
                        this->capacity() );

                // flag to false to avoid repeated issue of warning
                mPushBackMemWarn = false;
            }
#else
            mVector.push_back( value );
#endif
        }

        //------------------------------------------------------------------

        /**
         * @brief Appends the given element value to the end of the container
         *
         * @param[in] value    The value of the element to append
         * @param[in] tSizeInc The size by which the cell capacity is increased when not sufficient
         *
         */

        void
        push_back(
                T const &             value,
                moris::size_t const & tSizeInc )
        {
            if ( mVector.size() == mVector.capacity() )
            {
                MORIS_CHECK_MEMORY( mNumResizeCalls++ != MORIS_CELL_RESIZE_CALL_LIMIT,
                        "Vector::push_back: number of resize calls exceeds limit.\n" );

                MORIS_CHECK_MEMORY( mNumReserveCalls++ < MORIS_CELL_RESERVE_CALL_LIMIT,
                        "Vector::reserve: number of reserve calls exceeds limit.\n" );

                mVector.reserve( mVector.size() + tSizeInc );
            }

            mVector.push_back( value );
        }

        //------------------------------------------------------------------

        /**
         * @brief removes the last element of the container
         */

        void
        pop_back()
        {
            mVector.pop_back();
        }

        //------------------------------------------------------------------

        /**
         * @brief Resizes the container to contain count elements
         *
         * @param[in] aCount The new size of the Vector
         * @param[in] aValue The value to initialize the new elements with
         */

        void
        resize(
                moris::size_t const & aCount,
                T                     aValue = T() )
        {
            MORIS_CHECK_MEMORY( sizeof( T ) * aCount < MORIS_MAX_CELL_CAPACITY,
                    "Vector::resize: Maximum allowable capacity exceeded: %f MB.\n",
                    sizeof( T ) * aCount / 1e6 );

            // check that resize on large cells does not indicate overly conservative initial size allocation
            MORIS_CHECK_MEMORY(
                    sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? aCount > MORIS_CELL_RESIZE_FRACTION_LIMIT * this->capacity() : true,
                    "Vector::resize: Resize to less than 10 percent of capacity of large matrix: size %zu capacity %zu\n",
                    this->size(),
                    this->capacity() );

            // check that number of resize + shrink_to_fit calls does not exceed limit
            MORIS_CHECK_MEMORY( aCount != this->size() ? mNumResizeCalls++ != MORIS_CELL_RESIZE_CALL_LIMIT : true,
                    "Vector::shrink_to_fit: number of resize + shrink_to_fit calls exceeds limit.\n" );

            mVector.resize( aCount, aValue );

            MORIS_CHECK_MEMORY(
                    sizeof( T ) * this->capacity() > MORIS_CELL_RESIZE_CHECK_LIMIT * MORIS_MAX_CELL_CAPACITY ? this->size() > MORIS_CELL_UTILIZATION_FRACTION_LIMIT * this->capacity() : true,
                    "Vector::resize: After resize memory used is less than 75 percent of capacity of large matrix: size %zu capacity %zu\n",
                    this->size(),
                    this->capacity() );
        }

        //------------------------------------------------------------------

        /**
         * @brief Returns an iterator to the first element
         */
        auto
        begin()
                -> decltype( mVector.begin() )
        {
            return mVector.begin();
        }

        //------------------------------------------------------------------

        /**
         * @brief Returns an iterator to the first element
         */

        auto
        end()
                -> decltype( mVector.end() )
        {
            return mVector.end();
        }


        //------------------------------------------------------------------

        /**
         * @brief Returns a const iterator to the first element
         */
        auto
        begin() const
                -> decltype( mVector.begin() )
        {
            return mVector.begin();
        }

        //------------------------------------------------------------------

        /**
         * @brief Returns a const iterator to the last element
         */

        auto
        end() const
                -> decltype( mVector.end() )
        {
            return mVector.end();
        }

        //------------------------------------------------------------------

        /**
         * @brief Returns a const iterator to the first element
         */

        auto
        cbegin() const
                -> decltype( mVector.cbegin() )
        {
            return mVector.cbegin();
        }

        //------------------------------------------------------------------

        /**
         * @brief Returns an iterator to the first element
         */

        auto
        cend() const
                -> decltype( mVector.cend() )
        {
            return mVector.cend();
        }

        //------------------------------------------------------------------

        /**
         * @brief insert function warpper for std::insert
         *
         * @tparam InputIterator
         * @param pos posistion we want to insert the new enrties
         * @param aFirst templated iterator/pointer begining of the inupt
         * @param aLast  templated iterator/pointer end of the inupt
         */
        template< class InputIterator >
        void
        insert( moris::size_t const & pos, InputIterator aFirst, InputIterator aLast )
        {
            // We use the insert just to assemble data from smaller containers, this check is to prevent chaning the cell size
            // Technically insert has the capability to resize the vector (it never rewrites data)
            MORIS_ASSERT( this->size() == pos, "Vector is being reallocated such that it can take in the new data" );

            // defer to call to the std::vecotr
            mVector.insert( mVector.begin() + pos, aFirst, aLast );
        }

        //------------------------------------------------------------------

        /*!
         * @brief empalce back method warpper for std::emplace_back
         *
         * @param value r values being passed
         */

        template< typename... _Args >
        void
        emplace_back( _Args&&... __args )
        {
            mVector.emplace_back( std::forward< _Args >( __args )... );
        }

        //------------------------------------------------------------------

    };    // class Vector

    //------------------------------------------------------------------

    // Free functions
    template< typename T >
    void
    unique( Vector< T >& aVector )
    {
        // get ref to data
        std::vector< T >& tVec = aVector.data();

        // sort data
        std::sort( tVec.begin(), tVec.end() );

        // trim vector
        tVec.erase( std::unique( tVec.begin(), tVec.end() ), tVec.end() );
    }

    //------------------------------------------------------------------

    // https://stackoverflow.com/questions/25921706/creating-a-vector-of-indices-of-a-sorted-vector
    //  extended to create a unique with the first index of unique values
    template< typename T >
    Vector< moris::moris_index >
    unique_index( Vector< T >& aVector )
    {
        std::vector< T > x = aVector.data();

        std::vector< int > y( x.size() );
        std::size_t        n( 0 );
        std::generate( std::begin( y ), std::end( y ), [ & ] { return n++; } );

        std::sort( std::begin( y ),
                std::end( y ),
                [ & ]( int i1, int i2 ) { return x[ i1 ] < x[ i2 ]; } );

        Vector< moris::moris_index > tUniqueInds;
        for ( moris::uint i = 0; i < aVector.size(); i++ )
        {
            if ( i == 0 )
            {
                tUniqueInds.push_back( y[ i ] );
            }

            else if ( aVector( y[ i - 1 ] ) != aVector( y[ i ] ) )
            {
                tUniqueInds.push_back( y[ i ] );
            }
        }

        return tUniqueInds;
    }

    //------------------------------------------------------------------

    /*!
     * Iterates through cell and prints each cell.
     * Will only work on data types that allow std::cout calls
     */
    template< typename T >
    void
    print(
            Vector< T > const & aVector,
            std::string       aStr = "Vector" )
    {
        std::cout << "Vector Name: " << aStr << "\n";
        std::cout << "Number of entries = " << aVector.size() << "\n";
        for ( moris::uint i = 0; i < aVector.size(); i++ )
        {
            std::cout << aVector( i ) << "\n";
        }

        std::cout << std::endl;
    }

    //------------------------------------------------------------------

    // General is_container trait with default value as false
    template< typename T, typename = void >
    struct is_moris_cell : std::false_type
    {
    };

    // Specialization for moris cell type
    template< typename T >
    struct is_moris_cell< Vector< T > > : std::true_type
    {
    };

    template< typename T >
    std::string
    print_nested_cells( const T& aVector )
    {
        // initialize the return string
        std::string aReturnStr = "";

        // Check if the container type is a nested container
        if constexpr ( is_moris_cell< T >::value )
        {
            aReturnStr += "[";

            // Keep track of the first element to avoid printing an extra comma
            bool first = true;
            for ( const auto& iElement : aVector )
            {
                if ( !first )
                {
                    aReturnStr += ", ";
                }

                // Recursively print each nested element
                aReturnStr += print_nested_cells( iElement );
                first = false;
            }
            aReturnStr += "]";
        }
        else
        {

            aReturnStr += std::to_string( aVector );
        }

        return aReturnStr;
    }

    //------------------------------------------------------------------
    /**
     * @brief prints out a row vector for infinite nested cells
     *
     * @tparam T
     * @param aVector
     * @param aStr
     */

    template< typename T >
    void
    print_as_row_vector(
            Vector< T > const & aVector,
            std::string       aStr = "Vector" )
    {
        std::cout << aStr << " = " << print_nested_cells( aVector ) << std::endl;
    }

    //------------------------------------------------------------------

    inline Vector< char >
    string_to_char( Vector< std::string >& strings );

    //------------------------------------------------------------------

    template< typename T >
    void
    shrink_to_fit_all( Vector< T >& aVector )
    {
        aVector.shrink_to_fit();
    }

    //------------------------------------------------------------------

    template< typename T >
    void
    shrink_to_fit_all( Vector< Vector< T > >& aVector )
    {
        // trim inner cells
        for ( uint iI = 0; iI < aVector.size(); ++iI )
        {
            aVector( iI ).shrink_to_fit_all();
        }

        // trim outer cell
        aVector.shrink_to_fit();
    }

    //------------------------------------------------------------------

    template< typename T >
    void
    write_to_txt_file(
            Vector< T > const & aVector,
            std::string       aFileName )
    {
        // Open a file stream for writing
        std::ofstream tOutFile( aFileName );

        // Check if the file was opened successfully
        MORIS_ASSERT( tOutFile.is_open(), "Failed to open output file" );

        // Write the vector's elements to the file in a vertical format
        std::copy( aVector.begin(), aVector.end(), std::ostream_iterator< T >( tOutFile, "\n" ) );

        // Close the file stream
        tOutFile.close();
    }
}    // namespace moris

#endif /* MORIS_CONTAINERS_CL_Vector_HPP_ */
