#ifndef MORIS_CONTAINERS_CL_CELL_HPP_
#define MORIS_CONTAINERS_CL_CELL_HPP_

// C++ header files.
#include <vector>
#include <algorithm> // for unique
#include <iostream>


// MORIS library header files.
#include "typedefs.hpp" // COR/src

namespace moris
{

    template< typename T >
    class Cell
    {
    private:

        /**
         * MORIS cell
         */
        std::vector< T > mCell;

    public:

        /**
         * moris::Cell constructor
         */
        Cell() = default;

        /**
         * moris::Cell constructor
         *
         * @param[in] aCell A Cell
         */
        Cell(
                std::initializer_list< T > const & aCell )
        : mCell( aCell.begin(), aCell.end() )
        {
        }

        /**
         * moris::Cell constructor
         *
         * @param[in] aSize Size of the Cell to be initialized
         * @param[in] aValue Value of the elements of the initialized Cell
         */
        template< typename A >
        Cell(
                moris::uint const aSize,
                A           const aValue )
        : mCell( aSize, aValue )
        {
        }

        Cell(
                moris::uint const aSize )
        : mCell( aSize )
        {
        }

        /**
         * moris::Cell destructor
         */
        ~Cell() = default; // 'default' tells the compiler to automatically
                             // delete the underlying Cell

        /**
         * @brief Operator to replace the contents of the container
         *
         * @param[in] val Contents to be copied into the container
         */
        const moris::Cell< T > &
        operator=(moris::Cell<T> const & val)
        {
            mCell = val.data();

            return *this;
        }

        /**
         * @brief Operator to replace the contents of the container
         *
         * @param[in] val Replacement Cell
         *
         * @return Assignment.
         */
        template< typename A >
        const moris::Cell< T > &
        operator=(
                const A & val )
        {
            mCell = val;

            return *this;
        }

        /**
         * Replaces the contents of the container
         *
         * @param[in] aCount number of entries
         * @param[in] aval Value to set as contents of container
         */
        void
        assign(
                moris::size_t const aCount,
                T             const & aval)
        {
            mCell.assign( aCount, aval );
        }

        /**
         * @brief Operator for accessing a specified Cell element.
         *
         * @param[in] i_index Position of an element in the Cell.The function
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
#ifndef NDEBUG
                ( mCell.at( i_index ) )
#else
                ( mCell[ i_index ] )
#endif
        )
        {
#ifndef NDEBUG
            return( mCell.at( i_index ) );
#else
            return( mCell[ i_index ] );
#endif
        }

        /**
         * const version of above
         */
        auto
        operator()(
                moris::size_t const i_index ) const
        -> decltype(
#ifndef NDEBUG
                ( mCell.at( i_index ) )
#else
                ( mCell[ i_index ] )
#endif
        )
        {
#ifndef NDEBUG
            return( mCell.at( i_index ) );
#else
            return( mCell[ i_index ] );
#endif
        }

        /**
         * Direct access to the underlying array.
         *
         * @return Direct access to the underlying array.
         */
        std::vector< T > const &
        data() const
        {
            return mCell;
        }

        std::vector< T > &
        data()
        {
            return mCell;
        }

        /**
         * Checks whether the Cell container is empty.
         *
         * @return bool value indicating whether the Cell container is empty.
         */
        bool
        empty()
        {
            return mCell.empty();
        }

        /**
         * Returns the last element of the Cell.
         *
         * @return the last element of the Cell.
         */
        auto
        back()
        ->decltype( mCell.back() )
        {
            return mCell.back();
        }

        /**
         * Returns the number of elements in the Cell.
         *
         * @return Number of elements in the Cell.
         */
        // moris::size_t
        auto
        size()
        ->decltype (mCell.size() )
        {
            return mCell.size();
        }

        /**
         * const version of above
         */
        auto
        size() const
        ->decltype (mCell.size() )
        {
            return mCell.size();
        }

        /**
         * Returns the maximum possible number of elements in the Cell.
         *
         * @return Maximum possible number of elements in the Cell.
         */
        auto
        max_size()
        ->decltype (mCell.max_size() )
        {
            return mCell.max_size();
        }

        /**
         * @brief Increase the capacity of the container to a value that's
         *        greater or equal to new_cap
         *
         * @param[in] new_cap The new capacity of the Cell
         */
        void
        reserve(
                moris::size_t const & new_cap)
        {
            mCell.reserve( new_cap );
        }

        /**
         * @brief Returns the number elements that can be held in currently
         *        allocated storage.
         *
         * @return Number elements that can be held in currently allocated storage.
         */
        auto
        capacity()
        ->decltype (mCell.capacity() )
        {
            return mCell.capacity();
        }

        /**
         * @brief Requests the removal of unused capacity
         */
        void
        shrink_to_fit()
        {
            mCell.shrink_to_fit();
        }

        /**
         * Clears the contents of the Cell.
         */
        void
        clear()
        {
            mCell.clear();
        }

        /**
         * @brief Inserts elements at the specified location in the container
         *
         * @param[in] pos Iterator at which the content will be inserted.
         * @param[in] value Element value to insert
         */
        void
        insert(
                moris::size_t const & pos,
                T             const & value)
        {
            mCell.insert( mCell.begin() + pos, value );
        }

        /**
         * Appends a moris::Cell, similar to how push_back appends a single value.
         *
         * @param[in] aCell Cell to be appended
         */
        void
        append(
                moris::Cell < T > const & aCell )
        {
            mCell.insert( mCell.end(), aCell.data().begin(), aCell.data().end() );
        }

        /**
         * @brief Removes specified element from the container
         *
         * @param[in] pos Position (not index) of element to be remove
         */
        void
        erase(
                moris::size_t const & pos )
        {
            mCell.erase( mCell.begin() + pos );
        }

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
            mCell.erase( mCell.begin() + pos1, mCell.begin() + pos2 );
        }

        /**
         * @brief Appends the given element value to the end of the container
         *
         * @param[in] value The value of the element to append
         */
        void
        push_back(
                T const & value )
        {
            mCell.push_back( value );
        }

        /**
         * @brief removes the last element of the container
         */
        void
        pop_back()
        {
            mCell.pop_back();
        }

        /**
         * @brief Resizes the container to contain count elements
         *
         * @param[in] aCount The new size of the Cell
         * @param[in] aValue The value to initialize the new elements with
         */
        void
        resize(
                moris::size_t const & aCount,
                T                     aValue = T() )
        {
            mCell.resize( aCount, aValue  );
        }

        /**
         * @brief Returns an iterator to the first element
         */
        auto
        begin()
        -> decltype( mCell.begin() )
        {
            return mCell.begin();
        }

        /**
         * @brief Returns an iterator to the first element
         */
        auto
        end()
        -> decltype( mCell.end() )
        {
            return mCell.end();
        }
    };


    // Free functions
    template< typename T >
    void
    unique( Cell< T > & aCell )
    {
        // get ref to data
        std::vector< T > & tVec = aCell.data();

        // sort data
        std::sort( tVec.begin(), tVec.end() );

        // trim vector
        tVec.erase( std::unique( tVec.begin(), tVec.end() ), tVec.end() );
    }


    /*!
     * Iterates through cell and prints each cell.
     * Will only work on data types that allow std::cout calls
     */
    template< typename T >
    void
    print(Cell< T > & aCell,
          std::string aStr = "Cell")
    {
        std::cout<<"Cell Name: "<<aStr<<"\n";
        std::cout<<"Number of entries = "<<aCell.size()<<"\n";
        for(auto aEntry: aCell)
        {
            std::cout<<aEntry<<"\n";
        }

        std::cout<<std::endl;
    }
}

#endif /* MORIS_CONTAINERS_CL_Cell_HPP_ */
