/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell.hpp
 *
 */

//#ifndef SRC_CONTAINERS_CL_XTK_CELL_HPP_
//#define SRC_CONTAINERS_CL_XTK_CELL_HPP_
//
//// C++ header files.
//#include <vector>
//
//// XTK library header files.
//#include "xtk_typedefs.hpp"
//
//namespace xtk
//{
//
//template<typename T>
//class Cell
//{
//private:
//
//    /**
//     * XTK cell
//     */
//    std::vector<T> mCell;
//
//public:
//
//    /**
//     * Cell constructor
//     */
//    Cell() = default;
//
//    /**
//     * Cell constructor
//     *
//     * @param[in] aCell A Cell
//     */
//    Cell(std::initializer_list<T> const & aCell) :
//            mCell(aCell.begin(), aCell.end())
//    {
//    }
//
//    /**
//     * Cell constructor
//     *
//     * @param[in] aSize Size of the Cell to be initialized
//     * @param[in] aValue Value of the elements of the initialized Cell
//     */
//    template<typename A>
//    Cell(uint const aSize, A const aValue) :
//            mCell(aSize, aValue)
//    {
//    }
//
//    Cell(uint const aSize) :
//            mCell(aSize)
//    {
//    }
//
//    /**
//     * Cell destructor
//     */
//    ~Cell() = default; // 'default' tells the compiler to automatically
//                       // delete the underlying Cell
//
//    /**
//     * @brief Operator to replace the contents of the container
//     *
//     * @param[in] val Contents to be copied into the container
//     */
//    const Vector<T> &
//    operator=(Vector<T> const & val)
//    {
//        mCell = val.data();
//
//        return *this;
//    }
//
//    /**
//     * @brief Operator to replace the contents of the container
//     *
//     * @param[in] val Replacement Cell
//     *
//     * @return Assignment.
//     */
//    template<typename A>
//    const Vector<T> &
//    operator=(const A & val)
//    {
//        mCell = val;
//
//        return *this;
//    }
//
//    /**
//     * Replaces the contents of the container
//     *
//     * @param[in] aCount number of entries
//     * @param[in] aval Value to set as contents of container
//     */
//    void assign(size_t const aCount, T const & aval)
//    {
//        mCell.assign(aCount, aval);
//    }
//
//    /**
//     * @brief Operator for accessing a specified Cell element.
//     *
//     * @param[in] i_index Position of an element in the Cell.The function
//     *            at() automatically checks whether an element is within the
//     *            bounds of valid elements in the container, throwing an
//     *            out_of_range exception if it is not. This is in contrast
//     *            with member operator[], that does not check against bounds
//     *            and hence is faster.
//     */
//    T & operator()(size_t const i_index)
//
//    {
//#ifdef MORIS_HAVE_DEBUG
//        return (mCell.at(i_index));
//#else
//        return( mCell[ i_index ] );
//#endif
//    }
//
//    /**
//     * const version of above
//     */
//    T const & operator()(size_t const i_index) const
//    {
//#ifdef MORIS_HAVE_DEBUG
//        return (mCell.at(i_index));
//#else
//        return( mCell[ i_index ] );
//#endif
//    }
//
//    /**
//     * Direct access to the underlying array.
//     *
//     * @return Direct access to the underlying array.
//     */
//    std::vector< T > const &
//    data() const
//    {
//        return mCell;
//    }
//
//    std::vector< T > &
//    data()
//    {
//        return mCell;
//    }
//
//    /**
//     * Checks whether the Cell container is empty.
//     *
//     * @return bool value indicating whether the Cell container is empty.
//     */
//    bool empty()
//    {
//        return mCell.empty();
//    }
//   // const bool empty() const
//   //  {
//   //    return mCell.empty();
//   // }
//
//    /**
//     * Returns the last element of the Cell.
//     *
//     * @return the last element of the Cell.
//     */
//    T const & back()
//    {
//        return mCell.back();
//    }
//
//    /**
//     * Returns the number of elements in the Cell.
//     *
//     * @return Number of elements in the Cell.
//     */
//    // size_t
//    size_t size()
//    {
//        return mCell.size();
//    }
//
//    /**
//     * const version of above
//     */
//    size_t size() const
//    {
//        return mCell.size();
//    }
//
//    /**
//     * Returns the maximum possible number of elements in the Cell.
//     *
//     * @return Maximum possible number of elements in the Cell.
//     */
//    size_t max_size()
//    {
//        return mCell.max_size();
//    }
//
//    /**
//     * @brief Increase the capacity of the container to a value that's
//     *        greater or equal to new_cap
//     *
//     * @param[in] new_cap The new capacity of the Cell
//     */
//    void reserve(size_t const & new_cap)
//    {
//        mCell.reserve(new_cap);
//    }
//
//    /**
//     * @brief Returns the number elements that can be held in currently
//     *        allocated storage.
//     *
//     * @return Number elements that can be held in currently allocated storage.
//     */
//    size_t capacity() const
//    {
//        return mCell.capacity();
//    }
//
//    /**
//     * @brief Requests the removal of unused capacity
//     */
//    void shrink_to_fit()
//    {
//        mCell.shrink_to_fit();
//    }
//
//    /**
//     * Clears the contents of the Cell.
//     */
//    void clear()
//    {
//        mCell.clear();
//    }
//
//    /**
//     * @brief Inserts elements at the specified location in the container
//     *
//     * @param[in] pos Iterator at which the content will be inserted.
//     * @param[in] value Element value to insert
//     */
//    void insert(size_t const & pos, T const & value)
//    {
//        mCell.insert(mCell.begin() + pos, value);
//    }
//
//    /**
//     * Appends a Cell, similar to how push_back appends a single value.
//     *
//     * @param[in] aCell Cell to be appended
//     */
//    void append(Vector<T> const & aCell)
//    {
//        mCell.insert(mCell.end(), aCell.data().begin(), aCell.data().end());
//    }
//
//    /**
//     * @brief Removes specified elements from the container
//     *
//     * @param[in] pos Iterator to the element to remove
//     */
//    void erase(size_t const & pos)
//    {
//        mCell.erase(mCell.begin() + pos);
//    }
//
//    /**
//     * @brief Appends the given element value to the end of the container
//     *
//     * @param[in] value The value of the element to append
//     */
//    void push_back(T const & value)
//    {
//        mCell.push_back(value);
//    }
//
//    /**
//     * @brief removes the last element of the container
//     */
//    void pop_back()
//    {
//        mCell.pop_back();
//    }
//
//    /**
//     * @brief Resizes the container to contain count elements
//     *
//     * @param[in] aCount The new size of the Cell
//     * @param[in] aValue The value to initialize the new elements with
//     */
//    void resize(size_t const & aCount, T const & aValue)
//    {
//        mCell.resize(aCount, aValue);
//    }
//
//    auto begin()
//    ->decltype(mCell.begin())
//    {
//        return mCell.begin();
//    }
//
//    auto end()
//    ->decltype(mCell.end())
//    {
//        return mCell.end();
//    }
//};
//}
//
//
//#endif /* SRC_CONTAINERS_CL_XTK_CELL_HPP_ */

