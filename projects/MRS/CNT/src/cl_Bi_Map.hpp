/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bi_Map.hpp
 *
 */

//#ifndef SRC_CONTAINERS_CL_BI_MAP_HPP_
//#define SRC_CONTAINERS_CL_BI_MAP_HPP_
//// C++ header files.
//#include <numeric>
//#include <utility>
//
//// Third-party header files.
//#include <boost/bimap.hpp>
//
//// moris library header files.
//#include "moris_typedefs.hpp" // COR/src
//
//namespace moris
//{
//    /**
//     * wrapper class for boost::bimap
//     *
//     * A Bi_Map can be thought of the combination of two std::map's.
//     * A Bi_Map has "left" and "right" sides and can map items from one side
//     * to the other. A Bi_Map should be used over a Map when two-way mapping
//     * will be used frequently (i.e. a Map is used primarily for one-way
//     * mapping and only occasional searching).
//     */
//    template<typename T1, typename T2 >
//    class Bi_Map
//    {
//    public:
//
//        /**
//         * default Bi_Map constructor
//         */
//        Bi_Map() = default;
//
//        /**
//         * Bi_Map destructor
//         */
//        ~Bi_Map() = default;
//
//        /**
//         * @brief Insert elements
//         *
//         * Extends the container by inserting new elements,
//         * effectively increasing the container size
//         * by the number of elements inserted.
//         *
//         * Because element keys in a map are unique, the insertion operation checks
//         * whether each inserted element has a key equivalent to the one of an element
//         * already in the container, and if so, the element is not inserted, returning
//         * an iterator to this existing element (if the function returns a value).
//         *
//         * @param[in] aLeftItem  Left-side item to be inserted
//         * @param[in] aRightItem Right-side item to be inserted
//         */
//        void
//        insert(
//                T1 aLeftItem,
//                T2 aRightItem )
//        {
//            mBiMap.insert(typename boost::bimap<T1,T2>::value_type(aLeftItem, aRightItem));
//        }
//
//        void
//        insert(
//                std::pair<T1, T2> pair)
//        {
//            this->insert(pair.first, pair.second);
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * right-to-left map access
//         *
//         * @param[in] aRightItem Item from the right side of map
//         *
//         * @return Corresponding item on the left side of map
//         */
//        T1
//        getLeft(
//                T2 aRightItem )
//        {
//            return mBiMap.right.at(aRightItem);
//        }
//
//        /**
//         * const version of above
//         */
//        const T1
//        getLeft(
//                T2 aRightItem ) const
//        {
//            return mBiMap.right.at(aRightItem);
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * left-to-right map access
//         *
//         * @param[in] aLeftItem Item from the left side of map
//         *
//         * @return Corresponding item on the rightt side of map
//         */
//        T2
//        getRight(
//                T1 aLeftItem )
//        {
//            return mBiMap.left.at(aLeftItem);
//        }
//
//        /**
//         * const version of above
//         */
//        const T2
//        getRight(
//                T1 aLeftItem ) const
//        {
//            return mBiMap.left.at(aLeftItem);
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * size of the bi-direction map
//         */
//        moris::size_t
//        size()
//        {
//            return mBiMap.size();
//        }
//
//        /**
//         * const version of above
//         */
//        const moris::size_t
//        size() const
//        {
//            return mBiMap.size();
//        }
//
//        //----------------------------------------------------------------------
//
//    protected:
//        /**
//         * Underlying boost::bimap
//         */
//        boost::bimap<T1,T2> mBiMap;
//    };
//}
//
//#endif /* SRC_CONTAINERS_CL_BI_MAP_HPP_ */

