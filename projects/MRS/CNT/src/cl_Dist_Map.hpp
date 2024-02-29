/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Dist_Map.hpp
 *
 */

//#ifndef SRC_CONTAINERS_CL_DIST_MAP_HPP_
//#define SRC_CONTAINERS_CL_DIST_MAP_HPP_
//
//
//// C++ header files.
//#include <numeric>
//#include <utility>
//
//// Third-party header files.
//#include <boost/bimap.hpp>
//
//// moris library header files.
//#include "moris_typedefs.hpp" // COR/src
//#include "cl_Bi_Map.hpp" // CON/src
//#include "fn_zip.hpp" // CON/src
//
//namespace moris
//{
//    /**
//     * @brief Distributed map for LIDs (local indices) and GIDs (global identifiers).
//     *
//     * Wrapper around boost::bimap.
//     */
//    class Dist_Map : public moris::Bi_Map< moris::size_t, moris::lint>
//    {
//    public:
//
//        /*
//         * typedefs used to abbreviate code
//         */
//        typedef boost::bimap<moris::size_t, moris::lint> bimap_t;
//        typedef bimap_t::value_type pair_t;
//
//        //----------------------------------------------------------------------
//
//        /**
//         * Dist_Map default constructor
//         */
//        Dist_Map() = default;
//
//        //----------------------------------------------------------------------
//
//        /**
//         * Dist_Map default constructor
//         *
//         * In this constructor, the local indices will automatically be set
//         * to 0 through n-1, where n is the number of global ids passed in.
//         *
//         * @param[in] aGlbIds Global identifiers
//         */
//        template<typename T>
//        Dist_Map(
//                T const & aGlbIds)
//        : moris::Bi_Map<moris::size_t, moris::lint >()
//        {
//            T locInds(aGlbIds.size());
//            std::iota(locInds.begin(), locInds.end(), 0 );
//
//            for (auto tup: moris::zip(locInds, aGlbIds))
//            {
//                moris::size_t locInds;
//                moris::lint   aGlbIds;
//
//                boost::tie(locInds, aGlbIds) = tup;
//
//                this->insert(locInds, aGlbIds);
//            }
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * Dist_Map destructor
//         */
//        ~Dist_Map() = default;
//
//        //----------------------------------------------------------------------
//
//        /**
//         * Get the local index from the given global identifier
//         *
//         * @param[in] aGlbId Global identifier
//         *
//         * @return Local index
//         */
//        moris::size_t
//        locId(
//                moris::lint const & aGlbId)
//        {
//            return this->getLeft(aGlbId);
//        }
//
//        /*
//         * const version of above
//         */
//        const moris::size_t
//        locId(
//                moris::lint const & aGlbId) const
//        {
//            return this->getLeft(aGlbId);
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * Get the global identifier from the given local index
//         *
//         * @param[in] aLocId Local index
//         *
//         * @return Global identifier
//         */
//        moris::lint
//        glbId(
//                moris::size_t const & aLocId)
//        {
//            return this->getRight(aLocId);
//        }
//
//        /*
//         * const version of above
//         */
//        const moris::lint
//        glbId(
//                moris::size_t const & aLocId) const
//        {
//            return this->getRight(aLocId);
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * @brief Insert elements.
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
//         * @param[in] aLocId Local index of new pair
//         * @param[in] aGlbId Global identifier of new pair
//         */
//        void
//        insert(
//                moris::size_t aLocId,
//                moris::lint   aGlbId )
//        {
//            this->mBiMap.insert(pair_t(aLocId, aGlbId));
//        }
//
//        void
//        insert(
//                std::pair<moris::size_t, moris::lint> pair)
//        {
//            this->insert(pair.first, pair.second);
//        }
//
//        //----------------------------------------------------------------------
//
//        /**
//         * @brief Access data
//         *
//         * Returns the memory used internally by the container.
//         *
//         * @return mBiMap
//         */
//        const bimap_t &
//        data() const
//        {
//            return this->mBiMap;
//        }
//
//        /**
//         * Non-const version of data
//         */
//        bimap_t &
//        data()
//        {
//            return const_cast<bimap_t &>(static_cast<moris::Dist_Map const *>(this)->data());
//        }
//    };
//
//}// namespace moris
//
//#endif /* SRC_CONTAINERS_CL_DIST_MAP_HPP_ */

