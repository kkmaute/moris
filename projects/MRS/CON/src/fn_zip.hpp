/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_zip.hpp
 *
 */

#ifndef SRC_CONTAINERS_FN_ZIP_HPP_
#define SRC_CONTAINERS_FN_ZIP_HPP_

// Third-party header files.
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

namespace moris
{
    /**
     * @brief Variadic zip function.
     *
     * @see [Sequence-zip function for c++11?]
     * (http://stackoverflow.com/questions/8511035/sequence-zip-function-for-c11)
     *
     * @param[in] containers Variadic number of containers.
     *
     * @return Boost zip.
     */
    template<typename... T>
    auto
    zip(
            const T&... containers)
    -> boost::iterator_range<boost::zip_iterator<decltype(boost::make_tuple(std::begin(containers)...))>>
    {
        auto zip_begin = boost::make_zip_iterator(boost::make_tuple(std::begin(containers)...));
        auto zip_end = boost::make_zip_iterator(boost::make_tuple(std::end(containers)...));
        return boost::make_iterator_range(zip_begin, zip_end);
    }

}// namespace moris

#endif /* SRC_CONTAINERS_FN_ZIP_HPP_ */

