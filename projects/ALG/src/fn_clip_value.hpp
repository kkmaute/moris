#ifndef SRC_ALGORITHMS_FN_CLIP_VALUE_HPP_
#define SRC_ALGORITHMS_FN_CLIP_VALUE_HPP_

#include <iostream>

// MORIS header files.
#include "typedefs.hpp"    // COR/src
#include "linalg_typedefs.hpp"

namespace moris
{
    /**
     * clip_value function
     *
     * clips value based on threshold, i.e., small value close to zero maintaining sign
     *
     *    if abs(value) < threshold && value >= 0 :  threshold
     *    if abs(value) < threshold && value <  0 : -threshold
     *    otherwise:                                 aValue
     *
     * @param[in] aValue     value to be clipped
     * @param[in] aThreshold threshold, default: MORIS_REAL_EPS
     *
     * @return clipped value
     */
    inline real
    clip_value(
            const real&       aValue,
            const real&       aThreshold = MORIS_REAL_EPS,
            const real&       aNominator = MORIS_REAL_MAX,
            const std::string aCaller    = __builtin_FUNCTION() )
    {
        if ( std::abs( aValue ) < aThreshold )
        {
            std::ostringstream streamObj1;
            std::ostringstream streamObj2;

            streamObj1 << std::scientific << std::setprecision( 15 ) << aValue;
            streamObj2 << std::scientific << std::setprecision( 15 ) << aThreshold;

            std::string str = "Clipping " + streamObj1.str() + " with " + streamObj2.str() + " Caller: " + aCaller;

            real tRetValue = aValue < 0.0 ? -aThreshold : aThreshold;

            if ( aNominator < 0.99 * MORIS_REAL_MAX )
            {
                real tExact = aNominator / aValue;
                real tApprx = aNominator / tRetValue;

                if ( std::abs( tExact - tApprx ) > 0.1 * std::abs( tExact ) )
                {
                    //                    std::ostringstream streamObj1;
                    //                    std::ostringstream streamObj2;
                    //                    std::ostringstream streamObj3;
                    //
                    //                    streamObj1 << std::scientific << std::setprecision( 15 ) << aValue;
                    //                    streamObj2 << std::scientific << std::setprecision( 15 ) << tRetValue;
                    //                    streamObj3 << std::scientific << std::setprecision( 15 ) << aNominator;

                    //                    std::cout << "Clipping error exceeds 10 perc : exact value " << streamObj1.str() << " clipped value " << streamObj2.str() << " nominator " << streamObj3.str() << std::endl;

                    std::cout << str << " E! " << std::endl;
                }
                else
                {
                    std::cout << str << std::endl;
                }
            }

            return tRetValue;
        }

        return aValue;
    }
}    // namespace moris

#endif /* SRC_ALGORITHMS_FN_CLIP_VALUE_HPP_ */
