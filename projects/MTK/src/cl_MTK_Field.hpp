/*
 * cl_MTK_Field.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_
#define PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_

#include <string>
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Block;
//------------------------------------------------------------------------------

        class Field
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! A short description of this field
            std::string    mLabel;

            //! an id that defines this field
            moris_id       mID;

            //! pointer to mesh or block object this field refers to
            const Block   * mBlock = nullptr;

            //! B-Spline coefficients of this field
            Mat< real >  * mCoefficients = nullptr;

            //! Node values of this field
            Mat< real >  * mNodeValues = nullptr;

            const bool mOwnNodeValues;

            //! Dimensionality of the field
            const uint     mNumberOfDimensions = 1;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Field(
                    const std::string & aLabel,
                    const moris_id      aID,
                    const Block *       aBlock ) :
                        mLabel( aLabel ),
                        mID( aID ),
                        mBlock( aBlock ),
                        mOwnNodeValues( true )
            {
                mCoefficients = new Mat<real>;
                mNodeValues = new Mat<real>;
            }
//------------------------------------------------------------------------------

            Field(
                    const std::string & aLabel,
                    const moris_id      aID,
                    const Block *       aBlock,
                    Mat<real>   *       aNodeValues ) :
                        mLabel( aLabel ),
                        mID( aID ),
                        mBlock( aBlock ),
                        mNodeValues( aNodeValues ),
                        mOwnNodeValues( false )
            {
                mCoefficients = new Mat<real>;
            }

//------------------------------------------------------------------------------

            virtual ~Field()
            {
                delete mCoefficients;
                if( mOwnNodeValues )
                {
                    delete mNodeValues;
                }
            };

//------------------------------------------------------------------------------

            Mat< real > *
            get_coefficients()
            {
                return mCoefficients;
            }

//------------------------------------------------------------------------------

            Mat< real > *
            get_node_values()
            {
                return mNodeValues;
            }

//------------------------------------------------------------------------------

            void
            evaluate_node_values( const Mat< real > & aCoefficients );

//------------------------------------------------------------------------------

            void
            evaluate_node_values();

//------------------------------------------------------------------------------

            uint
            get_interpolation_order() const;

//------------------------------------------------------------------------------

            uint
            get_number_of_dimensions() const
            {
                return mNumberOfDimensions;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_ */
