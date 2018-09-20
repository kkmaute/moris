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

            //! pointer to mesh or block object this field refers to
            const Block   * mBlock = nullptr;

            //! B-Spline coefficients of this field
            Mat< real >  mCoefficients;

            //! Node values of this field
            Mat< real >  mNodeValues;

            //! Dimensionality of the field
            const uint     mNumberOfDimensions = 1;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Field(
                    const std::string & aLabel,
                    const Block *       aBlock ) :
                        mLabel( aLabel ),
                        mBlock( aBlock )
            {

            }

//------------------------------------------------------------------------------

            virtual ~Field()
            {

            };

//------------------------------------------------------------------------------

            Mat< real > &
            get_coefficients()
            {
                return mCoefficients;
            }

//------------------------------------------------------------------------------

            const Mat< real > &
            get_coefficients() const
            {
                return mCoefficients;
            }

//------------------------------------------------------------------------------

            virtual Mat< real > &
            get_node_values()
            {
                return mNodeValues;
            }

//------------------------------------------------------------------------------

            virtual const Mat< real > &
            get_node_values() const
            {
                return mNodeValues;
            }

//------------------------------------------------------------------------------

            const std::string &
            get_label() const
            {
                return mLabel;
            }

//------------------------------------------------------------------------------

            void
            evaluate_scalar_function(
                    real (*aFunction)( const Mat< real > & aPoint ) );

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

            virtual const Block *
            get_block() const
            {
                return mBlock;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_ */
