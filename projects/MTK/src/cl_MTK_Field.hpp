/*
 * cl_MTK_Field.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_
#define PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_

#include <string>
#include <memory>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Mesh;

        class Field
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! pointer to mesh or block object this field refers to
            //const Block   * mBlock = nullptr;
            const std::shared_ptr< Mesh > mMesh;

            //! Dimensionality of the field
            const uint     mNumberOfDimensions = 1;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
             //TODO: Add a notion of entity rank attached to field (not all fields are node based)
             //! A short description of this field
            std::string    mLabel;

            //! B-Spline coefficients of this field
            Matrix< DDRMat >  mCoefficients;

            //! Node values of this field
            Matrix< DDRMat >  mNodeValues;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Field(  const  std::string             & aLabel,
                    const  std::shared_ptr< Mesh >   aMesh ) :
                        mMesh( aMesh )
            {
                this->set_label( aLabel );
            }

//------------------------------------------------------------------------------

            virtual ~Field()
            {

            };

//------------------------------------------------------------------------------

            virtual Matrix< DDRMat > &
            get_coefficients()
            {
                return mCoefficients;
            }

//------------------------------------------------------------------------------

            virtual const Matrix< DDRMat > &
            get_coefficients() const
            {
                return mCoefficients;
            }

//------------------------------------------------------------------------------

            virtual Matrix< DDRMat > &
            get_node_values()
            {
                return mNodeValues;
            }

//------------------------------------------------------------------------------

            virtual const Matrix< DDRMat > &
            get_node_values() const
            {
                return mNodeValues;
            }

//------------------------------------------------------------------------------

            virtual const std::string &
            get_label() const
            {
                return mLabel;
            }

//------------------------------------------------------------------------------

            virtual void
            set_label( const std::string & aLabel )
            {
                mLabel = aLabel;
            }

//------------------------------------------------------------------------------

            void
            evaluate_scalar_function(
                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) );

//------------------------------------------------------------------------------

            void
            evaluate_node_values( const Matrix< DDRMat > & aCoefficients );

//------------------------------------------------------------------------------

            void
            evaluate_node_values();

//------------------------------------------------------------------------------

            uint
            get_number_of_dimensions() const
            {
                return mNumberOfDimensions;
            }

//------------------------------------------------------------------------------

            Interpolation_Order
            get_interpolation_order() const;

//------------------------------------------------------------------------------

            uint
            get_num_nodes() const
            {
                return mNodeValues.length();
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_ */
