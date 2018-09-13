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
        class Mesh;
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

            //! pointer to mesh object this field refers to
            const Mesh   * mMesh;

            //! B-Spline coefficients of this field
            Mat< real >  * mCoefficients = nullptr;

            //! Node values of this field
            Mat< real >  * mNodeValues = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Field(
                    const std::string & aLabel,
                    const moris_id      aID,
                    const Mesh *        aMesh ) :
                        mLabel( aLabel ),
                        mID( aID ),
                        mMesh( aMesh )
            {
                mCoefficients = new Mat< real >;
                mNodeValues   = new Mat< real >;
            }

//------------------------------------------------------------------------------

            virtual ~Field()
            {
                delete mCoefficients;
                delete mNodeValues;
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
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_FIELD_CPP_ */
