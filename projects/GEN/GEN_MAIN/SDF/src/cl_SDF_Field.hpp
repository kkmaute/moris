/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Field.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_FIELD_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_FIELD_HPP_

#include <string>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Mesh.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        class Field
        {
            std::string        mLabel;

            // pointer to SDF Mesh Wrapper
            Mesh             & mMesh;

            // pointer to node data
            Matrix< DDRMat > & mData;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            // constructor
            Field(
                    const std::string & aLabel,
                    Mesh              & aMesh,
                    Matrix< DDRMat >  & aData );

//-------------------------------------------------------------------------------

            // destructor
            ~Field(){};

//-------------------------------------------------------------------------------

            void
            save_field_to_hdf5( const std::string & aFilePath );

//-------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_FIELD_HPP_ */

