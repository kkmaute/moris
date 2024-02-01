/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_STK.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_STK_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_STK_HPP_

#include <string>
#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Mesh.hpp"
#include "cl_MTK_Fields_Info.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        class STK
        {
            sdf::Mesh           & mMesh;

            //! struc required by MTK
            mtk::MtkMeshData      mMeshData;

            //! Container of all fields
            mtk::MtkFieldsInfo    mFieldsInfo;

            // Fields the field info contain references
            moris::Cell<mtk::Scalar_Field_Info<DDRMat>> mFields;

            //! set this to 3D
            uint                  mNumberOfDimensions = 3;

            //! connectivity passed to MTK
            Matrix< IdMat >           mElementTopology;

            //! element IDs passed to MTK
            Matrix< IdMat >           mElementLocalToGlobal;

            //! node IDs passed to MTK
            Matrix< IdMat >           mNodeLocalToGlobal;

            //! node coordinates passed to MTK
            Matrix< DDRMat >          mNodeCoords;

            //! node owners passed to MTK
            Matrix< IdMat >           mNodeOwner;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * constructor
             */
            STK( sdf::Mesh & aMesh ) : mMesh( aMesh ) {};

//-------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~STK(){};

//-------------------------------------------------------------------------------

            void
            create_mesh_data(
                    moris::Cell< Matrix<DDRMat> > & aFields,
                    moris::Cell< std::string >    & aFieldLabels,
                    const double                    aTimeStep=0.0 );

//-------------------------------------------------------------------------------

            void
            save_to_file( const std::string & aFilePath );

//-------------------------------------------------------------------------------
        };
//-------------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_STK_HPP_ */

