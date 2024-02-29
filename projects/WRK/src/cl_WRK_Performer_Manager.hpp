/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Performer_Manager.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_

#include "moris_typedefs.hpp"                       //MRS/COR/src
#include "cl_Vector.hpp"                            //MRS/CNT/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris
{
    class Library_IO;
    //------------------------------------------------------------------------------
    namespace hmr
    {
        class HMR;
    }
    namespace mtk
    {
        class Mesh_Manager;
    }
    namespace gen
    {
        class Geometry_Engine;
    }
    namespace xtk
    {
        class Model;
    }
    namespace mdl
    {
        class Model;
    }
    namespace opt
    {
        class Manager;
    }

    namespace wrk
    {
        class Workflow;
        class Workflow_HMR_XTK;
        class Workflow_STK_XTK;
        class Workflow_STK_FEM;
        class Remeshing_Mini_Performer;
        class Reinitialize_Performer;
        class DataBase_Performer;

        //------------------------------------------------------------------------------

        class Performer_Manager
        {
                std::shared_ptr< Library_IO > mLibrary = nullptr;

                Vector< std::shared_ptr< mtk::Mesh_Manager > >        mMTKPerformer;
                Vector< std::shared_ptr< hmr::HMR > >                 mHMRPerformer;
                Vector< std::shared_ptr< gen::Geometry_Engine > >      mGENPerformer;
                Vector< std::shared_ptr< xtk::Model > >               mXTKPerformer;
                Vector< std::shared_ptr< mdl::Model > >               mMDLPerformer;
                Vector< std::shared_ptr< opt::Manager > >             mOPTPerformer;

                Vector< std::shared_ptr< wrk::Remeshing_Mini_Performer > > mRemeshingMiniPerformer;
                Vector< std::shared_ptr< wrk::Reinitialize_Performer > >   mReinitializePerformer;
                Vector< std::shared_ptr< wrk::DataBase_Performer > >       mDataBasePerformer;

                friend class wrk::Workflow;
                friend class wrk::Workflow_HMR_XTK;
                friend class wrk::Workflow_STK_XTK;
                friend class wrk::Workflow_STK_FEM;

            public:

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 * @param[ in ] aMesh          mesh for this problem
                 */
                Performer_Manager( std::shared_ptr< Library_IO > aLibrary );

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Performer_Manager();
                //------------------------------------------------------------------------------
        };
    } /* namespace mdl */
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_WRK_PERFORMER_MANAGER_HPP_ */

