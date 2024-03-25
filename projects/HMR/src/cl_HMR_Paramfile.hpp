/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Paramfile.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_PARAMFILE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_PARAMFILE_HPP_

#include <string>

#include "cl_HMR_Field_Param.hpp"
#include "cl_HMR_State.hpp"
#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"

#include "cl_Communication_Tools.hpp"
#include "moris_typedefs.hpp"

#include "cl_Vector.hpp"

//#include "cl_Matrix.hpp"
//#include "linalg_typedefs.hpp"

#include "cl_Map.hpp"

namespace moris
{
    class XML_Parser;
}
namespace moris::hmr
{
// -----------------------------------------------------------------------------

    class Parameters;

// -----------------------------------------------------------------------------
    struct Mesh_Param
    {
        moris_id    mID     = gNoID;
        uint        mOrder  = 0;
        std::string mPath;
    };

// -----------------------------------------------------------------------------

    class Paramfile
    {
        //! HMR runtime mode
        const enum State    mState;

        // the parser object
        XML_Parser *        mParser = nullptr;

        Parameter_List mParameterList;

        //! container with mesh settings
        Vector< Mesh_Param >  mMeshParams;

        //! container with field setting
        Vector< Field_Param > mFieldParams;

        //! list of target files for mappint
        //Vector< std::string > mTargets;

        map< moris_id, uint > mMeshMap;
        map< moris_id, uint > mFieldMap;

        std::string mInputDatabase;
        std::string mOutputDatabase;
        std::string mCoefficients;

        Matrix< IdMat > mFieldIDs;
        Matrix< IdMat > mMeshIDs;

        std::string  mLibraryPath;
        std::string  mUserFunction;
        std::string  mUnionMesh;

        sint mInitialBSplineRefinement     = -1;
        sint mAdditionalLagrangeRefinement = -1;

// -----------------------------------------------------------------------------
    public:
// -----------------------------------------------------------------------------

        Paramfile( const std::string & aPath, const enum State aState );

// -----------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Paramfile();

// -----------------------------------------------------------------------------

        Parameter_List & get_parameter_list();

// -----------------------------------------------------------------------------

        const std::string & get_input_db_path() const
        {
            return mInputDatabase;
        }

// -----------------------------------------------------------------------------

        const std::string & get_output_db_path() const
        {
            return mOutputDatabase;
        }

// -----------------------------------------------------------------------------

        const std::string & get_coefficient_db_path() const
        {
            return mCoefficients;
        }

// -----------------------------------------------------------------------------

        const std::string & get_library_path() const
        {
            return mLibraryPath;
        }

// -----------------------------------------------------------------------------

        const std::string & get_user_function_name() const
        {
            return mUserFunction;
        }

// -----------------------------------------------------------------------------

        uint get_number_of_meshes() const
        {
            return mMeshIDs.length();
        }

// -----------------------------------------------------------------------------

        uint get_mesh_order( uint aIndex ) const
        {
            return mMeshParams( mMeshMap.find( mMeshIDs( aIndex ) ) ).mOrder;
        }

// -----------------------------------------------------------------------------

        const std::string & get_mesh_path( uint aIndex ) const
        {
            return mMeshParams( mMeshMap.find( mMeshIDs( aIndex ) ) ).mPath;
        }

// -----------------------------------------------------------------------------

        uint get_number_of_fields() const
        {
            return mFieldIDs.length();
        }

// -----------------------------------------------------------------------------

        const Field_Param & get_field_params( uint aIndex ) const
        {
            return mFieldParams( mFieldMap.find( mFieldIDs( aIndex ) ) );
        }

// -----------------------------------------------------------------------------

        const std::string & get_union_mesh_path() const
        {
            return mUnionMesh;
        }

// -----------------------------------------------------------------------------
    private:
// -----------------------------------------------------------------------------

        void load_mesh_params();

// -----------------------------------------------------------------------------

        void
        load_field_params();

// -----------------------------------------------------------------------------

        void load_state_params();

// -----------------------------------------------------------------------------

        void load_parameter_list();

// -----------------------------------------------------------------------------

        void update_parameter_list();

// -----------------------------------------------------------------------------

        void load_user_refinement_parameters();

// -----------------------------------------------------------------------------
    };

// -----------------------------------------------------------------------------
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_PARAMFILE_HPP_ */

