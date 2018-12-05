/*
 * cl_HMR_Paramfile.hpp
 *
 *  Created on: Nov 27, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_PARAMFILE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_PARAMFILE_HPP_

#include <string>
//#include "assert.hpp"

#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

#include "cl_Cell.hpp"

//#include "cl_Matrix.hpp"
//#include "linalg_typedefs.hpp"



#include "HMR_Globals.hpp"
#include "HMR_Tools.hpp"

#include "cl_HMR_State.hpp"
#include "cl_Map.hpp"

namespace moris
{
    class XML_Parser;

    namespace hmr
    {
// -----------------------------------------------------------------------------

        class Parameters;

// -----------------------------------------------------------------------------
        struct Mesh_Param
        {
            moris_id    mID     = gNoID;
            uint        mOrder  = 0;
            std::string mPath   = "";

            Mesh_Param(){};
            ~Mesh_Param(){};
        };
// -----------------------------------------------------------------------------

        struct Field_Param
        {
            std::string mLabel;
            moris_id    mID = gNoID;
            uint        mInputBSplineOrder = 0;
            uint        mOutputBSplineOrder = 0;
            std::string mSource;
            std::string mTarget;

            Field_Param(){};
            ~Field_Param(){};
        };

// -----------------------------------------------------------------------------

        class Paramfile
        {
            //! HMR runtime mode
            const enum State    mState;

            // the parser object
            XML_Parser *        mParser = nullptr;

            ParameterList       mParameterList;

            //! container with mesh settings
            Cell< Mesh_Param >  mMeshParams;

            //! container with field setting
            Cell< Field_Param > mFieldParams;

            //! list of target files for mappint
            //Cell< std::string > mTargets;

            map< moris_id, uint > mMeshMap;
            map< moris_id, uint > mFieldMap;

            std::string mInputDatabase  = "";
            std::string mOutputDatabase = "";
            std::string mCoefficients   = "";

            Matrix< IdMat > mFieldIDs;
            Matrix< IdMat > mMeshIDs;

            std::string  mLibraryPath = "";
            std::string  mUserFunction = "";
            Cell< real > mUserParameters;

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

            ParameterList &
            get_parameter_list();

// -----------------------------------------------------------------------------

            const std::string &
            get_input_db_path() const
            {
                return mInputDatabase;
            }

// -----------------------------------------------------------------------------

            const std::string &
            get_output_db_path() const
            {
                return mOutputDatabase;
            }

// -----------------------------------------------------------------------------

            const std::string &
            get_coefficient_db_path() const
            {
                return mCoefficients;
            }

// -----------------------------------------------------------------------------

            const uint
            get_number_of_meshes() const
            {
                return mMeshIDs.length();
            }

// -----------------------------------------------------------------------------

            const uint &
            get_mesh_order( const uint & aIndex ) const
            {
                return mMeshParams( mMeshMap.find( mMeshIDs( aIndex ) ) ).mOrder;
            }

// -----------------------------------------------------------------------------

            const std::string &
            get_mesh_path( const uint & aIndex ) const
            {
                return mMeshParams( mMeshMap.find( mMeshIDs( aIndex ) ) ).mPath;
            }

// -----------------------------------------------------------------------------
        private:
// -----------------------------------------------------------------------------

            void
            load_mesh_params();

// -----------------------------------------------------------------------------

            void
            load_field_params();

// -----------------------------------------------------------------------------

            void
            load_state_params();

// -----------------------------------------------------------------------------

            void
            load_parameter_list();

// -----------------------------------------------------------------------------

            void
            update_parameter_list();

// -----------------------------------------------------------------------------
        };

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_CL_HMR_PARAMFILE_HPP_ */
