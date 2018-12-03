
#include "assert.hpp"

#include "HMR_Tools.hpp"
#include "cl_HMR_Paramfile.hpp"

#include "cl_XML_Parser.hpp"

#include "cl_HMR_Parameters.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_unique.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

    Paramfile::Paramfile( const std::string & aPath, const enum State aState ) :
            mState( aState )
    {
        // create parser
        mParser = new XML_Parser( aPath );

        // load meshes
        this->load_mesh_params();

        // load fields
        this->load_field_params();

        // load parameters
        this->load_parameter_list();

        // make parameters consistent with shipped fields and meshes
        this->update_parameter_list();

    }

// -----------------------------------------------------------------------------

    Paramfile::~Paramfile()
    {
        delete mParser;
    }
// -----------------------------------------------------------------------------

    Parameters *
    Paramfile::create_parameters()
    {
        return nullptr;
    }

// -----------------------------------------------------------------------------

    void
    Paramfile::load_mesh_params()
    {
        // count meshes
        uint tNumberOfMeshes
            =  mParser->count_keys_in_subtree( "moris.hmr", "mesh" );

        // allocate cell
        mMeshParams.resize( tNumberOfMeshes, Mesh_Param() );

        // clear map
        mMeshMap.clear();

        // list with IDs
        Matrix< IdMat > tIDs( tNumberOfMeshes, 1 );

        for( uint m=0; m<tNumberOfMeshes; ++m )
        {
            Cell< std::string > tFirst;
            Cell< std::string > tSecond;
            mParser->get_keys_from_subtree( "moris.hmr", "mesh", m, tFirst, tSecond );

            Mesh_Param & tMesh = mMeshParams( m );

            // copy key to settings struct
            for( uint k=0; k<tFirst.size(); ++k )
            {
                std::string tKey = tFirst( k );

                if ( tKey == "id" )
                {
                    tMesh.mID = stoi( tSecond( k ) );
                }
                else if ( tKey == "order" )
                {
                    tMesh.mOrder = stoi( tSecond( k ) );
                }
                else if ( tKey == "path" )
                {
                    tMesh.mPath = tSecond( k );
                }
            }

            tIDs( m ) = tMesh.mID;

            // add field to map
            mMeshMap[ tMesh.mID ] = m;

            // check for consistency
            MORIS_ERROR( tMesh.mID != gNoID, "Mesh ID is not set." );
            MORIS_ERROR( tMesh.mPath.size() > 0, "Mesh path must not be empty." );
            MORIS_ERROR( tMesh.mOrder > 0 && tMesh.mOrder <= 3, "Mesh order must be between 0 and 3." );

        }

        // make IDs unique
        Matrix< IdMat > tUniqueIDs;
        unique( tIDs, tUniqueIDs );

        // make sure that each ID exists only once
        MORIS_ERROR( tUniqueIDs.length() == tNumberOfMeshes,
                "Duplicate mesh IDs in input file" );
    }

// -----------------------------------------------------------------------------

        void
        Paramfile::load_field_params()
        {
            // count fields
            uint tNumberOfFields
                =  mParser->count_keys_in_subtree( "moris.hmr", "field" );

            // allocate cell
            mFieldParams.resize( tNumberOfFields, Field_Param() );

            // clear map
            mFieldMap.clear();

            // list with IDs
            Matrix< IdMat > tIDs( tNumberOfFields, 1 );

            for( uint f=0; f<tNumberOfFields; ++f )
            {
                Cell< std::string > tFirst;
                Cell< std::string > tSecond;
                mParser->get_keys_from_subtree( "moris.hmr", "field", f, tFirst, tSecond );

                Field_Param & tField = mFieldParams( f );

                // copy key to settings struct
                for( uint k=0; k<tFirst.size(); ++k )
                {
                    std::string tKey = tFirst( k );
                    if ( tKey == "label" )
                    {
                        tField.mLabel = tSecond( k );
                    }
                    else if ( tKey == "id" )
                    {
                        tField.mID = stoi( tSecond( k ) );
                    }
                    else if ( tKey == "input_bspline_order" )
                    {
                        tField.mInputBSplineOrder = stoi( tSecond( k ) );
                    }
                    else if ( tKey == "output_bspline_order" )
                    {
                        tField.mOutputBSplineOrder = stoi( tSecond( k ) );
                    }
                    else if ( tKey == "source" )
                    {
                        tField.mSource = tSecond( k );
                    }
                    else if ( tKey == "target" )
                    {
                        tField.mTarget = tSecond( k );
                    }
                }

                // copy ID into array
                tIDs( f ) = tField.mID;

                // add field to map
                mFieldMap[ tField.mID ] = f;

                // check for consistency
                MORIS_ERROR( tField.mID != gNoID, "Field ID is not set." );
                MORIS_ERROR( tField.mLabel.size() > 0, "Field label must not be empty." );
                MORIS_ERROR( tField.mSource.size() > 0, "Field source path must not be empty." );
                MORIS_ERROR( tField.mInputBSplineOrder <= gMaxBSplineOrder, "BSpline order of field too big." );
                MORIS_ERROR( tField.mOutputBSplineOrder <= gMaxBSplineOrder, "BSpline order of field too big." );
            }

            // make IDs unique
            Matrix< IdMat > tUniqueIDs;
            unique( tIDs, tUniqueIDs );

            // make sure that each ID exists only once
            MORIS_ERROR( tUniqueIDs.length() == tNumberOfFields,
                    "Duplicate field IDs in input file" );
        }

// -----------------------------------------------------------------------------

        void
        Paramfile::load_state_params()
        {

            Cell< std::string > tFirst;
            Cell< std::string > tSecond;

            // string for fields
            //std::string tFields;
            switch( mState )
            {
                case( State::INITIALIZE_MESH ) :
                {
                    mParser->get_keys_from_subtree( "moris.hmr",
                            "initialize", 0, tFirst, tSecond );

                    break;
                }
                case( State::REFINE_MESH ) :
                {
                    mParser->get_keys_from_subtree( "moris.hmr",
                            "refine", 0, tFirst, tSecond );
                    break;
                }
                case( State::MAP_FIELDS ) :
                {
                    mParser->get_keys_from_subtree( "moris.hmr",
                            "map", 0, tFirst, tSecond );
                    break;
                }
                default :
                {
                    break;
                }
            }

            // copy key to settings struct
            for( uint k=0; k<tFirst.size(); ++k )
            {
                std::string tKey = tFirst( k );
                if ( tKey == "input_database" )
                {
                    mOutputDatabase = tSecond( k );
                }
                else if ( tKey == "output_database" )
                {
                    mOutputDatabase = tSecond( k );
                }
                else if( tKey == "coefficients" )
                {
                    mCoefficients = tSecond( k );
                }
                else if( tKey == "meshes" )
                {
                    string_to_mat( tSecond( k ), mMeshIDs );
                }
                else if( tKey == "fields" )
                {
                    string_to_mat( tSecond( k ), mFieldIDs );
                }
                else if ( tKey == "library" )
                {
                    mLibraryPath = tSecond( k );
                }
                else if ( tKey == "function" )
                {
                    mUserFunction = tSecond( k );
                }
            }
        }

// -----------------------------------------------------------------------------

        void
        Paramfile::load_parameter_list()
        {
            // call parameter list initializer from cl_Parameters.hpp
            mParameterList = create_hmr_parameter_list();

            for( uint k=0; k<tFirst.size(); ++k )
            {
                std::string tKey = tFirst( k );

                if( tKey == "number_of_elements_per_dimension" )
                {
                    mParameterList.set("number_of_elements_per_dimension", tSecond( k ) );
                }
                else if( tKey == "domain_dimensions" )
                {
                    mParameterList.set( "domain_dimensions", tSecond( k )  );
                }
                else if( tKey == "domain_offset" )
                {
                    mParameterList.set("domain_offset", tSecond( k ) );
                }
                else if( tKey == "domain_sidesets" )
                {
                    mParameterList.set("domain_sidesets", tSecond( k ) );
                }
                else if( tKey == "refinement_buffer" )
                {
                    mParameterList.set( "refinement_buffer", ( sint ) std::stoi( tSecond( k ) ) );
                }
                else if( tKey == "staircase_buffer" )
                {
                    mParameterList.set( "staircase_buffer", ( sint ) std::stoi( tSecond( k ) ) );
                }
                else if( tKey == "bspline_orders" )
                {
                    mParameterList.set( "bspline_orders", tSecond( k ) );
                }
                else if( tKey == "lagrange_orders" )
                {
                    mParameterList.set( "lagrange_orders", tSecond( k ) );
                }
                else if( tKey == "initial_bspline_refinement" )
                {
                    mParameterList.set( "initial_bspline_refinement", ( sint ) std::stoi( tSecond( k ) ) );
                }
                else if ( tKey == "verbose" )
                {
                    mParameterList.set( "verbose", ( sint ) string_to_bool( tSecond( k ) ) );
                }
                else if ( tKey == "truncate_bsplines" )
                {
                    mParameterList.set( "truncate_bsplines", ( sint ) string_to_bool( tSecond( k ) ) );
                }
                else if ( tKey == "additional_lagrange_refinement" )
                {
                    mParameterList.set(  "additional_lagrange_refinement", ( sint ) std::stoi( tSecond( k ) ) );
                }
                else if ( tKey == "use_multigrid" )
                {
                    mParameterList.set(  "use_multigrid", ( sint ) string_to_bool( tSecond( k ) ) );
                }
                else if ( tKey == "max_refinement_level" )
                {
                    sint tValue = ( sint ) stoi( tSecond( k ) );
                    if( tValue > 0 )
                    {
                        mParameterList.set(  "max_refinement_level", tValue );
                    }
                    else
                    {
                        mParameterList.set(  "max_refinement_level", ( sint ) gMaxNumberOfLevels-1 );
                    }
                }
            }
        }

// -----------------------------------------------------------------------------

        // this function merges parameters in the parameters tag with those
        // from other tags to make everythig consistent
        void
        Paramfile::update_parameter_list()
        {
            // - - - - - - - - - - - - - - - - - - - - - - - -
            // Part 1: Make Lagrange orders consistent
            // - - - - - - - - - - - - - - - - - - - - - - - -

            Matrix< DDUMat > tOrdersFromXML;
            string_to_mat( mParameterList.get< std::string>("lagrange_orders"),
                    tOrdersFromXML );

            // allocate vector
            Matrix< DDUMat > tOrders(
                    tOrdersFromXML.length()
                    + mMeshParams.size(), 1 );

            // initialize counter
            uint tCount = 0;

            // copy orders from XML
            for( uint k=0; k<tOrdersFromXML.length(); ++k )
            {
                tOrders( tCount++ ) = tOrdersFromXML( k );
            }

            // add Lagrange orders from meshes
            for( uint k=0; k< mMeshParams.size(); ++k )
            {
                tOrders( tCount++ ) = mMeshParams( k ).mOrder;
            }

            // make orders unique
            Matrix< DDUMat > tOrdersUnique;
            unique( tOrders, tOrdersUnique );

            // convert matrix to string
            std::string tString;
            mat_to_string( tOrdersUnique, tString );

            // update entry in parameter list
            mParameterList.set( "lagrange_orders", tString );

            // - - - - - - - - - - - - - - - - - - - - - - - -
            // Part 2: Make B-Spline orders consistent
            // - - - - - - - - - - - - - - - - - - - - - - - -

            string_to_mat( mParameterList.get< std::string>("bspline_orders"),
                               tOrdersFromXML );

            // allocate vector
            tOrders.set_size(
                    tOrdersFromXML.length()
                    + 2*mFieldParams.size(), 1 );

            // reset counter
            tCount = 0;

            // copy orders from XML
            for( uint k=0; k<tOrdersFromXML.length(); ++k )
            {
                tOrders( tCount++ ) = tOrdersFromXML( k );
            }

            // add B-Spline orders from fields
            for( uint k=0; k< mFieldParams.size(); ++k )
            {
                tOrders( tCount++ ) = mFieldParams( k ).mInputBSplineOrder;
                tOrders( tCount++ ) = mFieldParams( k ).mOutputBSplineOrder;
            }

            // make orders unique
            unique( tOrders, tOrdersUnique );

            // convert matrix to string
            mat_to_string( tOrdersUnique, tString );

            // update entry in parameter list
            mParameterList.set( "bspline_orders", tString );
        }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
