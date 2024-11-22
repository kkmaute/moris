/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Paramfile.cpp
 *
 */

#include "cl_HMR_Paramfile.hpp"

#include "cl_HMR_Parameters.hpp"
#include "HMR_Tools.hpp"
#include "assert.hpp"

#include "fn_Parsing_Tools.hpp"
#include "cl_XML_Parser.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_unique.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris::hmr
{
// -----------------------------------------------------------------------------

Paramfile::Paramfile( const std::string & aPath,
                      const enum State aState ) : mState( aState ) , mParameterList( "" )
{
    // create parser
    mParser = new XML_Parser( aPath );

    // load meshes
    this->load_mesh_params();

    // load fields
    this->load_field_params();

    // load parameters
    this->load_parameter_list();

    // state dependent parameters
    this->load_state_params();

    // make parameters consistent with shipped fields and meshes
    this->update_parameter_list();

    switch( aState )
    {
        case ( State::REFINE_MESH  ) :
        {
            // load user parameters for refinement
            this->load_user_refinement_parameters();
        }
        default :
        {
            // do nothing
        }
    }
}

// -----------------------------------------------------------------------------

Paramfile::~Paramfile()
{
    delete mParser;
}
// -----------------------------------------------------------------------------

Parameter_List & Paramfile::get_parameter_list()
{
    return mParameterList;
}

// -----------------------------------------------------------------------------

void Paramfile::load_mesh_params()
{
    // count meshes
    uint tNumberOfMeshes
        =  mParser->count_keys_in_subtree( "moris.hmr", "mesh" );

    // make sure that number of meshes is at least one
    MORIS_ERROR( tNumberOfMeshes>0,
                          "Paramfile::load_mesh_params: at least one Lagrange meshes needs to be defined." );

    // allocate cell
    mMeshParams.resize( tNumberOfMeshes, Mesh_Param() );

    // clear map
    mMeshMap.clear();

    // list with IDs
    Matrix< IdMat > tIDs( tNumberOfMeshes, 1 );

    for( uint m=0; m<tNumberOfMeshes; ++m )
    {
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;
        mParser->get_keys_from_subtree( "moris.hmr", "mesh", m, tFirst, tSecond );

        Mesh_Param & tMesh = mMeshParams( m );

        // copy key to settings struct
        for( uint k=0; k<tFirst.size(); ++k )
        {
            std::string & tKey = tFirst( k );

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

    void Paramfile::load_field_params()
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
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;
            mParser->get_keys_from_subtree( "moris.hmr", "field", f, tFirst, tSecond );

            Field_Param & tField = mFieldParams( f );

            // copy key to settings struct
            for( uint k=0; k<tFirst.size(); ++k )
            {
                std::string & tKey = tFirst( k );
                if ( tKey == "label" )
                {
                    tField.mLabel = tSecond( k );
                }
                else if ( tKey == "id" )
                {
                    tField.mID = stoi( tSecond( k ) );
                }
                else if ( tKey == "input_lagrange_order" )
                {
                    tField.mInputLagrangeOrder = stoi( tSecond( k ) );
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
                else if ( tKey == "l2_alpha" )
                {
                    tField.mL2alpha = stod( tSecond( k ) );
                }
            }

            // copy ID into array
            tIDs( f ) = tField.mID;

            // add field to map
            mFieldMap[ tField.mID ] = f;

            // check for consistency
            MORIS_ERROR( tField.mID != gNoID,                            "Field ID is not set." );
            MORIS_ERROR( tField.mLabel.size() > 0,                       "Field label must not be empty." );
            MORIS_ERROR( tField.mSource.size() > 0,                      "Field source path must not be empty." );
            MORIS_ERROR( tField.mInputBSplineOrder <= gMaxBSplineOrder,  "BSpline input order of field too big." );
            MORIS_ERROR( tField.mInputBSplineOrder > 0,                  "BSpline input order needs to be larger than 0." );
            MORIS_ERROR( tField.mOutputBSplineOrder <= gMaxBSplineOrder, "BSpline order of field too big." );
            MORIS_ERROR( tField.mOutputBSplineOrder >  0,                "BSpline output order needs to be larger than 0." );
        }

        // make IDs unique
        Matrix< IdMat > tUniqueIDs;
        unique( tIDs, tUniqueIDs );

        // make sure that each ID exists only once
        MORIS_ERROR( tUniqueIDs.length() == tNumberOfFields,
                "Duplicate field IDs in input file" );
    }

// -----------------------------------------------------------------------------

    void Paramfile::load_state_params()
    {
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

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
            std::string & tKey = tFirst( k );
            if ( tKey == "input_database" )
            {
                mInputDatabase = tSecond( k );
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
                string_to_matrix( tSecond( k ), mMeshIDs );
            }
            else if( tKey == "fields" )
            {
                string_to_matrix( tSecond( k ), mFieldIDs );
            }
            else if ( tKey == "library" )
            {
                mLibraryPath = tSecond( k );
            }
            else if ( tKey == "function" )
            {
                mUserFunction = tSecond( k );
            }
//                else if ( tKey == "initial_refinement" )
//                {
//                    mInitialBSplineRefinement = stoi( tSecond( k ) );
//                }
            else if ( tKey == "additional_lagrange_refinement" )
            {
                mAdditionalLagrangeRefinement = stoi( tSecond( k ) );
            }
            else if ( tKey == "union_mesh" )
            {
                mUnionMesh = tSecond( k );
            }
        }
    }

// -----------------------------------------------------------------------------

    void Paramfile::load_parameter_list()
    {
        // call parameter list initializer from cl_Parameters.hpp
        mParameterList = prm::create_hmr_parameter_list();
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;
        mParser->get_keys_from_subtree( "moris.hmr", "parameters", 0, tFirst, tSecond );

        for( uint k=0; k<tFirst.size(); ++k )
        {
            std::string & tKey = tFirst( k );

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
//                else if( tKey == "initial_refinement" )
//                {
//                    mParameterList.set( "initial_refinement", ( sint ) std::stoi( tSecond( k ) ) );
//                }
            else if ( tKey == "severity_level" )
            {
                mParameterList.set( "severity_level", ( sint ) std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "truncate_bsplines" )
            {
                mParameterList.set( "truncate_bsplines", ( sint ) string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "additional_lagrange_refinement" )
            {
                mParameterList.set(  "additional_lagrange_refinement", ( sint ) std::stoi( tSecond( k ) ) );
            }
            else if ( tKey == "use_number_aura" )
            {
                mParameterList.set(  "use_number_aura", ( sint ) string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "use_multigrid" )
            {
                mParameterList.set(  "use_multigrid", ( sint ) string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "use_refinement_interrelation" )
            {
                mParameterList.set(  "use_refinement_interrelation", ( sint ) string_to_bool( tSecond( k ) ) );
            }
            else if ( tKey == "renumber_lagrange_nodes" )
            {
                mParameterList.set(  "renumber_lagrange_nodes", ( sint ) string_to_bool( tSecond( k ) ) );
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
    void Paramfile::update_parameter_list()
    {
        // - - - - - - - - - - - - - - - - - - - - - - - -
        // Part 1: Make Lagrange orders consistent
        // - - - - - - - - - - - - - - - - - - - - - - - -

        Matrix< DDUMat > tOrdersFromXML;
        string_to_matrix( mParameterList.get< std::string >( "lagrange_orders" ),
                tOrdersFromXML );

        // allocate vector
        Matrix< DDUMat > tOrders(
                tOrdersFromXML.length()
                + mMeshIDs.length(), 1 );

        // initialize counter
        uint tCount = 0;

        // copy orders from XML
        for( uint k=0; k<tOrdersFromXML.length(); ++k )
        {
            tOrders( tCount++ ) = tOrdersFromXML( k );
        }

        // add Lagrange orders from meshes
        for( uint k=0; k< mMeshIDs.length(); ++k )
        {
            tOrders( tCount++ ) = mMeshParams( mMeshMap.find( mMeshIDs( k ) ) ).mOrder;
        }

        // make orders unique
        Matrix< DDUMat > tOrdersUnique;
        unique( tOrders, tOrdersUnique );

        // make sure that Lagrange orders are withing bounds
        MORIS_ERROR( tOrdersUnique.min()>0,
                              "Paramfile::update_parameter_list: minimum order of Lagrange meshes needs to be larger 0" );

        // convert matrix to string
        std::string tString;
        matrix_to_string( tOrdersUnique, tString );

        // update entry in parameter list
        mParameterList.set( "lagrange_orders", tString );

        // - - - - - - - - - - - - - - - - - - - - - - - -
        // Part 2: Make B-Spline orders consistent
        // - - - - - - - - - - - - - - - - - - - - - - - -

        string_to_matrix( mParameterList.get< std::string >( "bspline_orders" ),
                tOrdersFromXML );

        // allocate vector
        tOrders.set_size(
                tOrdersFromXML.length()
                + 2*mFieldIDs.length(), 1 );

        // reset counter
        tCount = 0;

        // copy orders from XML
        for( uint k=0; k<tOrdersFromXML.length(); ++k )
        {
            tOrders( tCount++ ) = tOrdersFromXML( k );
        }

        // add B-Spline orders from fields
        for( uint k=0; k< mFieldIDs.length(); ++k )
        {
            tOrders( tCount++ ) = mFieldParams( mFieldMap.find( mFieldIDs( k ) ) ).mInputBSplineOrder;
            tOrders( tCount++ ) = mFieldParams( mFieldMap.find( mFieldIDs( k ) ) ).mOutputBSplineOrder;
        }

        // make orders unique
        unique( tOrders, tOrdersUnique );

        // make sure that Lagrange orders are withing bounds
        MORIS_ERROR( tOrdersUnique.min()>0,
                              "Paramfile::update_parameter_list: minimum order of Bsplines needs to be larger 0" );

        // convert matrix to string
        matrix_to_string( tOrdersUnique, tString );

        // update entry in parameter list
        mParameterList.set( "bspline_orders", tString );

        // update B-SPline and Lagrange refinement levels
//            if( mInitialBSplineRefinement >=0 )
//            {
//                mParameterList.set( "initial_refinement", mInitialBSplineRefinement );
//            }
        if( mAdditionalLagrangeRefinement >= 0 )
        {
            mParameterList.set( "additional_lagrange_refinement", mAdditionalLagrangeRefinement );
        }
    }

// -----------------------------------------------------------------------------

    void Paramfile::load_user_refinement_parameters()
    {
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;
        mParser->get_keys_from_subtree( "moris.hmr", "refine", 0, tFirst, tSecond );

        // reset counter
        for( uint k=0; k<tFirst.size(); ++k )
        {
            std::string & tKey = tFirst( k );

            if( tKey == "library" )
            {
                mLibraryPath = tSecond( k );
            }
            else if ( tKey == "function" )
            {
                mUserFunction = tSecond( k );
            }
            else if ( tKey == "max_refinement_level" )
            {
                // this overwrites the current setting
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
            else
            {
                std::string tType  = tKey.substr( 0, tKey.find_first_of( '_' ) );
                std::string tLabel = tKey.substr( tKey.find_first_of( '_' ) + 1, tKey.size() );
                if( tType == "real" )
                {
                    mParameterList.insert( tLabel, ( moris::real ) stod( tSecond( k ) ) );
                }
                else if( tType == "int" )
                {
                    mParameterList.insert( tLabel, ( moris::sint ) stoi( tSecond( k ) ) );
                }
                else if ( tType == "bool" )
                {
                    mParameterList.insert( tLabel, ( moris::sint ) string_to_bool( tSecond( k )  ) );
                }
                else if ( tType == "string" )
                {
                    mParameterList.insert( tLabel, tSecond( k ) );
                }
            }
        }
    }

// -----------------------------------------------------------------------------
} /* namespace moris */
