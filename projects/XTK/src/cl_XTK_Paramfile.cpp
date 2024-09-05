/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Paramfile.cpp
 *
 */

#include "cl_XTK_Paramfile.hpp"

namespace moris::xtk
{
    Paramfile::Paramfile( const std::string &aPath )
    {
        // create parser
        mParser = new XML_Parser( aPath );

        // load the parsed mesh information
        this->load_xtk_problems();
    }

    Paramfile::~Paramfile()
    {
        delete mParser;
    }

    void
    Paramfile::load_xtk_problems()
    {
        // count number of problems
        uint tNumberXTKProblems = mParser->count_keys_in_subtree( "moris", "xtk_problem" );

        // make sure that number of problems is at least one
        MORIS_ERROR( tNumberXTKProblems > 0,
                "Paramfile::load_xtk_problems: no xtk problems defined in XML file" );

        mXTKProblems.resize( tNumberXTKProblems );

        // set up xtk problems
        for ( uint m = 0; m < tNumberXTKProblems; ++m )
        {

            // set up mesh parameters
            this->parse_xtk_problem_input_mesh( m );

            // iterate through geometry and setup parameters
            this->parse_xtk_problem_geometry( m );

            // parse decomposition methods
            this->parse_xtk_problem_decomp( m );

            // parse operators
            this->parse_xtk_problem_operators( m );

            // parse output file
            this->parse_xtk_problem_output( m );

            // parse stl
            this->parse_xtk_problem_obj( m );
        }
    }

    void
    Paramfile::parse_xtk_problem_input_mesh( moris::uint aProblemIndex )
    {
        // first stores the keys, second stores the vals
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

        // here I need to be only looking at the first xtk_problem
        mParser->get_keys_from_subtree( "moris.xtk_problem", "mesh_input", 0, tFirst, tSecond );

        MORIS_ERROR( tFirst.size() != 0, "No input mesh parsed in XML file using element name mesh_input" );
        // iterate through keys related to the mesh input
        for ( uint k = 0; k < tFirst.size(); k++ )
        {
            std::string &tKey = tFirst( k );
            if ( tKey == "path" )
            {
                mXTKProblems( aProblemIndex ).mInputMeshFile = tSecond( k );
            }
            else if ( tKey == "field_type" )
            {
                mXTKProblems( aProblemIndex ).mMeshType = get_mesh_type_enum( tSecond( k ) );
            }
            else
            {
                MORIS_ERROR( 0, "Invalid mesh parameter parsed: %s", tFirst( k ).c_str() );
            }
        }
    }

    void
    Paramfile::parse_xtk_problem_geometry( moris::uint aProblemIndex )
    {
        // first stores the keys, second stores the vals
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

        // here I need to be only looking at the first xtk_problem
        mParser->get_keys_from_subtree( "moris.xtk_problem", "geometry", 0, tFirst, tSecond );

        MORIS_ERROR( tFirst.size() != 0, "No input mesh parsed in XML file. Make sure to use geometry" );

        // iterate through keys related to the geometry
        for ( uint k = 0; k < tFirst.size(); k++ )
        {
            std::string &tKey = tFirst( k );
            if ( tKey == "geometry_type" )
            {
                mXTKProblems( aProblemIndex ).mGeometryName = tSecond( k );
            }
            else if ( tKey == "real_parameters" )
            {
                mXTKProblems( aProblemIndex ).mRealGeomParams = convert_str_to_cell_real( tSecond( k ) );
            }

            else if ( tKey == "real_labels" )
            {
                mXTKProblems( aProblemIndex ).mRealGeomLabels = convert_str_to_cell_str( tSecond( k ) );
            }
            else if ( tKey == "int_parameters" )
            {
                MORIS_ERROR( 0, "Integer parameters not hooked up" );
            }
            else if ( tKey == "int_labels" )
            {
                MORIS_ERROR( 0, "Integer parameters not hooked up" );
            }
            else
            {
                MORIS_ERROR( 0, "Unrecognized geometry parameter parsed: %s. ", tKey.c_str() );
            }
        }
    }

    void
    Paramfile::parse_xtk_problem_decomp( moris::uint aProblemIndex )
    {
        // first stores the keys, second stores the vals
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

        // here I need to be only looking at the first xtk_problem
        mParser->get_keys_from_subtree( "moris.xtk_problem", "refine", 0, tFirst, tSecond );

        MORIS_ERROR( tFirst.size() != 0, "No refinement parsed in XML file. Make sure to use refine" );
        MORIS_ERROR( tFirst.size() == 1 || tFirst.size() == 2, "There are only two possible refinements possible, more than two have been passed in" );

        // iterate through keys related to the geometry
        for ( uint k = 0; k < tFirst.size(); k++ )
        {
            std::string &tKey = tFirst( k );

            if ( tKey == "subdivision_method" )
            {
                if ( tSecond( k ) == "Hex8 Regular Subdivision" )
                {
                    MORIS_ERROR( k == 0, " Hex8 Regular Subdivision should show up first in XML file" );
                    mXTKProblems( aProblemIndex ).mSubdivisionMethods.push_back( get_decomp_enum( tSecond( k ) ) );
                    mXTKProblems( aProblemIndex ).mSubdivisionStrings.push_back( tSecond( k ) );
                }

                else if ( tSecond( k ) == "Tet4 Node Hierarchy" )
                {
                    mXTKProblems( aProblemIndex ).mSubdivisionMethods.push_back( get_decomp_enum( tSecond( k ) ) );
                    mXTKProblems( aProblemIndex ).mSubdivisionStrings.push_back( tSecond( k ) );
                }
                else
                {
                    MORIS_ERROR( 0, "Invalid decomposition parsed: %s", tKey.c_str() );
                }
            }
            else
            {
                MORIS_ERROR( 0, "Invalid key in refinement" );
            }
        }
    }

    void
    Paramfile::parse_xtk_problem_operators( moris::uint aProblemIndex )
    {
        // first stores the keys, second stores the vals
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

        // here I need to be only looking at the first xtk_problem
        mParser->get_keys_from_subtree( "moris.xtk_problem", "operators", 0, tFirst, tSecond );

        // iterate through keys related to the geometry
        for ( uint k = 0; k < tFirst.size(); k++ )
        {
            std::string &tKey = tFirst( k );

            if ( tKey == "compute_sens" )
            {
                mXTKProblems( aProblemIndex ).mComputeSens = mParser->to_bool( tSecond( k ) );
            }
            else if ( tKey == "unzip" )
            {
                mXTKProblems( aProblemIndex ).mUnzip = mParser->to_bool( tSecond( k ) );
            }
            else if ( tKey == "enrich" )
            {
                mXTKProblems( aProblemIndex ).mEnrich = mParser->to_bool( tSecond( k ) );
            }
            else if ( tKey == "ghost_stab" )
            {
                mXTKProblems( aProblemIndex ).mGhost = mParser->to_bool( tSecond( k ) );
            }
            else
            {
                MORIS_ERROR( 0, "Invalid key in operators" );
            }
        }
    }

    void
    Paramfile::parse_xtk_problem_output( moris::uint aProblemIndex )
    {
        // first stores the keys, second stores the vals
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

        // here I need to be only looking at the first xtk_problem
        mParser->get_keys_from_subtree( "moris.xtk_problem", "output", 0, tFirst, tSecond );

        // iterate through keys related to the geometry
        for ( uint k = 0; k < tFirst.size(); k++ )
        {
            std::string &tKey = tFirst( k );

            if ( tKey == "output_mesh" )
            {
                mXTKProblems( aProblemIndex ).mExport         = true;
                mXTKProblems( aProblemIndex ).mOutputMeshFile = tSecond( k );
            }
            else if ( tKey == "output_file" )
            {
                mXTKProblems( aProblemIndex ).mOutputData = true;
                mXTKProblems( aProblemIndex ).mDataFile   = tSecond( k );
            }
            else
            {
                MORIS_ERROR( 0, "Parsing error: output parameter not recognized: %s", tKey.c_str() );
            }
        }
    }

    void
    Paramfile::parse_xtk_problem_obj( moris::uint aProblemIndex )
    {
        // first stores the keys, second stores the vals
        Vector< std::string > tFirst;
        Vector< std::string > tSecond;

        // here I need to be only looking at the first xtk_problem
        mParser->get_keys_from_subtree( "moris.xtk_problem", "obj", 0, tFirst, tSecond );

        // iterate through keys related to the geometry
        for ( uint k = 0; k < tFirst.size(); k++ )
        {
            std::string &tKey = tFirst( k );

            std::cout << "key = " << tKey << '\n';

            if ( tKey == "dump_obj" )
            {
                mXTKProblems( aProblemIndex ).mWriteobj = mParser->to_bool( tSecond( k ) );
            }
            else if ( tKey == "phase" )
            {
                mXTKProblems( aProblemIndex ).mPhaseForobj = 0;
            }

            else if ( tKey == "obj_file" )
            {
                std::cout << "here" << '\n';
                mXTKProblems( aProblemIndex ).mobjOutputFile = tSecond( k );
            }

            else
            {
                MORIS_ERROR( 0, "Parsing error: obj parameter not recognized: %s", tKey.c_str() );
            }
        }
    }

}    // namespace moris::xtk
