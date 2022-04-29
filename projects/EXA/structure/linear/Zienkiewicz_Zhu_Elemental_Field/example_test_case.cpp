/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example specific interface to moris
 *
 */
#include <catch.hpp>

#include "paths.hpp"

#include "cl_Logger.hpp"    // MRS/IOS/src
#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

TEST_CASE( "Zienkiewicz_Zhu_Elemental_Field",
        "[moris],[example],[Zienkiewicz_Zhu_Elemental_Field]" )
{
    if ( par_size() == 1 )
    {
        // define command line call
        int argc = 2;

        char tString1[] = "";
        char tString2[] = "Zienkiewicz_Zhu_Elemental_Field.so";

        char* argv[ 2 ] = { tString1, tString2 };

        // call to performance manager main interface
        int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

        // catch test statements should follow
        REQUIRE( tRet == 0 );

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tFieldRefPath = tPrefix + "/projects/EXA/structure/linear/Zienkiewicz_Zhu_Elemental_Field/Field_Zienkiewicz_Zhu_Ref.hdf5";
        std::string tLabel        = "FieldZienkiewiczZhu";

        Matrix< DDRMat > tRefValues;
        Matrix< DDRMat > tValues;

        hid_t  tFileRef = open_hdf5_file( tFieldRefPath );
        herr_t tStatus  = 0;
        load_matrix_from_hdf5_file(
                tFileRef,
                tLabel,
                tRefValues,
                tStatus );

        tStatus = close_hdf5_file( tFileRef );

        std::string tFieldPath = "./Field_Zienkiewicz_Zhu.hdf5";

        hid_t tFile = open_hdf5_file( tFieldPath );
        tStatus     = 0;
        load_matrix_from_hdf5_file(
                tFile,
                tLabel,
                tValues,
                tStatus );

        tStatus = close_hdf5_file( tFile );

        MORIS_ERROR( tStatus == 0, "zienkiewicz_zhu_elemental: Status returned != 0, Error in reading values" );

        CHECK( tRefValues.numel() == tValues.numel() );

        for ( uint Ik = 0; Ik < tRefValues.n_cols(); Ik++ )
        {
            for ( uint Ii = 0; Ii < tRefValues.n_rows(); Ii++ )
            {
                CHECK( std::abs( tRefValues( Ii, Ik ) - tValues( Ii, Ik ) ) < 1e-9 );
            }
        }
    }
}
