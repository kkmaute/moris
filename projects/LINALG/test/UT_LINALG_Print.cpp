/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include <catch.hpp>
#include "fn_print.hpp"
#include <cstdio>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32
#include <io.h>
#define MORIS_TEST_DUP _dup
#define MORIS_TEST_DUP2 _dup2
#define MORIS_TEST_CLOSE _close
#define MORIS_TEST_FILENO _fileno
#else
#include <unistd.h>
#define MORIS_TEST_DUP dup
#define MORIS_TEST_DUP2 dup2
#define MORIS_TEST_CLOSE close
#define MORIS_TEST_FILENO fileno
#endif

namespace moris
{
    namespace
    {
        template< typename Action >
        std::string capture_cout( Action aAction )
        {
            std::ostringstream tOutput;
            std::ios           tFormatting( nullptr );
            tFormatting.copyfmt( std::cout );
            auto* tOriginal = std::cout.rdbuf( tOutput.rdbuf() );
            try
            {
                aAction();
            } catch ( ... )
            {
                std::cout.rdbuf( tOriginal );
                std::cout.copyfmt( tFormatting );
                throw;
            }
            std::cout.rdbuf( tOriginal );
            std::cout.copyfmt( tFormatting );
            return tOutput.str();
        }

        // The specialized real-matrix printer uses fprintf rather than cout.
        // Restore the descriptor before assertions so failures remain visible.
        class Stdout_Capture
        {
            FILE* mFile       = std::tmpfile();
            int   mDescriptor = -1;

          public:
            Stdout_Capture()
            {
                if ( !mFile )
                {
                    throw std::runtime_error( "Cannot create stdout capture file" );
                }
                std::fflush( stdout );
                mDescriptor = MORIS_TEST_DUP( MORIS_TEST_FILENO( stdout ) );
                if ( mDescriptor < 0 || MORIS_TEST_DUP2( MORIS_TEST_FILENO( mFile ), MORIS_TEST_FILENO( stdout ) ) < 0 )
                {
                    if ( mDescriptor >= 0 )
                    {
                        MORIS_TEST_CLOSE( mDescriptor );
                    }
                    std::fclose( mFile );
                    throw std::runtime_error( "Cannot redirect stdout" );
                }
            }
            ~Stdout_Capture()
            {
                std::fflush( stdout );
                MORIS_TEST_DUP2( mDescriptor, MORIS_TEST_FILENO( stdout ) );
                MORIS_TEST_CLOSE( mDescriptor );
                std::fclose( mFile );
            }
            std::string str()
            {
                std::fflush( stdout );
                std::rewind( mFile );
                std::string tResult;
                char        tBuffer[ 512 ];
                size_t      tCount = 0;
                while ( ( tCount = std::fread( tBuffer, 1, sizeof( tBuffer ), mFile ) ) > 0 )
                {
                    tResult.append( tBuffer, tCount );
                }
                return tResult;
            }
            Stdout_Capture( const Stdout_Capture& )            = delete;
            Stdout_Capture& operator=( const Stdout_Capture& ) = delete;
        };
    }    // namespace

    TEST_CASE( "LINALG print formats include values and dimensions", "[linalgebra],[print]" )
    {
        const Matrix< DDUMat > tA     = { { 1, 2 }, { 3, 4 } };
        const auto             tPlain = capture_cout( [ & ]() { print( tA, "plain" ); } );
        CHECK( tPlain.find( "plain" ) != std::string::npos );
        CHECK( tPlain.find( "Num Rows: 2 | Num Cols: 2" ) != std::string::npos );
        CHECK( tPlain.find( "4" ) != std::string::npos );
        const auto tFancy = capture_cout( [ & ]() { print_fancy( tA, "fancy" ); } );
        CHECK( tFancy.find( "fancy" ) != std::string::npos );
        CHECK( tFancy.find( "1 |" ) != std::string::npos );
        const auto tList = capture_cout( [ & ]() { print_std_initializer_list( tA, "list" ); } );
        CHECK( tList.find( "{{1, 2}, \n{3, 4}}" ) != std::string::npos );
        CHECK( tList.find( "list" ) != std::string::npos );
        const Matrix< DDUMat > tRow     = { { 5, 6 } };
        const auto             tRowText = capture_cout( [ & ]() { print_as_row_vector( tRow, "row" ); } );
        CHECK( tRowText.find( "row = [" ) != std::string::npos );
        CHECK( tRowText.find( "5,\n6" ) != std::string::npos );
        const auto tMacros = capture_cout( [ & ]() { PRINT( tA ); PRINTFANCY( tA ); PRINTLIST( tA ); PRINTROW( tRow ); } );
        CHECK( tMacros.find( "tA" ) != std::string::npos );
        CHECK( tMacros.find( "tRow" ) != std::string::npos );
        CHECK_THROWS( capture_cout( [ & ]() { print_as_row_vector( tA, "invalid" ); } ) );
        const auto tEmpty = capture_cout( [ & ]() { print_std_initializer_list( Matrix< DDUMat >(), "empty" ); } );
        CHECK( tEmpty.find( "{}" ) != std::string::npos );
    }

    TEST_CASE( "LINALG real printer supports named and default titles", "[linalgebra],[print]" )
    {
        const Matrix< DDRMat > tA = { { 1.25, -2.0 } };
        std::string            tNamed, tDefault;
        {
            Stdout_Capture tCapture;
            print( tA, "real_values" );
            tNamed = tCapture.str();
        }
        {
            Stdout_Capture tCapture;
            print( tA, "" );
            tDefault = tCapture.str();
        }
        CHECK( tNamed.find( "real_values = [" ) != std::string::npos );
        CHECK( tNamed.find( "+1.250000000000000e+00" ) != std::string::npos );
        CHECK( tNamed.find( "-2.000000000000000e+00" ) != std::string::npos );
        CHECK( tDefault.find( "morisMat = [" ) != std::string::npos );
    }

    TEST_CASE( "LINALG printers traverse matrices and borrowed matrix pointers", "[linalgebra],[print]" )
    {
        Matrix< DDUMat >                  tA        = { { 7 } };
        const Vector< Matrix< DDUMat > >  tMatrices = { tA, tA };
        const Vector< Matrix< DDUMat >* > tPointers = { &tA };
        const auto                        tValues   = capture_cout( [ & ]() { print( tMatrices ); } );
        CHECK( tValues.find( "Number of entries = 2" ) != std::string::npos );
        CHECK( tValues.find( "Cell 1" ) != std::string::npos );
        const auto tBorrowed = capture_cout( [ & ]() { print( tPointers, "borrowed" ); } );
        CHECK( tBorrowed.find( "Cell Name: borrowed" ) != std::string::npos );
        CHECK( tBorrowed.find( "Number of entries = 1" ) != std::string::npos );
        CHECK( tBorrowed.find( "7" ) != std::string::npos );
    }
}    // namespace moris
