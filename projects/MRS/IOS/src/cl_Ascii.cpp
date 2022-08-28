/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Ascii.cpp
 *
 */

#include <fstream>
#include <iostream>

#include "cl_Ascii.hpp"
#include "assert.hpp"

namespace moris
{

    //------------------------------------------------------------------------------

    bool
    file_exists( const std::string & aPath )
    {
        // test if file exists
        std::ifstream tFile( aPath );

        // save result into output variable
        bool aFileExists;

        if( tFile )
        {
            // close file
            tFile.close();
            aFileExists = true;
        }
        else
        {
            aFileExists = false;
        }

        return aFileExists;
    }

    //------------------------------------------------------------------------------

    Ascii::Ascii( const std::string & aPath, const FileMode & aMode ) :
                mMode( aMode )
    {
        // test if path is absolute
        if( aPath.substr( 0,1 ) == "/" )
        {
            mPath = aPath;
        }
        // test if path is relative
        else if( aPath.substr( 0,2 ) == "./" )
        {
            mPath = aPath;
        }
        else
        {
            MORIS_ERROR( false, " correct Ascii file path");
            //mPath = std::sprint( "%s/%s", std::getenv( "PWD" ), aPath.c_str() );
        }

        switch ( aMode )
        {
            case( FileMode::OPEN_RDONLY ) :
            {
                this->load_buffer();
                break;
            }
            case( FileMode::NEW ) :
            {
                mBuffer.clear();
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Unknown file mode for ASCII file" );
            }
        }
    }

    //------------------------------------------------------------------------------

    Ascii::~Ascii()
    {
        MORIS_ERROR( ! mChangedSinceLastSave,
                "File %s was changed but never saved.",
                mPath.c_str() );

        mBuffer.clear();
    }
    //------------------------------------------------------------------------------

    bool Ascii::save()
    {
        MORIS_ERROR( mMode != FileMode::OPEN_RDONLY,
                "File %s can't be saved since it is opened in write protected mode.",
                mPath.c_str() );

        // open file
        std::ofstream tFile( mPath.c_str(),  std::ofstream::trunc );

        if( tFile )
        {
            // save buffer to file
            for( std::string & tLine : mBuffer )
            {
                tFile << tLine << std::endl;
            }

            tFile.close();
        }
        else
        {
            MORIS_ERROR( false,
                    "Something went wrong while trying to save %s.",
                    mPath.c_str() );
        }

        mChangedSinceLastSave = false;

        return mChangedSinceLastSave;
    }

    //------------------------------------------------------------------------------

    moris::uint Ascii::length() const
    {
        return mBuffer.size();
    }

    //------------------------------------------------------------------------------

    std::string & Ascii::line( const moris::uint aLineNumber )
    {
        return mBuffer( aLineNumber );
    }

    //------------------------------------------------------------------------------

    const std::string & Ascii::line( const moris::uint aLineNumber ) const
    {
        return mBuffer( aLineNumber );
    }

    //------------------------------------------------------------------------------

    void Ascii::print( const std::string & aLine )
    {
        mBuffer.push_back( aLine );

        mChangedSinceLastSave = true;
    }

    //------------------------------------------------------------------------------

    void Ascii::load_buffer()
    {
        // tidy up buffer
        mBuffer.clear();

        // make sure that file exists
        MORIS_ERROR( file_exists( mPath ),
                "File %s does not exist.",
                mPath.c_str() );

        // open file
        std::ifstream tFile( mPath );

        // test if file can be opened
        if( tFile )
        {
            // temporary container for string
            std::string tLine;

            while ( std::getline( tFile, tLine ) )
            {
                mBuffer.push_back( tLine );
            }

            // close file
            tFile.close();
        }
        else
        {
            MORIS_ERROR( false, "Someting went wrong while opening file\n %s",
                    mPath.c_str() );
        }
    }

    //------------------------------------------------------------------------------

}

