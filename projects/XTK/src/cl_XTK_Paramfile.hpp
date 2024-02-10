/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Paramfile.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_PARAMFILE_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_PARAMFILE_HPP_

#include "cl_XML_Parser.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_XTK_Enums.hpp"

#include <string>
#include <sstream>
#include <iostream>

using namespace moris;

namespace moris::xtk
{

    class XTK_Problem_Params
    {
      public:
        XTK_Problem_Params()
                : mInputMeshFile( "" )
                , mMeshType( mtk::MeshType::UNDEFINED )
                , mRealGeomParams( 0 )
                , mIntGeomParams( 0 )
                , mSubdivisionMethods( 0 )
                , mComputeSens( false )
                , mUnzip( false )
                , mEnrich( false )
                , mGhost( false )
                , mExport( false )
                , mOutputMeshFile( "" )
                , mWriteobj( false ){};

        // input mesh
        std::string   mInputMeshFile;
        mtk::MeshType mMeshType;

        // geometry
        std::string           mGeometryName;
        Vector< real >        mRealGeomParams;
        Vector< std::string > mRealGeomLabels;
        Vector< moris_index > mIntGeomParams;
        Vector< std::string > mIntGeomLabels;

        // decomposition
        Vector< enum Subdivision_Method > mSubdivisionMethods;
        Vector< std::string >             mSubdivisionStrings;

        // compute sens
        bool mComputeSens;

        // unzip
        bool mUnzip;

        // enrich
        bool mEnrich;

        // ghost
        bool mGhost;

        // export mesh file
        bool        mExport;
        std::string mOutputMeshFile;

        // write obj file
        bool        mWriteobj;
        moris_index mPhaseForobj;
        std::string mobjOutputFile;

        // dump data to hdf5
        bool        mOutputData;
        std::string mDataFile; /*HDF5*/

        // print geometry parameters
        void print_geom()
        {
            // TODO: support discrete with a field name
            std::cout << "Geometry:  " << mGeometryName << std::endl;
            for ( moris::uint i = 0; i < mRealGeomParams.size(); i++ )
            {
                std::cout << "    " << mRealGeomLabels( i ) << " = " << mRealGeomParams( i ) << std::endl;
            }
        }

        void print_operations()
        {
            std::cout << "Operations to Perform:" << std::endl;

            for ( moris::uint i = 0; i < mSubdivisionStrings.size(); i++ )
            {
                std::cout << "    " << mSubdivisionStrings( i ) << std::endl;
            }

            if ( mComputeSens )
            {
                std::cout << "    Sensitivity Computation" << std::endl;
            }

            if ( mUnzip )
            {
                std::cout << "    Unzip Interface" << std::endl;
            }

            if ( mEnrich )
            {
                std::cout << "    Basis Enrichment" << std::endl;
            }

            if ( mGhost )
            {
                std::cout << "    Ghost Penalization" << std::endl;
            }
        }
    };

    class Paramfile
    {
      public:
        Paramfile( const std::string &aPath );

        ~Paramfile();

        Vector< XTK_Problem_Params > &
        get_xtk_problem_params()
        {
            return mXTKProblems;
        }

      private:
        // the xml parser
        XML_Parser *mParser = nullptr;

        // XTK problems to run
        Vector< XTK_Problem_Params > mXTKProblems;

        /*!
         * Load mesh parameters from XML file
         */
        void
        load_xtk_problems();

        /*!
         * Parse the problems mesh
         */
        void
        parse_xtk_problem_input_mesh( moris::uint aProblemIndex );

        /*!
         * Parse the problems geometry
         */
        void
        parse_xtk_problem_geometry( moris::uint aProblemIndex );

        /*!
         * Parse the problems geometry
         */
        void
        parse_xtk_problem_decomp( moris::uint aProblemIndex );

        /*!
         * Parse the problems geometry
         */
        void
        parse_xtk_problem_operators( moris::uint aProblemIndex );

        /*!
         * Parse the problems geometry
         */
        void
        parse_xtk_problem_output( moris::uint aProblemIndex );

        /*!
         * Parse the obj output information
         */
        void
        parse_xtk_problem_obj( moris::uint aProblemIndex );

        mtk::MeshType
        get_mesh_type_enum( std::string const &aMeshStr )
        {
            if ( aMeshStr == "STK" )
            {
                return mtk::MeshType::STK;
            }

            else
            {
                MORIS_ERROR( 0, "Mesh str not recognized." );
                return mtk::MeshType::STK;
            }
        }

        enum Subdivision_Method
        get_decomp_enum( std::string const &aDecompStr )
        {
            if ( aDecompStr == "Hex8 Regular Subdivision" )
            {
                return Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;
            }
            else if ( aDecompStr == "Tet4 Node Hierarchy" )
            {
                return Subdivision_Method::C_HIERARCHY_TET4;
            }
            else
            {
                MORIS_ERROR( 0, "Decomposition str not recognized: %s. Please use Hex8 Regular Subdivision or Tet4 Node Hierarchy.", aDecompStr.c_str() );
                return Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;
            }
        }

        Vector< real >
        convert_str_to_cell_real( std::string const &aStr )
        {
            std::stringstream ss( aStr );
            Vector< real >    result;

            while ( ss.good() )
            {
                std::string substr;
                std::getline( ss, substr, ',' );
                result.push_back( std::stod( substr ) );
            }

            return result;
        }

        Vector< std::string >
        convert_str_to_cell_str( std::string const &aStr )
        {
            std::stringstream     ss( aStr );
            Vector< std::string > result;

            while ( ss.good() )
            {
                std::string substr;
                std::getline( ss, substr, ',' );
                result.push_back( substr );
            }

            return result;
        }
    };
}    // namespace moris::xtk
#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_PARAMFILE_HPP_ */
