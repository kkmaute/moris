/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_File.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_FILE_HPP_
#define SRC_HMR_CL_HMR_FILE_HPP_

#include <string>

#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
// HD5 c-interface
#include "hdf5.h"

#include "moris_typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "linalg_typedefs.hpp"
#include "fn_sum.hpp" //LINALG/src

namespace moris::hmr
{

//-------------------------------------------------------------------------------
    class File
    {

        //! file ID for HDF5
        hid_t  mFileID;

        //! error handler
        herr_t mStatus;

//-------------------------------------------------------------------------------
    public :
//-------------------------------------------------------------------------------

        // default constructor
        File() = default;

//-------------------------------------------------------------------------------
        // default destructor
        ~File() = default;

//-------------------------------------------------------------------------------

        /** creates a new HDF5 file
         *
         * @param[ in ] aPath   path to file
         */
        void create( const std::string & aPath );

//-------------------------------------------------------------------------------

        /**
         * opens an existing HDF5 file
         *
         * @param[ in ] aPath   path to file
         */
        void open( const std::string & aPath );

//-------------------------------------------------------------------------------

        /**
         * stores the contents of aParameters into the file
         *
         * @param[ in ] aParameters   settings object
         */
        void save_settings( const Parameters * aParameters );

//-------------------------------------------------------------------------------

        /**
         * reads the settings from the file and writes them into the object
         */
        void load_settings( Parameters * aParameters );

//-------------------------------------------------------------------------------

        /**
         * stores the refinement pattern of the lagrange mesh and all related BSpline meshes
         * of the current proc into the file
         *
         * @param[ in ]  aMesh     pointer to Lagrange Mesh
         */
        void save_refinement_pattern( Lagrange_Mesh_Base * aLagrangeMesh );

        //-------------------------------------------------------------------------------

        void save_refinement_pattern(
                Background_Mesh_Base          * aBackgroundMesh,
                const moris::Matrix< DDUMat > & tPatternToSave );

        //-------------------------------------------------------------------------------

        static void save_refinement_pattern(
                Background_Mesh_Base             * aBackgroundMesh,
                const moris::Matrix< DDUMat >    & tPatternToSave,
                Matrix< DDLUMat >                & aElementCounterPerLevelAndPattern,
                Vector< Matrix< DDLUMat > > & aElementPerPattern);

//-------------------------------------------------------------------------------

        /**
         * laods the refinement pattern of the current proc
         * and initializes a new mesh object
         *
         * @param[ inout ] aMesh           aMesh   pointer to background mesh
         * @param[ in ]    aMode           false: input, true: output
         *
         */
        void load_refinement_pattern( Background_Mesh_Base * aMesh );

//-------------------------------------------------------------------------------

        /**
         * closes the interface to the file
         */
        void close();

    };

//-------------------------------------------------------------------------------

    /**
     * free function needed by loading constructor
     */
    Parameters * create_hmr_parameters_from_hdf5_file( const std::string & aPath );

//-------------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_FILE_HPP_ */

