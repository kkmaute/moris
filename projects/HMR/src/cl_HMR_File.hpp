/*
 * cl_HMR_File.hpp
 *
 *  Created on: Jun 29, 2018
 *      Author: messe
 */

//https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html

#ifndef SRC_HMR_CL_HMR_FILE_HPP_
#define SRC_HMR_CL_HMR_FILE_HPP_

#include <string>

// HD5 c-interface
#include "hdf5.h"

#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src
#include "fn_sum.hpp" //LNA/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace hmr
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
            File(){}

//-------------------------------------------------------------------------------
            // default destructor
            ~File(){}

//-------------------------------------------------------------------------------

            /** creates a new HDF5 file
             *
             * @param[ in ] aPath   path to file
             */
            void
            create( const std::string & aPath );

//-------------------------------------------------------------------------------

            /**
             * opens an existing HDF5 file
             *
             * @param[ in ] aPath   path to file
             */
            void
            open( const std::string & aPath );

//-------------------------------------------------------------------------------

            /**
             * stores the contents of aParameters into the file
             *
             * @param[ in ] aParameters   settings object
             */

            void
            save_settings( const Parameters * aParameters );

//-------------------------------------------------------------------------------

            /**
             * reads the settings from the file and writes them into the object
             */
            void
            load_settings( Parameters * aParameters );

//-------------------------------------------------------------------------------

            /**
             * stores the refinement pattern of the current proc into the file
             * @param[ in ] aMesh   pointer to background mesh
             */
            void
            save_refinement_pattern( Background_Mesh_Base * aMesh );

//-------------------------------------------------------------------------------

            /**
             * laods the refinement pattern of the current proc
             * and initializes a new mesh object
             *
             * @param[ in ] aMesh              settings object
             *
             * @return Background_Mesh_Base *  pointer to new background mesh
             */
            Background_Mesh_Base *
            load_refinement_pattern( const Parameters       * aParameters );

//-------------------------------------------------------------------------------

            /**
             * closes the interface to the file
             */
            void
            close();

//-------------------------------------------------------------------------------
        private:
//-------------------------------------------------------------------------------

            /**
             * adds the proc number to the filename
             */
            std::string
            parralize_filename( const std::string & aPath );

//-------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------

    } /* namespace hmr */

} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_FILE_HPP_ */
