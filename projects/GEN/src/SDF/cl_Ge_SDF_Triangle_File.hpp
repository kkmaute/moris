//
// Created by messe on 1/26/18.
//

#ifndef MORIS_CL_GETRIANGLEFILE_HPP
#define MORIS_CL_GETRIANGLEFILE_HPP

#include <fstream>
#include <iostream>

#include "cl_Cell.hpp" // CON/src
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"
#include "GeUtilities.hpp"
/*! The TriangleFile contains a Matrix of node coordinates
*   and a Matrix cof the node numbers for each triangle. */
namespace ge {
    class SDF_Triangle_File
    {
        struct Settings{
            moris::real mMeshHighPass = 1e-9; /**< If the absolute of an x-,y-, or z-coordinate is less than
                                                *   this value, it is treated as zero.*/

            moris::uint mQuadSplitMode = 1; /**< If a GMSH-File is read, eventual quadrangles
                                          * are automatically split into triangles.
                                           * Valid values for mQuadSplitMode are:
                                           *
                                           * 1: split quadrangle into two triangles from upper left to lower right
                                           * 2: split quadrangle into two triangles  from upper right to lower left
                                           * 4: split quadrangle into four triangles. An adittional node is created per each quadrangle.
                                           *
                                           */

            moris::bool_t mSaveVTKasASCII = true; /**< Switch that tells if VTK is saved
                                                  * as ASCII or binary */

        };

        Settings mSettings;

        moris::uint mNumberOfTriangles;
        moris::uint mNumberOfNodes;

        moris::Matrix< moris::DDRMat > mNodeCoords;
        moris::Matrix< moris::DDUMat > mTriangles;

// =============================================================================
    public:
// =============================================================================

        /**
        * @brief Loads the triangle from a given file.
        *
        * @param[in] aFilePath       - Path to the ascii file containing the mesh.
        *                              Must be of type *.obj or *.msh.
        *                              If a *.msh file is read, quadrangles are
        *                              automatically split into triangles.
        *
        */
        SDF_Triangle_File(const std::string& aFilePath);
// -----------------------------------------------------------------------------

        ~SDF_Triangle_File() = default;

// -----------------------------------------------------------------------------

        /**
         * @brief Saves the node coordinates and triangles to a file.
         *
         * @param[in] aFilePath  destination path for the file that is to be saved.
         *
         */
        void save_to_file(const std::string& aFilePath);

// -----------------------------------------------------------------------------
        /**
       * @brief returns a moris::Mat containing the node coordinates.
       *
       */
        moris::Matrix< moris::DDRMat >
        const& get_node_coords() const
        {
            return mNodeCoords;
        }

// -----------------------------------------------------------------------------

        /**
        * @brief returns a moris::Mat containing the triangle nodes.
        *
        */
        moris::Matrix< moris::DDUMat >
        const & get_triangles()  const
        {
            return mTriangles;
        }
// -----------------------------------------------------------------------------
        /**
        * @brief returns the number of triangles
        *
        */
        moris::uint
        get_number_of_triangles() const
        {
            return  mNumberOfTriangles;
        }
// =============================================================================
    private:
// =============================================================================
        /**
        * @brief private subroutine for loading file into buffer
        *
        */
        void load_ascii_to_buffer(const std::string& aFilePath,
                                  moris::Cell<std::string>& aBuffer);

// -----------------------------------------------------------------------------

        /**
          * @brief private subroutine for loading an obj file
          *
          * @param[in] aFilePath       - Path to the ascii file containing the mesh.
          */
        void load_from_obj_file(const std::string& aFilePath);

// -----------------------------------------------------------------------------
        /**
          * @brief private subroutine for loading a msh file
          *
          * @param[in] aFilePath       - Path to the ascii file containing the mesh.
          *
          */
        void load_from_msh_file(const std::string& aFilePath);

// -----------------------------------------------------------------------------
        /**
          * @brief private subroutine for saving a msh file
          *
          * @param[in] aFilePath       - Path to the ascii file to contain the mesh.
          *
          */
        void save_to_msh(const std::string& aFilePath);

// -----------------------------------------------------------------------------
        /**
          * @brief private subroutine for saving a vtk file
          *
          * @param[in] aFilePath       - Path to the ascii file to contain the mesh.
          *
          */
        void save_to_vtk(const std::string& aFilePath);

// -----------------------------------------------------------------------------
        /**
          * @brief private subroutine for saving an obj file
          *
          * @param[in] aFilePath       - Path to the ascii file to contain the mesh.
          *
          */
        void save_to_obj(const std::string& aFilePath);

// -----------------------------------------------------------------------------
    };
}

#endif //MORIS_CL_GETRIANGLEFILE_HPP
