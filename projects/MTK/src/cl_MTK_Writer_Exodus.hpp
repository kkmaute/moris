/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Writer_Exodus.hpp
 *
 */

#ifndef MORIS_CL_MTK_WRITER_EXODUS_HPP
#define MORIS_CL_MTK_WRITER_EXODUS_HPP

#include <exodusII.h>
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

// define maximum number of entities that can be written
// note: not clear why there are issues when exceeding these limits
#define EXODUS_MAX_NUM_SIDESET 1000

namespace moris
{
    namespace mtk
    {
        class Writer_Exodus
        {
            //------------------------------------------------------------------------------

          public:
            std::string mPermFileName;

            // indicate write status of exodus file
            int mExoID = -1;

            //------------------------------------------------------------------------------

          private:
            
            //------------------------------------------------------------------------------
            
            // pointer to the mtk mesh which is written to EXODUS
            moris::mtk::Mesh* mMesh;

            // name maps
            moris::map< std::string, int > mBlockNamesMap;
            moris::map< std::string, int > mSideSetNamesMap;
            moris::map< std::string, int > mNodalFieldNamesMap;
            moris::map< std::string, int > mElementalFieldNamesMap;
            moris::map< std::string, int > mSideSetFieldNamesMap;
            moris::map< std::string, int > mGlobalVariableNamesMap;

            // indices of non-empty sets across all procs
            Vector< uint > mElementBlockIndices;
            Vector< uint > mSideSetIndices;
            Vector< uint > mNodeSetIndices;

            // map mtk element indices to exodus element indices
            Matrix< IndexMat > mMtkExodusElementIndexMap;

            // map exodus element indices to position in exodus file
            Matrix< IndexMat > mExodusElementIndexOrderMap;

            // punch card marking which facets on a given side set are actually present in the final exodus mesh
            // note: facets within the aura may not be part of it
            // input: (1) ordinal of side set in exo mesh, (2) index of facet within side set of enr. IG mesh || output: true for facets in exo mesh
            Vector< Vector< bool > > mFacetUsedInExodus;

            // name of temporary file
            std::string mTempFileName;

            // time step used for writing exodus data
            uint mTimeStep = 0;

            // number of nodes (same in MTK and exodus mesh)
            uint mNumNodes;

            // number of elements in MTK mesh
            uint mNumMtkElements;

            // number of unique elements in exodus mesh
            uint mNumUniqueExodusElements;

            // number of total elements in exodus mesh
            uint mNumTotalExodusElements;

            // flag for using MTK node and element ID maps versus ad-hod maps
            bool mMtkIndexMap = true;

            //------------------------------------------------------------------------------

          public:

            //------------------------------------------------------------------------------
            /**
             * Constructor
             *
             * @param  aMeshPointer Pointer to an MTK mesh
             */
            explicit Writer_Exodus( moris::mtk::Mesh* aMeshPointer );

            //------------------------------------------------------------------------------
            
            /**
             * Constructor
             */
            Writer_Exodus();

            //------------------------------------------------------------------------------
            
            /** Destructor */
            ~Writer_Exodus();

            //------------------------------------------------------------------------------
            
            /**
             * Changes how Exodus handles errors
             *
             * @param abort Causes fatal errors to force program exit.
             * @param debug Causes certain messages to print for debugging use.
             * @param verbose Causes all error messages to print when true, otherwise no error messages will print.
             */
            void set_error_options(
                    bool abort,
                    bool debug,
                    bool verbose );

            //------------------------------------------------------------------------------
            
            /**
             *  Opens an Exodus file and stores the ID for future operations
             *
             *  @param aExodusFileName Name of the Exodus file.
             *  @param aVersion Version of the database. Current version is 4.72 as of programming.
             */
            void open_file(
                    std::string& aExodusFileName,
                    bool         aReadOnly = true,
                    float        aVersion  = 4.72 );

            //------------------------------------------------------------------------------
            
            /**
             * Closes the open Exodus database *and* renames it to the permanent file name stored under mPermFileName. This
             * must be called in order for the Exodus file to be able to be read properly.
             */
            void close_file( bool aRename = true );

            //------------------------------------------------------------------------------
            
            /**
             * Creates an Exodus file and writes everything MTK provides about the given mesh.
             *
             * @param aFilePath The path of the final file
             * @param aFileName The name of the final file
             * @param aTempPath The path of the temporary file
             * @param aTempName The name of a temporary file
             */
            void write_mesh(
                    std::string        aFilePath,
                    const std::string& aFileName,
                    std::string        aTempPath,
                    const std::string& aTempName );

            //------------------------------------------------------------------------------
            
            /**
             * Save temporary to permanent Exodus file.
             */
            void save_mesh();

            //------------------------------------------------------------------------------
            
            /**
             * Creates an Exodus file and writes everything MTK provides about the mesh.
             *
             * @param aFilePath The path of the final file
             * @param aFileName The name of the final file
             * @param aTempPath The path of the temporary file
             * @param aTempName The name of a temporary file
             * @param aCoordinates The coordinates of the points to be written
             */
            void write_points(
                    std::string        aFilePath,
                    const std::string& aFileName,
                    std::string        aTempPath,
                    const std::string& aTempName,
                    Matrix< DDRMat >   aCoordinates );

            //------------------------------------------------------------------------------
            
            /**
             * Sets the number of variables to be written for point data (no mesh)
             *
             * @param aFieldNames The names of the fields that can be written
             */
            void set_point_fields( Vector< std::string > aFieldNames );

            //------------------------------------------------------------------------------
            
            /**
             * Sets the number of variables to be written for nodal data
             *
             * @param aFieldNames The names of the fields that can be written
             */
            void set_nodal_fields( Vector< std::string > aFieldNames );

            //------------------------------------------------------------------------------
            
            /**
             * Sets the number of variables to be written for elemental data
             *
             * @param aFieldNames The names of the fields that can be written
             */
            void set_elemental_fields( Vector< std::string > aFieldNames );

            //------------------------------------------------------------------------------

            /**
             * @brief Sets the number of variables to be written for side set data
             * 
             * @param aFieldNames The names of the fields that can be written
             */
            void
            set_side_set_fields( Vector< std::string > aFieldNames );

            //------------------------------------------------------------------------------
            
            /**
             * Sets the number of variables to be written globally
             *
             * @param aFieldNames The names of the fields that can be written
             */
            void set_global_variables( Vector< std::string > aFieldNames );

            //------------------------------------------------------------------------------
            
            /**
             *  Writes a time to be used for subsequent fields
             *
             *  @param aTimeValue the time for the next time index
             */
            void set_time( moris::real aTimeValue );

            //------------------------------------------------------------------------------
            
            /**
             *  Writes a point field at the current time step.
             *
             *  @param aFieldName The name of the field being written
             *  @param aFieldValues Matrix of values to write for this field.
             */
            void write_point_field(
                    std::string             aFieldName,
                    const Matrix< DDRMat >& aFieldValues );

            //------------------------------------------------------------------------------
            
            /**
             *  Writes a field to the mesh nodes at the current time step.
             *
             *  @param aFieldName The name of the field being written
             *  @param aFieldValues Matrix of values to write for this field.
             */
            void write_nodal_field(
                    std::string             aFieldName,
                    const Matrix< DDRMat >& aFieldValues );

            //------------------------------------------------------------------------------
            
            /**
             *  Writes a field to the mesh elements at the current time step.
             *
             *  @param aBlockName The name of the block that will receive the field
             *  @param aFieldName The name of the field being written
             *  @param aFieldValues Matrix of values to write
             */
            void write_elemental_field(
                    std::string             aBlockName,
                    std::string             aFieldName,
                    const Matrix< DDRMat >& aFieldValues );

            //------------------------------------------------------------------------------

            /**
             * @brief Writes a field to a specified side set at the current time step.
             *
             * @param aSideSetName The name of the side set that will receive the field
             * @param aFieldName The name of the field being written
             * @param aFieldValues Matrix of values to write
             */
            void
            write_side_set_field(
                    std::string             aSideSetName,
                    std::string             aFieldName,
                    const Matrix< DDRMat >& aFieldValues );

            //------------------------------------------------------------------------------
            
            /**
             *  Writes all global variables at the current time step.
             *
             *  @param aVariableNames  cell of names of the variable being written
             *  @param aVariableValues vector of values of the global variables
             */
            void write_global_variables(
                    Vector< std::string >& aVariableNames,
                    const Matrix< DDRMat >&     aVariableValues );

            //------------------------------------------------------------------------------

          private:
            
            //------------------------------------------------------------------------------
            /**
             * Creates an Exodus database at the given file path and string
             *
             * @param aFilePath The path of the final file
             * @param aFileName The name of the final file
             * @param aTempPath The path of the temporary file
             * @param aTempName The name of a temporary file
             */
            void create_file(
                    std::string        aFilePath,
                    const std::string& aFileName,
                    std::string        aTempPath,
                    const std::string& aTempName );

            //------------------------------------------------------------------------------
            
            /**
             * Creates an Exodus database and initializes it at the given file path and string using an MTK mesh
             *
             * @param aFilePath The path of the final file
             * @param aFileName The name of the final file
             * @param aTempPath The path of the temporary file
             * @param aTempName The name of a temporary file
             */
            void create_init_mesh_file(
                    std::string        aFilePath,
                    const std::string& aFileName,
                    std::string        aTempPath,
                    const std::string& aTempName );

            //------------------------------------------------------------------------------
            
            /**
             * Determine number of non-empty node sets across all procs.
             */
            void get_node_sets();

            //------------------------------------------------------------------------------
            
            /**
             *  Determine number of non-empty side sets across all procs.
             */
            void get_side_sets();

            //------------------------------------------------------------------------------
            
            /**
             * Determine number of non-empty blocks across all procs and number of local elements.
             */
            void get_block_sets();

            //------------------------------------------------------------------------------
            
            /**
             * Writes the coordinates of the nodes in the MTK mesh to Exodus.
             */
            void write_nodes();

            //------------------------------------------------------------------------------
            
            /**
             * Writes the node sets in the MTK mesh.
             *
             * @warning This will probably not work, it hasn't been tested yet (I need a mesh with node sets)
             */
            void write_node_sets();

            //------------------------------------------------------------------------------
            
            /**
             * Writes the element blocks in the MTK mesh. Currently supports element and face blocks.
             */
            void write_blocks();

            //------------------------------------------------------------------------------
            
            /**
             * Writes the side sets in the MTK mesh.
             */
            void write_side_sets();

            //------------------------------------------------------------------------------
            
            /**
             * Gets a more detailed description of the elements in the block for exodus from the MTK CellTopology.
             *
             * @param aCellTopology The type of element in MTK.
             * @return Character string describing the Exodus element block.
             */
            const char* get_exodus_block_topology( CellTopology aCellTopology );

            //------------------------------------------------------------------------------
            
            /**
             * Gets the number of nodes in a given element type.
             *
             * @param aCellTopology The type of element in MTK.
             * @return The number of nodes per element of this topology.
             */
            int get_nodes_per_element( CellTopology aCellTopology );

            //------------------------------------------------------------------------------

        };    // end: class Writer_Exodus

        //------------------------------------------------------------------------------

    }    // namespace mtk
}    // namespace moris

#endif    // MORIS_CL_MTK_WRITER_EXODUS_HPP
