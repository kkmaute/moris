//
// Created by christopherson on 9/13/19.
//

#ifndef MORIS_CL_MTK_WRITER_EXODUS_HPP
#define MORIS_CL_MTK_WRITER_EXODUS_HPP

#include <exodusII.h>
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

class Writer_Exodus
{
public:
    std::string                     mPermFileName;

private:
    int                             mExoid;
    moris::mtk::Mesh*               mMesh;
    std::string                     mTempFileName;
    moris::uint                     mTimeStep = 0;
    moris::Matrix<moris::IndexMat>  mMtkExodusElementIndexMap;
    moris::map<std::string, int>    mNodalFieldNames;
    moris::map<std::string, int>    mElementalFieldNames;
    moris::map<std::string, int>    mGlobalVariableNames;

public:
    /**
    * Constructor
     *
    * @param  aMeshPointer Pointer to an MTK mesh
    * @param  aFilePath File path where temporary and permanent files are saved
    * @param  aFileName Name of the final file to be saved
    */
    explicit Writer_Exodus(moris::mtk::Mesh* aMeshPointer);

    /** Destructor */
    ~Writer_Exodus();

    /**
     * Changes how Exodus handles errors
     *
     * @param abort Causes fatal errors to force program exit.
     * @param debug Causes certain messages to print for debugging use.
     * @param verbose Causes all error messages to print when true, otherwise no error messages will print.
     */
    void set_error_options(bool abort, bool debug, bool verbose);

    /**
     * Creates an Exodus file and writes everything MTK provides about the mesh.
     */
    void write_mesh(std::string aFilePath, const std::string& aFileName);

    /**
     * Sets the number of variables to be written for nodal data
     *
     * @param aFieldNames The names of the fields that can be written
     */
    void set_nodal_fields(moris::Cell<std::string> aFieldNames);

    /**
     * Sets the number of variables to be written for elemental data
     *
     * @param aFieldNames The names of the fields that can be written
     */
    void set_elemental_fields(moris::Cell<std::string> aFieldNames);

    /**
     * Sets the number of variables to be written globally
     *
     * @param aFieldNames The names of the fields that can be written
     */
    void set_global_variables(moris::Cell<std::string> aFieldNames);

    /**
     *  Writes a time to be used for subsequent fields
     *
     *  @param aTimeValue the time for the next time index
     */
    void set_time(moris::real aTimeValue);

    /**
     *  Writes a field to the mesh nodes at the current time step.
     *
     *  @param aFieldName The name of the field being written
     *  @param aFieldValues Matrix of values to write for this field.
     */
    void write_nodal_field(std::string aFieldName, moris::Matrix<moris::DDRMat> aFieldValues);

    /**
     *  Writes a field to the mesh elements at the current time step.
     *
     *  @param aBlockIndex The index of the block that will receive the field
     *  @param aFieldName The name of the field being written
     *  @param aFieldValues Matrix of values to write
     */
    void write_elemental_field(moris::uint aBlockIndex, std::string aFieldName, moris::Matrix<moris::DDRMat> aFieldValues);

    /**
     *  Writes a global variable at the current time step.
     *
     *  @param aVariableName The name of the variable being written
     *  @param aVariableValue The value of the global variable
     */
    void write_global_variable(std::string aVariableName, moris::real aVariableValue);

    /**
     *  Opens an Exodus file and stores the ID for future operations
     *
     *  @param aExodusFileName Name of the Exodus file.
     *  @param aVersion Version of the database. Current version is 4.72 as of programming.
     */
    void open_file(std::string aExodusFileName, bool aReadOnly = true, float aVersion = 4.72);

    /**
     * Closes the open Exodus database *and* renames it to the permanent file name stored under mPermFileName. This
     * must be called in order for the Exodus file to be able to be read properly.
     */
    void close_file(bool aRename = true);


private:
    /**
     * Creates an Exodus database and initializes it based on what is provided by the associated MTK mesh.
     *
     * @note This will not write anything in the MTK mesh to the database, only the parameters that define it
     * (e.g. number of nodes).
     */
    void create_file(std::string aFilePath, const std::string& aFileName);

    /**
     * Writes the coordinates of the nodes in the MTK mesh to Exodus.
     */
    void write_nodes();

    /**
     * Writes the node sets in the MTK mesh.
     *
     * @warning This will probably not work, it hasn't been tested yet (I need a mesh with node sets)
     */
    void write_node_sets();

    /**
     * Writes the element blocks in the MTK mesh. Currently supports element and face blocks.
     */
    void write_blocks();

    /**
     * Writes the side sets in the MTK mesh.
     */
    void write_side_sets();

    /**
     * Gets an exodus element type from an MTK geometry type.
     *
     * @param aCellTopology The type of element in MTK.
     * @return Character string of an element type that Exodus can recognize.
     */
    const char* get_exodus_element_type(CellTopology aCellTopology);

    /**
     * Gets a more detailed description of the elements in the block for exodus from the MTK CellTopology.
     *
     * @param aCellTopology The type of element in MTK.
     * @return Character string describing the Exodus element block.
     */
    const char* get_exodus_block_description(CellTopology aCellTopology);

    /**
     * Gets the number of nodes in a given element type.
     *
     * @param aCellTopology The type of element in MTK.
     * @return The number of nodes per element of this topology.
     */
    int get_nodes_per_element(CellTopology aCellTopology);

    /**
     * Converts a moris::Cell of std::string's to char**
     *
     * @param aStringCell Cell of strings to be converted.
     * @return The corresponding character array.
     */
    char** string_cell_to_char_array(moris::Cell<std::string> aStringCell);

protected:
};

#endif //MORIS_CL_MTK_WRITER_EXODUS_HPP
