//
// Created by christopherson on 9/13/19.
//

#ifndef MORIS_CL_MTK_WRITER_EXODUS_HPP
#define MORIS_CL_MTK_WRITER_EXODUS_HPP

#include <exodusII.h>
#include "cl_MTK_Mesh_Core.hpp"

class Writer_Exodus //TODO wrap with : public Writer
{
public:
    moris::mtk::Mesh*   mMesh;
    std::string         mTempFileName;
    std::string         mPermFileName;
    int                 mExoid;
    /**
    *  Constructor
    *  @param  aMeshPointer Pointer to an MTK mesh
    *  @param  aFilePath File path where temporary and permanent files are saved
    *  @param  aFileName Name of the final file to be saved
    *  @post   don't know what this is
    *  @throw  also not sure what this is
    */
    Writer_Exodus(moris::mtk::Mesh* aMeshPointer, std::string aFilePath, const std::string& aFileName);

    /** Destructor */
    ~Writer_Exodus();

    /**
     *  Changes how Exodus handles errors
     *  @param abort Causes fatal errors to force program exit.
     *  @param debug Causes certain messages to print for debugging use.
     *  @param verbose Causes all error messages to print when true, otherwise no error messages will print.
     */
    void set_error_options(bool abort, bool debug, bool verbose);

    /**
     *  Creates an Exodus file and writes everything MTK provides about the mesh.
     */
    void write_mesh();

private:
    /**
     *  Creates an Exodus database and initializes it based on what is provided by the associated MTK mesh.
     *
     *  @note This will not write anything in the MTK mesh to the database, only the parameters that define it
     *  (e.g. number of nodes).
     */
    void create_file();

    /**
     *  Opens an Exodus file and stores the ID for future operations
     *  @param aExodusFileName Name of the Exodus file.
     *  @param aVersion Version of the database. Current version is 4.72 as of programming.
     */
    void open_file(const char* aExodusFileName, float aVersion = 4.72);

    /**
     *  Closes the open Exodus database *and* renames it to the permanent file name stored under mPermFileName. This
     *  must be called in order for the Exodus file to be able to be read properly.
     */
    void close_file();

    /**
     *  Writes the coordinates of the nodes in the MTK mesh.
     */
    void write_nodes();

    /**
     *  Writes the node sets in the MTK mesh.
     *  @warning This will probably not work, I need a way to get the node ids in a set first.
     */
    void write_node_sets();

    /**
     *  Writes the element blocks in the MTK mesh. Currently supports element and face blocks.
     */
    void write_blocks();

    /**
     *  Gets an exodus element type from an MTK geometry type.
     *  @param aCellTopology The type of element in MTK.
     *  @return Character string of an element type that Exodus can recognize.
     */
    const char* get_exodus_element_type(CellTopology aCellTopology);

    /**
     *  Gets a more detailed description of the elements in the block for exodus from the MTK CellTopology.
     *  @param aCellTopology The type of element in MTK.
     *  @return Character string describing the Exodus element block.
     */
    const char* get_exodus_block_description(CellTopology aCellTopology);

    /**
     * Gets the number of nodes in a given element type.
     * @param aCellTopology The type of element in MTK.
     * @return The number of nodes per element of this topology.
     */
    int64_t get_nodes_per_element(CellTopology aCellTopology);

    /**
     * Gets the number of edges in a given element type.
     * @param aCellTopology The type of element in MTK.
     * @return The number of edges per element of this topology.
     */
    int64_t get_edges_per_element(CellTopology aCellTopology);

    /**
     * Gets the number of faces in a given element type.
     * @param aCellTopology The type of element in MTK.
     * @return The number of faces per element of this topology.
     */
    int64_t get_faces_per_element(CellTopology aCellTopology);

    /**
     * Converts a moris::Cell of std::string's to char**
     * @param aStringCell Cell of strings to be converted.
     * @return The corresponding character array.
     */
    char** string_cell_to_char_array(moris::Cell<std::string> aStringCell);

protected:
};

#endif //MORIS_CL_MTK_WRITER_EXODUS_HPP
