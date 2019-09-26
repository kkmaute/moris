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
    moris::mtk::Mesh* mMesh;
    std::string mTempFileName;
    std::string mPermFileName;
    int mExoid;
    moris::real mStartTime;
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
    void open_file(const char* aExodusFileName, moris::real aVersion = 4.72);

    /**
     *  Closes the open Exodus database *and* renames it to the permanent file name stored under mPermFileName. This
     *  must be called in order for the Exodus file to be able to be read properly.
     */
    void close_file();

    /**
     *  Writes the coordinates of the nodes in the MTK mesh.
     */
     void write_coordinates();

private:


protected:
};

#endif //MORIS_CL_MTK_WRITER_EXODUS_HPP
