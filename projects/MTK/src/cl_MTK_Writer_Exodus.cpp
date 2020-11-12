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

#include <iostream>

namespace moris
{
    namespace mtk
    {
        //--------------------------------------------------------------------------------------------------------------
        // Public
        //--------------------------------------------------------------------------------------------------------------

        Writer_Exodus::Writer_Exodus(mtk::Mesh* aMeshPointer)
        : mMesh(aMeshPointer)
        {
            this->set_error_options(true, true, true);
        }

        //--------------------------------------------------------------------------------------------------------------

        Writer_Exodus::Writer_Exodus()
        {
            mMesh = nullptr;
            this->set_error_options(true, true, true);
        }

        //--------------------------------------------------------------------------------------------------------------

        Writer_Exodus::~Writer_Exodus()
        = default;

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::set_error_options(
                bool abort,
                bool debug,
                bool verbose)
        {
            ex_opts(abort * EX_ABORT | debug * EX_DEBUG | verbose * EX_VERBOSE);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::open_file(
                std::string  & aExodusFileName,
                bool           aReadOnly,
                float          aVersion)
        {
            MORIS_ERROR(mExoid == -1, "Exodus file is currently open, call close_file() before opening a new one.");

            int tCPUWordSize = sizeof(real), tIOWordSize = 0;

            if (aReadOnly)
            {
                mExoid = ex_open(aExodusFileName.c_str(), EX_READ, &tCPUWordSize, &tIOWordSize, &aVersion);
            }
            else
            {
                mExoid = ex_open(aExodusFileName.c_str(), EX_WRITE, &tCPUWordSize, &tIOWordSize, &aVersion);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::close_file(bool aRename)
        {
            // check that mesh is open
            MORIS_ERROR(mExoid > 0, "Exodus cannot be saved as it is not open\n.");

            ex_close(mExoid);
            mExoid = -1;

            if (aRename)
            {
                MORIS_ERROR( std::rename(mTempFileName.c_str(), mPermFileName.c_str()) == 0 ,"Cannot save exodus file");
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_mesh(
                std::string         aFilePath,
                const std::string & aFileName,
                std::string         aTempPath,
                const std::string & aTempName,
                const uint          aParSize,
                const uint          aParRank)
        {
            MORIS_ERROR(mMesh != nullptr, "No mesh has been given to the Exodus Writer!");

            this->create_init_mesh_file(
                    aFilePath,
                    aFileName,
                    aTempPath,
                    aTempName,
                    aParSize,
                    aParRank);

            this->write_nodes();
            this->write_node_sets();
            this->write_blocks();
            this->write_side_sets();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_points(
                std::string         aFilePath,
                const std::string & aFileName,
                std::string         aTempPath,
                const std::string & aTempName,
                Matrix<DDRMat>      aCoordinates,
                const uint          aParSize,
                const uint          aParRank)
        {
            // Create the actual file
            this->create_file(
                    aFilePath,
                    aFileName,
                    aTempPath,
                    aTempName,
                    aParSize,
                    aParRank);

            // Initialize database
            int tNumDimensions = aCoordinates.n_cols();
            int tNumPoints = aCoordinates.n_rows();

            ex_put_init(mExoid, "MTK", tNumDimensions, tNumPoints, int(1), int(1), int(0), int(0));

            // Set the point coordinates
            int tSpatialDim = aCoordinates.n_cols();
            bool tYDim = tSpatialDim >= 2;
            bool tZDim = tSpatialDim >= 3;

            // Set up coordinate and node map arrays based on the number of vertices
            MORIS_ERROR(aCoordinates.n_rows() > 0, "Points need to be given to create a point field");

            Matrix<IdMat> tNodeMap(aCoordinates.n_rows(), 1, 0);

            // Coordinate arrays
            Matrix<DDRMat> tXCoordinates(aCoordinates.n_rows(), 1, 0.0);
            Matrix<DDRMat> tYCoordinates(aCoordinates.n_rows(), 1, 0.0);
            Matrix<DDRMat> tZCoordinates(aCoordinates.n_rows(), 1, 0.0);

            for (uint tNodeIndex = 0; tNodeIndex < aCoordinates.n_rows(); tNodeIndex++)
            {
                // Place in coordinate arrays
                tXCoordinates(tNodeIndex) = aCoordinates(tNodeIndex, 0);
                tYCoordinates(tNodeIndex) = aCoordinates(tNodeIndex, 1 * tYDim) * tYDim;
                tZCoordinates(tNodeIndex) = aCoordinates(tNodeIndex, 2 * tZDim) * tZDim;

                // Get global ids for id map
                tNodeMap(tNodeIndex) = tNodeIndex + 1;
            }

            // Write coordinates
            ex_put_coord(mExoid, tXCoordinates.data(), tYCoordinates.data(), tZCoordinates.data());

            // Write node id map
            ex_put_id_map(mExoid, EX_NODE_MAP, tNodeMap.data());

            // Create single block
            if (tNumDimensions <= 2)
            {
                ex_put_block(mExoid, EX_ELEM_BLOCK, 1, "CIRCLE", tNumPoints, 1, 0, 0, 0);
            }
            else
            {
                ex_put_block(mExoid, EX_ELEM_BLOCK, 1, "SPHERE", tNumPoints, 1, 0, 0, 0);
            }
            std::string tBlockName("Points");
            ex_put_name(mExoid, EX_ELEM_BLOCK, 1, tBlockName.c_str());

            // Create point elements
            Matrix<IndexMat> tConnectivityArray(tNumPoints, 1);
            Matrix<IndexMat> tElementIdMap(tNumPoints, 1);

            for (int tNodeIndex = 0; tNodeIndex < tNumPoints; tNodeIndex++)
            {
                tConnectivityArray(tNodeIndex) = tNodeIndex + 1;
                tElementIdMap(tNodeIndex) = tNodeIndex + 1;
            }
            ex_put_conn(mExoid, EX_ELEM_BLOCK, 1, tConnectivityArray.data(), nullptr, nullptr);

            // Write the element map
            ex_put_id_map(mExoid, EX_ELEM_MAP, tElementIdMap.data());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::set_point_fields(moris::Cell<std::string> aFieldNames)
        {
            // Set the field names
            if (aFieldNames.size() > 0)
            {
                // Write the number of nodal fields
                ex_put_variable_param(mExoid, EX_NODAL, aFieldNames.size());

                // Write the nodal field names and store as a map
                for (uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++)
                {
                    ex_put_variable_name(mExoid, EX_NODAL, tFieldIndex + 1, aFieldNames(tFieldIndex).c_str());
                    mNodalFieldNamesMap[aFieldNames(tFieldIndex)] = tFieldIndex;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::set_nodal_fields(moris::Cell<std::string> aFieldNames)
        {
            if (aFieldNames.size() > 0)
            {
                // Write the number of nodal fields
                ex_put_variable_param(mExoid, EX_NODAL, aFieldNames.size());

                // Write the nodal field names and store as a map
                for (uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++)
                {
                    ex_put_variable_name(mExoid, EX_NODAL, tFieldIndex + 1, aFieldNames(tFieldIndex).c_str());
                    mNodalFieldNamesMap[aFieldNames(tFieldIndex)] = tFieldIndex;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::set_elemental_fields(moris::Cell<std::string> aFieldNames)
        {
            if (aFieldNames.size() > 0)
            {
                // Write the number of elemental fields
                ex_put_variable_param(mExoid, EX_ELEM_BLOCK, aFieldNames.size());

                // Write the elemental field names and store as a map
                for (uint tFieldIndex = 0; tFieldIndex < aFieldNames.size(); tFieldIndex++)
                {
                    ex_put_variable_name(mExoid, EX_ELEM_BLOCK, tFieldIndex + 1, aFieldNames(tFieldIndex).c_str());
                    mElementalFieldNamesMap[aFieldNames(tFieldIndex)] = tFieldIndex;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::set_global_variables(moris::Cell<std::string> aVariableNames)
        {
            if (aVariableNames.size() > 0)
            {
                // Write the number of global fields
                ex_put_variable_param(mExoid, EX_GLOBAL, aVariableNames.size());

                // Write the global field names and store as a map
                for (uint tVariableIndex = 0; tVariableIndex < aVariableNames.size(); tVariableIndex++)
                {
                    ex_put_variable_name(mExoid, EX_GLOBAL, tVariableIndex + 1, aVariableNames(tVariableIndex).c_str());
                    mGlobalVariableNamesMap[aVariableNames(tVariableIndex)] = tVariableIndex;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::save_mesh()
        {
            // check that mesh is open
            MORIS_ERROR(mExoid > 0, "Exodus cannot be saved as it is not open\n.");

            // close mesh
            ex_close(mExoid);

            // write log information
            std::string tMessage = "Copying " + mTempFileName + " to " + mPermFileName + ".";
            MORIS_LOG( tMessage.c_str() );

            // copy temporary file on permanent file
            std::ifstream src(mTempFileName.c_str(),  std::ios::binary);
            std::ofstream dest(mPermFileName.c_str(), std::ios::binary);
            dest << src.rdbuf();

            // open mesh file again
            int tCPUWordSize = sizeof(real), tIOWordSize = 0;
            float tVersion;

            mExoid = ex_open(
                    mTempFileName.c_str(),
                    EX_WRITE,
                    &tCPUWordSize,
                    &tIOWordSize,
                    &tVersion);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::set_time(real aTimeValue)
        {
            ex_put_time(mExoid, ++mTimeStep, &aTimeValue);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_point_field(
                std::string            aFieldName,
                const Matrix<DDRMat> & aFieldValues)
        {
            // skip if no nodal values exist
            if ( aFieldValues.numel() == 0 )
            {
                return;
            }

            // Field name to index
            int tFieldIndex = mNodalFieldNamesMap[aFieldName];

            MORIS_ASSERT(mNodalFieldNamesMap.size() == mNodalFieldNamesMap.size(), aFieldName.append(
                    " is not a point field name on this mesh!").c_str());

            // Write the field
            int tErrMsg = ex_put_var(mExoid, mTimeStep, EX_NODAL, tFieldIndex + 1, 1, aFieldValues.numel(), aFieldValues.data());

            // Check for error
            MORIS_ERROR(tErrMsg==0,"Point field %s could not be written to exodus file.",aFieldName.c_str());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_nodal_field(
                std::string            aFieldName,
                const Matrix<DDRMat> & aFieldValues)
        {
            // skip if no nodal values exist
            if ( aFieldValues.numel() == 0 )
            {
                return;
            }

            // Field name to index
            int tFieldIndex = mNodalFieldNamesMap[aFieldName];

            MORIS_ASSERT(mNodalFieldNamesMap.size() == mNodalFieldNamesMap.size(), aFieldName.append(
                    " is not a nodal field name on this mesh!").c_str());

            // Check number of field values = number of nodes
            MORIS_ERROR(aFieldValues.numel() == mMesh->get_num_nodes(), aFieldName.append(
                    " field was attempted to be written with " + std::to_string(aFieldValues.numel()) + " values, but there are " +
                    std::to_string(mMesh->get_num_nodes()) + " nodes in this mesh.").c_str());

            // Ensure that time step is larger than or equal 1
            int tTimeStep = mTimeStep > 1 ? mTimeStep : 1;

            // Write the field
            int tErrMsg = ex_put_var(mExoid, tTimeStep, EX_NODAL, tFieldIndex + 1, 1, aFieldValues.numel(), aFieldValues.data());

            // Check for error
            MORIS_ERROR(tErrMsg==0,"Nodal field %s could not be written to exodus file.",aFieldName.c_str());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_elemental_field(       
                std::string            aBlockName,
                std::string            aFieldName,
                const Matrix<DDRMat> & aFieldValues)
        {
            // skip if no elemental values exist
            if ( aFieldValues.numel() == 0 )
            {
                return;
            }

            MORIS_ERROR(mBlockNamesMap.key_exists(aBlockName), aBlockName.append(
                    " is not a block name on this mesh!").c_str());

            // Block name to index
            int tBlockIndex = mBlockNamesMap[aBlockName];

            // Check that block index is valid
            int tNumBlocks = ex_inquire_int(mExoid, EX_INQ_ELEM_BLK);

            MORIS_ERROR(tNumBlocks > tBlockIndex,"Index of block set is larger than number of blocks");

            // Field name to index
            int tFieldIndex = mElementalFieldNamesMap[aFieldName];

            MORIS_ERROR(mElementalFieldNamesMap.size() == mElementalFieldNamesMap.size(), aFieldName.append(
                    " is not an elemental field name on this mesh!").c_str());

            // Check number of field values = number of elements
            MORIS_ERROR(aFieldValues.numel() == mMesh->get_set_cells(aBlockName).size(), aFieldName.append(
                    " field was attempted to be written with " + std::to_string(aFieldValues.numel()) + " values, but there are " +
                    std::to_string(mMesh->get_set_cells(aBlockName).size()) + " elements in block " + aBlockName +
                    ".").c_str());

            // Check that number of field values = number of element stored in mesh
            ex_block tBlockInfo;
            tBlockInfo.id   = tBlockIndex + 1;
            tBlockInfo.type = EX_ELEM_BLOCK;

            ex_get_block_param(
                    mExoid,
                    &tBlockInfo);

            MORIS_ERROR(tBlockInfo.num_entry == (int) aFieldValues.numel(),
                    "Number of entries in field does not match number of elements stored in mesh for current block.");

            // Ensure that time step is larger than or equal 1
            int tTimeStep = mTimeStep > 1 ? mTimeStep : 1;

            // Write the field
            int tErrMsg = ex_put_var(mExoid, tTimeStep , EX_ELEM_BLOCK, tFieldIndex + 1, tBlockIndex + 1,
                    aFieldValues.numel(), aFieldValues.data());

            // Check for error
            MORIS_ERROR(tErrMsg==0,"Elemental field %s could not be written for element block %s to exodus file.",
                    aFieldName.c_str(),aBlockName.c_str());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_global_variables(
                moris::Cell<std::string>    & aVariableNames,
                const Matrix<DDRMat> & aVariableValues)
        {
            // number of global variables
            uint tNumVarialbes = aVariableNames.size();

            // check that number of variables names is not larger than number of values
            MORIS_ASSERT(tNumVarialbes <= aVariableValues.numel(),
                    "Number of global variables names larger than number of values.");

            // allocated vector of sorted global variables
            Matrix<DDRMat> tSortedValues(tNumVarialbes,1);

            for (uint ig=0;ig<tNumVarialbes;++ig)
            {
                // Variable name to index
                uint tVariableIndex = mGlobalVariableNamesMap[ aVariableNames(ig) ];

                // check that index is valid
                MORIS_ASSERT(tVariableIndex < tNumVarialbes,"Global variable index too large.");

                // store global variable in sorted list
                tSortedValues(tVariableIndex) =  aVariableValues(ig);
            }

            // Write the variables
            int tErrMsg = ex_put_var(mExoid, mTimeStep, EX_GLOBAL, 1, 0, tNumVarialbes, tSortedValues.data());

            // check for error
            MORIS_ERROR(tErrMsg==0,"Global variables could not be written to exodus file.");
        }

        //--------------------------------------------------------------------------------------------------------------
        // Private
        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::create_file(
                std::string         aFilePath,
                const std::string & aFileName,
                std::string         aTempPath,
                const std::string & aTempName,
                const uint          aParSize,
                const uint          aParRank)
        {
            MORIS_ERROR(mExoid == -1, "Exodus file is currently open, call close_file() before creating a new one.");

            // Add temporary and permanent file names to file paths
            if (!aFilePath.empty())
            {
                aFilePath += "/";
            }

            if (!aTempPath.empty())
            {
                aTempPath += "/";
            }

            mTempFileName = aTempPath + aTempName;
            mPermFileName = aFilePath + aFileName;

            // determine number of processors and rank of this processor
            uint tParSize = aParSize > 0 ? aParSize : par_size();
            uint tParRank = aParSize > 0 ? aParRank : par_rank();

            // Make file name parallel, if necessary
            if (tParSize > 1)
            {

                // Get par size and rank as strings
                std::string tParSizeStr     =  std::to_string(tParSize);
                std::string tParRankBaseStr =  std::to_string(tParRank);

                // Make sure all par rank strings have the same number of leading zeros
                std::string tParRankStr = std::string( tParSizeStr.length() - tParRankBaseStr.length(), '0').append( tParRankBaseStr );

                // Append to temporary and permament file names
                mTempFileName += "." + tParSizeStr + "." + tParRankStr;
                mPermFileName += "." + tParSizeStr + "." + tParRankStr;
            }

            // Create the database
            int cpu_ws = sizeof( real );         // word size in bytes of the floating point variables used in moris
            int io_ws  = sizeof( real );         // word size as stored in exodus

            mExoid = ex_create(mTempFileName.c_str(), EX_CLOBBER, &cpu_ws, &io_ws);

            MORIS_ERROR(mExoid > -1,"Exodus file cannot be created: %s",mTempFileName.c_str());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::create_init_mesh_file(
                std::string         aFilePath,
                const std::string & aFileName,
                std::string         aTempPath,
                const std::string & aTempName,
                const uint          aParSize,
                const uint          aParRank)
        {
            // Create the actual file
            this->create_file(
                    aFilePath,
                    aFileName,
                    aTempPath,
                    aTempName,
                    aParSize,
                    aParRank);

            // Number of dimensions
            int tNumDimensions = mMesh->get_spatial_dim();

            // Number of nodes
            int tNumNodes = mMesh->get_num_nodes();

            // Determine number of non-empty node sets
            moris::Cell<std::string> tNodeSetNames = mMesh->get_set_names(EntityRank::NODE);
            mNumNodeSets = 0;
            for (uint tNodeSetIndex = 0; tNodeSetIndex < tNodeSetNames.size(); tNodeSetIndex++)
            {
                // FIXME: should be replaced by call to MTK to just get number of nodes in node set (local or global?)
                Matrix<IndexMat> tNodeIndices = mMesh->get_set_entity_loc_inds(
                        EntityRank::NODE,
                        tNodeSetNames(tNodeSetIndex));

                if ( tNodeIndices.numel() > 0 )
                {
                    mNumNodeSets++;
                }
            }

            // Determine number of non-empty side sets
            moris::Cell<std::string> tSideSetNames = mMesh->get_set_names(mMesh->get_facet_rank());
            mNumSideSets = 0;
            for (uint tSideSetIndex = 0; tSideSetIndex < tSideSetNames.size(); tSideSetIndex++)
            {
                // FIXME: should be replaced by call to MTK to just get number of faces in side set (local or global?)
                Matrix<IndexMat> tSideSetElements;
                Matrix<IndexMat> tSideSetOrdinals;

                mMesh->get_sideset_elems_loc_inds_and_ords(
                        tSideSetNames(tSideSetIndex),
                        tSideSetElements,
                        tSideSetOrdinals);

                if (tSideSetOrdinals.numel() > 0 )
                {
                    mNumSideSets++;
                }
            }

            // Determine number of non-empty blocks and elements
            int tNumElements  = 0;
            mNumElementBlocks = 0;

            moris::Cell<std::string> tBlockNames = mMesh->get_set_names(EntityRank::ELEMENT);

            for (uint tBlockIndex = 0; tBlockIndex < tBlockNames.size(); tBlockIndex++)
            {
                uint tNumElementsInBlock = mMesh->get_element_indices_in_block_set(tBlockIndex).length();

                if ( mMesh->get_element_indices_in_block_set(tBlockIndex).length() >0 )
                {
                    mNumElementBlocks ++;
                    tNumElements += tNumElementsInBlock;
                }
            }

            // Initialize database
            ex_put_init(mExoid, "MTK", tNumDimensions, tNumNodes, tNumElements, mNumElementBlocks, mNumNodeSets, mNumSideSets);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_nodes()
        {
            // spatial dimension
            int tSpatialDim = mMesh->get_spatial_dim();
            bool tYDim = tSpatialDim >= 2;
            bool tZDim = tSpatialDim >= 3;

            // Set up coordinate and node map arrays based on the number of vertices
            Matrix<IdMat> tNodeMap(mMesh->get_num_nodes(), 1, 0);

            // Coordinate arrays
            Matrix<DDRMat> tXCoordinates(mMesh->get_num_nodes(), 1, 0.0);
            Matrix<DDRMat> tYCoordinates(mMesh->get_num_nodes(), 1, 0.0);
            Matrix<DDRMat> tZCoordinates(mMesh->get_num_nodes(), 1, 0.0);

            for (uint tNodeIndex = 0; tNodeIndex < mMesh->get_num_nodes(); tNodeIndex++)
            {
                // Get coordinates
                Matrix<DDRMat> tNodeCoordinates = mMesh->get_node_coordinate(tNodeIndex);

                // Place in coordinate arrays
                tXCoordinates(tNodeIndex) = tNodeCoordinates(0);
                tYCoordinates(tNodeIndex) = tNodeCoordinates(1 * tYDim) * tYDim;
                tZCoordinates(tNodeIndex) = tNodeCoordinates(2 * tZDim) * tZDim;

                // Get global ids for id map
                tNodeMap(tNodeIndex) = mMesh->get_glb_entity_id_from_entity_loc_index(tNodeIndex, EntityRank::NODE);
            }

            // Write coordinates
            ex_put_coord(mExoid, tXCoordinates.data(), tYCoordinates.data(), tZCoordinates.data());

            // Write node id map
            ex_put_id_map(mExoid, EX_NODE_MAP, tNodeMap.data());
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_node_sets()
        {
            // Get the number of node sets and their names
            moris::Cell<std::string> tNodeSetNames = mMesh->get_set_names(EntityRank::NODE);

            // Write each node set
            int tNumNodeSets = 0;

            for (uint tNodeSetIndex = 0; tNodeSetIndex <  tNodeSetNames.size(); tNodeSetIndex++)
            {
                Matrix<IndexMat> tNodeIndices = mMesh->get_set_entity_loc_inds(
                        EntityRank::NODE,
                        tNodeSetNames(tNodeSetIndex));

                if ( tNodeIndices.numel() > 0 )
                {
                    ex_put_set_param(mExoid, EX_NODE_SET, tNumNodeSets + 1, tNodeIndices.numel(), 0);
                    ex_put_set      (mExoid, EX_NODE_SET, tNumNodeSets + 1, tNodeIndices.data(), nullptr);

                    // Increase counter for non-empty node sets
                    tNumNodeSets++;
                }
            }

            // Check for correct non-empty node set count
            MORIS_ERROR(tNumNodeSets == mNumNodeSets, "Incorrect number of non-empty node sets written to file." );

        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_blocks()
        {
            // Start the element maps
            mMtkExodusElementIndexMap.set_size(mMesh->get_num_elems()+1, 1, 0);
            Matrix<IdMat> tElementIdMap(0, 1, 0);
            uint tElementIdMapStartIndex = 0;

            // All of the block names
            moris::Cell<std::string> tBlockNames = mMesh->get_set_names(EntityRank::ELEMENT);

            // Loop through each block
            int tNumElementBlocks = 0;

            for (uint tBlockIndex = 0; tBlockIndex < tBlockNames.size(); tBlockIndex++)
            {
                // Get the block element identifiers
                Matrix<IndexMat> tElementIndices = mMesh->get_element_indices_in_block_set(tBlockIndex);
                Matrix<IdMat> tElementIDs        = mMesh->get_element_ids_in_block_set(tBlockIndex);

                MORIS_ASSERT(tElementIndices.length() == tElementIDs.length(),
                        "Vector of element IDs must match vector of element indices for Exodus file writing.");

                // Number of elements in the block
                uint tNumElementsInBlock = tElementIndices.length();

                // Ignore empty blocks
                if (tNumElementsInBlock > 0)
                {
                    // Add name to map
                    mBlockNamesMap[tBlockNames(tBlockIndex)] = tNumElementBlocks;

                    // Resize element map
                    tElementIdMap.resize(tElementIdMapStartIndex + tNumElementsInBlock, 1);

                    // Get the CellTopology of this block
                    enum CellTopology tMorisBlockTopology = mMesh->get_blockset_topology(tBlockNames(tBlockIndex));
                    const char* tExodusBlockTopology = this->get_exodus_block_topology(tMorisBlockTopology);

                    // Get the number of nodes/edges/faces/attributes per element
                    int tNumNodesPerElement      = this->get_nodes_per_element(tMorisBlockTopology);
                    int tNumEdgesPerElement      = 0;
                    int tNumFacesPerElement      = 0;
                    int tNumAttributesPerElement = 0;

                    // Make a block and name it
                    ex_put_block(mExoid, EX_ELEM_BLOCK, tNumElementBlocks + 1, tExodusBlockTopology, tNumElementsInBlock,
                            tNumNodesPerElement, tNumEdgesPerElement, tNumFacesPerElement, tNumAttributesPerElement);
                    ex_put_name(mExoid, EX_ELEM_BLOCK, tNumElementBlocks + 1, tBlockNames(tBlockIndex).c_str());

                    // Construct matrix of node indices per element
                    Matrix<IndexMat> tConnectivityArray(tNumNodesPerElement * tNumElementsInBlock, 1, 0);

                    // Loop through the elements in this block
                    uint tConnectivityIndex = 0;
                    for (uint tElementNumber = 0; tElementNumber < tNumElementsInBlock; tElementNumber++)
                    {
                        // Get the vertex indices of this element
                        Matrix<IndexMat> tNodeIndices = mMesh->get_nodes_connected_to_element_loc_inds(tElementIndices(tElementNumber));

                        // Assign each vertex individually
                        for (int tNodeNum = 0; tNodeNum < tNumNodesPerElement; tNodeNum++)
                        {
                            tConnectivityArray(tConnectivityIndex) = tNodeIndices(tNodeNum) + 1;
                            tConnectivityIndex++;
                        }

                        // Get the global id and index of this element and add to maps
                        tElementIdMap(tElementIdMapStartIndex + tElementNumber) = tElementIDs(tElementNumber);
                        mMtkExodusElementIndexMap(tElementIndices(tElementNumber)) = tElementIdMapStartIndex + tElementNumber;
                    }

                    // Update location in element map
                    tElementIdMapStartIndex += tNumElementsInBlock;

                    // Write connectivity
                    ex_put_conn(mExoid, EX_ELEM_BLOCK, tNumElementBlocks + 1, tConnectivityArray.data(), nullptr, nullptr);

                    // Increase counter for non-empty element blocks
                    tNumElementBlocks++;
                }
            }

            // Write the element map
            ex_put_id_map(mExoid, EX_ELEM_MAP, tElementIdMap.data());

            // Check for correct non-empty side set count
            MORIS_ERROR(tNumElementBlocks == mNumElementBlocks, "Incorrect number of non-empty element blocks written to file." );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Writer_Exodus::write_side_sets()
        {
            // Get side set names
            moris::Cell<std::string> tSideSetNames = mMesh->get_set_names(mMesh->get_facet_rank());

            // Write side sets
            int tNumSideSets = 0;

            for (uint tSideSetIndex = 0; tSideSetIndex < tSideSetNames.size(); tSideSetIndex++)
            {
                // Get the side set element ids
                Matrix<IndexMat> tSideSetElements;
                Matrix<IndexMat> tSideSetOrdinals;

                mMesh->get_sideset_elems_loc_inds_and_ords(tSideSetNames(tSideSetIndex), tSideSetElements, tSideSetOrdinals);

                // Skip if side set is empty
                if (tSideSetElements.numel() > 0 )
                {
                    // Change element and ordinal to what Exodus wants
                    for (uint tElementNum = 0; tElementNum < tSideSetOrdinals.numel(); tElementNum++)
                    {
                        tSideSetElements(tElementNum) = mMtkExodusElementIndexMap(tSideSetElements(tElementNum)) + 1;
                        tSideSetOrdinals(tElementNum)++;
                    }

                    // Write the side set
                    ex_put_set_param(mExoid, EX_SIDE_SET, tNumSideSets + 1, tSideSetElements.numel(), 0);
                    ex_put_set      (mExoid, EX_SIDE_SET, tNumSideSets + 1, tSideSetElements.data(), tSideSetOrdinals.data());
                    ex_put_name     (mExoid, EX_SIDE_SET, tNumSideSets + 1, tSideSetNames(tSideSetIndex).c_str());

                    // Increase counter for non-empty side sets
                    tNumSideSets++;
                }
            }

            // Check for correct non-empty side set count
            MORIS_ERROR(tNumSideSets == mNumSideSets, "Incorrect number of non-empty side sets written to file." );
        }

        //--------------------------------------------------------------------------------------------------------------
        // Static (for now)
        //--------------------------------------------------------------------------------------------------------------

        const char* Writer_Exodus::get_exodus_block_topology(CellTopology aCellTopology)
        {
            switch (aCellTopology)
            {
                case CellTopology::TRI3:
                    return "TRI3";
                case CellTopology::QUAD4:
                    return "QUAD4";
                case CellTopology::TET4:
                    return "TET4";
                case CellTopology::TET10:
                    return "TET10";
                case CellTopology::HEX8:
                    return "HEX8";
                case CellTopology::PRISM6:
                    return "PRISM6";
                default:
                    MORIS_ERROR(false, "This element is invalid or it hasn't been implemented yet!");
                    return "";
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        int Writer_Exodus::get_nodes_per_element(CellTopology aCellTopology)
        {
            switch (aCellTopology)
            {
                case CellTopology::TRI3:
                    return 3;
                case CellTopology::QUAD4:
                    return 4;
                case CellTopology::TET4:
                    return 4;
                case CellTopology::TET10:
                    return 10;
                case CellTopology::HEX8:
                    return 8;
                case CellTopology::PRISM6:
                    return 6;
                default:
                    MORIS_ERROR(false, "This element is invalid or it hasn't been implemented yet!");
                    return 0;
            }
        }
    }
}

