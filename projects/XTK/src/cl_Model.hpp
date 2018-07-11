
#ifndef MORIS_XTK_MODEL_HPP_
#define MORIS_XTK_MODEL_HPP_

#include "core.hpp"
#include "assert.hpp"

// MORIS project header files.
#include "cl_Mesh_XTK.hpp"                // For use of the simple mesh class internal to xtk
#include "cl_XTK_Enums.hpp"               // For use of the xtk specific enums
#include "cl_Entity_Tracker.hpp"
#include "cl_Mat.hpp"              // For use of Mat // LNA/src
#include "fn_sort.hpp"             // For use of sort function // LNA/src
#include "cl_Cell.hpp"         // For use of Cell // CON/src
#include "cl_Mesh.hpp" // MTK/src               // For use of the Mesh wrapper class around STK // MTK/src
#include "cl_GeometryEngine.hpp"  // For use of geometry engine class
#include "cl_Pairing.hpp" // TOL/src
#include "mpi.h"

// Keenan Doble 10/01/16

namespace moris
{
    namespace xtk
    {
        class Model
        {
        public:
            /*
             * Primary constructor (for full execution of program this is used)
             */
            Model(mesh*                aMesh,
                  GeometryEngine*      aGeometryEngine,
                  uint                 aModelDimension);

            /*
             * Simple constructor (minimum amount of information required in model)
             */
            Model(uint    aModelDimension);

            /* Test case constructor
             * Used for tests: formulate_node_request();
             * Only used when you do not need to use the geometry engine's functionality
             */
            Model(mesh*    aMesh,
                  uint     aModelDimension);

            ~Model();


            /*
             * Request bucket stores all of the requests for creating entities of same
             * parent entity rank and child entity rank
             */

            class RequestHandler
            {
            public:
                RequestHandler(uint              aNumExpectedRequests,
                               uint              aNumChildrenAllowed,
                               enum EntityRank   aParentEntityRank,
                               enum EntityRank   aChildEntityRank,
                               mesh            * aParentMesh):
                                   mEntityTracker(aParentEntityRank,aChildEntityRank,aParentMesh->get_num_entities_universal(aParentEntityRank),aNumChildrenAllowed)
                {
                    mParentEntityRank   = aParentEntityRank;
                    mChildEntityRank    = aChildEntityRank;
                    mRequestCounter     = 0;
                    mMeshPtr            = aParentMesh;
                    mNumChildrenAllowed = aNumChildrenAllowed;
                    mRequestInfo.set_size(aNumExpectedRequests,2,-1);
                }

                ~RequestHandler(){}

                /*
                 * When you set a request, you get a pointer to the location where
                 * the index is going to be created
                 */
                uint*
                set_request_info(uint              aParentEntityIndex,
                                   Mat<real>  aChildCoords);

                /*
                 * Same as above but for when a secondary entity index is needed
                 * Returns a pointer to the entity index
                 */
                uint*
                set_request_info(uint              aParentEntityIndex,
                                 uint              aSecondaryEntityIndex,
                                 Mat<real>  aChildCoords);


                void handle_requests();

                // -------------------DEBUG FUNCTIONS --------------------------
                void
                print()
                {
                    MORIS_LOG_INFO<<"Printing the Request Manager Data Structure";
                    uint temp = mRequestInfo.n_rows();
                    MORIS_LOG_INFO<<"Number of Request: "<< temp;

//                    for(uint i = 0; i<temp; i++)
//                    {
//                            MORIS_LOG_INFO<<mRequestInfo.row(i)<<"  "<< mCoords.row(i);
//                    }
                }

                uint
                get_num_requests()
                {
                    return mRequestCounter;
                }
            private:
                uint                      mRequestCounter;
                uint                      mNumChildrenAllowed;
                enum EntityRank           mParentEntityRank;
                enum EntityRank           mChildEntityRank;
                Mat<uint>                 mRequestInfo;
                xtk::EntityTracker        mEntityTracker;
                mesh                    * mMeshPtr;
            };

            /* @brief Decomposes a mesh using a geometry engine (split up across processors)
             * @param[in] aMesh - constant pointer to a constant mesh
             * @param[in] aGeometryEngine - constant pointer to a geometry engine which contains geometry information
             * @param[in] aMethod - Decomposition method to use. (ie. REGULAR_HIER specifies regular subdivision followed by hierarchy subdivision)
             */

            Cell<MeshXTK>
            decompose(enum Decomposition    aMethod);

            //FUNCTIONS BELOW THIS LINE ARE FOR UNIT TEST USE ONLY-------------------------

            /* Do not use this function outside of unit tests
             * Purpose: Serve as an interface to private handle node request function for testing purpose
             *
             */
            void
            handle_node_request_test_interface(Cell<MeshXTK>    & aXTKMeshList);
            /* Do not use this function outside of unit tests
             * Purpose: Serve as an interface to private regular subdivision function for testing purpose
             *
             */
            void
            regular_subdivision_test_interface(Cell<MeshXTK>    & aXTKMeshList);

            Cell<Mat<uint>>
            pack_XTK_mesh(Cell<MeshXTK>  & aXTKMesh,
                          uint             aPackType);


        private:
            mesh                    * mParentMesh;       // Pointer to a constant parent mesh (STK) to use for information
            GeometryEngine          * mGeometryEngines;  // Will change to a cell of pointers to geometry engines
            uint                      mModelDimension;   // Dimension of the model (1 - 1D, 2- 2D, 3-3D)

            // Decomposition Functions------------------------------------------------------

            /* regular_subdivision uses template topologies to subdivide a rectangular type element into triangular
             * type elements. For 2D, the rectangular element is subdivided into 4 triangular elements. For 3D, a
             * hexahedral element is subdivided into 24 tetrahedral elements.
             *
             * @param[in]  aGeometryObjectList         cell containing reference to geometry objects
             * @param[out] aSimpleMesh                 a simple mesh (only nodes are unique)
             */
            void
            regular_subdivision(Cell<MeshXTK>    & aXTKMeshList);

            /*
             * @brief Handles a request to create a new entity. Makes sure that a new entity is not created multiple time, also computes coordinates of new nodes created
             *
             * @param[in] aRequestEntityDim - Dimension of the entity being created
             * i.e. aRequestDim = 0 signifies a new node request
             *
             * @param[in] aGID          First global entity Id which will be used by function
             * @param[in] aRequest      The request for a new entity. Needs to contain the following:
             *                                   Column 0.)   Parent Entity Dimension
             *                                   Column 1.)   Parent Entity Id
             *                                   Column 2-4.) Coordinates
             *                                   Column 5.)   Last column is reserved for the new entities Id
             *
             * @param[in] aRequestTracker....reference to a list of used Parent entities
             *
             * @param[out] A matrix of new entity ids
             */

//            void
//            conformal_subdivision(Cell<geomeng::GeometryObject> const    & aGeometryObjectList,
//                                  Cell<MeshXTK>                          & aSimpleMesh);



            /*formulates node requests in the geometry objects. Dependent on the type of decomposition
             * @param[in] aReqType- specifies which template mesh is going to be used later on
             *
             *see cl_model.cpp for template on how to code a new request type
             */
            void
            subdivide(enum RequestType                 aReqType,
                      Cell<MeshXTK>    & aXTKMesh);


            /*
             * Handles the node requests using an entity tracking matrix and fulfil_node_request function in
             * geometry objects. Requests are formulated in the geometry engine intersection_check() and stored
             * in the geometry objects.
             *
             * @param[in]  aGeometryObjectList - Cell of geometry objects which contain node requests.
             * @param[out] aXTKMeshList - A list of xtk meshes.
             */
            void
            handle_node_request(Cell<MeshXTK>                    & aXTKMeshList);



            /*
             * Processor manager
             * Dictates whether a processor is active and whether or not the current processor sends or receives from it
             * Makes sure a processor is not added to list more than once
             * Places on send or receive depending on processor rank
             *
             * @param[in] aCurrentProcRank   - Current processor rank
             * @param[in] aOtherProcRank     - Other processor rank
             * @param[in] aActiveProcTracker - row corresponds to proc rank
             *                                 column 0 - active flag  (0 - not active, 1- active sender, 2-active receiver)
             *                                 column 1 - active index (index in the active lists)
             * @param[in] aActSendProcList   - active processors to send to
             * @param[in] aActRecProcList    - active processors to receive from
             * @param[in] aSendIds           - information to send to processors (does not populate)
             * @param[in] aRecIds            - information to received from processors (does not populate)
             *
             * Returns a matrix containing which type of processor communication and the index of the active processor
             */
//            uint active_processor_manager(uint                 aCurrentProcRank,
//                                                 uint                 aOtherProcRank,
//                                                 Mat<uint>   & aActiveProcTracker,
//                                                 Mat<uint>   & aActiveProcList);


            // debugging tool to print the tSendIds and aRecIds to loggers
            void print_send_communication_Ids(Cell<Cell<Cell<uint>>>  & aSendIds,
                                              Cell<uint>                            & aActiveSendProc);

            void print_recv_communication_Ids(Cell<Cell<Cell<uint>>>  & aRecIds,
                                              Cell<uint>                            & aActiveRecProc);

            /*
             * Contains the template used for regular subdivision.
             */

//            void
//            regular_subdivision_template(Mat<uint>   & aElementNodes,
//                                         MeshXTK            & aSimpleMesh);



            //Enrichment Functions----------------------------------------------------------

            //determine_enrichment_strategy(enum Enrichment       aMethod);
            //Post Processing Functions----------------------------------------------------

            //Generate_mesh();
            void
            create_regular_subdivided_mesh_output(Cell<MeshXTK>  & aXTKMesh);

        };
    }
}


#endif  /* MORIS_XTK_MODEL_HPP_ */
