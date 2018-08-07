/*
 * cl_HMR_MTK.hpp
 *
 *  Created on: May 16, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_MTK_HPP_
#define SRC_HMR_CL_HMR_MTK_HPP_

#include <string>


#include "typedefs.hpp" //COR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_HMR_Basis.hpp"
#include "cl_Database.hpp" //MTK/src
#include "cl_Mesh_Enums.hpp" //MTK/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src


namespace moris
{
    namespace hmr
    {

// ----------------------------------------------------------------------------
    /**
     * \brief   This class provides an interface to the moris::MTK module.
     *
     * MTK does not support luint, therefore, all luint values must be casted
     * to uint. Node Coordinates and the Element connectivity need to be
     * transposed for MTK to understand.
     */
    class MTK
    {
        //! pointer to container of user defined settings
        const Parameters      * mParameters;

        //! struc required by MTK
        MtkMeshData           mMeshData;

        //! struc required by MTK
        MtkFieldsInfo         mFieldsInfo;

        //! node and element data passed to MTK
        Cell< Mat< real > >   mFieldData;

        //! 1D, 2D or 3D
        uint                  mNumberOfDimensions;

        //! connectivity passed to MTK
        Mat< uint >           mElementTopology;

        //! element IDs passed to MTK
        Mat< uint >           mElementLocalToGlobal;

        //! node IDs passed to MTK
        Mat< uint >           mNodeLocalToGlobal;

        //! node coordinates passed to MTK
        Mat< real >           mNodeCoords;

        //! node owners passed to MTK
        Mat< uint >           mNodeOwner;

        //! container for all elements on proc
        Cell< Element* > & mAllElementsOnProc;

        //! container for all nodes on proc
        Cell< Basis* > &      mAllBasisOnProc;

        //! switch telling if debug information is to be printed
        bool                  mVerboseFlag = true;

// ----------------------------------------------------------------------------
    public:
// ----------------------------------------------------------------------------
            /**
             * Default constructor of MTK object
             *
             * @param[in] aParameters ref to container of user defined settings
             */
            MTK( const Parameters * aParameters,
                    Cell< Element* > & aAllElementsOnProc,
                    Cell< Basis* >    & aAllNodesOnProc ) :
                    mParameters( aParameters ),
                    mAllElementsOnProc( aAllElementsOnProc ),
                    mAllBasisOnProc( aAllNodesOnProc ),
                    mVerboseFlag( aParameters->is_verbose() )
             {};

// ----------------------------------------------------------------------------

            /**
             * Default MTK destructor. Does nothing
             */
            ~MTK(){};

// ----------------------------------------------------------------------------

            /**
             * Converts the data passed from Lagrange_Mesh to a format
             * MTK can understand.
             *
             *
             * @return void
             */
            void
            create_mesh_data(
                    const luint                    & aNumberOfElements,
                    const uint                     & aNumberOfNodesPerElement,
                    const luint                    & aNumberOfNodes) ;

// ----------------------------------------------------------------------------

            void
            add_node_data(
                    const std::string & aLabel,
                    const Mat< real > & aData );

// ----------------------------------------------------------------------------

            void
            add_element_data(
                    const std::string & aLabel,
                    const Mat< real > & aData );

// ----------------------------------------------------------------------------

            /**
             * Tells MTK to dump the calculated mesh data into an output file
             *
             * @return void
             */
            void
            save_to_file( const std::string & aFilePath );

    };
// ----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_MTK_HPP_ */
