/*
 * cl_HMR_MTK.hpp
 *
 *  Created on: May 16, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_STK_HPP_
#define SRC_HMR_CL_HMR_STK_HPP_

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

    // forward declaration of Lagrange mesh
    class Lagrange_Mesh_Base;

// ----------------------------------------------------------------------------
    /**
     * \brief   This class provides an interface to the STK writer module.
     *
     * MTK does not support luint, therefore, all luint values must be casted
     * to uint. Node Coordinates and the Element connectivity need to be
     * transposed for MTK to understand.
     */
    class STK
    {
        //! pointer to container of user defined settings
        const Parameters      * mParameters;

        Lagrange_Mesh_Base *  mMesh;

        //! struc required by MTK
        MtkMeshData           mMeshData;

        //! struc required by MTK
        MtkFieldsInfo         mFieldsInfo;

        //! 2D or 3D
        uint                  mNumberOfDimensions;

        //! node and element data passed to MTK
        Cell< Mat< real > >   mFieldData;

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

// ----------------------------------------------------------------------------
    public:
// ----------------------------------------------------------------------------
            /**
             * Default constructor of MTK object
             *
             * @param[in] aMesh Lagrange Mesh object this MTK refers to
             */
             STK( Lagrange_Mesh_Base * aMesh );

// ----------------------------------------------------------------------------

            /**
             * Default MTK destructor. Does nothing
             */
            ~STK(){};

// ----------------------------------------------------------------------------

            /**
             * Converts the data passed from Lagrange_Mesh to a format
             * MTK can understand.
             *
             *
             * @return void
             */
            void
            create_mesh_data() ;

// ----------------------------------------------------------------------------

         /*   void
            add_node_data(
                    const std::string & aLabel,
                    const Mat< real > & aData );

// ----------------------------------------------------------------------------

            void
            add_element_data(
                    const std::string & aLabel,
                    const Mat< real > & aData ); */

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

#endif /* SRC_HMR_CL_HMR_STK_HPP_ */
