/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_STK.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_STK_HPP_
#define SRC_HMR_CL_HMR_STK_HPP_

#include <string>

#include "cl_HMR_Basis.hpp"
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Fields_Info.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

namespace moris::hmr
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
        mtk::MtkMeshData      mMeshData;

        //! struc required by MTK
        mtk::MtkFieldsInfo    mFieldsInfo;

        // Fields the field info contain references
        moris::Cell< mtk::Scalar_Field_Info<DDRMat> > mRealScalarFields;
        moris::Cell< mtk::Scalar_Field_Info<DDSMat> > mSintScalarFields;

        //! struc for sidesets
        mtk::MtkSetsInfo      mSetsInfo;

        //! 2D or 3D
        uint                  mNumberOfDimensions;

        //! node and element data passed to MTK
        Cell< Matrix< DDRMat > >   mRealScalarFieldData;

        //! connectivity passed to MTK
        Matrix< IdMat >           mElementTopology;

        //! element IDs passed to MTK
        Matrix< IdMat >           mElementLocalToGlobal;

        //! node IDs passed to MTK
        Matrix< IdMat >           mNodeLocalToGlobal;

        //! node coordinates passed to MTK
        Matrix< DDRMat >          mNodeCoords;

        //! node owners passed to MTK
        Matrix< IdMat >           mNodeOwner;

        //! node sharing passed to MTK
        Matrix< IdMat >           mNodeSharing;

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
            void create_mesh_data( const double aTimeStep = 0.0 ) ;

    // ----------------------------------------------------------------------------

         /*   void
            add_node_data(
                    const std::string & aLabel,
                    const Matrix< DDRMat > & aData );

    // ----------------------------------------------------------------------------

            void
            add_element_data(
                    const std::string & aLabel,
                    const Matrix< DDRMat > & aData ); */

    // ----------------------------------------------------------------------------

            /**
             * Tells MTK to dump the calculated mesh data into an output file
             *
             * @return void
             */
            void save_to_file( const std::string & aFilePath );

    // ----------------------------------------------------------------------------

            /**
             * Return node sharing data. Needed for unit test
             */
            Matrix< IdMat > const & get_node_sharing_data()
            {
                return mNodeSharing;
            }

    // ----------------------------------------------------------------------------
    private:
            void flag_old_and_new_elements();

    };
    // ----------------------------------------------------------------------------

} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_STK_HPP_ */

