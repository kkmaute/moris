/*
 * cl_Equation_Object.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_OBJECT_HPP_
#define SRC_FEM_CL_EQUATION_OBJECT_HPP_

#include "linalg.hpp"
#include "cl_Pdof_Host.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_IWG.hpp"

#include <memory>

namespace moris
{
//FIXME will be removed soon
class Linear_Solver;
    namespace MSI
    {
    class Pdof_Host;
    class Equation_Object
    {
    private:
    moris::Cell< mtk::Vertex* >             mNodeObj;
    moris::Cell< Pdof_Host * >              mMyPdofHosts;             // Pointer to the pdof hosts of this equation object

    moris::Cell< enum Dof_Type >            mEqnObjDofTypeList;       // List of dof types of this equation obj
    moris::Mat< moris::uint >               mTimeSteps;               // List of time levels  for each dof type
    moris::Cell< Pdof* >                    mFreePdofs;               // List of the pdof pointers of this equation obj

    moris::Mat< moris::sint >               mUniqueAdofList; // Unique adof list for this equation object
    moris::map < moris::uint, moris::uint > mUniqueAdofMap;  // FIXME replace this map with an MAT. is basically used like a map right now

    // FIXME rest will be replaced
    moris::Mat< moris::real > mResidual;
    moris::Mat< moris::real > mJacobian;

    moris::Mat< moris::real > mPdofValues;

    moris::fem::Element* mElement = nullptr;

    // Integrationorder for dof types

    // dof types eg Temp

    //FIXME will be deleted soon. just for testing
    std::shared_ptr< Linear_Solver > mLin;

    //-------------------------------------------------------------------------------------------------

    public:
        Equation_Object()
        {
//            mDofType1.resize( 2, Dof_Type::TEMP );
//            mDofType1( 1 ) = Dof_Type::UX;
        };

    //-------------------------------------------------------------------------------------------------
        Equation_Object( const moris::Cell< mtk::Vertex* > & aNodeObjs ) : mNodeObj( aNodeObjs )
        {
            mTimeSteps.resize( 1, 1 );
            mTimeSteps( 0, 0 ) = 0;
        };

        Equation_Object( mtk::Cell * aCell, fem::IWG *aIWG )
        {
            // FIXME just temporary
            mElement = new fem::Element( aCell, aIWG );

            mNodeObj = mElement->get_vertex_pointers();

            mTimeSteps.set_size( 1, 1, 0 );

            mEqnObjDofTypeList.resize( 1, Dof_Type::TEMP );

            // FIXME: just temporary
            mPdofValues.set_size( mNodeObj.size(), 1, 0.0 );

            mElement->compute_jacobian_and_residual( mJacobian, mResidual, mPdofValues );

            /*
            // get number of nodes of element
            auto tNumberOfNodes = mElement->get_number_of_nodes();

            // create matrix for field
            Mat<real> tPhi_Hat( tNumberOfNodes, 1 );

            // get nodes of ekement
            Mat<real> tNodes = mElement->get_node_coords();

            // evaluate field
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
            	tPhi_Hat( k ) = std::sqrt( tNodes( k, 0 )*tNodes( k, 0 ) +  tNodes( k, 1 )*tNodes( k, 1 ) );
            } */


            //mResidual = mJacobian * tPhi_Hat;

            //mJacobian.print("J");

            //mResidual.print("r");

        }
    //-------------------------------------------------------------------------------------------------

        ~Equation_Object()
        {
            // delete element pointer if it was created
            if ( mElement != NULL )
            {
                delete mElement;
            }
        };

    //-------------------------------------------------------------------------------------------------
        /**
         * @brief Get function to get the dof types used by this equation object. This function is tested by the test [Dof_Mgn_create_unique_dof_type_list]
         * [Dof_Mgn_create_unique_dof_type_map_matrix]
         *
         * @param[in] aDofType   List of dof types.
         *
         */
        void get_dof_types( moris::Cell< enum Dof_Type > & aDofType ) { aDofType = mEqnObjDofTypeList; }

        /**
         * @brief Returns the number of nodes, elements and ghosts related to this equation object.
         *
         */
        // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
        const moris::uint get_num_pdof_hosts() { return mNodeObj.size(); }

        /**
         * @brief Returns the maximal pdof host (node) index of this equation object
         *
         */
        const moris::uint get_max_pdof_hosts_ind();

        /**
         * @brief Creates the pdof hosts of this equation object, if not created earlier, and puts them into the local pdof host list. This function is tested by the test [Eqn_Obj_create_pdof_host]
         *
         * @param[in] aNumUsedDofTypes   Number of globally used dof types
         * @param[in] aPdofTypeMap       Map which maps the dof type enum values to a consecutive list of dof type indices.
         * @param[in] aPdofHostList      List of pdof hosts.
         *
         */
        void create_my_pdof_hosts( const moris::uint                    aNumUsedDofTypes,
                                   const moris::Mat< moris::sint >    & aPdofTypeMap,
                                         moris::Cell< Pdof_Host * >   & aPdofHostList );

        /**
         * @brief This function creates a list of pdof pointers related to this equation object. This function is tested by the test [Eqn_Obj_create_my_pdof_list]
         * [Dof_Mgn_create_unique_dof_type_map_matrix]
         *
         */
        void create_my_pdof_list();

        /**
         * @brief This function creates a unique list of adofs Ids corresponding to this equation object. This function is tested by the test [Eqn_Obj_create_my_list_of_adof_ids]
         *
         */
        void create_my_list_of_adof_ids();

        /**
         * @brief This function creates a map relating the adof ids to the positions for this equation object . This function is tested by the test [Eqn_Obj_create_adof_map]
         *
         */
        void set_unique_adof_map();

        /**
         * @brief This function creates a PADofMap witch can be used to for a calculation from pdofs to adofs . This function is tested by the test [Eqn_Obj_PADofMap]
         *
         */
        void build_PADofMap( moris::Mat< moris::real > & aPADofMap );

        //-------------------------------------------------------------------------------------------------
        void get_egn_obj_jacobian( moris::Mat< moris::real > & aEqnObjMatrix )
        {
            moris::Mat< moris::real> tTMatrix;

            this->build_PADofMap( tTMatrix );

            aEqnObjMatrix = trans( tTMatrix )* mJacobian *  tTMatrix ;
        };

        //-------------------------------------------------------------------------------------------------
        void get_equation_obj_residual( moris::Mat< moris::real > & aEqnObjRHS )
        {
            moris::Mat< moris::real> tTMatrix;

            this->build_PADofMap( tTMatrix );

            aEqnObjRHS = trans( tTMatrix ) * mResidual;
        };

        void get_equation_obj_dof_ids( moris::Mat< int > & aEqnObjAdofId )
        {
            aEqnObjAdofId = mUniqueAdofList;

        };

        void get_pdof_values( Mat < real > & aValues );

        //FIXME will be deleted soon
        void set_solver( std::shared_ptr< Linear_Solver > aLin);

    };
    }
}

#endif /* SRC_FEM_CL_EQUATION_OBJECT_HPP_ */
