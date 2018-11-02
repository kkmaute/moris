/*
 * cl_Equation_Object.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_OBJECT_HPP_
#define SRC_FEM_CL_EQUATION_OBJECT_HPP_

#include <memory>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"

#include "fn_trans.hpp"
#include "op_times.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
class Dist_Vector;
    namespace fem
    {
        class Node_Base;
    }
//FIXME will be removed soon
    class Linear_Solver;
    namespace MSI
    {
    class Pdof;
    class Pdof_Host;
    class Equation_Object
    {

    protected:
    moris::Cell< fem::Node_Base * >         mNodeObj;
    moris::Cell< Pdof_Host * >              mMyPdofHosts;             // Pointer to the pdof hosts of this equation object

    moris::Cell< enum Dof_Type >            mEqnObjDofTypeList;       // List of dof types of this equation obj
    Matrix< DDUMat >                        mTimeSteps;               // List of time levels  for each dof type
    moris::Cell< Pdof* >                    mFreePdofs;               // List of the pdof pointers of this equation obj

    Matrix< DDSMat >                        mUniqueAdofList; // Unique adof list for this equation object
    moris::map < moris::uint, moris::uint > mUniqueAdofMap;  // FIXME replace this map with an MAT. is basically used like a map right now

    // FIXME rest will be replaced

    //! weak BCs of element
    Matrix< DDRMat >   mNodalWeakBCs;

    Matrix< DDRMat > mResidual;
    Matrix< DDRMat > mJacobian;

    Matrix< DDRMat > mPdofValues;

    Dist_Vector * mSolVec;

    moris::uint mEqnObjInd;

    public:
        Equation_Object() {};

//-------------------------------------------------------------------------------------------------

        Equation_Object( const moris::Cell< fem::Node_Base * > & aNodeObjs );

//-------------------------------------------------------------------------------------------------

        virtual ~Equation_Object(){};

//-------------------------------------------------------------------------------------------------
        /**
         * @brief Get function to get the dof types used by this equation object. This function is tested by the test [Dof_Mgn_create_unique_dof_type_list]
         * [Dof_Mgn_create_unique_dof_type_map_matrix]
         *
         * @param[in] aDofType   List of dof types.
         *
         */
        void get_dof_types( moris::Cell< enum Dof_Type > & aDofType ) { aDofType = mEqnObjDofTypeList; }
//-------------------------------------------------------------------------------------------------
        /**
         * @brief Returns the number of nodes, elements and ghosts related to this equation object.
         *
         */
        // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
        moris::uint get_num_pdof_hosts() { return mNodeObj.size(); }

//------------------------------------------------------------------------------------------------
        /**
         * @brief Returns the maximal pdof host (node) index of this equation object
         *
         */
        moris::uint get_max_pdof_hosts_ind();

//-------------------------------------------------------------------------------------------------
        /**
         * @brief Creates the pdof hosts of this equation object, if not created earlier, and puts them into the local pdof host list. This function is tested by the test [Eqn_Obj_create_pdof_host]
         *
         * @param[in] aNumUsedDofTypes   Number of globally used dof types
         * @param[in] aPdofTypeMap       Map which maps the dof type enum values to a consecutive list of dof type indices.
         * @param[in] aPdofHostList      List of pdof hosts.
         *
         */
        void create_my_pdof_hosts( const moris::uint                  aNumUsedDofTypes,
                                   const Matrix< DDSMat >           & aPdofTypeMap,
                                         moris::Cell< Pdof_Host * > & aPdofHostList );

//-------------------------------------------------------------------------------------------------
        /**
         * @brief This function creates a list of pdof pointers related to this equation object. This function is tested by the test [Eqn_Obj_create_my_pdof_list]
         * [Dof_Mgn_create_unique_dof_type_map_matrix]
         *
         */
        void create_my_pdof_list();

//-------------------------------------------------------------------------------------------------
        /**
         * @brief This function creates a unique list of adofs Ids corresponding to this equation object. This function is tested by the test [Eqn_Obj_create_my_list_of_adof_ids]
         *
         */
        void create_my_list_of_adof_ids();

//-------------------------------------------------------------------------------------------------
        /**
         * @brief This function creates a map relating the adof ids to the positions for this equation object . This function is tested by the test [Eqn_Obj_create_adof_map]
         *
         */
        void set_unique_adof_map();

//-------------------------------------------------------------------------------------------------
        /**
         * @brief This function creates a PADofMap witch can be used to for a calculation from pdofs to adofs . This function is tested by the test [Eqn_Obj_PADofMap]
         *
         */
        void build_PADofMap( Matrix< DDRMat > & aPADofMap );

//-------------------------------------------------------------------------------------------------
        void get_egn_obj_jacobian( Matrix< DDRMat > & aEqnObjMatrix )
        {
            Matrix< DDRMat > tTMatrix;

            this->build_PADofMap( tTMatrix );

            aEqnObjMatrix = trans( tTMatrix ) * mJacobian *  tTMatrix ;
        };

//-------------------------------------------------------------------------------------------------
        void get_equation_obj_residual( Matrix< DDRMat > & aEqnObjRHS, Dist_Vector * aSolutionVector )
        {
            mSolVec = aSolutionVector;

            this->compute_jacobian_and_residual();

            Matrix< DDRMat > tTMatrix;

            this->build_PADofMap( tTMatrix );

            aEqnObjRHS = trans( tTMatrix ) * mResidual;
        };

//-------------------------------------------------------------------------------------------------
        void get_equation_obj_dof_ids( Matrix< DDSMat > & aEqnObjAdofId )
        {
            aEqnObjAdofId = mUniqueAdofList;
        };
//-------------------------------------------------------------------------------------------------

        /**
         * returns a moris::Mat with indices of vertices that are connected to this element
         */
        moris_index
        get_node_index( const moris_index aElementLocalNodeIndex ) const ;

//-------------------------------------------------------------------------------------------------
        // void get_pdof_values( Mat < real > & aValues );
//        void
//        extract_values( std::shared_ptr< Linear_Solver > aLin );

        //void get_pdof_values( std::shared_ptr< Linear_Solver > aLin );

//-------------------------------------------------------------------------------------------------

        //void get_adof_values( Mat < real > & aValues );

//-------------------------------------------------------------------------------------------------
        virtual Matrix< DDSMat > get_adof_indices()
        {
            MORIS_ERROR( false, "this function does nothing");

            return Matrix< DDSMat >(0,0);
        }

//-------------------------------------------------------------------------------------------------
        //FIXME will be deleted soon
        void set_solver( std::shared_ptr< Linear_Solver > aLin);

//-------------------------------------------------------------------------------------------------
        virtual void compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, "this function does nothing");
        }

//-------------------------------------------------------------------------------------------------

        virtual moris::real compute_integration_error( moris::real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
        {
            MORIS_ERROR( false, "this function does nothing");
            return 0.0;
        }

//-------------------------------------------------------------------------------------------------

        /**
         * retrun Neumann boundary conditions, writable version
         */
        Matrix< DDRMat > &
        get_weak_bcs()
        {
            return mNodalWeakBCs;
        }

//-------------------------------------------------------------------------------------------------

        /**
         * retrun Neumann boundary conditions, const version
         */
        const Matrix< DDRMat > &
        get_weak_bcs() const
        {
            return mNodalWeakBCs;
        }

//-------------------------------------------------------------------------------------------------

        /**
         * how many nodes are connected to this element
         */
        uint
        get_num_nodes() const
        {
            return mNodeObj.size();
        }

//-------------------------------------------------------------------------------------------------
    };
    }
}

#endif /* SRC_FEM_CL_EQUATION_OBJECT_HPP_ */
