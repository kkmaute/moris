/*
 * cl_Equation_Object.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_OBJECT_HPP_
#define SRC_FEM_CL_EQUATION_OBJECT_HPP_

#include <memory>

#include "linalg.hpp"
#include "cl_Pdof_Host.hpp"
//#include "cl_FEM_Element.hpp"
//#include "cl_FEM_IWG.hpp"
#include "cl_MSI_Node.hpp"

namespace moris
{
//FIXME will be removed soon
class Linear_Solver;
    namespace MSI
    {
    class Pdof_Host;
    class Equation_Object
    {

    protected:
    moris::Cell< MSI::Node * >   mNodeObj;
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

        //moris::fem::Element* mElement = nullptr;

        // Integrationorder for dof types

        // dof types eg Temp

        //FIXME will be deleted soon. just for testing
        std::shared_ptr< Linear_Solver > mLin;

//-------------------------------------------------------------------------------------------------
    public:
//-------------------------------------------------------------------------------------------------
        Equation_Object()
        {
//            mDofType1.resize( 2, Dof_Type::TEMP );
//            mDofType1( 1 ) = Dof_Type::UX;
        };

//-------------------------------------------------------------------------------------------------
        Equation_Object( const moris::Cell< MSI::Node* > & aNodeObjs ) : mNodeObj( aNodeObjs )
        {
            mTimeSteps.resize( 1, 1 );
            mTimeSteps( 0, 0 ) = 0;
        };

    //-------------------------------------------------------------------------------------------------

        ~Equation_Object()
        {
            /*// delete element pointer if it was created
            if ( mElement != NULL )
            {
                delete mElement;
            } */
        };

    //-------------------------------------------------------------------------------------------------
        void get_dof_types( moris::Cell< enum Dof_Type > & aDofType ) { aDofType = mEqnObjDofTypeList; }

        // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
        const moris::uint get_num_pdof_hosts() { return mNodeObj.size(); }

        const moris::uint get_max_pdof_hosts_ind();

        void create_my_pdof_hosts( const moris::uint                    aNumUsedDofTypes,
                                   const moris::Mat< moris::sint >    & aPdofTypeMap,
                                         moris::Cell< Pdof_Host * >   & aPdofHostList);

        void create_my_pdof_list();

        void create_my_list_of_adof_ids();

        void set_unique_adof_map();

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
