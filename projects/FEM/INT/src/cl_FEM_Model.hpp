/*
 * cl_FEM_Model.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_

#include "typedefs.hpp"                       //MRS/COR/src
#include "cl_Cell.hpp"                        //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_PRM_FEM_Parameters.hpp" //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_GEN_Dv_Enums.hpp"

#include "cl_MSI_Equation_Model.hpp"
#include "cl_FEM_Set_User_Info.hpp"


namespace moris
{

//------------------------------------------------------------------------------
    namespace mtk
    {
       class Mesh_Manager;
    }

    namespace fem
    {
        class IWG;
        class Node_Base;
        class Cell;
        class Set;
        class Field_Interpolator;
        class Property;
        class Constitutive_Model;
        class Stabilization_Parameter;
        class IWG;
        class IQI;
    }

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Design_Variable_Interface;
        enum class Dof_Type;
    }

    namespace fem
    {
//------------------------------------------------------------------------------

        class FEM_Model : public  MSI::Equation_Model
        {
            // pointer to reference mesh
            mtk::Mesh_Manager* mMeshManager = nullptr;
            moris_index        mMeshPairIndex;

            // list of IP node pointers
            moris::Cell< fem::Node_Base* > mIPNodes;

            // list of QI values
            moris::Cell< moris::real >       mQi;

            // parameter list to build the fem model
            moris::Cell< moris::Cell< ParameterList > > mParameterList;

            // unpacked fem inputs
            moris::Cell< fem::Set_User_Info > mSetInfo;

            // file path for property functions input
            std::string mFilePath;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aSetInfo       cell of set user info
             */
            FEM_Model
            (       mtk::Mesh_Manager                 * aMeshManager,
              const moris_index                       & aMeshPairIndex,
                    moris::Cell< fem::Set_User_Info > & aSetInfo );

//------------------------------------------------------------------------------
            /**
             * constructor with fem input
             * @param[ in ] aMesh          mesh for this problem
             * @param[ in ] aMeshPairIndex mesh pair index
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aInputFilePath a file path for property functions
             */
            FEM_Model
            (       mtk::Mesh_Manager                           * aMeshManager,
              const moris_index                                 & aMeshPairIndex,
                    moris::Cell< moris::Cell< ParameterList > > & aParameterList,
                    std::string                                   aInputFilePath );

//------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            FEM_Model(){};

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~FEM_Model();

//------------------------------------------------------------------------------
            /**
             * initialize the FEM model from parameter lists
             * @param[ in ] aParameterList a list of list of parameter lists
             * @param[ in ] aInputFilePath a file path for property functions
             */
            void initialize
            ( moris::Cell< moris::Cell< ParameterList > > & aParameterList,
              std::string                                   aInputFilePath );

//------------------------------------------------------------------------------
            /**
             * set file path
             * @param[ in ] aFilePath string for file path
             */
            void set_file_path( std::string aFilePath )
            {
                mFilePath = aFilePath;
            };

//------------------------------------------------------------------------------
            /**
             * get equation sets for test
             */
            moris::Cell< MSI::Equation_Set * > & get_equation_sets( )
            {
                return mFemSets;
            };

//------------------------------------------------------------------------------
            /**
             * get equation objects
             */
            moris::Cell< MSI::Equation_Object * > & get_equation_objects( )
            {
                return mFemClusters;
            };

//------------------------------------------------------------------------------
            /**
             * MTK set to fem set index map
             */
            map< moris_index, moris_index > & get_mesh_set_to_fem_set_index_map()
            {
                return mMeshSetToFemSetMap;
            };

//------------------------------------------------------------------------------
            /**
             * create a list of property pointers
             * param[ in ] aProperties    a list of property pointers to fill
             * param[ in ] aParameterList a list of parameter lists for properties
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_properties
            ( moris::Cell< std::shared_ptr< fem::Property > > & aProperties,
              moris::map< std::string, uint >                 & aPropertyMap,
              moris::Cell< ParameterList >                    & aParameterList,
              moris::map< std::string, MSI::Dof_Type >        & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >               & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create a list of constitutive model pointers
             * param[ in ] aCMs           a list of CM pointers to fill
             * param[ in ] aParameterList a list of parameter lists for CMs
             * param[ in ] aProperties    a list of property pointers
             *                            to assign to CMs
             * param[ in ] aPropertyMap   a map from property name to property index
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_constitutive_models
            ( moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & aCMs,
              moris::map< std::string, uint >                           & aCMMap,
              moris::Cell< ParameterList >                              & aParameterList,
              moris::Cell< std::shared_ptr< fem::Property > >           & aProperties,
              moris::map< std::string, uint >                           & aPropertyMap,
              moris::map< std::string, MSI::Dof_Type >                  & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                         & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create a list of stabilization parameter pointers
             * param[ in ] aSPs           a list of SP pointers to fill
             * param[ in ] aParameterList a list of parameter lists for SPs
             * param[ in ] aProperties    a list of property pointers
             *                            to assign to SPs
             * param[ in ] aPropertyMap   a map from property name to property
             *                            index in aProperties
             * param[ in ] aCMs           a list of CM pointers
             *                            to assign to SPs
             * param[ in ] aCMMap         a map from CM name to CM index
             *                            in aCMs
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_stabilization_parameters
            ( moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & aSPs,
              moris::map< std::string, uint >                                & aSPMap,
              moris::Cell< ParameterList >                                   & aParameterList,
              moris::Cell< std::shared_ptr< fem::Property > >                & aProperties,
              moris::map< std::string, uint >                                & aPropertyMap,
              moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      & aCMs,
              moris::map< std::string, uint >                                & aCMMap,
              moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                              & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create a list of IWG pointers
             * param[ in ] aIWGs          a list of IWG pointers to fill
             * param[ in ] aParameterList a list of parameter lists
             * param[ in ] aProperties    a list of property pointers
             *                            to assign to IWGs
             * param[ in ] aPropertyMap   a map from property name to property
             *                            index in aProperties
             * param[ in ] aCMs           a list of CM pointers
             *                            to assign to IWGs
             * param[ in ] aCMMap         a map from CM name to CM index
             *                            in aCMs
             * param[ in ] aSPs           a list of SP pointers
             *                            to assign to IWGs
             * param[ in ] aSPMap         a map from SP name to SP index
             *                            in aSPs
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_IWGs
            ( moris::Cell< std::shared_ptr< fem::IWG > >                     & aIWGs,
              moris::Cell< ParameterList >                                   & aParameterList,
              moris::Cell< std::shared_ptr< fem::Property > >                & aProperties,
              moris::map< std::string, uint >                                & aPropertyMap,
              moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      & aCMs,
              moris::map< std::string, uint >                                & aCMMap,
              moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & aSPs,
              moris::map< std::string, uint >                                & aSPMap,
              moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                              & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create an IQI
             * param[ in ] aIQIs          a list of IQI pointers to fill
             * param[ in ] aParameterList a list of parameter lists
             * param[ in ] aProperties    a list of property pointers
             *                            to assign to IQIs
             * param[ in ] aPropertyMap   a map from property name to property
             *                            index in aProperties
             * param[ in ] aCMs           a list of CM pointers
             *                            to assign to IQIs
             * param[ in ] aCMMap         a map from CM name to CM index
             *                            in aCMs
             * param[ in ] aSPs           a list of SP pointers
             *                            to assign to IQIs
             * param[ in ] aSPMap         a map from SP name to SP index
             *                            in aSPs
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_IQIs
            ( moris::Cell< std::shared_ptr< fem::IQI > >                     & aIQIs,
              moris::Cell< ParameterList >                                   & aParameterList,
              moris::Cell< std::shared_ptr< fem::Property > >                & aProperties,
              moris::map< std::string, uint >                                & aPropertyMap,
              moris::Cell< std::shared_ptr< fem::Constitutive_Model > >      & aCMs,
              moris::map< std::string, uint >                                & aCMMap,
              moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & aSPs,
              moris::map< std::string, uint >                                & aSPMap,
              moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                              & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create fem set info
             * param[ in ] aSetInfo          a list of set user info to fill
             * param[ in ] aIWGParameterList a list of IWG parameter lists
             * param[ in ] aIWGs             a list of IWG pointers
             * param[ in ] aIQIParameterList a list of IQI parameter lists
             * param[ in ] aIQIs             a list of IQI pointers
            */
            void create_fem_set_info
            ( moris::Cell< Set_User_Info >               & aSetInfo,
              moris::Cell< ParameterList >               & aIWGParameterList,
              moris::Cell< std::shared_ptr< fem::IWG > > & aIWGs,
              moris::Cell< ParameterList >               & aIQIParameterList,
              moris::Cell< std::shared_ptr< fem::IQI > > & aIQIs );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_ */
