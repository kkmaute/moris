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
#include "fn_Exec_load_user_library.hpp"


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

//            // file path for property functions input
//            std::string mFilePath;

            // fixme remove ?
            moris::Cell< std::shared_ptr< fem::Property > > mProperties;
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mCMs;
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > mSPs;
            moris::Cell< std::shared_ptr< fem::IWG > > mIWGs;
            moris::Cell< std::shared_ptr< fem::IQI > > mIQIs;

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
             * @param[ in ] aLibrary       a file path for property functions
             */
            FEM_Model
            (       mtk::Mesh_Manager                           * aMeshManager,
              const moris_index                                 & aMeshPairIndex,
                    moris::Cell< moris::Cell< ParameterList > >   aParameterList,
                    std::shared_ptr< Library_IO >                 aLibrary );


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
             * @param[ in ] aLibrary       a file path for property functions
             */
            void initialize( std::shared_ptr< Library_IO > aLibrary );

//------------------------------------------------------------------------------
            /**
             * set parameter list
             * @param[ in ] aParameterList a list of parameetr for the FEM model
             */
            void set_parameter_list( moris::Cell< moris::Cell< ParameterList > > aParameterList )
            {
                mParameterList = aParameterList;
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
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             * @param[ in ] aLibrary       a file path for property functions
             */
            void create_properties
            ( moris::map< std::string, uint >                 & aPropertyMap,
              moris::map< std::string, MSI::Dof_Type >        & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >               & aDvTypeMap,
              std::shared_ptr< Library_IO >                     aLibrary );

//------------------------------------------------------------------------------
            /**
             * create a list of constitutive model pointers
             * param[ in ] aCMs           a list of CM pointers to fill
             * param[ in ] aProperties    a list of property pointers
             *                            to assign to CMs
             * param[ in ] aPropertyMap   a map from property name to property index
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_constitutive_models
            ( moris::map< std::string, uint >                           & aCMMap,
              moris::map< std::string, uint >                           & aPropertyMap,
              moris::map< std::string, MSI::Dof_Type >                  & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                         & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create a list of stabilization parameter pointers
             * param[ in ] aPropertyMap   a map from property name to property
             *                            index in aProperties
             * param[ in ] aCMMap         a map from CM name to CM index
             *                            in aCMs
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_stabilization_parameters
            ( moris::map< std::string, uint >                                & aSPMap,
              moris::map< std::string, uint >                                & aPropertyMap,
              moris::map< std::string, uint >                                & aCMMap,
              moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                              & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create a list of IWG pointers
             * param[ in ] aPropertyMap   a map from property name to property
             *                            index in aProperties
             * param[ in ] aCMMap         a map from CM name to CM index
             *                            in aCMs
             * param[ in ] aSPMap         a map from SP name to SP index
             *                            in aSPs
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_IWGs
            ( moris::map< std::string, uint >                                & aPropertyMap,
              moris::map< std::string, uint >                                & aCMMap,
              moris::map< std::string, uint >                                & aSPMap,
              moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                              & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create an IQI
             * param[ in ] aPropertyMap   a map from property name to property
             *                            index in aProperties
             * param[ in ] aCMMap         a map from CM name to CM index
             *                            in aCMs
             * param[ in ] aSPMap         a map from SP name to SP index
             *                            in aSPs
             * param[ in ] aMSIDofTypeMap a map from std::string to MSI::Dof_Type
             * param[ in ] aDvTypeMap     a map from std::string to GEN_DV
             */
            void create_IQIs
            ( moris::map< std::string, uint >                                & aPropertyMap,
              moris::map< std::string, uint >                                & aCMMap,
              moris::map< std::string, uint >                                & aSPMap,
              moris::map< std::string, MSI::Dof_Type >                       & aMSIDofTypeMap,
              moris::map< std::string, GEN_DV >                              & aDvTypeMap );

//------------------------------------------------------------------------------
            /**
             * create fem set info
            */
            void create_fem_set_info();

//------------------------------------------------------------------------------

        };
//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */


#endif /* PROJECTS_FEM_MDL_SRC_CL_FEM_MODEL_HPP_ */
