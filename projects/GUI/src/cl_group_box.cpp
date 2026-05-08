#include "cl_group_box.hpp"

namespace moris
{

    Moris_Group_Box::Moris_Group_Box( QWidget *a_parent, Parameter &a_param,const QStringList &a_options )
            : QWidget( a_parent )
            , mParameter( &a_param )
            , mOptions( a_options )
            , mFormLayout( new QFormLayout() )
    {
        setLayout( mFormLayout );
        refresh_data_parameter();
        }
        

        Parameter &Moris_Group_Box::get_parameter()
        {
            MORIS_ERROR(mParameter, "Moris_Group_Box::get_parameter() called with null mParameter.");
            return *mParameter;
        }
        

        void Moris_Group_Box::setParameter( Parameter &a_parameter )
        {
            mParameter = &a_parameter;
            refresh_data_parameter();
        }
        

        void Moris_Group_Box::set_property_list(const QStringList &a_options)
        {
            mOptions = a_options;
            refresh_data_parameter();
        }
        

        void Moris_Group_Box::clear_rows()
        {
            while ( mFormLayout && mFormLayout->rowCount() > 0 )
            {
                mFormLayout->removeRow( 0 );
            }
            mWidget.clear();
        }
        

        void Moris_Group_Box::add_row( const std::string &a_key, const QString &a_selectedText)
        {
            QComboBox *tComboBox = new QComboBox( this );
            {
                QSignalBlocker blocker( tComboBox );
                tComboBox->addItems( mOptions );

                if (!a_selectedText.isEmpty() )
                {
                    int tIndex = tComboBox->findText( a_selectedText );
                    if ( tIndex >= 0 )
                    {
                        tComboBox->setCurrentIndex( tIndex );
                    }
                    else{
                        tComboBox->addItem( a_selectedText );
                        tComboBox->setCurrentIndex( tComboBox->count() - 1 );
                    }
                }
                else if(tComboBox->count() > 0)
                {
                    tComboBox->setCurrentIndex( 0 );
                }
            }

            if ( mParameter &&mParameter->is_locked() )
            {
                tComboBox->setDisabled( true );
            }

            connect(
                    tComboBox, 
                    QOverload< int>::of(&QComboBox::currentIndexChanged), 
                    this, 
                    &Moris_Group_Box::on_property_selection_changed );
            mWidget[ a_key ] = tComboBox;
            mFormLayout->addRow( QString::fromStdString( a_key ), tComboBox );
        }

        
        void Moris_Group_Box::build_default_rows()
        {
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCM  = 
                    tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );

            std::map<std::string, uint> tMap = tCM->get_property_map();

            for ( const auto &iMap : tMap )
            {
                add_row (iMap.first);
            }
        }
        

        void Moris_Group_Box::build_rows_from_serialized( const std::string &a_serialized )
        {
            std::stringstream tStream( a_serialized );
            std::string       tPair;
            while ( std::getline( tStream, tPair, ';' ) )
            {
                std::stringstream tPairStream( tPair );
                std::string       tKey, tValue;

                // stored format is value,key
                if ( std::getline( tPairStream, tValue, ',' ) && std::getline( tPairStream, tKey ) )
                {
                    add_row( tKey, QString::fromStdString( tValue ) );
                }
            }
        }
        

        std::string Moris_Group_Box::serialize_current_state() const
        {
            std::string tResult;

            for ( auto it = mWidget.begin(); it != mWidget.end(); ++it )
            {
                // Append the value
                tResult += it->second->currentText().toStdString();
                tResult += ",";

                // Append the text from the key
                tResult += it->first;

                // If it's not the last element, add a semicolon
                if ( std::next( it ) != mWidget.end() )
                {
                    tResult += ";";
                }
            }

            return tResult;
        }
        

        void Moris_Group_Box::refresh_data_parameter()
        {
            clear_rows();

            if ( !mParameter )
            {
                setDisabled( true );
                return;
            }

            setDisabled( false );

            const std::string tSerialized = mParameter->get_string();

            if ( tSerialized.empty() )
            {
                build_default_rows();
            }
            else
            {
                build_rows_from_serialized( tSerialized );
            }

            if ( mParameter->is_locked() )
            {
                setDisabled( true );
            }
        }
        

        void Moris_Group_Box::on_property_selection_changed(int a_index )
        {
            Q_UNUSED (a_index);

            if ( !mParameter || mParameter->is_locked() )
            {
                return;
            }

            mParameter->set_value( objectName().toStdString(), 
                                    serialize_current_state(), 
                                    false );

        }
        
        
        void Moris_Group_Box::on_combo_box_selection_changed( int a_index )
        {
            Q_UNUSED (a_index);

            if ( !mParameter || mParameter->is_locked() )
            {
                return;
            }

            clear_rows();

            fem::CM_Factory tCMFactory;
            std::shared_ptr< fem::Constitutive_Model > tCM  = 
                    tCMFactory.create_CM( (fem::Constitutive_Type)a_index);

            std::map<std::string, uint> tMap = tCM->get_property_map();

            for ( const auto &iMap : tMap )
            {
                add_row (iMap.first);
            }

            mParameter->set_value( objectName().toStdString(), 
                                    serialize_current_state(), 
                                    false );
        }
    }
