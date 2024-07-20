#ifndef MORIS_PAIR_BOX_HPP
#define MORIS_PAIR_BOX_HPP

#include "cl_Parameter.hpp"

#include <QWidget>
#include <QComboBox>
#include <QLineEdit>
#include <QStringList>
#include <QHBoxLayout>

namespace moris {
    class Parameter;
}

class Moris_Pair_Box : public QWidget {
    Q_OBJECT

public:
    // Constructor
    explicit Moris_Pair_Box(QWidget *parent, moris::Parameter &param, const QStringList &options = {});

    QComboBox *morisPairComboBox; 
    QLineEdit *morisPairLineEdit; 

signals:
    // Public signals
    void comboBoxTextChanged(const QString &name, const QString &text);
    void lineEditTextChanged(const QString &name, const QString &text);

public slots:
    // Private slots for handling changes
    void onComboBoxTextChanged(const QString &text);
    void onLineEditTextChanged(const QString &text);

private:
    moris::Parameter &parameter; // Parameter object
};

#endif // MORIS_PAIR_BOX_HPP
