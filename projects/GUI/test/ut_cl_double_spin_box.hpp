#pragma once

#include <QObject>
#include <QtTest>

namespace moris
{
    class Parameter;
    class Moris_Double_Spin_Box;
}

class ut_cl_double_spin_box : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void testInitialSigned();
    void testInitialUnsigned();
    void testvaluechanged();
    void testReadOnlyParameter();
};
