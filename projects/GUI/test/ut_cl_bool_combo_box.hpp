#pragma once

#include <QObject>
#include <QtTest>

namespace moris
{
    class Parameter;
    class Moris_Combo_Box;
}

class ut_cl_bool_combo_box : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void testInitialTrue();
    void testInitialFalse();
    void testIndexChanged();
    void testReadOnlyParameter();

};