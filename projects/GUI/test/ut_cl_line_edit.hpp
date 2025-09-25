#pragma once

#include <QObject>
#include <QtTest>

namespace moris
{
    class Parameter;
    class Moris_Line_Edit;
}




class ut_cl_line_edit : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void testTextChanged();
};