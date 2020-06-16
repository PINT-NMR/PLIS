#ifndef TABLE_H
#define TABLE_H

#include "qcustomplot.h"
#include "data.h"

class Table
{
public:
    Table(QTableWidget *_table);
    QTableWidget *table;

    void addRow(QVector<Data*> &dataSets, int current_datSet, QString mode);
    void write_In_Table(QVector<Data*> &dataSets, int current_dataSet, QString mode);
    void removeRow(QVector<Data*> &dataSets, int current_dataSet, QString mode);
    void cellChanged(QVector<Data*> &dataSets, int current_dataSet, QString mode, int row, int column);
    void paste_data_to_table(QVector<Data*> &dataSets, int current_dataSet);

};

#endif // TABLE_H
