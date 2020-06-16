#include "table.h"
#include <iostream>
#include <sstream>
#include <string>

// 200601
// - Creates this class and moved content from Plis to it.


Table::Table(QTableWidget *_table) : table(_table) {}


void Table::addRow(QVector<Data*> &dataSets, int current_dataSet, QString mode) {
    if (table->currentRow()==-1)
        table->setCurrentCell(0,0);
    if (mode=="standard" || mode=="chemshift") {
        table->insertRow(table->currentRow());
        dataSets[current_dataSet]->insertData(table->currentRow(), 0., 0., 0., 0., 1.e-4);
        dataSets[current_dataSet]->update_protein_conc();
    }
    else {
        dataSets[current_dataSet]->insertCPMGdata(table->currentRow(), 0., 0., 0.);
        table->insertRow(table->currentRow());
    }
    write_In_Table(dataSets, current_dataSet, mode);
}

void Table::write_In_Table(QVector<Data*> &dataSets, int current_dataSet, QString mode) {
    if(mode=="standard" || mode=="chemshift"){
        table->setRowCount(dataSets[current_dataSet]->responceVector.size());
        for (int i=0;i<3;i++){
            for (int j=0; j<dataSets[current_dataSet]->responceVector.size(); j++){
                QTableWidgetItem *theItem = new QTableWidgetItem();
                if (i==0)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->responceVector.at(j));
                else if(i==1)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->concVector.at(j));
                else
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->volumeVector.at(j));
                table->setItem(j,i,theItem);
            }
        }
    }
    else if(mode=="cpmg"){
        table->setRowCount(dataSets[current_dataSet]->n_cpmgVector.size());
        for (int i=0;i<3;i++){
            for (int j=0; j<dataSets[current_dataSet]->n_cpmgVector.size(); j++){
                QTableWidgetItem *theItem = new QTableWidgetItem();
                if (i==0)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->n_cpmgVector.at(j));
                else if(i==1)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->R2effVector.at(j));
                else
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->dyVector.at(j));
                table->setItem(j,i,theItem);
            }
        }
    }
}

void Table::removeRow(QVector<Data*> &dataSets, int current_dataSet, QString mode)
{
    if (table->currentRow()==-1)
        table->setCurrentCell(0,0);
    int row{table->currentRow()},col{table->currentColumn()};
    if (dataSets[current_dataSet]->responceVector.size()!=1 && (mode=="standard" || mode=="chemshift")){
        dataSets[current_dataSet]->removeData(table->currentRow());
        dataSets[current_dataSet]->update_protein_conc();
        table->removeRow(table->currentRow());
    }
    else if(dataSets[current_dataSet]->n_cpmgVector.size()!=1 && mode=="cpmg"){
        dataSets[current_dataSet]->removeCPMGdata(table->currentRow());
        table->removeRow(table->currentRow());
    }
    write_In_Table(dataSets, current_dataSet, mode);
    table->setCurrentCell(row,col);
}

void Table::cellChanged(QVector<Data*> &dataSets, int current_dataSet, QString mode, int row, int column)
{
    QString text{table->item(row,column)->text()};
    double changed_value{text.toDouble()};
    if(mode=="standard" || mode=="chemshift") {
        if (column==0)
            dataSets[current_dataSet]->responceVector[row]=changed_value;
        else if (column==1)
            dataSets[current_dataSet]->concVector[row]=changed_value;
        else {
            dataSets[current_dataSet]->volumeVector[row]=changed_value;
            dataSets[current_dataSet]->update_protein_conc();
        }
    }
    else if(mode=="cpmg"){
        if (column==0)
            dataSets[current_dataSet]->n_cpmgVector[row]=changed_value;
        else if (column==1)
            dataSets[current_dataSet]->R2effVector[row]=changed_value;
        else
            dataSets[current_dataSet]->dyVector[row]=changed_value;
    }
}

void Table::paste_data_to_table(QVector<Data*> &dataSets, int current_dataSet)
{
    std::istringstream iss{QApplication::clipboard()->text().toStdString()};
    double test;
    std::string line;
    int row{table->currentRow()};
    int start_column{table->currentColumn()};
    int column=start_column;
    while (!iss.eof()) {
        std::getline(iss,line);
        std::replace(line.begin(), line.end(), ',','.');
        std::istringstream iss2{line};
        while (!iss2.eof()) {
            if (column>2)
                break;
            QTableWidgetItem *theItem = new QTableWidgetItem();
            iss2 >> test;
            if (row>dataSets[current_dataSet]->responceVector.size()-1) {
                dataSets[current_dataSet]->pushBackData(0., 0., 0., 0., 1.e-4);
                if (column==0)
                    dataSets[current_dataSet]->responceVector[row]=test;
                else if(column==1)
                    dataSets[current_dataSet]->concVector[row]=test;
                else if(column==2)
                    dataSets[current_dataSet]->volumeVector[row]=test;
                table->setRowCount(dataSets[current_dataSet]->responceVector.size());
            }
            theItem->setData(Qt::EditRole, test);
            table->setItem(row,column, theItem);
            column++;
        }
        row++;
        column=start_column;
    }
    dataSets[current_dataSet]->update_protein_conc();
}

