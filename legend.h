#ifndef LEGEND_H
#define LEGEND_H

#include "plis.h"

void legend_clicked(QCPLegend*,QCPAbstractLegendItem* ai,QMouseEvent*, QVector<Data*> &dataSets, QVector<QCPGraph*> &graphSets);

QString legend_double_clicked(QCPLegend*,QCPAbstractLegendItem* ai,QMouseEvent*, QVector<Data*> &dataSets, QString new_name);


#endif // LEGEND_H
