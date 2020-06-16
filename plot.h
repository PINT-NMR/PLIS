#ifndef PLOT_H
#define PLOT_H

#include "data.h"
#include "plis.h"

void initPlotArea(QCustomPlot *plot_area, QString xlabel, QString ylabel, QCPTextElement *title);
void updatePlotArea(QCustomPlot *plot_area, QString xlabel, QString ylabel, QFont xaxis_font, QFont yaxis_font);
void setGraphStyle(QCustomPlot *plot_area, Data *dataSet, QString mode);
void clearAndReplot(QCustomPlot *plot_area);

#endif // PLOT_H
