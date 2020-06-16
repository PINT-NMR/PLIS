#include "plot.h"


void initPlotArea(QCustomPlot *plot_area, QString xlabel, QString ylabel, QCPTextElement *title)
{
    plot_area->xAxis->setLabel(xlabel);
    plot_area->yAxis->setLabel(ylabel);
    plot_area->xAxis->setRange(0,100);
    plot_area->yAxis->setRange(0,500);
    plot_area->legend->setVisible(true);
    plot_area->legend->setFont(QFont("Helvetica",9));
    plot_area->legend->setRowSpacing(-3);
    plot_area->legend->setSelectableParts(QCPLegend::spItems);
    plot_area->plotLayout()->insertRow(0);
    plot_area->plotLayout()->addElement(0,0, title);
}

void updatePlotArea(QCustomPlot *plot_area, QString xlabel, QString ylabel, QFont xaxis_font, QFont yaxis_font)
{
    plot_area->rescaleAxes();
    plot_area->xAxis->setLabel(xlabel);
    plot_area->xAxis->setLabelFont(xaxis_font);
    plot_area->xAxis->setSelectedLabelFont(xaxis_font);
    plot_area->yAxis->setLabel(ylabel);
    plot_area->yAxis->setLabelFont(yaxis_font);
    plot_area->yAxis->setSelectedLabelFont(yaxis_font);
    plot_area->replot();
}

void setGraphStyle(QCustomPlot *plot_area, Data *dataSet, QString mode)
{
    plot_area->addGraph();
    plot_area->graph()->setName(dataSet->name);
    plot_area->graph()->setPen(dataSet->pen);
    plot_area->graph()->setScatterStyle(dataSet->style);
    plot_area->graph()->setLineStyle(QCPGraph::lsNone);
    if (mode=="standard" || mode=="chemshift")
        plot_area->graph()->setData(dataSet->concVector,dataSet->responceVector);
    else if(mode=="cpmg")
        plot_area->graph()->setData(dataSet->n_cpmgVector,dataSet->R2effVector);
}

void clearAndReplot(QCustomPlot *plot_area) {
    plot_area->clearGraphs();
    plot_area->replot();
}
