#include "legend.h"

void legend_clicked(QCPLegend*,QCPAbstractLegendItem* ai,QMouseEvent*, QVector<Data*> &dataSets, QVector<QCPGraph*> &graphSets){
    if(ai!=nullptr){
        for (int i=0;i< dataSets.size();i++){
            if(dataSets[i]->name== static_cast<QCPPlottableLegendItem*>(ai)->plottable()->name()){
                if(ai->textColor()=="Grey"){
                    ai->setTextColor("black");
                    static_cast<QCPPlottableLegendItem*>(ai)->plottable()->setVisible(true);
                    dataSets[i]->data_visible=true;
                }
                else{
                    ai->setTextColor("Grey");
                    static_cast<QCPPlottableLegendItem*>(ai)->plottable()->setVisible(false);
                    dataSets[i]->data_visible=false;
                }
            }
        }
        for (int i=0;i<graphSets.size();i++){
            if (graphSets[i]!=nullptr)
                if(graphSets[i]->name() == static_cast<QCPPlottableLegendItem*>(ai)->plottable()->name()){
                    if(ai->textColor()=="Grey"){
                        ai->setTextColor("black");
                        static_cast<QCPPlottableLegendItem*>(ai)->plottable()->setVisible(true);
                        dataSets[i]->curve_visible=true;
                    }
                    else{
                        ai->setTextColor("Grey");
                        static_cast<QCPPlottableLegendItem*>(ai)->plottable()->setVisible(false);
                        dataSets[i]->curve_visible=false;
                    }
                }
        }
    }
}

QString legend_double_clicked(QCPLegend*,QCPAbstractLegendItem* ai,QMouseEvent*, QVector<Data*> &dataSets, QString new_name){
    for (int i=0;i< dataSets.size();i++){
        if(dataSets[i]->name== static_cast<QCPPlottableLegendItem*>(ai)->plottable()->name()){
            QString old_name=dataSets[i]->name;
            dataSets[i]->name=new_name;
            return old_name;
        }
    }
    return "";
}
