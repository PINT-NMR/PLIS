#include "optimize.h"
#include "ui_optimize.h"

Optimize::Optimize(QWidget *parent, Data* &d) :
    QDialog(parent), data(d),
    ui(new Ui::Optimize)
{
    ui->setupUi(this);
    ui->plot_area->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
                                    QCP::iSelectLegend | QCP::iSelectPlottables | QCP::iSelectItems);

    ui->kd2->hide();
    ui->kd2_label->hide();
    ui->kdc->hide();
    ui->kdc_label->hide();
    ui->kd3->hide();
    ui->kd3_label->hide();
    ui->intmed1->hide();
    ui->intmed_label->hide();
    ui->kd2_check->hide();
    ui->kd3_check->hide();
    ui->intmed_check->hide();
    ui->Kdc_check->hide();

    ui->chi2->setText(QString::number(data->chi2));

    ui->plot_area->addGraph();
    if(data->num_bind_site==4)
        ui->plot_area->graph()->setData(data->n_cpmgVector,data->R2effVector);
    else
        ui->plot_area->graph()->setData(data->concVector,data->responceVector);
    ui->plot_area->graph()->setScatterStyle(data->style);
    ui->plot_area->graph()->setLineStyle(QCPGraph::lsNone);
    ui->plot_area->graph()->setPen(data->pen);
    ui->plot_area->addGraph();
    ui->plot_area->graph()->setData(data->calc_data[0],data->calc_data[1]);
    ui->plot_area->graph()->setPen(data->pen_curve);
    ui->plot_area->graph()->setLineStyle(QCPGraph::lsLine);
    ui->plot_area->replot();
    ui->plot_area->rescaleAxes();

    ui->start->setValue(data->a[0]);
    if(data->num_bind_site==1){
        ui->kd1->setValue(data->a[2]);
        ui->end->setValue(data->a[1]);
    }
    else if(data->num_bind_site==2){
        ui->kd2->show();
        ui->kd2_label->show();
        ui->kd2->setValue(data->a[4]);
        ui->kd1->setValue(data->a[3]);
        ui->intmed1->show();
        ui->intmed_label->show();
        ui->intmed1->setValue(data->a[1]);
        ui->end->setValue(data->a[2]);
        ui->kd2_check->show();
        ui->intmed_check->show();
    }
    else if(data->num_bind_site==3){
        ui->kd2->show();
        ui->kd2_label->show();
        ui->kd2->setValue(data->a[3]);
        ui->kd1->setValue(data->a[2]);
        ui->end->setValue(data->a[1]);
        ui->kdc->setValue(data->a[4]);
        ui->kdc->show();
        ui->kdc_label->show();
        ui->Kdc_check->show();
        ui->kd2_check->show();
    }
    else if(data->num_bind_site==4){
        ui->kd2_label->show();
        ui->kd2_check->show();
        ui->kd2->show();
        ui->kd3->show();
        ui->kd3_label->show();
        ui->kd3->show();
        ui->kd3_check->show();
        ui->kdc->show();
        ui->kdc_label->show();
        ui->Kdc_check->show();
        ui->kd1_label->setText("R20:");
        ui->kd2_label->setText("k_ex:");
        ui->kd3_label->setText("Population B (PB):");
        ui->kdc_label->setText("Delta W:");
        ui->kd1->setMaximum(30);
        ui->kd1->setValue(data->a[0]);
        ui->kd2->setValue(data->a[1]);
        ui->kd3->setValue(data->a[2]);
        ui->kdc->setValue(data->a[3]);
        ui->start->hide();
        ui->start_check->hide();
        ui->start_label->hide();
        ui->end->hide();
        ui->end_check->hide();
        ui->end_label->hide();
    }
}

Optimize::~Optimize()
{
    delete ui;
    delete last;
}

void Optimize::on_pushButton_3_clicked()
{
    Data temp;
    temp=*data;
    if(temp.num_bind_site==1){
        temp.a[0]=ui->start->value();
        temp.a[1]=ui->end->value();
        temp.a[2]=ui->kd1->value();
        Fitmrq fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, fit_one_site_dilution, fit_one_site_dilution};
        if(ui->start_check->isChecked())
            fit.ia[0]=false;
        if(ui->end_check->isChecked())
            fit.ia[1]=false;
        if(ui->kd1_check->isChecked())
            fit.ia[2]=false;
        fit.ia[3]=false;
        fit.ia[4]=false;
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.chi2=fit.chisq;
            std::vector<std::vector<double>> temp_vector{};
            temp_vector=fit.plot();
            temp.calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
            temp.calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
            ui->kd1->setValue(temp.a[2]);
            ui->chi2->setText(QString::number(temp.chi2));
            ui->start->setValue(temp.a[0]);
            ui->end->setValue(temp.a[1]);
            ui->plot_area->removeGraph(1);
            ui->plot_area->addGraph();
            ui->plot_area->graph()->setData(temp.calc_data[0],temp.calc_data[1]);
            ui->plot_area->graph()->setPen(temp.pen_curve);
            ui->plot_area->replot();
            ui->plot_area->rescaleAxes();
            if(last!=nullptr)
                delete last;
            last=nullptr;
            last=new Data(temp);

        }
    }
    else if(temp.num_bind_site==2){
        temp.a[0]=ui->start->value();
        temp.a[2]=ui->end->value();
        temp.a[3]=ui->kd1->value();
        temp.a[1]=ui->intmed1->value();
        temp.a[4]=ui->kd2->value();
        Fitmrq fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, fit_two_sites_dilution, fit_two_sites_dilution};
        if(ui->start_check->isChecked())
            fit.ia[0]=false;
        if(ui->intmed_check->isChecked())
            fit.ia[1]=false;
        if(ui->end_check->isChecked())
            fit.ia[2]=false;
        if(ui->kd1_check->isChecked())
            fit.ia[3]=false;
        if(ui->kd2_check->isChecked())
            fit.ia[4]=false;
        fit.ia[5]=false;
        fit.ia[6]=false;
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.chi2=fit.chisq;
            std::vector<std::vector<double>> temp_vector{};
            temp_vector=fit.plot();
            temp.calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
            temp.calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
            ui->kd1->setValue(temp.a[3]);
            ui->chi2->setText(QString::number(temp.chi2));
            ui->start->setValue(temp.a[0]);
            ui->end->setValue(temp.a[2]);
            ui->intmed1->setValue(temp.a[1]);
            ui->kd2->setValue(temp.a[4]);
            ui->plot_area->removeGraph(1);
            ui->plot_area->addGraph();
            ui->plot_area->graph()->setData(temp.calc_data[0],temp.calc_data[1]);
            ui->plot_area->graph()->setPen(temp.pen_curve);
            ui->plot_area->replot();
            ui->plot_area->rescaleAxes();
            if(last!=nullptr)
                delete last;
            last=nullptr;
            last=new Data(temp);
        }
    }
    else if(temp.num_bind_site==3){
        temp.a[0]=ui->start->value();
        temp.a[1]=ui->end->value();
        temp.a[2]=ui->kd1->value();
        temp.a[3]=ui->kd2->value();
        temp.a[4]=ui->kdc->value();
        Fitmrq2 fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.comp_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, fit_comp, fit_comp};
        if(ui->start_check->isChecked())
            fit.ia[0]=false;
        if(ui->Kdc_check->isChecked())
            fit.ia[4]=false;
        if(ui->end_check->isChecked())
            fit.ia[1]=false;
        if(ui->kd1_check->isChecked())
            fit.ia[2]=false;
        if(ui->kd2_check->isChecked())
            fit.ia[3]=false;
        fit.ia[5]=false;
        fit.ia[6]=false;
        fit.ia[7]=false;
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.chi2=fit.chisq;
            std::vector<std::vector<double>> temp_vector{};
            temp_vector=fit.plot();
            temp.calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
            temp.calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
            ui->kd1->setValue(temp.a[2]);
            ui->chi2->setText(QString::number(temp.chi2));
            ui->start->setValue(temp.a[0]);
            ui->end->setValue(temp.a[1]);
            ui->kd2->setValue(temp.a[3]);
            ui->kdc->setValue(temp.a[4]);
            ui->plot_area->removeGraph(1);
            ui->plot_area->addGraph();
            ui->plot_area->graph()->setData(temp.calc_data[0],temp.calc_data[1]);
            ui->plot_area->graph()->setPen(temp.pen_curve);
            ui->plot_area->replot();
            ui->plot_area->rescaleAxes();
            if(last!=nullptr)
                delete last;
            last=nullptr;
            last=new Data(temp);
        }
    }
    else if(temp.num_bind_site==4){
        temp.a[0]=ui->kd1->value();
        temp.a[2]=ui->kd3->value();
        temp.a[3]=ui->kdc->value();
        temp.a[1]=ui->kd2->value();
        Fitmrq3 fit{temp.n_cpmgVector.toStdVector(),temp.R2effVector.toStdVector(),
                    temp.dyVector.toStdVector(), temp.a,
                    carverrichards,carverrichards};
        if(ui->kd1_check->isChecked())
            fit.ia[0]=false;
        if(ui->kd2_check->isChecked())
            fit.ia[1]=false;
        if(ui->kd3_check->isChecked())
            fit.ia[2]=false;
        if(ui->Kdc_check->isChecked())
            fit.ia[3]=false;
        fit.ia[4]=false;
        int did_fit_work{};
        did_fit_work=fit.fit();

        if (did_fit_work==1)
            QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        else if(did_fit_work==2)
            QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        else{
            temp.model="CPMG-model";
            temp.num_bind_site=4;
            temp.a=fit.a;
            temp.chi2=fit.chisq;
            temp.kd_cpmg=CPMG_kd(temp.a);
            std::vector<std::vector<double>> temp_vector{};
            temp_vector=fit.plot();
            temp. calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
            temp. calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
            ui->kd1->setValue(temp.a[0]);
            ui->chi2->setText(QString::number(temp.chi2));
            ui->kd2->setValue(temp.a[1]);
            ui->kd3->setValue(temp.a[2]);
            ui->kdc->setValue(temp.a[3]);
            ui->plot_area->removeGraph(1);
            ui->plot_area->addGraph();
            ui->plot_area->graph()->setData(temp.calc_data[0],temp.calc_data[1]);
            ui->plot_area->graph()->setPen(temp.pen_curve);
            ui->plot_area->replot();
            ui->plot_area->rescaleAxes();
            if(last!=nullptr)
                delete last;
            last=nullptr;
            last=new Data(temp);
        }
    }
    if(last!=nullptr){
    Jackknife jack{last};
    jack.compute();
    }
}

void Optimize::on_apply_clicked()
{
    if (last!=nullptr){
        data->a=last->a;
        data->chi2=last->chi2;
        if(data->num_bind_site==1)
            data->kd_error=last->kd_error;
        else if(data->num_bind_site==2){
            data->kd1_error=last->kd1_error;
            data->kd2_error=last->kd2_error;
        }
        else if(data->num_bind_site==3){
            data->kd1_error=last->kd1_error;
            data->kd2_error=last->kd2_error;
            data->kdc_error=last->kdc_error;
        }
        else if(data->num_bind_site==4){
            data->kd_cpmg=last->kd_cpmg;
            data->kd_error=last->kd_error;
        }
        data->calc_data=last->calc_data;
    }
    hide();
}

void Optimize::on_pushButton_2_clicked()
{
    hide();
}
