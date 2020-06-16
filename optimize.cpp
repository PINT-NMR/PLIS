#include "optimize.h"
#include "ui_optimize.h"

// 200615
// - absolute values are now taken of calculatd Kd:s

void showLabel(QLabel *label, QString text) {
    label->show();
    label->setText(text);
}

void showValue(QDoubleSpinBox *spinbox, double value) {
    spinbox->show();
    spinbox->setValue(value);
}

void showGroup(QLabel *label, QDoubleSpinBox *spinbox, QCheckBox *check, QString text, double value)
{
    showLabel(label, text);
    showValue(spinbox, value);
    check->show();
}

void hideGroup(QLabel *label, QDoubleSpinBox *spinbox, QCheckBox *check)
{
    label->hide();
    spinbox->hide();
    check->hide();
}

Optimize::Optimize(QWidget *parent, Data* &d, QString _mode) :
    QDialog(parent),
    ui(new Ui::Optimize), data(d), mode(_mode)
{
    ui->setupUi(this);
    ui->plot_area->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
                                    QCP::iSelectLegend | QCP::iSelectPlottables | QCP::iSelectItems);

    hideGroup(ui->kd2_label, ui->kd2, ui->kd2_check);
    hideGroup(ui->kdc_label, ui->kdc, ui->Kdc_check);
    hideGroup(ui->kd4_label, ui->kd4, ui->kd4_check);
    hideGroup(ui->intmed_label, ui->intmed1, ui->intmed_check);
    hideGroup(ui->end_label, ui->intmed1, ui->intmed_check);
    hideGroup(ui->end_label2, ui->end2, ui->end2_check);
    hideGroup(ui->end_label3, ui->end3, ui->end3_check);

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
    if (data->num_bind_site==1) {
        showGroup(ui->kd1_label, ui->kd1, ui->kd1_check, "Kd:", data->a[2]);
        showGroup(ui->start_label, ui->start, ui->start_check, "s(P):", data->a[0]);
        showGroup(ui->end_label, ui->end, ui->end_check, "s(PL):", data->a[1]);

    }
    else if(data->num_bind_site==2) {
        showGroup(ui->kd1_label, ui->kd1, ui->kd1_check, "Kd1:", data->a[3]);
        showGroup(ui->kd2_label, ui->kd2, ui->kd2_check, "Kd2:", data->a[4]);
        showGroup(ui->start_label, ui->start, ui->start_check, "s(P):", data->a[0]);
        showGroup(ui->intmed_label, ui->intmed1, ui->intmed_check, "s(PL):", data->a[1]);
        showGroup(ui->end_label, ui->end, ui->end_check, "s(PL2):", data->a[2]);
    }
    else if (data->num_bind_site==5) {
        showGroup(ui->kd1_label, ui->kd1, ui->kd1_check, "Kd1:", data->a[5]);
        showGroup(ui->kd2_label, ui->kd2, ui->kd2_check, "Kd2:", data->a[6]);
        showGroup(ui->kdc_label, ui->kdc, ui->Kdc_check, "Kd3:", data->a[7]);
        showGroup(ui->kd4_label, ui->kd4, ui->kd4_check, "Kd4:", data->a[8]);
        showGroup(ui->start_label, ui->start, ui->start_check, "s(P):", data->a[0]);
        showGroup(ui->intmed_label, ui->intmed1, ui->intmed_check, "s(PL):", data->a[1]);
        showGroup(ui->end_label, ui->end, ui->end_check, "s(PL2):", data->a[2]);
        showGroup(ui->end_label2, ui->end2, ui->end2_check, "s(PL3):", data->a[3]);
        showGroup(ui->end_label3, ui->end3, ui->end3_check, "s(PL4):", data->a[4]);
    }
    else if (data->num_bind_site==3) {
        showGroup(ui->kd1_label, ui->kd1, ui->kd1_check, "Kd1:", data->a[2]);
        showGroup(ui->kd2_label, ui->kd2, ui->kd2_check, "Kd2:", data->a[3]);
        showGroup(ui->kdc_label, ui->kdc, ui->Kdc_check, "Kd3:", data->a[4]);
        showGroup(ui->start_label, ui->start, ui->start_check, "s(Comp):", data->a[0]);
        showGroup(ui->end_label, ui->end, ui->end_check, "s(CompL):", data->a[1]);
    }
    else if (data->num_bind_site==4) {
        showGroup(ui->kd1_label, ui->kd1, ui->kd1_check, "R<sub>2,0</sub>:", data->a[0]);
        showGroup(ui->kd2_label, ui->kd2, ui->kd2_check, "k<sub>ex</sub>:", data->a[1]);
        showGroup(ui->kdc_label, ui->kdc, ui->Kdc_check, "p<sub>B</sub>:", data->a[2]);
        showGroup(ui->kd4_label, ui->kd4, ui->kd4_check, "∆𝜔:", data->a[3]);
        hideGroup(ui->start_label, ui->start, ui->start_check);
        hideGroup(ui->end_label, ui->end, ui->end_check);
        ui->kd1->setMaximum(100.);
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
    std::vector<bool> ia{};

    if(temp.num_bind_site==1){
        temp.a[0]=ui->start->value();
        temp.a[1]=ui->end->value();
        temp.a[2]=ui->kd1->value();
        ia = {!ui->start_check->isChecked(), !ui->end_check->isChecked(), !ui->kd1_check->isChecked(), false, false};
        void (*voidptr)(const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_one_site_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_one_site_dilution};
        if (mode=="chemshift") {
            voidptr = fit_one_site;
            dblptr = fit_one_site;
        }
        Fitmrq fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, ia, voidptr, dblptr};
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Singular matrix.");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Too many iterations.\nThe best fit could not be found.");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.a[2] = fabs(fit.a[2]);
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
    else if (temp.num_bind_site==2) {
        temp.a[0]=ui->start->value();
        temp.a[1]=ui->intmed1->value();
        temp.a[2]=ui->end->value();
        temp.a[3]=ui->kd1->value();
        temp.a[4]=ui->kd2->value();

        ia = { !ui->start_check->isChecked(), !ui->intmed_check->isChecked(), !ui->end_check->isChecked(),
               !ui->kd1_check->isChecked(), !ui->kd2_check->isChecked(), false, false};
        void (*voidptr)(const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_two_sites_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_two_sites_dilution};
        if (mode=="chemshift") {
            voidptr = fit_two_sites;
            dblptr = fit_two_sites;
        }

        Fitmrq fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, ia, voidptr, dblptr};
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Singular matrix.");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Too many iterations.\nThe best fit could not be found");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.a[3] = fabs(fit.a[3]);
            temp.a[4] = fabs(fit.a[4]);
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
    else if (temp.num_bind_site==5) {
        temp.a[0]=ui->start->value();
        temp.a[1]=ui->intmed1->value();
        temp.a[2]=ui->end->value();
        temp.a[3]=ui->end2->value();
        temp.a[4]=ui->end3->value();
        temp.a[5]=ui->kd1->value();
        temp.a[6]=ui->kd2->value();
        temp.a[7]=ui->kdc->value();
        temp.a[8]=ui->kd4->value();

        ia = { !ui->start_check->isChecked(), !ui->intmed_check->isChecked(), !ui->end_check->isChecked(),
               !ui->end2_check->isChecked(), !ui->end3_check->isChecked(), !ui->kd1_check->isChecked(),
               !ui->kd2_check->isChecked(), !ui->Kdc_check->isChecked(), !ui->kd4_check->isChecked(), false, false};
        void (*voidptr)(const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_four_sites_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_four_sites_dilution};
        if (mode=="chemshift") {
            voidptr = fit_four_sites;
            dblptr = fit_four_sites;
        }
        Fitmrq fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, ia, voidptr, dblptr};
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Singular matrix.");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Too many iterations.\nThe best fit could not be found");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.a[5] = fabs(fit.a[5]);
            temp.a[6] = fabs(fit.a[6]);
            temp.a[7] = fabs(fit.a[7]);
            temp.a[8] = fabs(fit.a[8]);
            temp.chi2=fit.chisq;
            std::vector<std::vector<double>> temp_vector{};
            temp_vector=fit.plot();
            temp.calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
            temp.calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
            ui->chi2->setText(QString::number(temp.chi2));
            ui->kd1->setValue(temp.a[5]);
            ui->kd2->setValue(temp.a[6]);
            ui->kdc->setValue(temp.a[7]);
            ui->kd4->setValue(temp.a[8]);
            ui->start->setValue(temp.a[0]);
            ui->intmed1->setValue(temp.a[1]);
            ui->end->setValue(temp.a[2]);
            ui->end2->setValue(temp.a[3]);
            ui->end3->setValue(temp.a[4]);
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
    else if (temp.num_bind_site==3) {
        temp.a[0]=ui->start->value();
        temp.a[1]=ui->end->value();
        temp.a[2]=ui->kd1->value();
        temp.a[3]=ui->kd2->value();
        temp.a[4]=ui->kdc->value();
        ia = {!ui->start_check->isChecked(), !ui->end_check->isChecked(), !ui->kd1_check->isChecked(),
                  !ui->kd2_check->isChecked(), !ui->Kdc_check->isChecked(), false, false, false};
        void (*voidptr)(const double, const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_comp_dilution};
        double (*dblptr)(const double, const double, const double, const std::vector<double> &){fit_comp_dilution};
        if (mode=="chemshift") {
            voidptr = fit_comp;
            dblptr = fit_comp;
        }
        Fitmrq2 fit{temp.concVector.toStdVector(),
                    temp.protein_conc_vector.toStdVector(),
                    temp.comp_vector.toStdVector(),
                    temp.responceVector.toStdVector(),
                    temp.volumeVector.toStdVector(),
                    temp.error_vector.toStdVector(),
                    temp.a, ia, voidptr, dblptr};
        int did_fit_work{};
        did_fit_work=fit.fit();
        if (did_fit_work==1){ //singular matrix!
            QMessageBox::information(this, "Warning: Something went wrong!", "Singular matrix.");
        }
        else if(did_fit_work==2){
            QMessageBox::information(this, "Warning: Something went wrong!", "Too many iterations.\nThe best fit could not be found.");
        }
        else {
            QMessageBox::information(this, "Optimization done", "Success!\nOptimization done.");
            temp.a=fit.a;
            temp.a[2] = fabs(fit.a[2]);
            temp.a[3] = fabs(fit.a[3]);
            temp.a[4] = fabs(fit.a[4]);
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
    else if(temp.num_bind_site==4) {
        temp.a[0]=ui->kd1->value();
        temp.a[1]=ui->kd2->value();
        temp.a[2]=ui->kd4->value();
        temp.a[3]=ui->kdc->value();
        ia = {!ui->kd1_check->isChecked(), !ui->kd2_check->isChecked(), !ui->kd4_check->isChecked(), !ui->Kdc_check->isChecked(), false};
        Fitmrq3 fit{temp.n_cpmgVector.toStdVector(),temp.R2effVector.toStdVector(),
                    temp.dyVector.toStdVector(), temp.a, ia,
                    carverrichards,carverrichards};
        int did_fit_work{};
        did_fit_work=fit.fit();

        if (did_fit_work==1)
            QMessageBox::warning(this, "Warning: Something went wrong!", "Singular matrix.");
        else if(did_fit_work==2)
            QMessageBox::warning(this, "Warning: Something went wrong!", "Too many iterations.\nThe best fit could not be found.");
        else {
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
            ui->kd4->setValue(temp.a[2]);
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
    if (last!=nullptr) {
    Jackknife jack{last, ia};
    jack.compute();
    }
}

void Optimize::on_apply_clicked()
{
    if (last!=nullptr){
        data->a=last->a;
        data->dyda = last->dyda;
        data->chi2=last->chi2;
        if(data->num_bind_site==1)
            data->kd_error=last->kd_error;
        else if(data->num_bind_site==2){
            data->kd1_error=last->kd1_error;
            data->kd2_error=last->kd2_error;
        }
        else if (data->num_bind_site==5) {
            data->kd1_error = last->kd1_error;
            data->kd2_error = last->kd2_error;
            data->kd3_error = last->kd3_error;
            data->kd4_error = last->kd4_error;
        }
        else if (data->num_bind_site==3) {
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
