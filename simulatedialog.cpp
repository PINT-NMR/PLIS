#include "simulatedialog.h"
#include "ui_simulatedialog.h"

SimulateDialog::SimulateDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SimulateDialog)
{
    ui->setupUi(this);
    ui->kd->setValidator(v);
    ui->kd2->setValidator(v);
    ui->kd1->setValidator(v);
    ui->protein_conc->setValidator(v);
    ui->start_comp->setValidator(v);
    ui->protein_conc_2->setValidator(v);
    ui->protein_conc_comp->setValidator(v);
    ui->kd1_comp->setValidator(v);
    ui->kd2_comp->setValidator(v);
    ui->kdc->setValidator(v);
    mode="standard";

    ui->tabWidget->removeTab(3);

#ifdef Q_OS_MAC
    ui->tabWidget->tabBar()->setFont(QFont("Helvetica", 12));
#endif

}

SimulateDialog::SimulateDialog(QWidget *parent, QString m) :
    QDialog(parent),
    ui(new Ui::SimulateDialog)
{
    mode=m;
    ui->setupUi(this);
    ui->ligand_conc->setValidator(v);
    ui->dw->setValidator(v);
    ui->kex->setValidator(v);
    ui->pb->setValidator(v);
    ui->r20->setValidator(v);
    for (int i=0;i<3;i++)
        ui->tabWidget->removeTab(0);

    if(mode=="cpmg"){
        ui->label->setText("N_cpmg range:");
        ui->ligand->hide();
        ui->maxLigandBox->hide();
        ui->minLigandBox->hide();
        ui->label_42->hide();
        ui->minVolume->setMaximum(100);
        ui->maxVolume->setMaximum(2000);
        ui->minVolume->setValue(25);
        ui->maxVolume->setValue(1000);
        ui->datapointBox->setValue(40);
    }

#ifdef Q_OS_MAC
    ui->tabWidget->tabBar()->setFont(QFont("Helvetica", 12));
#endif

}


SimulateDialog::~SimulateDialog()
{
    delete ui;
    delete v;
}

void SimulateDialog::setup(QVector<Data*> &datasets, QString mode){
    dataSets_local=datasets;
    Data* temp=new Data;
    if(mode=="standard"){
        temp->a[0]=0;
        temp->a[1]=1;
        temp->a[2]=10;
        temp->concVector.clear();
        temp->volumeVector.clear();
        temp->protein_conc_vector.clear();
        temp->error_vector.clear();
        temp->protein_conc=10;
    }
    else if(mode=="cpmg"){
        temp->a.clear();
        temp->a.push_back(5);
        temp->a.push_back(0.005);
        temp->a.push_back(1000);
        temp->a.push_back(1000);
        temp->a.push_back(0);
    }
    dataSets_local.push_back(temp);
    current_dataSet=dataSets_local.size()-1;
    dataSets_local[current_dataSet]->name="Simulated dataset";
}

void SimulateDialog::on_simulateButton_clicked()
{
    if(mode=="standard"){
        double volume_step, ligand_stock;
        volume_step=(ui->maxVolume->value() - ui->minVolume->value()) / (ui->datapointBox->value()-1);
        ligand_stock=(ui->maxLigandBox->value()*ui->maxVolume->value()-ui->minLigandBox->value()*ui->minVolume->value())
                / (ui->maxVolume->value()-ui->minVolume->value());
        for (int i=0; i<ui->datapointBox->value(); i++){
            dataSets_local[current_dataSet]->concVector.push_back((ui->minLigandBox->value()*ui->minVolume->value()+i*ligand_stock*volume_step)
                                                                  / (ui->minVolume->value()+volume_step*(i)));
            dataSets_local[current_dataSet]->volumeVector.push_back(ui->minVolume->value()+volume_step*(i));
            dataSets_local[current_dataSet]->protein_conc_vector.push_back(0);
            dataSets_local[current_dataSet]->error_vector.push_back(0.001);
        }
        dataSets_local[current_dataSet]->responceVector.clear();
        dataSets_local[current_dataSet]->data_visible=true;
        dataSets_local[current_dataSet]->curve_visible=false;

        if (ui->tabWidget->currentIndex()==0)
            one_site_sim();
        else if(ui->tabWidget->currentIndex()==1)
            two_site_sim();
        else if(ui->tabWidget->currentIndex()==2)
            comp_sim();
    }
    else if(mode=="cpmg"){
        QString s=ui->r20->text();
        s.replace(',','.');
        dataSets_local[current_dataSet]->a[0]=s.toDouble();
        s=ui->dw->text();
        s.replace(',','.');
        dataSets_local[current_dataSet]->a[3]=s.toDouble();
        s=ui->kex->text();
        s.replace(',','.');
        dataSets_local[current_dataSet]->a[1]=s.toDouble();
        s=ui->pb->text();
        s.replace(',','.');
        dataSets_local[current_dataSet]->a[2]=s.toDouble();
        s=ui->ligand_conc->text();
        s.replace(',','.');
        dataSets_local[current_dataSet]->a[4]=s.toDouble();
        dataSets_local[current_dataSet]->ligand_cpmg=s.toDouble();
        double n=((ui->maxVolume->value()-ui->minVolume->value())/(ui->datapointBox->value()-1));
        for (int i=0;i<ui->datapointBox->value();i++){
            dataSets_local[current_dataSet]->n_cpmgVector.push_back(n*i+ui->minVolume->value());
            dataSets_local[current_dataSet]->dyVector.push_back(1.e-3);
            dataSets_local[current_dataSet]->R2effVector.push_back(
                        carverrichards(dataSets_local[current_dataSet]->n_cpmgVector[i],
                                       dataSets_local[current_dataSet]->a));
        }
        dataSets_local[current_dataSet]->data_visible=true;
        dataSets_local[current_dataSet]->curve_visible=false;
        simulateGenerateData(dataSets_local[current_dataSet]->R2effVector,ui->noiseBox->value());
        hide();
    }
    setResult(1);

}

void SimulateDialog::one_site_sim(){
    dataSets_local[current_dataSet]->a[0]=ui->responce_start->value();
    dataSets_local[current_dataSet]->a[1]=ui->responce_end->value();
    QString s=ui->protein_conc->text();
    s=ui->kd->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a[2]=s.toDouble();
    s=ui->protein_conc->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->protein_conc=s.toDouble();
    dataSets_local[current_dataSet]->a.push_back(dataSets_local[current_dataSet]->concVector[0]);
    dataSets_local[current_dataSet]->a.push_back(dataSets_local[current_dataSet]->protein_conc);
    dataSets_local[current_dataSet]->update_protein_conc();

    for (int i=0; i<dataSets_local[current_dataSet]->concVector.size(); i++){
        dataSets_local[current_dataSet]->responceVector.push_back(
                    fit_one_site_dilution(dataSets_local[current_dataSet]->concVector[i], dataSets_local[current_dataSet]->protein_conc_vector[i],
                                 dataSets_local[current_dataSet]->a));
    }
    simulateGenerateData(dataSets_local[current_dataSet]->responceVector,ui->noiseBox->value());
    hide();
}

void SimulateDialog::two_site_sim(){
    dataSets_local[current_dataSet]->a[0]=ui->responce_s->value();
    dataSets_local[current_dataSet]->a[1]=ui->responce_I->value();
    dataSets_local[current_dataSet]->a[2]=ui->responce_L->value();
    QString s=ui->kd1->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    s=ui->kd2->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    dataSets_local[current_dataSet]->a.push_back(dataSets_local[current_dataSet]->concVector[0]); //since we need start conc for the calculations of zbrent
    s=ui->protein_conc_2->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    dataSets_local[current_dataSet]->protein_conc=s.toDouble();
    dataSets_local[current_dataSet]->update_protein_conc();
    for (int i=0;i<3;i++)
        dataSets_local[current_dataSet]->dyda.push_back(0);
    for (int i=0; i<dataSets_local[current_dataSet]->concVector.size(); i++){
        dataSets_local[current_dataSet]->responceVector.push_back(
                    fit_two_sites_dilution(dataSets_local[current_dataSet]->concVector[i], dataSets_local[current_dataSet]->protein_conc_vector[i],
                                 dataSets_local[current_dataSet]->a));
    }
    simulateGenerateData(dataSets_local[current_dataSet]->responceVector,ui->noiseBox->value());
    hide();

}

void SimulateDialog::comp_sim(){
    dataSets_local[current_dataSet]->a[0]=ui->responce_s_comp->value();
    dataSets_local[current_dataSet]->a[1]=ui->responce_end_comp->value();
    QString s=ui->kd1_comp->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a[2]=s.toDouble();
    s=ui->kd2_comp->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    s=ui->kdc->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    dataSets_local[current_dataSet]->a.push_back(dataSets_local[current_dataSet]->concVector[0]); //since we need start conc for the calculations of zbrent
    s=ui->protein_conc_comp->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    dataSets_local[current_dataSet]->protein_conc=s.toDouble();
    s=ui->start_comp->text();
    s.replace(',','.');
    dataSets_local[current_dataSet]->a.push_back(s.toDouble());
    dataSets_local[current_dataSet]->update_protein_conc();
    dataSets_local[current_dataSet]->update_comp_vector(dataSets_local[current_dataSet]->a.back());
    for (unsigned int i=0;i<dataSets_local[current_dataSet]->a.size();i++)
        dataSets_local[current_dataSet]->dyda.push_back(0);
    for (int i=0; i<dataSets_local[current_dataSet]->concVector.size(); i++){
        dataSets_local[current_dataSet]->responceVector.push_back(
                    fit_comp(dataSets_local[current_dataSet]->concVector[i], dataSets_local[current_dataSet]->protein_conc_vector[i],
                             dataSets_local[current_dataSet]->comp_vector[i], dataSets_local[current_dataSet]->a));
    }
    simulateGenerateData(dataSets_local[current_dataSet]->responceVector,ui->noiseBox->value());
    hide();
}

void SimulateDialog::on_responce_s_valueChanged(double arg1)
{
    if (ui->responce_s->value()>=ui->responce_I->value())
        ui->responce_I->setValue(arg1+0.01);
}

void SimulateDialog::on_responce_I_valueChanged(double arg1)
{
    if(ui->responce_I->value()<=ui->responce_s->value())
        ui->responce_s->setValue(arg1-0.01);

    if (ui->responce_I->value()>=ui->responce_L->value())
        ui->responce_L->setValue(arg1+0.01);

}

void SimulateDialog::on_responce_L_valueChanged(double arg1)
{
    if(ui->responce_L->value()<=ui->responce_I->value())
       ui->responce_I->setValue(arg1-0.01);
}

void SimulateDialog::on_responce_start_valueChanged(double arg1)
{
    if (ui->responce_start->value()>=ui->responce_end->value())
        ui->responce_end->setValue(arg1+0.01);
}

void SimulateDialog::on_responce_end_valueChanged(double arg1)
{
    if (ui->responce_start->value()>=ui->responce_end->value())
        ui->responce_start->setValue(arg1-0.01);
}

void SimulateDialog::on_cancelButton_clicked()
{
    reject();
}

void SimulateDialog::simulateGenerateData(QVector<double> &x, double noise)
{
     const double mean = 0.0;
     std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
     std::normal_distribution<double> dist(mean,noise);
     for (int i=0;i<x.size();i++)
         x[i]+=dist(generator);
}
