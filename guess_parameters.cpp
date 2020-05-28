#include "guess_parameters.h"
#include "ui_guess_parameters.h"

Guess_parameters::Guess_parameters(QWidget *parent, Data *&d) :
    QDialog(parent),
    ui(new Ui::Guess_parameters)
{
    ui->setupUi(this);
#ifdef Q_OS_MAC
    ui->label_8->setFont(QFont("Helvetica",12));
    ui->label_9->setFont(QFont("Helvetica",12));
    ui->label_6->setFont(QFont("Helvetica",12));
    ui->label_7->setFont(QFont("Helvetica",12));
    ui->label_16->setFont(QFont("Helvetica",12));
    ui->label_13->setFont(QFont("Helvetica",12));
    ui->label_11->setFont(QFont("Helvetica",12));
    ui->label_12->setFont(QFont("Helvetica",12));
    ui->label_26->setFont(QFont("Helvetica",12));
    ui->label_23->setFont(QFont("Helvetica",12));
    ui->label_27->setFont(QFont("Helvetica",12));
    ui->label_25->setFont(QFont("Helvetica",12));
    ui->label_41->setFont(QFont("Helvetica",12));
    ui->label_43->setFont(QFont("Helvetica",12));
    ui->label_44->setFont(QFont("Helvetica",12));
    ui->label_45->setFont(QFont("Helvetica",12));
    ui->label_37->setFont(QFont("Helvetica",12));
    ui->label_36->setFont(QFont("Helvetica",12));
    ui->label_58->setFont(QFont("Helvetica",12));
    ui->label_52->setFont(QFont("Helvetica",12));
    ui->label_51->setFont(QFont("Helvetica",12));
    ui->label_62->setFont(QFont("Helvetica",12));
    ui->label_64->setFont(QFont("Helvetica",12));
    ui->label_65->setFont(QFont("Helvetica",12));
    ui->label_66->setFont(QFont("Helvetica",12));
#endif
    data=d;
    if(data->num_bind_site==1){
        ui->tabWidget->setCurrentIndex(0);
        ui->tabWidget->removeTab(1);
        ui->tabWidget->removeTab(1);
        ui->tabWidget->removeTab(1);
    }
    else if(data->num_bind_site==2){
        ui->tabWidget->removeTab(2);
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(1);
    }
    else if(data->num_bind_site==3){
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(1);
    }
    else if(data->num_bind_site==4){
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(0);
    }
}

Guess_parameters::~Guess_parameters()
{
    delete ui;
}

void Guess_parameters::on_pushButton_clicked()
{
    if(ui->tabWidget->tabText(ui->tabWidget->currentIndex())=="One site"){
        if(ui->tabWidget_one->currentIndex()==0){
            data->a[2]=ui->kd1_one->value();
            data->a[0]=ui->start_guess_one->value();
            data->a[1]=ui->end_guess_one->value();
        }
        else{
            QVector<double> kd{},start{},end{};
            if(ui->kd_check->isChecked())
                kd.push_back(ui->kd_min->value());
            else{
                double d=(ui->kd_max->value()-ui->kd_min->value())/(ui->kd_number->value()-1);
                for (int i=0;i<ui->kd_number->value();i++)
                    kd.push_back(i*d+ui->kd_min->value());
            }
            if(ui->start_check_range->isChecked())
                start.push_back(ui->start_min->value());
            else{
                double d=(ui->start_max->value()-ui->start_min->value())/(ui->start_number->value()-1);
                for (int i=0;i<ui->start_number->value();i++)
                    start.push_back(i*d+ui->start_min->value());
            }
            if(ui->end_check_range->isChecked())
                end.push_back(ui->end_min->value());
            else{
                double d=(ui->end_max->value()-ui->end_min->value())/(ui->end_number->value()-1);
                for (int i=0;i<ui->end_number->value();i++)
                    end.push_back(i*d+ui->end_min->value());
            }
            guess_parameters_range(data,kd,start,end);
        }
    }
    else if(ui->tabWidget->tabText(ui->tabWidget->currentIndex())=="Two site"){
        if(ui->tabWidget_two->currentIndex()==0){
            data->a[3]=ui->kd1_two->value();
            data->a[0]=ui->start_guess_two->value();
            data->a[2]=ui->end_guess_two->value();
            data->a[1]=ui->intermediate_guess_two->value();
            data->a[4]=ui->kd2_two->value();
        }
        else{
            QVector<double> kd{}, kd2{},start{},inter{},end{};
            if(ui->kd_check_range_2->isChecked())
                kd.push_back(ui->kd_min_2->value());
            else{
                double d=(ui->kd_max_2->value()-ui->kd_min_2->value())/(ui->kd_number_2->value()-1);
                for (int i=0;i<ui->kd_number_2->value();i++)
                    kd.push_back(i*d+ui->kd_min_2->value());
            }
            if(ui->start_check_range_2->isChecked())
                start.push_back(ui->start_min_2->value());
            else{
                double d=(ui->start_max_2->value()-ui->start_min_2->value())/(ui->start_number_2->value()-1);
                for (int i=0;i<ui->start_number_2->value();i++)
                    start.push_back(i*d+ui->kd_min_2->value());
            }
            if(ui->end_check_range_2->isChecked())
                end.push_back(ui->end_min_2->value());
            else{
                double d=(ui->end_max_2->value()-ui->end_min_2->value())/(ui->end_number_2->value()-1);
                for (int i=0;i<ui->end_number_2->value();i++)
                    end.push_back(i*d+ui->end_min_2->value());
            }
            if(ui->intermediate_check_range->isChecked())
                inter.push_back(ui->intermediate_min->value());
            else{
                double d=(ui->intermediate_max->value()-ui->intermediate_min->value())/(ui->intermediate_number->value()-1);
                for (int i=0;i<ui->intermediate_number->value();i++)
                    inter.push_back(i*d+ui->intermediate_min->value());
            }
            if(ui->kd2_check_range->isChecked())
                kd2.push_back(ui->kd2_min->value());
            else{
                double d=(ui->kd2_max->value()-ui->kd2_min->value())/(ui->kd2_number_2->value()-1);
                for (int i=0;i<ui->kd2_number_2->value();i++)
                    kd2.push_back(i*d+ui->kd2_min->value());
            }
            guess_parameters_range_2_sites(data,kd,kd2,start,inter,end);
        }
    }
    else if(ui->tabWidget->tabText(ui->tabWidget->currentIndex())=="Competitive binding"){
        if(ui->tabWidget_comp->currentIndex()==0){
            data->a[2]=ui->kd1_comp->value();
            data->a[0]=ui->start_guess_comp->value();
            data->a[1]=ui->end_guess_comp->value();
            data->a[3]=ui->kd2_comp->value();
            data->a[4]=ui->kdc->value();
        }
        else{
            QVector<double> kd{}, kd2{},start{},kdc{},end{};
            if(ui->kd_check_range_comp->isChecked())
                kd.push_back(ui->kd_min_comp->value());
            else{
                double d=(ui->kd_max_comp->value()-ui->kd_min_comp->value())/(ui->kd_num_comp->value()-1);
                for (int i=0;i<ui->kd_num_comp->value();i++)
                    kd.push_back(i*d+ui->kd_min_comp->value());
            }
            if(ui->start_check_range_comp->isChecked())
                start.push_back(ui->start_min_comp->value());
            else{
                double d=(ui->start_max_comp->value()-ui->start_min_comp->value())/(ui->start_number_comp->value()-1);
                for (int i=0;i<ui->start_number_comp->value();i++)
                    start.push_back(i*d+ui->start_min_comp->value());
            }
            if(ui->end_check_range_comp->isChecked())
                end.push_back(ui->end_min_comp->value());
            else{
                double d=(ui->end_max_comp->value()-ui->end_min_comp->value())/(ui->end_number_comp->value()-1);
                for (int i=0;i<ui->end_number_comp->value();i++)
                    end.push_back(i*d+ui->end_min_comp->value());
            }
            if(ui->kdc_check_range->isChecked())
                kdc.push_back(ui->kdc_min->value());
            else{
                double d=(ui->kdc_max->value()-ui->kdc_min->value())/(ui->kdc_num->value()-1);
                for (int i=0;i<ui->kdc_num->value();i++)
                    kdc.push_back(i*d+ui->kdc_min->value());
            }
            if(ui->kd2_check_range_comp->isChecked())
                kd2.push_back(ui->kd2_min_comp->value());
            else{
                double d=(ui->kd2_max_comp->value()-ui->kd2_min_comp->value())/(ui->kd2_num_comp->value()-1);
                for (int i=0;i<ui->kd2_num_comp->value();i++)
                    kd2.push_back(i*d+ui->kd2_min_comp->value());
            }
            guess_parameters_comp_range(data,kd,kd2,kdc,start,end);
        }
    }
    else if(ui->tabWidget->tabText(ui->tabWidget->currentIndex())=="CPMG"){
        if(ui->tabWidget_cpmg->currentIndex()==0){
            data->a[2]=ui->pb->value();
            data->a[0]=ui->r20->value();
            data->a[1]=ui->kex->value();
            data->a[3]=ui->dw->value();
        }
        else{
            QVector<double> r20{}, pb{},kex{},dw{};
            if(ui->r20_check_r->isChecked())
                r20.push_back(ui->r20_min->value());
            else{
                double d=(ui->r20_max->value()-ui->r20_min->value())/(ui->r20_num_r->value()-1);
                for (int i=0;i<ui->r20_num_r->value();i++)
                    r20.push_back(i*d+ui->r20_min->value());
            }
            if(ui->pb_check_r->isChecked())
                pb.push_back(ui->pb_min->value());
            else{
                double d=(ui->pb_max->value()-ui->pb_min->value())/(ui->pb_num_r->value()-1);
                for (int i=0;i<ui->pb_num_r->value();i++)
                    pb.push_back(i*d+ui->pb_min->value());
            }
            if(ui->kex_check_r->isChecked())
                kex.push_back(ui->kex_min->value());
            else{
                double d=(ui->kex_max->value()-ui->kex_min->value())/(ui->kex_num_r->value()-1);
                for (int i=0;i<ui->kex_num_r->value();i++)
                    kex.push_back(i*d+ui->kex_min->value());
            }
            if(ui->dw_check_r->isChecked())
                dw.push_back(ui->dw_min->value());
            else{
                double d=(ui->dw_max->value()-ui->dw_min->value())/(ui->dw_num_r->value()-1);
                for (int i=0;i<ui->dw_num_r->value();i++)
                    dw.push_back(i*d+ui->dw_min->value());
            }
            guess_parameters_cpmg_range(data,r20,kex,pb,dw);
        }
    }
    setResult(1);
    hide();
}

void Guess_parameters::on_kd1_check_one_stateChanged(int arg1)
{
    a2=arg1;
}

void Guess_parameters::on_start_check_one_stateChanged(int arg1)
{
    a0=arg1;
}

void Guess_parameters::on_end_check_one_stateChanged(int arg1)
{
    a1=arg1;
}

void Guess_parameters::on_kd1_check_two_stateChanged(int arg1)
{
    a3=arg1;
}

void Guess_parameters::on_kd2_check_two_stateChanged(int arg1)
{
    a4=arg1;
}

void Guess_parameters::on_start_check_two_stateChanged(int arg1)
{
    a0=arg1;
}

void Guess_parameters::on_intermediate_check_two_stateChanged(int arg1)
{
    a1=arg1;
}

void Guess_parameters::on_end_check_two_stateChanged(int arg1)
{
    a2=arg1;
}

void Guess_parameters::on_kd_min_valueChanged(double arg1)
{
    if(ui->kd_check->isChecked())
       ui->kd_max->setValue(arg1);
}

void Guess_parameters::on_kd_max_valueChanged(double arg1)
{
    if(ui->kd_check->isChecked())
        ui->kd_min->setValue(arg1);
}

void Guess_parameters::on_start_min_valueChanged(double arg1)
{
    if(ui->start_check_range->isChecked())
        ui->start_max->setValue(arg1);
}

void Guess_parameters::on_start_max_valueChanged(double arg1)
{
    if(ui->start_check_range->isChecked())
        ui->start_min->setValue(arg1);
}

void Guess_parameters::on_end_min_valueChanged(double arg1)
{
    if(ui->end_check_range->isChecked())
        ui->end_max->setValue(arg1);
}

void Guess_parameters::on_end_max_valueChanged(double arg1)
{
    if(ui->end_check_range->isChecked())
        ui->end_min->setValue(arg1);
}

void Guess_parameters::on_kd_check_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kd_max->setValue(ui->kd_min->value());
        a2=true;
    }
    else
        a2=false;
}

void Guess_parameters::on_start_check_range_stateChanged(int arg1)
{
    if(arg1==2){
        ui->start_max->setValue(ui->start_min->value());
        a0=true;
    }
    else
        a0=false;
}

void Guess_parameters::on_end_check_range_stateChanged(int arg1)
{
    if(arg1==2){
        ui->end_max->setValue(ui->end_min->value());
        a1=true;
    }
    else
        a1=false;
}

void Guess_parameters::on_kd_min_2_valueChanged(double arg1)
{
    if(ui->kd_check_range_2->isChecked())
        ui->kd_max_2->setValue(arg1);
}

void Guess_parameters::on_kd_max_2_valueChanged(double arg1)
{
    if(ui->kd_check_range_2->isChecked())
        ui->kd_min_2->setValue(arg1);
}

void Guess_parameters::on_kd2_min_valueChanged(double arg1)
{
    if(ui->kd2_check_range->isChecked())
        ui->kd2_max->setValue(arg1);
}

void Guess_parameters::on_kd2_max_valueChanged(double arg1)
{
    if(ui->kd2_check_range->isChecked())
        ui->kd2_min->setValue(arg1);
}

void Guess_parameters::on_start_min_2_valueChanged(double arg1)
{
    if(ui->start_check_range_2->isChecked())
        ui->start_max_2->setValue(arg1);
}

void Guess_parameters::on_start_max_2_valueChanged(double arg1)
{
    if(ui->start_check_range_2->isChecked())
        ui->start_min_2->setValue(arg1);
}

void Guess_parameters::on_intermediate_min_valueChanged(double arg1)
{
    if(ui->intermediate_check_range->isChecked())
        ui->intermediate_max->setValue(arg1);
}

void Guess_parameters::on_intermediate_max_valueChanged(double arg1)
{
    if(ui->intermediate_check_range->isChecked())
        ui->intermediate_min->setValue(arg1);
}

void Guess_parameters::on_end_min_2_valueChanged(double arg1)
{
    if(ui->end_check_range_2->isChecked())
        ui->end_max_2->setValue(arg1);
}

void Guess_parameters::on_end_max_2_valueChanged(double arg1)
{
    if(ui->end_check_range_2->isChecked())
        ui->end_min_2->setValue(arg1);
}

void Guess_parameters::on_kd_check_range_2_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kd_max_2->setValue(ui->kd_min_2->value());
        a3=true;
    }
    else
        a3=false;
}

void Guess_parameters::on_kd2_check_range_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kd2_max->setValue(ui->kd2_min->value());
        a4=true;
    }
    else
        a4=false;
}

void Guess_parameters::on_start_check_range_2_stateChanged(int arg1)
{
    if(arg1==2){
        ui->start_max_2->setValue(ui->start_min_2->value());
        a0=true;
    }
    else
        a0=false;
}

void Guess_parameters::on_intermediate_check_range_stateChanged(int arg1)
{
    if(arg1==2){
        ui->intermediate_max->setValue(ui->intermediate_min->value());
        a1=true;
    }
    else
        a1=false;
}

void Guess_parameters::on_end_check_range_2_stateChanged(int arg1)
{
    if(arg1==2){
        ui->end_max_2->setValue(ui->end_min_2->value());
        a2=true;
    }
    else
        a2=false;
}

void Guess_parameters::on_cancel_clicked()
{
    setResult(0);
    reject();
}

void Guess_parameters::on_kd1_check_two_2_stateChanged(int arg1)
{
    a2=arg1;
}

void Guess_parameters::on_kd2_check_two_2_stateChanged(int arg1)
{
    a3=arg1;
}

void Guess_parameters::on_kdc_check_stateChanged(int arg1)
{
    a4=arg1;
}

void Guess_parameters::on_start_check_two_2_stateChanged(int arg1)
{
    a0=arg1;
}

void Guess_parameters::on_end_check_two_2_stateChanged(int arg1)
{
    a1=arg1;
}

void Guess_parameters::on_kd_check_range_comp_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kd_max_comp->setValue(ui->kd_min_2->value());
        a2=true;
    }
    else
        a2=false;
}

void Guess_parameters::on_kd2_check_range_comp_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kd2_max_comp->setValue(ui->kd2_max_comp->value());
        a3=true;
    }
    else
        a3=false;
}

void Guess_parameters::on_kdc_check_range_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kdc_max->setValue(ui->kdc_max->value());
        a4=true;
    }
    else
        a4=false;
}

void Guess_parameters::on_start_check_range_comp_stateChanged(int arg1)
{
    if(arg1==2){
        ui->start_max_comp->setValue(ui->start_max_comp->value());
        a0=true;
    }
    else
        a0=false;
}

void Guess_parameters::on_end_check_range_comp_stateChanged(int arg1)
{
    if(arg1==2){
        ui->end_max_comp->setValue(ui->end_min_comp->value());
        a1=true;
    }
    else
        a1=false;
}

void Guess_parameters::on_kd_min_comp_valueChanged(double arg1)
{
    if(ui->kd_check_range_comp->isChecked())
        ui->kd_max_comp->setValue(arg1);
}

void Guess_parameters::on_kd_max_comp_valueChanged(double arg1)
{
    if(ui->kd_check_range_comp->isChecked())
        ui->kd_min_comp->setValue(arg1);
}

void Guess_parameters::on_kd2_min_comp_valueChanged(double arg1)
{
    if(ui->kd2_check_range_comp->isChecked())
        ui->kd2_max_comp->setValue(arg1);
}

void Guess_parameters::on_kd2_max_comp_valueChanged(double arg1)
{
    if(ui->kd2_check_range_comp->isChecked())
        ui->kd2_min_comp->setValue(arg1);
}

void Guess_parameters::on_kdc_min_valueChanged(double arg1)
{
    if(ui->kdc_check_range->isChecked())
        ui->kdc_max->setValue(arg1);
}

void Guess_parameters::on_kdc_max_valueChanged(double arg1)
{
    if(ui->kdc_check_range->isChecked())
        ui->kdc_min->setValue(arg1);
}

void Guess_parameters::on_start_min_comp_valueChanged(double arg1)
{
    if(ui->start_check_range_comp->isChecked())
        ui->start_max_comp->setValue(arg1);
}

void Guess_parameters::on_start_max_comp_valueChanged(double arg1)
{
    if(ui->start_check_range_comp->isChecked())
        ui->start_min_comp->setValue(arg1);
}

void Guess_parameters::on_end_min_comp_valueChanged(double arg1)
{
    if(ui->end_check_range_comp->isChecked())
        ui->end_max_comp->setValue(arg1);
}

void Guess_parameters::on_end_max_comp_valueChanged(double arg1)
{
    if(ui->end_check_range_comp->isChecked())
        ui->end_min_comp->setValue(arg1);
}

void Guess_parameters::on_r20_check_stateChanged(int arg1)
{
    a0=arg1;
}

void Guess_parameters::on_kex_check_stateChanged(int arg1)
{
    a1=arg1;
}

void Guess_parameters::on_pb_check_stateChanged(int arg1)
{
    a2=arg1;
}

void Guess_parameters::on_dw_check_stateChanged(int arg1)
{
    a3=arg1;
}

void Guess_parameters::on_r20_min_valueChanged(double arg1)
{
    if(ui->r20_check_r->isChecked())
        ui->r20_max->setValue(arg1);
}

void Guess_parameters::on_r20_max_valueChanged(double arg1)
{
    if(ui->r20_check_r->isChecked())
        ui->r20_min->setValue(arg1);
}

void Guess_parameters::on_kex_min_valueChanged(double arg1)
{
    if(ui->kex_check_r->isChecked())
        ui->kex_max->setValue(arg1);
}

void Guess_parameters::on_kex_max_valueChanged(double arg1)
{
    if(ui->kex_check_r->isChecked())
        ui->kex_min->setValue(arg1);
}

void Guess_parameters::on_pb_min_valueChanged(double arg1)
{
    if(ui->pb_check_r->isChecked())
        ui->pb_max->setValue(arg1);
}

void Guess_parameters::on_pb_max_valueChanged(double arg1)
{
    if(ui->pb_check_r->isChecked())
        ui->pb_min->setValue(arg1);
}

void Guess_parameters::on_dw_min_valueChanged(double arg1)
{
    if(ui->dw_check_r->isChecked())
        ui->dw_max->setValue(arg1);
}

void Guess_parameters::on_dw_max_valueChanged(double arg1)
{
    if(ui->dw_check_r->isChecked())
        ui->dw_min->setValue(arg1);
}

void Guess_parameters::on_dw_check_r_stateChanged(int arg1)
{
    if(arg1==2){
        ui->dw_max->setValue(ui->dw_min->value());
        a3=true;
    }
    else
        a3=false;
}

void Guess_parameters::on_kex_check_r_stateChanged(int arg1)
{
    if(arg1==2){
        ui->kex_max->setValue(ui->kex_min->value());
        a1=true;
    }
    else
        a1=false;
}

void Guess_parameters::on_pb_check_r_stateChanged(int arg1)
{
    if(arg1==2){
        ui->pb_max->setValue(ui->pb_min->value());
        a2=true;
    }
    else
        a2=false;
}

void Guess_parameters::on_r20_check_r_stateChanged(int arg1)
{
    if(arg1==2){
            ui->r20_max->setValue(ui->r20_min->value());
            a0=true;
        }
    else
        a0=false;
}

void Guess_parameters::on_tabWidget_one_currentChanged(int)
{
    ui->kd1_check_one->setChecked(false);
    ui->start_check_one->setChecked(false);
    ui->end_check_one->setChecked(false);
    ui->kd_check->setChecked(false);
    ui->start_check_range->setChecked(false);
    ui->end_check_range->setChecked(false);
    a0=false;
    a1=false;
    a2=false;
    a3=false;
    a4=false;
    a5=false;
}

void Guess_parameters::on_tabWidget_two_currentChanged(int)
{
    ui->kd1_check_two->setChecked(false);
    ui->kd2_check_two->setChecked(false);
    ui->start_check_two->setChecked(false);
    ui->end_check_two->setChecked(false);
    ui->intermediate_check_two->setChecked(false);
    ui->kd_check_range_2->setChecked(false);
    ui->kd2_check_range->setChecked(false);
    ui->start_check_range_2->setChecked(false);
    ui->intermediate_check_range->setChecked(false);
    ui->end_check_range_2->setChecked(false);
    a0=false;
    a1=false;
    a2=false;
    a3=false;
    a4=false;
    a5=false;
}

void Guess_parameters::on_tabWidget_comp_currentChanged(int)
{
    ui->kdc_check->setChecked(false);
    ui->kdc_check_range->setChecked(false);
    ui->start_check_range_comp->setChecked(false);
    ui->start_check_two_2->setChecked(false);
    ui->kd1_check_two_2->setChecked(false);
    ui->kd_check_range_comp->setChecked(false);
    ui->kd2_check_two_2->setChecked(false);
    ui->kd2_check_range_comp->setChecked(false);
    ui->end_check_range_comp->setChecked(false);
    ui->end_check_two_2->setChecked(false);
    a0=false;
    a1=false;
    a2=false;
    a3=false;
    a4=false;
    a5=false;
}

void Guess_parameters::on_tabWidget_cpmg_currentChanged(int)
{
    ui->r20_check->setChecked(false);
    ui->r20_check_r->setChecked(false);
    ui->kex_check->setChecked(false);
    ui->kex_check_r->setChecked(false);
    ui->dw_check->setChecked(false);
    ui->dw_check_r->setChecked(false);
    ui->pb_check->setChecked(false);
    ui->pb_check_r->setChecked(false);
    a0=false;
    a1=false;
    a2=false;
    a3=false;
    a4=false;
    a5=false;
}
