#include "plis.h"
#include "ui_plis.h"
#include "save_open.h"
#include "legend.h"
#include "calculations.h"

Plis::Plis(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Plis)
{
    ui->setupUi(this);
    //Hide ligand concentration and protein concentration unit
    ui->ligan_conc_unit->hide();
    ui->label_3->hide();
    ui->unit->removeItem(1);
    ui->unit->removeItem(1);

    Data* temp=new Data{};
    for (int i=0; i<9;i++)
    {
        temp->concVector.push_back(0.0);
        temp->responceVector.push_back(0.0);
        temp->volumeVector.push_back(0.0);
        temp->protein_conc_vector.push_back(0.0);
        temp->error_vector.push_back(1.e-4);
    }
    temp->curve_visible=false;
    temp->data_visible=true;
    temp->name="Dataset 1";
    dataSets.push_back(temp);
    graphSets.push_back(nullptr);

    write_In_Table();
    current_dataSet=0;
    ui->selected_dataSet->addItem("Dataset 1");

    set_scatterstyle();

    //Setup plot-area
    ui->actionStandard->setDisabled(true);
    s_action=ui->menuFit->actions();
    ui->frame_2->setStyleSheet("background-color: rgb(0, 0, 0);");

    ui->plot_area->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
                                    QCP::iSelectLegend | QCP::iSelectPlottables | QCP::iSelectItems);
    xlabel="Concentration of Ligand (µM)";
    ylabel="Responce";
    ui->protein_conc_square->setValidator(v);
    ui->table->setCurrentCell(0,0);
    ui->plot_area->xAxis->setLabel(xlabel);
    ui->plot_area->yAxis->setLabel(ylabel);
    ui->plot_area->xAxis->setRange(0,100);
    ui->plot_area->yAxis->setRange(0,500);
    ui->plot_area->legend->setVisible(true);
    ui->plot_area->legend->setFont(QFont("Helvetica",9));
    ui->plot_area->legend->setRowSpacing(-3);
    ui->plot_area->legend->setSelectableParts(QCPLegend::spItems);
    ui->plot_area->plotLayout()->insertRow(0);
    title=new QCPTextElement(ui->plot_area, "Title (double click to change)", QFont("Arial",14 ,QFont::Bold));
    ui->plot_area->plotLayout()->addElement(0,0, title);
    //ContextMenu, rightclicked
    ui->plot_area->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->plot_area, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));
    connect(ui->plot_area, SIGNAL(legendClick(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)),
            this, SLOT(selection_changed(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)));
    connect(ui->plot_area, SIGNAL(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)),
            this, SLOT(legend_dubbel_clicked(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)));
    connect(ui->plot_area, SIGNAL(plottableDoubleClick(QCPAbstractPlottable *, int, QMouseEvent *)), this, SLOT(graph_double_clicked(QCPAbstractPlottable*, int)));
    connect(ui->plot_area, SIGNAL(plottableClick(QCPAbstractPlottable *, int, QMouseEvent *)), this, SLOT(graph_clicked(QCPAbstractPlottable*, int)));
    connect(title, SIGNAL(doubleClicked(QMouseEvent*)), this, SLOT(change_title(QMouseEvent*)));
    connect(ui->plot_area, SIGNAL(axisDoubleClick(QCPAxis*,QCPAxis::SelectablePart,QMouseEvent*)),
            this,SLOT(axis_double_clicked(QCPAxis *, QCPAxis::SelectablePart, QMouseEvent *)));

#ifdef Q_OS_MAC
    title->setFont(QFont("Helvetica",24));
    xaxis_font=QFont("Helvetica",14);
    yaxis_font=QFont("Helvetica",14);
    ui->plot_area->legend->setFont(QFont("Helvetica", 14));
    ui->plot_area->legend->setSelectedFont(QFont("Helvetica", 14));
    QFont font{"Helvetica",16,QFont::Bold};
    font.setUnderline(true);
    ui->label_4->setFont(font);
    QFont f{"Arial Rounded MT Bold",12};
    f.setUnderline(true);
    ui->label->setFont(f);
    ui->label_2->setFont(f);
    ui->label_3->setFont(f);
#endif
    render_graphs();
}
Plis::~Plis()
{
    delete ui;
    delete v;
    delete sim;
    delete sim_cpmg;
    delete help;
}

//Toolbar clicked **************************************
void Plis::on_actionImport_data_triggered()
{
    QFileDialog dialog(this);
    QStringList namelist{};
    QMessageBox::StandardButton pop_up;

    dialog.setDirectory(QDir::current());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    dialog.setNameFilter("(*.txt)");
    if(dialog.exec()){
        if (dataSets.size()!=0){
            pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                        "Do you want to overrite the current dataset?\nThis step cannot be undone!",
                        QMessageBox::Yes | QMessageBox::No);
            if (pop_up==QMessageBox::Yes){
                dataSets.remove(current_dataSet);
                graphSets.remove(current_dataSet);
            }
        }
        for (int i=0; i<dialog.selectedFiles().size();i++){ //For-loop to add all dataset to a vector
            QString path{dialog.selectedFiles().at(i)};
            std::string f = path.toStdString();
            Data* temp=new Data;
            temp->read_Data(f);
            bool ok{};
            QString input{""};
            while (input=="")
            {
                input = QInputDialog::getText(this, "Please, name the dataset:", path, QLineEdit::Normal, QString(), &ok);
                if (!ok){
                    break;
                }
            }
            temp->curve_visible=false;
            temp->data_visible=true;
            temp->name=input;
            if (!ok)
                delete temp;
            else{
                namelist.push_back(input);
                dataSets.push_back(temp);
                int temp{current_dataSet};
                current_dataSet=dataSets.size()-1;
                set_scatterstyle();
                current_dataSet=temp;
                graphSets.push_back(nullptr);
            }
        }
    }
    current_dataSet=dataSets.size()-1;

    if(namelist.size()!=0)
        ui->selected_dataSet->addItems(namelist);

    if (pop_up==QMessageBox::Yes){
        if (ui->selected_dataSet->count()==1)
            ui->selected_dataSet->setItemText(0,"<None>");
        else
            ui->selected_dataSet->removeItem(current_dataSet);
    }

    if (current_dataSet!=-1){
        ui->selected_dataSet->setCurrentIndex(current_dataSet);
        write_In_Table();
    }
    else{
        ui->table->setRowCount(0);
        disconnect_when_no_dataset(true);
    }

    if (is_disconnected==true && current_dataSet!=-1)
        disconnect_when_no_dataset(false);
    render_graphs();
}

void Plis::on_actionExit_triggered()
{
    QMessageBox::StandardButton dialog;
    dialog = QMessageBox::warning(this, "PLIS",
                 "Do you want to quit?\nAny unsaved changes will be lost.",
                 QMessageBox::Ok | QMessageBox::Cancel);
    if(dialog == QMessageBox::Ok)
        QApplication::quit();
}

void Plis::closeEvent(QCloseEvent *event)
{
    QMessageBox::StandardButton dialog;
    dialog = QMessageBox::warning(this, "PLIS",
                 "Do you want to quit?\nAny unsaved changes will be lost.",
                 QMessageBox::Ok | QMessageBox::Cancel);
    if(dialog == QMessageBox::Ok)
        QApplication::quit();
    else
        event->ignore();
}

void Plis::on_actionAdd_dataset_triggered()
{
    if (is_disconnected==true)
        disconnect_when_no_dataset(false);
    Data* temp=new Data;
    bool ok{};
    QString input = QInputDialog::getText(this, "Would you like to change name of your dataset(s)?", "New dataset:", QLineEdit::Normal, QString(), &ok);
    if(!ok)
        delete temp;
    else{
        if(input!=""){
            input=get_name(input);
        }
        if(mode=="cpmg"){
            temp->n_cpmgVector.push_back(0);
            temp->R2effVector.push_back(0);
            temp->dyVector.push_back(0);
            temp->responceVector.clear();
            temp->concVector.clear();
            temp->volumeVector.clear();
        }
        dataSets.push_back(temp);
        graphSets.push_back(nullptr);
        temp->data_visible=true;
        temp->curve_visible=false;
        if (dataSets.size()!=1){
            if (input==""){
                QString new_name{};
                new_name=get_name("New data");
                ui->selected_dataSet->addItem(new_name);
                dataSets[dataSets.size()-1]->name=new_name;
            }
            else{
                ui->selected_dataSet->addItem(input);
                dataSets[dataSets.size()-1]->name=input;
            }
        }
        else{
            if (input==""){
                QString new_name{};
                new_name=get_name("New data");
                ui->selected_dataSet->addItem(new_name);
                dataSets[dataSets.size()-1]->name=new_name;
            }
            else{
                ui->selected_dataSet->addItem(input);
                dataSets[dataSets.size()-1]->name=input;
            }
            ui->selected_dataSet->removeItem(0);
        }
        write_In_Table();
        current_dataSet=dataSets.size()-1;
        set_scatterstyle();
        ui->selected_dataSet->setCurrentIndex(current_dataSet);
        ui->table->setCurrentCell(0,0);
    }
}

void Plis::on_actionSimulate_data_triggered()
{
    int done{-1};
    if(mode=="standard"){
        sim->setup(dataSets,mode);
        done=sim->exec();
    }
    else if(mode=="cpmg"){
        sim_cpmg->setup(dataSets,mode);
        done=sim_cpmg->exec();
    }
    if(done==1){
        if(mode=="standard"){
            dataSets.push_back(sim->dataSets_local.last());
            current_dataSet=sim->current_dataSet;
        }
        else{
            dataSets.push_back(sim_cpmg->dataSets_local.last());
            current_dataSet=sim_cpmg->current_dataSet;
        }
        set_scatterstyle();
        graphSets.push_back(nullptr);

        QString name=get_name("Simulated data");
        ui->selected_dataSet->addItem(name);
        dataSets[current_dataSet]->name=name;
        if (ui->selected_dataSet->currentText()=="<None>")
            ui->selected_dataSet->removeItem(0);
        ui->selected_dataSet->setCurrentIndex(current_dataSet);
        ui->table->setCurrentCell(0,0);
        if (is_disconnected==true){
            disconnect_when_no_dataset(false);
            is_disconnected=false;
        }
        write_In_Table();
        render_graphs();
    }
}


void Plis::on_actionRemove_dataset_triggered()
{
    QMessageBox::StandardButton pop_up;
    pop_up=QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies", "Warning: Are you sure that you want to remove the selected dataset?",
                                QMessageBox::Yes | QMessageBox::No);
    if (pop_up==QMessageBox::Yes){
        delete_dataSets(current_dataSet);
    }
    render_graphs();
}

void Plis::on_actionOpen_project_triggered(){
    QMessageBox::StandardButton pop_up;
    pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                "The current project will be closed. Any unsaved data will be lost",
                QMessageBox::Ok | QMessageBox::Cancel);
    int i{0};
    if(pop_up==QMessageBox::Ok){
        QString temp_mode{};
        i=open_project(this, xlabel, ylabel, dataSets, current_dataSet, is_disconnected, graphSets, shape, title, xaxis_font, yaxis_font,temp_mode);
        if(i==0){
            if(temp_mode=="cpmg" && mode=="standard"){
                change_mode();
                ui->protein_conc_square->setText(QString::number(dataSets[current_dataSet]->ligand_cpmg));
            }
            else if(temp_mode=="standard" && mode=="cpmg")
                change_mode();
            is_disconnected=false;
            disconnect_when_no_dataset(false);
            ui->selected_dataSet->disconnect();
            ui->selected_dataSet->clear();
            for (int i=0;i<dataSets.size();i++){
                ui->selected_dataSet->addItem(dataSets[i]->name);
                if (dataSets[i]->has_curve==true){
                    ui->plot_area->addGraph();
                    graphSets[i]=ui->plot_area->graph();
                }
            }
            connect(ui->selected_dataSet, SIGNAL(currentIndexChanged(int)), this, SLOT(on_selected_dataSet_currentIndexChanged(int)));
        }
    }
    if(i==0){
    render_graphs();
    write_In_Table();
    write_In_Result();
    if(dataSets[current_dataSet]->has_curve)
        ui->menuOptimize_fit->setDisabled(false);
    }
}

void Plis::on_actionSave_image_triggered()
{
    QString saveFilename = QFileDialog::getSaveFileName(this, "Save as", dataSets[current_dataSet]->name, "PNG(*.png);; TIFF(*.tiff *.tif);; "
                                                                                 "JPEG(*.jpg *.jpeg);; PDF(*.pdf)");
    if (saveFilename=="")
        return;
    QPixmap pixmap(ui->plot_area->size());
    ui->plot_area->render(&pixmap, QPoint(), QRegion());
    if (!saveFilename.contains(".pdf"))
        pixmap.save(saveFilename);
    else
        ui->plot_area->savePdf(saveFilename, true);
}

void Plis::on_actionSave_project_triggered()
{
    save_project(this, xlabel, ylabel, dataSets, current_dataSet, is_disconnected, shape, title, xaxis_font, yaxis_font,mode);
}

void Plis::on_actionPaste_data_to_table_triggered()
{
    istringstream iss{QApplication::clipboard()->text().toStdString()};
    double test;
    std::string line;
    int row{ui->table->currentRow()};
    int start_column{ui->table->currentColumn()};
    int column=start_column;
    while (!iss.eof()){
        std::getline(iss,line);
        std::replace(line.begin(), line.end(), ',','.');
        istringstream iss2{line};
        while (!iss2.eof()){
            if (column>2)
                break;
            QTableWidgetItem *theItem = new QTableWidgetItem();
            iss2 >> test;
            if (row>dataSets[current_dataSet]->responceVector.size()-1){
                dataSets[current_dataSet]->responceVector.push_back(0.0);
                dataSets[current_dataSet]->concVector.push_back(0.0);
                dataSets[current_dataSet]->volumeVector.push_back(0.0);
                dataSets[current_dataSet]->protein_conc_vector.push_back(0.0);
                dataSets[current_dataSet]->error_vector.push_back(1.e-4);
                if (column==0)
                    dataSets[current_dataSet]->responceVector[row]=test;
                else if(column==1)
                    dataSets[current_dataSet]->concVector[row]=test;
                else if(column==2)
                    dataSets[current_dataSet]->volumeVector[row]=test;
                ui->table->setRowCount(dataSets[current_dataSet]->responceVector.size());
            }
            theItem->setData(Qt::EditRole, test);
            ui->table->setItem(row,column, theItem);
            column++;
        }
        row++;
        column=start_column;
    }
    dataSets[current_dataSet]->update_protein_conc();
}

void Plis::on_actionRemove_all_datasets_triggered()
{
    QMessageBox::StandardButton pop_up;
    pop_up=QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies", "Warning: Are you sure that you want to remove all datasets?",
                                QMessageBox::Yes | QMessageBox::No);
    if (pop_up==QMessageBox::No)
        return;
    ui->plot_area->clearGraphs();
    ui->plot_area->replot();
    int temp{dataSets.size()};
    for (int i=0;i<temp ;i++)
        delete_dataSets(0);
}

void Plis::on_actionOne_bindning_site_triggered()
{
    if(!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "You already have a fit made to this dataset", QMessageBox::Ok);
        return;
    }
    if(std::find(dataSets[current_dataSet]->volumeVector.begin(), dataSets[current_dataSet]->volumeVector.end(), 0)
            != dataSets[current_dataSet]->volumeVector.end()){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "Make sure all volumes are larger than zero!", QMessageBox::Ok);
        return;
    }
    if (dataSets[current_dataSet]->protein_conc<=0.0){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                             "Don't forget to fill in a protein concentration.", QMessageBox::Ok);
        return;
    }
    dataSets[current_dataSet]->a.clear();
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->responceVector[0]);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->responceVector.last());
    dataSets[current_dataSet]->a.push_back(-1);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->concVector[0]);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->protein_conc);
    dataSets[current_dataSet]->model="One bind site";
    dataSets[current_dataSet]->num_bind_site=1;

    dataSets[current_dataSet]->dyda.clear();
    for (unsigned int i=0;i<dataSets[current_dataSet]->a.size();i++)
        dataSets[current_dataSet]->dyda.push_back(0);

    QMessageBox::StandardButton pop_up;
    int i{1};
    Guess_parameters guess{this,dataSets[current_dataSet]};
    while (1){
        pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                                      "Do you want to estimate the parameters?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (pop_up==QMessageBox::Yes){
            i=guess.exec();
            break;
        }
        else if(pop_up==QMessageBox::No){
            guess_parameters_one_site(dataSets[current_dataSet], kd_testVector);
            break;
        }
    }
    if(i==1){
        Fitmrq curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, fit_one_site_dilution, fit_one_site_dilution};
        if(guess.a0==true)
            curve_fit.ia[0]=false;
        if(guess.a1==true)
            curve_fit.ia[1]=false;
        if(guess.a2==true)
            curve_fit.ia[2]=false;
        curve_fit.ia[3]=false;
        curve_fit.ia[4]=false;

        make_fitmrq(curve_fit);

        Jackknife jack{dataSets[current_dataSet]};
        jack.compute();
        write_In_Result();
    }
    else {
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
}

void Plis::on_actionTwo_binding_sites_triggered()
{
    if(!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "You already have a fit made to this dataset",QMessageBox::Ok);
        return;
    }
    if(std::find(dataSets[current_dataSet]->volumeVector.begin(), dataSets[current_dataSet]->volumeVector.end(), 0)
            != dataSets[current_dataSet]->volumeVector.end()){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "Make sure all volumes are larger than zero!", QMessageBox::Ok);
        return;
    }
    if (dataSets[current_dataSet]->protein_conc<=0.0){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                             "Don't forget to fill in a protein concentration.",QMessageBox::Ok);
        return;
    }
    dataSets[current_dataSet]->a.clear();
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->responceVector[0]);
    for (int i=0;i<5;i++)
        dataSets[current_dataSet]->a.push_back(-1);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->concVector[0]);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->protein_conc);
    dataSets[current_dataSet]->num_bind_site=2;
    dataSets[current_dataSet]->model="Two bind site";

    dataSets[current_dataSet]->dyda.clear();
    for (unsigned int i=0;i<dataSets[current_dataSet]->a.size();i++)
        dataSets[current_dataSet]->dyda.push_back(0);

    QMessageBox::StandardButton pop_up;
    Guess_parameters guess{this, dataSets[current_dataSet]};
    int i{1};
    while (1){
        pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                                      "Do you want to estimate the parameters?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (pop_up==QMessageBox::Yes){
            i=guess.exec();
            break;
        }
        else if(pop_up==QMessageBox::No){
            QApplication::setOverrideCursor(Qt::WaitCursor);
            QProgressDialog progress("Optimizing...", "Abort Optimizing", 0, 12, this);
            progress.setCancelButton(nullptr);
            progress.setMinimumDuration(0);
            progress.setWindowModality(Qt::WindowModal);
            guess_parameters_two_sites(dataSets[current_dataSet], kd_testVector,progress);
            QApplication::setOverrideCursor(Qt::ArrowCursor);
            break;
        }
    }
    if(i==1){
        Fitmrq curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, fit_two_sites_dilution,fit_two_sites_dilution};
        if(guess.a0==true)
            curve_fit.ia[0]=false;
        if(guess.a1==true)
            curve_fit.ia[1]=false;
        if(guess.a2==true)
            curve_fit.ia[2]=false;
        if(guess.a3==true)
            curve_fit.ia[3]=false;
        if(guess.a4==true)
            curve_fit.ia[4]=false;
        curve_fit.ia[5]=false;
        curve_fit.ia[6]=false;
        make_fitmrq(curve_fit);
        Jackknife jack{dataSets[current_dataSet]};
        jack.compute();
        write_In_Result();
    }
    else{
        dataSets[current_dataSet]->num_bind_site=0;
        dataSets[current_dataSet]->model="No fit made";
    }
}

void Plis::on_actionCompetitive_binding_triggered()
{
    if(!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "You already have a fit made to this dataset", QMessageBox::Ok);
        return;
    }
    if(std::find(dataSets[current_dataSet]->volumeVector.begin(), dataSets[current_dataSet]->volumeVector.end(), 0)
            != dataSets[current_dataSet]->volumeVector.end()){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "Make sure all volumes are larger than zero!", QMessageBox::Ok);
        return;
    }
    if (dataSets[current_dataSet]->protein_conc<=0.0){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                             "Don't forget to fill in a protein concentration.", QMessageBox::Ok);
        return;
    }
    double comp_conc{};
    bool ok{};
    comp_conc=QInputDialog::getDouble(this, "Plis",
                                      "Please enter the concentration of the competitor (µM):",10,0.0001, 1000, 4,&ok);
    if (!ok){
        return;
    }
    dataSets[current_dataSet]->update_comp_vector(comp_conc);
    dataSets[current_dataSet]->a.clear();
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->responceVector[0]);
    dataSets[current_dataSet]->a.push_back(-1);
    dataSets[current_dataSet]->a.push_back(-1);
    dataSets[current_dataSet]->a.push_back(-1);
    dataSets[current_dataSet]->a.push_back(-1);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->concVector[0]);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->protein_conc);
    dataSets[current_dataSet]->a.push_back(comp_conc);
    dataSets[current_dataSet]->num_bind_site=3;
    dataSets[current_dataSet]->model="Competitive binding";
    dataSets[current_dataSet]->dyda.clear();
    for (unsigned int i=0;i<dataSets[current_dataSet]->a.size();i++)
        dataSets[current_dataSet]->dyda.push_back(0);

    QMessageBox::StandardButton pop_up;
    Guess_parameters guess{this, dataSets[current_dataSet]};
    int i{1};
    while (1){
        pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                                      "Do you want to estimate the parameters?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (pop_up==QMessageBox::Yes){
            i=guess.exec();
            break;
        }
        else if(pop_up==QMessageBox::No){
            QApplication::setOverrideCursor(Qt::WaitCursor);
            QProgressDialog progress("Optimizing...", "Abort Optimizing", 0, 12, this);
            progress.setCancelButton(nullptr);
            progress.setMinimumDuration(0);
            progress.setWindowModality(Qt::WindowModal);
            guess_parameters_comp(dataSets[current_dataSet], kd_testVector,progress);
            QApplication::setOverrideCursor(Qt::ArrowCursor);
            break;
        }
    }
    if(i==1){
        Fitmrq2 curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->comp_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, fit_comp, fit_comp};
        if(guess.a0==true)
            curve_fit.ia[0]=false;
        if(guess.a1==true)
            curve_fit.ia[1]=false;
        if(guess.a2==true)
            curve_fit.ia[2]=false;
        if(guess.a3==true)
            curve_fit.ia[3]=false;
        if(guess.a4==true)
            curve_fit.ia[4]=false;

        curve_fit.ia[5]=false;
        curve_fit.ia[6]=false;
        curve_fit.ia[7]=false;
        make_fitmrq2(curve_fit);
        Jackknife jack{dataSets[current_dataSet]};
        jack.compute();
        write_In_Result();
    }
    else{
        dataSets[current_dataSet]->num_bind_site=0;
        dataSets[current_dataSet]->model="No fit made";
    }
}

void Plis::CPMG(){
    if(!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "You already have a fit made to this dataset",
                             QMessageBox::Ok);
        return;
    }
    if (dataSets[current_dataSet]->ligand_cpmg<=0.0){
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                             "Don't forget to fill in a ligand concentration.", QMessageBox::Ok);
        return;
    }
    dataSets[current_dataSet]->a.clear();
    dataSets[current_dataSet]->a.push_back(0);
    dataSets[current_dataSet]->a.push_back(0);
    dataSets[current_dataSet]->a.push_back(0);
    dataSets[current_dataSet]->a.push_back(0);
    dataSets[current_dataSet]->a.push_back(dataSets[current_dataSet]->ligand_cpmg);
    dataSets[current_dataSet]->dyda.clear();
    for (unsigned int i=0;i<dataSets[current_dataSet]->a.size();i++)
        dataSets[current_dataSet]->dyda.push_back(0);

    dataSets[current_dataSet]->model="CPMG-model";
    dataSets[current_dataSet]->num_bind_site=4;
    QMessageBox::StandardButton pop_up;
    Guess_parameters guess{this, dataSets[current_dataSet]};
    int i{1};
    while (1){
        pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                                      "Do you want to estimate the parameters?",
                                      QMessageBox::Yes | QMessageBox::No);
        if (pop_up==QMessageBox::Yes){
            i=guess.exec();
            break;
        }
        else if(pop_up==QMessageBox::No){
            QApplication::setOverrideCursor(Qt::WaitCursor);
            QProgressDialog progress("Optimizing...", "Abort Optimizing", 0, 12, this);
            progress.setCancelButton(nullptr);
            progress.setMinimumDuration(0);
            progress.setWindowModality(Qt::WindowModal);
            guess_parameters_cpmg(dataSets[current_dataSet],progress);
            QApplication::setOverrideCursor(Qt::ArrowCursor);
            break;
        }
    }
    if(i==0){
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
        return;
    }

    Fitmrq3 fit{dataSets[current_dataSet]->n_cpmgVector.toStdVector(),dataSets[current_dataSet]->R2effVector.toStdVector(),
                dataSets[current_dataSet]->dyVector.toStdVector(), dataSets[current_dataSet]->a,
                carverrichards,carverrichards};
    if (guess.a0)
        fit.ia[0]=false;
    if (guess.a1)
        fit.ia[1]=false;
    if (guess.a2)
        fit.ia[2]=false;
    if (guess.a3)
        fit.ia[3]=false;
    fit.ia[4]=false;
    int did_fit_work=fit.fit();

    if (did_fit_work==1){ //singular matrix!
        QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
    else if(did_fit_work==2){
        QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
    else{
        dataSets[current_dataSet]->a=fit.a;
        dataSets[current_dataSet]->chi2=fit.chisq;
        dataSets[current_dataSet]->kd_cpmg=CPMG_kd(dataSets[current_dataSet]->a);
        std::vector<std::vector<double>> temp_vector{};
        temp_vector=fit.plot();
        dataSets[current_dataSet]-> calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
        dataSets[current_dataSet]-> calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
        ui->plot_area->addGraph();
        graphSets[current_dataSet]=ui->plot_area->graph();
        dataSets[current_dataSet]->curve_visible=true;
        render_graphs();
        Jackknife jack{dataSets[current_dataSet]};
        jack.compute_cpmg();
        dataSets[current_dataSet]->has_curve=true;
        write_In_Result();
        ui->menuOptimize_fit->setDisabled(false);
    }
}

void Plis::on_actionRemove_a_fitted_curve_triggered()
{
    if (dataSets[current_dataSet]->has_curve){
        QMessageBox::StandardButton pop_up;
        pop_up=QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies", "Warning: Are you sure that you want to remove the fitted curve?",
                                    QMessageBox::Yes | QMessageBox::No);
        if (pop_up==QMessageBox::No)
            return;
        QCPGraph * temp=graphSets[current_dataSet];
        graphSets[current_dataSet]=nullptr;
        ui->plot_area->removePlottable(temp);
        dataSets[current_dataSet]->has_curve=false;
        dataSets[current_dataSet]->curve_visible=false;
        ui->menuOptimize_fit->setDisabled(true);
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
        dataSets[current_dataSet]->a[2]=0;
        write_In_Result();
        render_graphs();
    }
}

void Plis::on_actionRemove_all_fitted_curves_triggered()
{
    QMessageBox::StandardButton pop_up;
    pop_up=QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies", "Warning: Are you sure that you want to remove all fitted curve?",
                                QMessageBox::Yes | QMessageBox::No);
    if (pop_up==QMessageBox::No)
        return;
    for (int i=0;i<graphSets.size();i++){
        if(graphSets[i]!=nullptr){
            ui->plot_area->removePlottable(graphSets[i]);
            graphSets[i]=nullptr;
            dataSets[i]->has_curve=false;
            dataSets[i]->curve_visible=false;
            dataSets[i]->model="No fit made";
            dataSets[i]->num_bind_site=0;
            dataSets[i]->a[2]=0;
        }
    }
    ui->menuOptimize_fit->setDisabled(true);
    write_In_Result();
    render_graphs();
}

//Something on screen clicked *****************************************
void Plis::on_add_row_button_clicked()
{
    if(ui->table->currentRow()==-1)
        ui->table->setCurrentCell(0,0);
    if(mode=="standard"){
        dataSets[current_dataSet]->responceVector.insert(ui->table->currentRow(),0.0);
        dataSets[current_dataSet]->volumeVector.insert(ui->table->currentRow(),0.0);
        dataSets[current_dataSet]->concVector.insert(ui->table->currentRow(),0.0);
        dataSets[current_dataSet]->error_vector.insert(ui->table->currentRow(),1.e-4);
        dataSets[current_dataSet]->protein_conc_vector.insert(ui->table->currentRow(),0.0);
        ui->table->insertRow(ui->table->currentRow());
        dataSets[current_dataSet]->update_protein_conc();
    }
    else{
        dataSets[current_dataSet]->n_cpmgVector.insert(ui->table->currentRow(),0.0);
        dataSets[current_dataSet]->R2effVector.insert(ui->table->currentRow(),0.0);
        dataSets[current_dataSet]->dyVector.insert(ui->table->currentRow(),0.0);
        ui->table->insertRow(ui->table->currentRow());
    }
    write_In_Table();
}

void Plis::on_remove_row_button_clicked()
{
    if (ui->table->currentRow()==-1)
        ui->table->setCurrentCell(0,0);
    int row{ui->table->currentRow()},col{ui->table->currentColumn()};
    if (dataSets[current_dataSet]->responceVector.size()!=1 && mode=="standard"){
        dataSets[current_dataSet]->responceVector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->volumeVector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->protein_conc_vector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->concVector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->error_vector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->update_protein_conc();
        ui->table->removeRow(ui->table->currentRow());
    }
    else if(dataSets[current_dataSet]->n_cpmgVector.size()!=1 && mode=="cpmg"){
        dataSets[current_dataSet]->n_cpmgVector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->R2effVector.remove(ui->table->currentRow());
        dataSets[current_dataSet]->dyVector.remove(ui->table->currentRow());
        ui->table->removeRow(ui->table->currentRow());
    }
    write_In_Table();
    ui->table->setCurrentCell(row,col);
}

void Plis::on_table_cellChanged(int row, int column)
{
    QString text{ui->table->item(row,column)->text()};
    double changed_value{text.toDouble()};
    if(mode=="standard"){
        if (column==0)
            dataSets[current_dataSet]->responceVector[row]=changed_value;
        else if (column==1)
            dataSets[current_dataSet]->concVector[row]=changed_value;
        else{
            dataSets[current_dataSet]->volumeVector[row]=changed_value;
            dataSets[current_dataSet]->update_protein_conc();
        }
    }
    else if(mode=="cpmg"){
        if (column==0)
            dataSets[current_dataSet]->n_cpmgVector[row]=changed_value;
        else if (column==1)
            dataSets[current_dataSet]->R2effVector[row]=changed_value;
        else{
            dataSets[current_dataSet]->dyVector[row]=changed_value;
        }
    }
}

void Plis::on_selected_dataSet_currentIndexChanged(int index)
{
    current_dataSet=index;
    write_In_Table();
    if(mode=="standard"){
        if (dataSets[current_dataSet]->ligand_unit=="uM")
            ui->table->horizontalHeaderItem(1)->setText("[Ligand] (µM)");
        else
            ui->table->horizontalHeaderItem(1)->setText("[Ligand] ("+dataSets[current_dataSet]->ligand_unit+")");
    }

    ui->plot_area->xAxis->setLabel(xlabel);
    if (dataSets[current_dataSet]->ligand_unit=="uM")
        ui->ligan_conc_unit->setCurrentIndex(0);
    else if(dataSets[current_dataSet]->ligand_unit=="mM")
        ui->ligan_conc_unit->setCurrentIndex(1);
    else if(dataSets[current_dataSet]->ligand_unit=="nM")
        ui->ligan_conc_unit->setCurrentIndex(2);

    if (dataSets[current_dataSet]->unit=="uM")
        ui->unit->setCurrentIndex(0);
    else if(dataSets[current_dataSet]->unit=="mM")
        ui->unit->setCurrentIndex(1);
    else if(dataSets[current_dataSet]->unit=="nM")
        ui->unit->setCurrentIndex(2);
    if(mode=="standard")
        ui->protein_conc_square->setText(QString::number(dataSets[current_dataSet]->protein_conc));
    else if(mode=="cpmg")
        ui->protein_conc_square->setText(QString::number(dataSets[current_dataSet]->ligand_cpmg));

    if(dataSets[current_dataSet]->has_curve)
        ui->menuOptimize_fit->setDisabled(false);
    else
        ui->menuOptimize_fit->setDisabled(true);

    write_In_Result();
}

void Plis::on_protein_conc_square_editingFinished()
{
    if(current_dataSet!=(-1)){
        QString s=ui->protein_conc_square->text();
        s.replace(',','.');
        if(mode=="standard"){
            dataSets[current_dataSet]->protein_conc=s.toDouble();
            dataSets[current_dataSet]->update_protein_conc();
        }
        else if(mode=="cpmg")
            dataSets[current_dataSet]->ligand_cpmg=s.toDouble();
    }
}

void Plis::on_unit_currentIndexChanged(const QString &arg1)
{
    if (arg1!=dataSets[current_dataSet]->unit){
        dataSets[current_dataSet]->unit=arg1;
        render_graphs();
    }
}

void Plis::on_ligan_conc_unit_currentIndexChanged(const QString &arg1)
{
    if (arg1!=dataSets[current_dataSet]->ligand_unit){
        dataSets[current_dataSet]->ligand_unit=arg1;
        ui->table->horizontalHeaderItem(1)->setText("[Ligand] ("+dataSets[current_dataSet]->ligand_unit+")");
        write_In_Table();
    }
    if (current_dataSet==0)
        ui->plot_area->xAxis->setLabel("Concentration of Ligand ("+dataSets[current_dataSet]->ligand_unit+")");
    render_graphs();
    write_In_Result();
}

void Plis::on_update_graph_button_clicked()
{
    render_graphs();
}

void Plis::selection_changed(QCPLegend* l,QCPAbstractLegendItem* ai,QMouseEvent* me){
    legend_clicked(l, ai, me, dataSets, graphSets);
}

void Plis::legend_dubbel_clicked(QCPLegend* l,QCPAbstractLegendItem* ai,QMouseEvent* me){
    QString old_name{}, new_name{};
    bool ok{};
    if (ai!=nullptr){
        new_name=QInputDialog::getText(this, "Plis",
                                       "Please enter the new name of the dataset:", QLineEdit::Normal, QString(), &ok);
        int namebox{};
        new_name=get_name(new_name);
        if (ok && !new_name.isEmpty()){
            old_name=legend_double_clicked(l, ai, me, dataSets, new_name);
            if (old_name!=""){
                namebox=ui->selected_dataSet->findText(old_name);
                ui->selected_dataSet->setItemText(namebox, new_name);
            }
            else
                QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                                                          "You can not change the name of a fitted curve!", QMessageBox::Ok);
        }
    }
    write_In_Result();
    render_graphs();
}

void Plis::contextMenuRequest(QPoint pos)
{
    QMenu *menu = new QMenu(this);
    menu->setAttribute(Qt::WA_DeleteOnClose);

    if (ui->plot_area->legend->selectTest(pos, false) >= 0) // context menu on legend requested
    {
      menu->addAction("Move to top left", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignTop|Qt::AlignLeft));
      menu->addAction("Move to top center", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignTop|Qt::AlignHCenter));
      menu->addAction("Move to top right", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignTop|Qt::AlignRight));
      menu->addAction("Move to bottom right", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignBottom|Qt::AlignRight));
      menu->addAction("Move to bottom left", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignBottom|Qt::AlignLeft));
    } else  // general context menu on graphs requested
    {
        if (ui->plot_area->graphCount() > 0){
            menu->addAction("Remove all datasets", this, SLOT(on_actionRemove_all_datasets_triggered()));
            menu->addAction("Remove all fitted graphs", this, SLOT(on_actionRemove_all_fitted_curves_triggered()));
        }
        if (ui->plot_area->selectedPlottables().size() > 0){
            menu->addSeparator();
            menu->addAction("Remove selected dataset",this,SLOT(on_actionRemove_dataset_triggered()));
            menu->addAction("Remove selected fitted graph", this, SLOT(on_actionRemove_a_fitted_curve_triggered()));
        }
    }
    menu->popup(ui->plot_area->mapToGlobal(pos));
}

void Plis::moveLegend()
{
  if (QAction* contextAction = qobject_cast<QAction*>(sender()))
  {
    bool ok;
    int dataInt = contextAction->data().toInt(&ok);
    if (ok)
    {
      ui->plot_area->axisRect()->insetLayout()->setInsetAlignment(0,static_cast<Qt::Alignment>(dataInt));
      ui->plot_area->replot();
    }
  }
}

void Plis::change_title(QMouseEvent *){
    if (QCPTextElement *t = qobject_cast<QCPTextElement*>(sender()))
    {
        Modify_text new_title{this, title};
        new_title.exec();
        ui->plot_area->replot();
    }
}

//Help-functions ***************************************
void Plis::render_graphs(){
    ui->plot_area->clearGraphs();
    int temp_current_dataSet{current_dataSet};
    for (int i=0; i<dataSets.size();i++)
    {
        current_dataSet=i;
        if (dataSets[i]->ligand_unit!=dataSets[0]->ligand_unit && mode=="standard"){
            QVector<double> temp{dataSets[i]->concVector};
            dataSets[i]->change_unit_of_data(dataSets[0]->ligand_unit, "conc");
            write_In_Graph();
            dataSets[i]->concVector=temp;
        }
        else{
            write_In_Graph();
        }

        if (graphSets[i]!=nullptr)
            draw_curve();
    }

    int curve{};
    for (int i=0;i<dataSets.size();i++){
        if (dataSets[i]->data_visible==false){
            ui->plot_area->legend->item(i+curve)->setTextColor("Grey");
            ui->plot_area->graph(i+curve)->setVisible(false);
        }
        if(dataSets[i]->has_curve)
            curve++;
    }
    curve=0;
    for (int i=0; i<dataSets.size();i++){
        if(dataSets[i]->has_curve)
            curve++;
        if (dataSets[i]->curve_visible==false && graphSets[i]!=nullptr){
            ui->plot_area->legend->item(i+curve)->setTextColor("Grey");
            ui->plot_area->graph(i+curve)->setVisible(false);
        }
    }
    current_dataSet=temp_current_dataSet;
    ui->plot_area->rescaleAxes();
    ui->plot_area->xAxis->setLabel(xlabel);
    ui->plot_area->xAxis->setLabelFont(xaxis_font);
    ui->plot_area->xAxis->setSelectedLabelFont(xaxis_font);
    ui->plot_area->yAxis->setLabel(ylabel);
    ui->plot_area->yAxis->setLabelFont(yaxis_font);
    ui->plot_area->yAxis->setSelectedLabelFont(yaxis_font);
    ui->plot_area->replot();
}

void Plis::change_penColor(int i){
    while (i>4)
        i-=5;
    switch(i)
    {
    case 0: dataSets[current_dataSet]->pen.setColor("black");
        break;
    case 1: dataSets[current_dataSet]->pen.setColor("red");
        break;
    case 2: dataSets[current_dataSet]->pen.setColor("green");
        break;
    case 3: dataSets[current_dataSet]->pen.setColor("blue");
        break;
    case 4: dataSets[current_dataSet]->pen.setColor("magenta");
        break;
    }
}

void Plis::disconnect_when_no_dataset(bool x){
    is_disconnected=x;
    ui->actionRemove_dataset->setDisabled(x);
    ui->add_row_button->setDisabled(x);
    ui->remove_row_button->setDisabled(x);
    ui->update_graph_button->setDisabled(x);
    ui->protein_conc_square->setDisabled(x);
    ui->unit->setDisabled(x);
    ui->ligan_conc_unit->setDisabled(x);
    ui->actionRemove_all_datasets->setDisabled(x);
    ui->actionRemove_a_fitted_curve->setDisabled(x);
    ui->actionRemove_all_fitted_curves->setDisabled(x);
    ui->actionPaste_data_to_table->setDisabled(x);
    ui->pop_result_button->setDisabled(x);
    ui->menuFit->setDisabled(x);
    ui->actionPaste_data_to_table->setDisabled(x);
    ui->actionCopy_data_from_table->setDisabled(x);
    ui->actionSave_project->setDisabled(x);
}

void Plis::write_In_Table(){
    if(mode=="standard"){
        ui->table->setRowCount(dataSets[current_dataSet]->responceVector.size());
        for (int i=0;i<3;i++){
            for (int j=0; j<dataSets[current_dataSet]->responceVector.size(); j++){
                QTableWidgetItem *theItem = new QTableWidgetItem();
                if (i==0)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->responceVector.at(j));
                else if(i==1)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->concVector.at(j));
                else
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->volumeVector.at(j));
                ui->table->setItem(j,i,theItem);
            }
        }
    }
    else if(mode=="cpmg"){
        ui->table->setRowCount(dataSets[current_dataSet]->n_cpmgVector.size());
        for (int i=0;i<3;i++){
            for (int j=0; j<dataSets[current_dataSet]->n_cpmgVector.size(); j++){
                QTableWidgetItem *theItem = new QTableWidgetItem();
                if (i==0)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->n_cpmgVector.at(j));
                else if(i==1)
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->R2effVector.at(j));
                else
                    theItem->setData(Qt::EditRole, dataSets[current_dataSet]->dyVector.at(j));
                ui->table->setItem(j,i,theItem);
            }
        }
    }
}

void Plis::write_In_Graph(){
    ui->plot_area->addGraph();
    ui->plot_area->graph()->setName(dataSets[current_dataSet]->name);
    ui->plot_area->graph()->setPen(dataSets[current_dataSet]->pen);
    ui->plot_area->graph()->setScatterStyle(dataSets[current_dataSet]->style);
    ui->plot_area->graph()->setLineStyle(QCPGraph::lsNone);
    if(mode=="standard")
        ui->plot_area->graph()->setData(dataSets[current_dataSet]->concVector,dataSets[current_dataSet]->responceVector);
    else if(mode=="cpmg")
        ui->plot_area->graph()->setData(dataSets[current_dataSet]->n_cpmgVector,dataSets[current_dataSet]->R2effVector);

}

void Plis::draw_curve(){
    ui->plot_area->addGraph();
    graphSets[current_dataSet]=ui->plot_area->graph();
    graphSets[current_dataSet]->setName("Curve - " + dataSets[current_dataSet]->name);
    graphSets[current_dataSet]->setPen(dataSets[current_dataSet]->pen_curve);
    graphSets[current_dataSet]->setLineStyle(QCPGraph::lsLine);
    /*
    if (dataSets[current_dataSet]->ligand_unit!=dataSets[0]->ligand_unit){
        QVector<double> temp{dataSets[current_dataSet]->calc_data[0]};
        dataSets[current_dataSet]->change_unit_of_data(dataSets[0]->ligand_unit, "calc");
        graphSets[current_dataSet]->setData(dataSets[current_dataSet]->calc_data[0],dataSets[current_dataSet]->calc_data[1]);
        dataSets[current_dataSet]->calc_data[0]=temp;
    }
    else
    */
        graphSets[current_dataSet]->setData(dataSets[current_dataSet]->calc_data[0],dataSets[current_dataSet]->calc_data[1]);

}

void Plis::write_In_Result(){
    if(current_dataSet!=-1){
        ui->name->setText(dataSets[current_dataSet]->name);
        ui->model->setText(dataSets[current_dataSet]->model);
        if(mode=="standard")
            ui->num_data->setText(QString::number(dataSets[current_dataSet]->responceVector.size()));
        ui->chi2_box->setText(QString::number(dataSets[current_dataSet]->chi2));
        ui->kd2_box->setVisible(false);
        ui->kd2_error->setVisible(false);
        ui->kd2_label->setVisible(false);
        ui->plus_minus2->setVisible(false);
        ui->kd3_box->setVisible(false);
        ui->kd3_error->setVisible(false);
        ui->Kd3_label->setVisible(false);
        ui->plus_minus3->setVisible(false);
        ui->kd4_box->setVisible(false);
        ui->kd4_error->setVisible(false);
        ui->Kd4_label->setVisible(false);
        ui->plus_minus4->setVisible(false);
        if (dataSets[current_dataSet]->num_bind_site==1){
            if(dataSets[current_dataSet]->a[2]>0.1){
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->a[2],'f',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
            else{
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->a[2],'e',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
        }
        else if(dataSets[current_dataSet]->num_bind_site==2){
            ui->kd2_box->setVisible(true);
            ui->kd2_error->setVisible(true);
            ui->kd2_label->setVisible(true);
            ui->plus_minus2->setVisible(true);
            if(dataSets[current_dataSet]->a[3]<0.1 || dataSets[current_dataSet]->a[4]<0.1){
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->a[3],'e',2));
                ui->kd2_box->setText(QString::number(dataSets[current_dataSet]->a[4],'e',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd1_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
                ui->kd2_error->setText(QString::number(dataSets[current_dataSet]->kd2_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
            else{
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->a[3],'f',2));
                ui->kd2_box->setText(QString::number(dataSets[current_dataSet]->a[4],'f',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd1_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
                ui->kd2_error->setText(QString::number(dataSets[current_dataSet]->kd2_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
        }
        else if(dataSets[current_dataSet]->num_bind_site==3){
            ui->kd2_box->setVisible(true);
            ui->kd2_error->setVisible(true);
            ui->kd2_label->setVisible(true);
            ui->plus_minus2->setVisible(true);
            ui->kd3_box->setVisible(true);
            ui->kd3_error->setVisible(true);
            ui->Kd3_label->setVisible(true);
            ui->plus_minus3->setVisible(true);
            if(dataSets[current_dataSet]->a[2]<0.1 || dataSets[current_dataSet]->a[3]<0.1 || dataSets[current_dataSet]->a[4]<0.1){
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->a[2],'e',2));
                ui->kd2_box->setText(QString::number(dataSets[current_dataSet]->a[3],'e',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd1_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
                ui->kd2_error->setText(QString::number(dataSets[current_dataSet]->kd2_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
                ui->kd3_box->setText(QString::number(dataSets[current_dataSet]->a[4],'e',2));
                ui->kd3_error->setText(QString::number(dataSets[current_dataSet]->kdc_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
            else{
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->a[2],'f',2));
                ui->kd2_box->setText(QString::number(dataSets[current_dataSet]->a[3],'f',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd1_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
                ui->kd2_error->setText(QString::number(dataSets[current_dataSet]->kd2_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
                ui->kd3_box->setText(QString::number(dataSets[current_dataSet]->a[4],'f',2));
                ui->kd3_error->setText(QString::number(dataSets[current_dataSet]->kdc_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
            ui->Kd3_label->setText("Kdc");
        }
        else if(dataSets[current_dataSet]->num_bind_site==4){
            ui->num_data->setText(QString::number(dataSets[current_dataSet]->n_cpmgVector.size()));
            if(dataSets[current_dataSet]->kd_cpmg<0.1){
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->kd_cpmg,'e',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd_error,'e',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
            else {
                ui->kd1_box->setText(QString::number(dataSets[current_dataSet]->kd_cpmg,'f',2));
                ui->kd1_error->setText(QString::number(dataSets[current_dataSet]->kd_error,'f',2) + " " + dataSets[current_dataSet]->ligand_unit);
            }
        }
        else if(dataSets[current_dataSet]->num_bind_site==0){
            ui->kd1_error->setText("");
            ui->kd1_box->setText("");
        }
    }
    else{
        ui->kd1_box->setText("");
        ui->name->setText("");
        ui->chi2_box->setText("");
        ui->kd1_error->setText("");
        ui->model->setText("");
        ui->num_data->setText("");
        ui->kd1_box->setText("");
        ui->kd2_box->setText("");
        ui->kd2_error->setText("");
    }
}

void Plis::delete_dataSets(int i){
    delete dataSets[i];
    dataSets[i]=nullptr;
    dataSets.remove(i);
    graphSets.remove(i); //kan bli minneslucka då jag inte gör delete på den ritade graphen
    if (dataSets.size()==0){
        ui->selected_dataSet->setItemText(i, "<None>");
        current_dataSet=-1; //if no dataset exist
        disconnect_when_no_dataset(true);
        ui->table->setRowCount(0);
        ui->menuOptimize_fit->setDisabled(true);
    }
    else{
        if (i!=0)
            current_dataSet-=1;
        ui->selected_dataSet->setCurrentIndex(0);
        ui->selected_dataSet->removeItem(i);
        write_In_Table();
    }
    write_In_Result();
}

QString Plis::get_name(QString name){
    int loop{0};
    QString new_name{};
    while (1){
        if(loop==0){
            int check=ui->selected_dataSet->findText(name);
            if (check==-1){
               new_name=name;
               break;
            }
        }
        else{
            int check=ui->selected_dataSet->findText(name+" (" + QString::number(loop)+")");
            if(check==-1){
                new_name=(name+" ("+QString::number(loop)+")");
                break;
            }
        }
        loop++;
    }
    return new_name;
}

void Plis::make_fitmrq(Fitmrq &curve_fit){
    int did_fit_work{};
    did_fit_work=curve_fit.fit();
    if (did_fit_work==1){ //singular matrix!
        QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
    else if(did_fit_work==2){
        QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
    else {
        dataSets[current_dataSet]->a=curve_fit.a;
        dataSets[current_dataSet]->chi2=curve_fit.chisq;

        std::vector<std::vector<double>> temp_vector{};
        temp_vector=curve_fit.plot();
        dataSets[current_dataSet]-> calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
        dataSets[current_dataSet]-> calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
        ui->plot_area->addGraph();
        graphSets[current_dataSet]=ui->plot_area->graph();
        dataSets[current_dataSet]->curve_visible=true;
        render_graphs();
        dataSets[current_dataSet]->has_curve=true;
        ui->menuOptimize_fit->setDisabled(false);
    }
}

void Plis::make_fitmrq2(Fitmrq2 &curve_fit){
    int did_fit_work{};
    did_fit_work=curve_fit.fit();
    if (did_fit_work==1){ //singular matrix!
        QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: The mathematical matrix could not be formed to create a fit");
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
    else if(did_fit_work==2){
        QMessageBox::warning(this, "Warning: Something went wrong!", "Warning: To many iterations.\nThe best fit could not be found");
        dataSets[current_dataSet]->model="No fit made";
        dataSets[current_dataSet]->num_bind_site=0;
    }
    else {
        dataSets[current_dataSet]->a=curve_fit.a;
        dataSets[current_dataSet]->chi2=curve_fit.chisq;

        std::vector<std::vector<double>> temp_vector{};
        temp_vector=curve_fit.plot();
        dataSets[current_dataSet]-> calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
        dataSets[current_dataSet]-> calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
        ui->plot_area->addGraph();
        graphSets[current_dataSet]=ui->plot_area->graph();
        dataSets[current_dataSet]->curve_visible=true;
        render_graphs();
        dataSets[current_dataSet]->has_curve=true;
        ui->menuOptimize_fit->setDisabled(false);
    }
}

void Plis::set_scatterstyle(){
    if (shapes.at(shape) != QCPScatterStyle::ssCustom){
        dataSets[current_dataSet]->style=shapes[shape];
        dataSets[current_dataSet]->style.setSize(10);
        change_penColor(shape);
    }
    else{
        QPainterPath customScatterPath;
        for (int j=0; j<3; ++j)
            customScatterPath.cubicTo(qCos(2*M_PI*j/3.0)*9, qSin(2*M_PI*j/3.0)*9, qCos(2*M_PI*(j+0.9)/3.0)*9, qSin(2*M_PI*(j+0.9)/3.0)*9, 0, 0);
        dataSets[current_dataSet]->style=QCPScatterStyle(customScatterPath, QPen(Qt::black, 0), QColor(40, 70, 255, 50), 10);
        dataSets[current_dataSet]->style.setSize(10);
    }
    shape+=1;
    if(shape>shapes.size()-1)
        shape=0;
}

void Plis::graph_double_clicked(QCPAbstractPlottable* a, int){
    for (int i=0; i<graphSets.size();i++){
        if (graphSets[i]!=nullptr)
            if (a->name()==graphSets[i]->name()){
                ui->selected_dataSet->setCurrentIndex(i);
                current_dataSet=i;
                break;
        }
    }
    for (int j=0; j<dataSets.size();j++)
        if (a->name()==dataSets[j]->name){
            ui->selected_dataSet->setCurrentIndex(j);
            current_dataSet=j;
            break;
        }

    Modify_data mod{this, shapes, dataSets[current_dataSet]};
    mod.exec();
    render_graphs();
}

void Plis::graph_clicked(QCPAbstractPlottable* a, int){
    for (int i=0; i<graphSets.size();i++){
        if (graphSets[i]!=nullptr)
            if (a->name()==graphSets[i]->name()){
                ui->selected_dataSet->setCurrentIndex(i);
                current_dataSet=i;
                break;
        }
    }
    for (int j=0; j<dataSets.size();j++)
        if (a->name()==dataSets[j]->name){
            ui->selected_dataSet->setCurrentIndex(j);
            current_dataSet=j;
            break;
        }
}

void Plis::axis_double_clicked(QCPAxis *axis, QCPAxis::SelectablePart sp, QMouseEvent *){
    if(sp==QCPAxis::spAxisLabel){
        Modify_text text{this, axis,xaxis_font,yaxis_font,xlabel,ylabel};
        text.exec();
    }
    render_graphs();
}

void Plis::on_pop_result_button_clicked()
{
    Show_result result{this,dataSets[current_dataSet]};
    result.exec();
}

void Plis::on_actionOptimize_fit_triggered()
{
    if(!ui->menuOptimize_fit->isEnabled())
        return;
    Optimize opt{this, dataSets[current_dataSet]};
    opt.exec();
    write_In_Result();
    render_graphs();

}

void Plis::on_actionStandard_triggered()
{
    QMessageBox::StandardButton input{};
    input=QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                               "Do you want to change state?\nYour current project will be discarded!",
                               QMessageBox::Yes | QMessageBox::No);

    if(input==QMessageBox::Yes){
        ui->actionStandard->setDisabled(true);
        ui->switch_to_CPGM->setDisabled(false);
        mode="standard";
        QStringList names{"Response", "[Ligand] (µM)","Volume"};
        ui->table->setHorizontalHeaderLabels(names);
        xlabel="Concentration of Ligand";
        ylabel="Response";
        ui->menuFit->clear();
        ui->menuFit->addActions(s_action);
        ui->menuFit->setTitle("Fit");
        ui->label_2->setText("Protein conc:");
        //ui->ligan_conc_unit->show();
        //ui->label_3->show();
        ui->plot_area->clearGraphs();
        ui->plot_area->replot();
        int temp{dataSets.size()};
        for (int i=0;i<temp ;i++)
            delete_dataSets(0);
        render_graphs();
    }
}

void Plis::on_switch_to_CPGM_triggered()
{
    QMessageBox::StandardButton input{};
    input=QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                               "Do you want to change state?\nYour current project will be discarded!",
                               QMessageBox::Yes | QMessageBox::No);

    if(input==QMessageBox::Yes){
        ui->actionStandard->setDisabled(false);
        ui->switch_to_CPGM->setDisabled(true);
        mode="cpmg";
        QStringList names{"N_CPMG (Hz)", "R2eff (/s)", "Error in R2eff"};
        ui->table->setHorizontalHeaderLabels(names);
        xlabel="N_CPMG (Hz)";
        ylabel="R2eff (/s)";
        ui->menuFit->setTitle("Fit");
        ui->menuFit->clear();
        ui->menuFit->addAction("Fit",this,SLOT(CPMG()));
        ui->label_2->setText("Ligand conc:");
        ui->label_3->hide();
        ui->ligan_conc_unit->hide();
        ui->plot_area->clearGraphs();
        ui->plot_area->replot();
        int temp{dataSets.size()};
        for (int i=0;i<temp ;i++)
            delete_dataSets(0);
        render_graphs();
    }
}

void Plis::change_mode(){

    if(mode=="cpmg"){
        ui->actionStandard->setDisabled(true);
        ui->switch_to_CPGM->setDisabled(false);
        mode="standard";
        QStringList names{"Response", "[Ligand]","Volume"};
        ui->table->setHorizontalHeaderLabels(names);
        xlabel="Concentration of Ligand";
        ylabel="Response";
        ui->menuFit->clear();
        ui->menuFit->addActions(s_action);
        ui->menuFit->setTitle("Fit");
        ui->label_2->setText("Protein\nConcentration");
        //ui->ligan_conc_unit->show();
        //ui->label_3->show();
    }

    else if(mode=="standard"){
        ui->actionStandard->setDisabled(false);
        ui->switch_to_CPGM->setDisabled(true);
        mode="cpmg";
        QStringList names{"N_CPMG (Hz)", "R2eff (/s)", "Error in R2eff"};
        ui->table->setHorizontalHeaderLabels(names);
        xlabel="N_CPMG (Hz)";
        ylabel="R2eff (/s)";
        ui->menuFit->setTitle("Fit");
        ui->menuFit->clear();
        ui->menuFit->addAction("Fit",this,SLOT(CPMG()));
        ui->label_2->setText("Ligand conc:");
        //ui->label_3->hide();
        //ui->ligan_conc_unit->hide();
    }
}

void Plis::on_actionCopy_data_from_table_triggered()
{
    QApplication::clipboard()->clear();
    QList<QTableWidgetItem*> list= ui->table->selectedItems();
    for (int i=0;i<list.size();i++){
        if (list[i]->column()==0 && QApplication::clipboard()->text()=="")
            QApplication::clipboard()->setText(list[i]->text());
        else if(list[i]->column()==0)
            QApplication::clipboard()->setText(QApplication::clipboard()->text() + "\n" + list[i]->text());
        else
            QApplication::clipboard()->setText(QApplication::clipboard()->text() + " " + list[i]->text());
    }
}

void Plis::on_actionHelp_triggered()
{
    help->show();
}

void Plis::on_actionInformation_about_PLIS_triggered()
{
    help->setup(1);
    help->show();
}

void Plis::keyReleaseEvent(QKeyEvent *event){
    if(event->key()==Qt::Key_Delete){
        if(ui->plot_area->selectedPlottables().size()>0){
            for (int i=0;i<dataSets.size();i++)
                if(dataSets[i]->name==ui->plot_area->selectedPlottables().at(0)->name()){
                    on_actionRemove_dataset_triggered();
                    return;
                }
            for (int i=0;i<graphSets.size();i++)
                if(graphSets[i]!=nullptr){
                    if(graphSets[i]->name()==ui->plot_area->selectedPlottables().at(0)->name()){
                        on_actionRemove_a_fitted_curve_triggered();
                        return;
}}}}}
