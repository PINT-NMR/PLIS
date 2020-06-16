#include "plis.h"
#include "ui_plis.h"
#include "save_open.h"
#include "legend.h"
#include "calculations.h"
#include "table.h"
#include "result.h"
#include "plot.h"

// 200530
// - the ia vector from fitmrq is now passed to jackknife
// 200601
// - created a new class Table and moved some contents to that
//   class
// 200607
// - created a new class Result and moved some contents there.
// 200608
// - created a new mode "chemshift". It is intended for titrations where
//   concentration does not affect relevant signal, e.g. chemical shifts

Plis::Plis(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Plis)
{
    ui->setupUi(this);
    ui->ligan_conc_unit->hide();
    ui->label_3->hide();
    ui->unit->removeItem(1);
    ui->unit->removeItem(1);

    Data* temp=new Data{};
    for (int i=0; i<9;i++)
        temp->pushBackData(0., 0., 0., 0., 1.e-4);
    temp->dataAndCurveVisible(true, false);
    temp->name="Dataset 1";
    dataSets.push_back(temp);
    graphSets.push_back(nullptr);

    Table table(ui->table);
    table.write_In_Table(dataSets, current_dataSet, mode);
    current_dataSet=0;
    ui->selected_dataSet->addItem("Dataset 1");

    set_scatterstyle();

    //Setup plot-area
    s_action=ui->menuFit->actions();
    ui->frame_2->setStyleSheet("background-color: rgb(0, 0, 0);");

    ui->plot_area->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectAxes |
                                    QCP::iSelectLegend | QCP::iSelectPlottables | QCP::iSelectItems);
    setAxisLabels("Concentration of Ligand (µM)", "Response (a.u.)");
    title=new QCPTextElement(ui->plot_area, "Title (double click to change)", QFont("Arial",14 ,QFont::Bold));
    ui->protein_conc_square->setValidator(v);
    ui->table->setCurrentCell(0,0);
    initPlotArea(ui->plot_area, xlabel, ylabel, title);

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
    QMessageBox::StandardButton pop_up{QMessageBox::No};

    dialog.setDirectory(QDir::current());
    dialog.setFileMode(QFileDialog::ExistingFiles);
    dialog.setNameFilter("(*.txt)");
    if (dialog.exec()) {
        if (dataSets.size()!=0) {
            if ((pop_up = QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                                     "Do you want to overrite the current dataset?\nThis step cannot be undone!",
                                     QMessageBox::Yes | QMessageBox::No)) == QMessageBox::Yes) {
                dataSets.remove(current_dataSet);
                graphSets.remove(current_dataSet);
            }
        }
        for (int i=0; i<dialog.selectedFiles().size();i++){ //For-loop to add all datasets to a vector
            QString path{dialog.selectedFiles().at(i)};
            std::string f = path.toStdString();
            Data* temp=new Data;
            temp->read_Data(f);
            bool ok{};
            QString input{""};
            while (input=="") {
                input = QInputDialog::getText(this, "Please, name the dataset:", path, QLineEdit::Normal, QString(), &ok);
                if (!ok)
                    break;
            }
            temp->dataAndCurveVisible(true, false);
            temp->name=input;
            if (!ok)
                delete temp;
            else {
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

    if (namelist.size()!=0)
        ui->selected_dataSet->addItems(namelist);

    if (pop_up==QMessageBox::Yes) {
        if (ui->selected_dataSet->count()==1)
            ui->selected_dataSet->setItemText(0,"<None>");
        else
            ui->selected_dataSet->removeItem(current_dataSet);
    }

    if (current_dataSet!=-1) {
        ui->selected_dataSet->setCurrentIndex(current_dataSet);
        Table table(ui->table);
        table.write_In_Table(dataSets, current_dataSet, mode);
    }
    else {
        ui->table->setRowCount(0);
        disconnect_when_no_dataset(true);
    }

    if (is_disconnected==true && current_dataSet!=-1)
        disconnect_when_no_dataset(false);
    render_graphs();
}

void Plis::on_actionExit_triggered()
{
    if (QMessageBox::warning(this, "PLIS",
                            "Do you want to quit?\nAny unsaved changes will be lost.",
                            QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok)
        QApplication::quit();
}

void Plis::closeEvent(QCloseEvent *event)
{
    if (QMessageBox::warning(this, "PLIS",
                            "Do you want to quit?\nAny unsaved changes will be lost.",
                            QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok)
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
    if (!ok)
        delete temp;
    else {
        if (input!="")
            input=get_name(input);

        if (mode=="cpmg") {
            temp->clearVectors();
            temp->pushBackCPMGdata(0., 0., 0.);
        }
        dataSets.push_back(temp);
        graphSets.push_back(nullptr);
        temp->dataAndCurveVisible(true, false);

        if (input=="") {
            QString new_name = get_name("New data");
            ui->selected_dataSet->addItem(new_name);
            dataSets[dataSets.size()-1]->name=new_name;
        }
        else {
            ui->selected_dataSet->addItem(input);
            dataSets[dataSets.size()-1]->name=input;
        }
        if (dataSets.size() == 1)
            ui->selected_dataSet->removeItem(0);

        Table table(ui->table);
        table.write_In_Table(dataSets, current_dataSet, mode);
        current_dataSet=dataSets.size()-1;
        set_scatterstyle();
        ui->selected_dataSet->setCurrentIndex(current_dataSet);
        ui->table->setCurrentCell(0,0);
    }
}

void Plis::on_actionSimulate_data_triggered()
{
    int done{-1};
    if(mode=="standard" || mode=="chemshift") {
        sim->setup(dataSets,mode);
        done=sim->exec();
    }
    else if(mode=="cpmg") {
        sim_cpmg->setup(dataSets,mode);
        done=sim_cpmg->exec();
    }
    if(done==1) {
        if (mode=="standard" || mode=="chemshift") {
            dataSets.push_back(sim->dataSets_local.last());
            current_dataSet=sim->current_dataSet;
        }
        else {
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
        if (is_disconnected==true) {
            disconnect_when_no_dataset(false);
            is_disconnected=false;
        }
        Table table(ui->table);
        table.write_In_Table(dataSets, current_dataSet, mode);
        render_graphs();
    }
}

void Plis::on_actionRemove_dataset_triggered()
{
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software", "Warning: Are you sure that you want to remove the selected dataset?",
                             QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        delete_dataSets(current_dataSet);
    render_graphs();
}

void Plis::on_actionOpen_project_triggered(){
    int i{0};
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "The current project will be closed. Any unsaved data will be lost",
                             QMessageBox::Ok | QMessageBox::Cancel) == QMessageBox::Ok) {
        QString temp_mode{};
        i=open_project(this, xlabel, ylabel, dataSets, current_dataSet, is_disconnected, graphSets, shape, title, xaxis_font, yaxis_font,temp_mode);
        if (i==0) {
            if(temp_mode=="cpmg" && (mode=="standard" || mode=="chemshift")){
                change_mode("cpmg");
                ui->protein_conc_square->setText(QString::number(dataSets[current_dataSet]->ligand_cpmg));
            }
            else if(temp_mode=="standard" && (mode=="cpmg" || mode=="chemshift"))
                change_mode("standard");
            else if(temp_mode=="chemshift" && (mode=="cpmg" || mode=="standard"))
                change_mode("chemshift");
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
    if (i==0) {
        render_graphs();
        Table table(ui->table);
        table.write_In_Table(dataSets, current_dataSet, mode);
        write_In_Result();
        if (dataSets[current_dataSet]->has_curve)
            ui->menuOptimize_fit->setDisabled(false);
    }
}

void Plis::on_actionSave_image_triggered()
{
    saveImage(this, ui->plot_area, dataSets, current_dataSet);
}

void Plis::on_actionSave_project_triggered()
{
    save_project(this, xlabel, ylabel, dataSets, current_dataSet, is_disconnected, shape, title, xaxis_font, yaxis_font,mode);
}

void Plis::on_actionPaste_data_to_table_triggered()
{
    Table table(ui->table);
    table.paste_data_to_table(dataSets, current_dataSet);
}

void Plis::on_actionRemove_all_datasets_triggered()
{
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software", "Warning: Are you sure that you want to remove all datasets?",
                             QMessageBox::Yes | QMessageBox::No) == QMessageBox::No)
        return;
    clearAndReplot(ui->plot_area);
    int temp{dataSets.size()};
    for (int i=0;i<temp ;i++)
        delete_dataSets(0);
}

void Plis::on_actionOne_bindning_site_triggered()
{
    if (!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true) {
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
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "Don't forget to fill in a protein concentration.", QMessageBox::Ok);
        return;
    }
    dataSets[current_dataSet]->setModel_a_dydaOneSite();

    int i{1};
    Guess_parameters guess{this,dataSets[current_dataSet], mode};
    while (1) {
        if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                                 "Do you want to estimate the parameters?",
                                 QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes){
            i=guess.exec();
            break;
        }
        else {
            guess_parameters_one_site(dataSets[current_dataSet], kd_testVector, mode);
            break;
        }
    }
    if(i==1) {
        void (*voidptr)(const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_one_site_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_one_site_dilution};
        if (mode=="chemshift") {
            voidptr = fit_one_site;
            dblptr = fit_one_site;
        }
        Fitmrq curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, voidptr, dblptr};
        curve_fit.ia = {!guess.a[0], !guess.a[1], !guess.a[2], false, false};
        if (!make_fitmrq(curve_fit)) {
            Jackknife jack{dataSets[current_dataSet], curve_fit.ia, mode};
            jack.compute();
        }
        write_In_Result();
    }
    else
        dataSets[current_dataSet]->setModelIfNoFit();
}

void Plis::on_actionTwo_binding_sites_triggered()
{
    if (!ui->menuFit->isEnabled())
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
    dataSets[current_dataSet]->setModel_a_dydaTwoSite();

    Guess_parameters guess{this, dataSets[current_dataSet], mode};
    int i{1};
    while (1) {
        if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                                 "Do you want to estimate the parameters?",
                                 QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
            i=guess.exec();
            break;
        }
        else {
            QApplication::setOverrideCursor(Qt::WaitCursor);
            QProgressDialog progress("Optimizing...", "Abort Optimizing", 0, 12, this);
            progress.setCancelButton(nullptr);
            progress.setMinimumDuration(0);
            progress.setWindowModality(Qt::WindowModal);
            guess_parameters_two_sites(dataSets[current_dataSet], kd_testVector,progress, mode);
            QApplication::setOverrideCursor(Qt::ArrowCursor);
            break;
        }
    }
    if (i==1) {
        void (*voidptr)(const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_two_sites_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_two_sites_dilution};
        if (mode=="chemshift") {
            voidptr = fit_two_sites;
            dblptr = fit_two_sites;
        }
        Fitmrq curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, voidptr, dblptr};
        curve_fit.ia = {!guess.a[0], !guess.a[1], !guess.a[2], !guess.a[3], !guess.a[4], false, false};
        if (!make_fitmrq(curve_fit)) {
            Jackknife jack{dataSets[current_dataSet], curve_fit.ia, mode};
            jack.compute();
        }
        write_In_Result();
    }
    else
        dataSets[current_dataSet]->setModelIfNoFit();
}

void Plis::on_actionFour_binding_sites_triggered()
{
    if (!ui->menuFit->isEnabled())
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
    if (dataSets[current_dataSet]->protein_conc<=0.0) {
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                             "Don't forget to fill in a protein concentration.",QMessageBox::Ok);
        return;
    }
    dataSets[current_dataSet]->setModel_a_dydaFourSite();

    Guess_parameters guess{this, dataSets[current_dataSet], mode};
    int i = 1;
    while (1) {
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "For four-site binding you must provide estimates "
                             "of the parameters.", QMessageBox::Ok);
        i = guess.exec();
        break;
    }
    if (i == 1) {
        void (*voidptr)(const double, const double, const std::vector<double> &, double &, std::vector<double> &){fit_four_sites_dilution};
        double (*dblptr)(const double, const double, const std::vector<double> &){fit_four_sites_dilution};
        if (mode=="chemshift") {
            voidptr = fit_four_sites;
            dblptr = fit_four_sites;
        }

        Fitmrq curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, voidptr, dblptr};
        curve_fit.ia = {1, 1, 1, 1, 1, 1, 1, 1, 1, false, false};
        if (!make_fitmrq(curve_fit)) {
            Jackknife jack{dataSets[current_dataSet], curve_fit.ia, mode};
            jack.compute();
        }
        write_In_Result();
    }
    else
        dataSets[current_dataSet]->setModelIfNoFit();
}


void Plis::on_actionCompetitive_binding_triggered()
{
    if(!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true) {
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
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "Don't forget to fill in a protein concentration.", QMessageBox::Ok);
        return;
    }
    double comp_conc{};
    bool ok{};
    comp_conc=QInputDialog::getDouble(this, "PLIS - Protein Ligand Interaction Software",
                                      "Please enter the concentration of the competitor (µM):",10,0.0001, 1000, 4,&ok);
    if (!ok) return;

    dataSets[current_dataSet]->update_comp_vector(comp_conc);
    dataSets[current_dataSet]->setModel_a_dydaCompTwoSite(comp_conc);

    Guess_parameters guess{this, dataSets[current_dataSet], mode};
    int i{1};
    while (1) {
        if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                                 "Do you want to estimate the parameters?",
                                 QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes){
            i=guess.exec();
            break;
        }
        else {
            QApplication::setOverrideCursor(Qt::WaitCursor);
            QProgressDialog progress("Optimizing...", "Abort Optimizing", 0, 12, this);
            progress.setCancelButton(nullptr);
            progress.setMinimumDuration(0);
            progress.setWindowModality(Qt::WindowModal);
            guess_parameters_comp(dataSets[current_dataSet], kd_testVector,progress, mode);
            QApplication::setOverrideCursor(Qt::ArrowCursor);
            break;
        }
    }
    if(i==1) {
        void (*voidptr)(const double, const double, const double, const std::vector<double> &,
                        double &, std::vector<double> &){fit_comp_dilution};
        double (*dblptr)(const double, const double, const double, const std::vector<double> &){fit_comp_dilution};
        if (mode=="chemshift") {
            voidptr = fit_comp;
            dblptr = fit_comp;
        }
        Fitmrq2 curve_fit{dataSets[current_dataSet]->concVector.toStdVector(),
                    dataSets[current_dataSet]->protein_conc_vector.toStdVector(),
                    dataSets[current_dataSet]->comp_vector.toStdVector(),
                    dataSets[current_dataSet]->responceVector.toStdVector(),
                    dataSets[current_dataSet]->volumeVector.toStdVector(),
                    dataSets[current_dataSet]->error_vector.toStdVector(),
                    dataSets[current_dataSet]->a, voidptr, dblptr};
        curve_fit.ia = {!guess.a[0], !guess.a[1],!guess.a[2],!guess.a[3],!guess.a[4],false,false,false};
        if (!make_fitmrq2(curve_fit)) {
            Jackknife jack{dataSets[current_dataSet], curve_fit.ia, mode};
            jack.compute();
        }
        write_In_Result();
    }
    else
        dataSets[current_dataSet]->setModelIfNoFit();
}

void Plis::handleFitDidNotWork(int did_fit_work, Data *dataSet) {
    QVector<QString> msg{"Singular matrix.", "Too many iterations"};
    QMessageBox::warning(this, "Warning: Something went wrong!", msg[did_fit_work-1]);
    dataSet->setModelIfNoFit();
}

void Plis::CPMG() {
    if(!ui->menuFit->isEnabled())
        return;
    if (dataSets[current_dataSet]->has_curve==true) {
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                             "You already have a fit made to this dataset",
                             QMessageBox::Ok);
        return;
    }
    if (dataSets[current_dataSet]->ligand_cpmg<=0.0) {
        QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Studies",
                             "Don't forget to enter ligand concentration.", QMessageBox::Ok);
        return;
    }
    dataSets[current_dataSet]->setModel_a_dydaCPMG();

    Guess_parameters guess{this, dataSets[current_dataSet], mode};
    int i{1};
    while (1) {
        if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                                 "Do you want to estimate the parameters?",
                                 QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
            i=guess.exec();
            break;
        }
        else {
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
    if(i==0) {
        dataSets[current_dataSet]->setModelIfNoFit();
        return;
    }

    Fitmrq3 fit{dataSets[current_dataSet]->n_cpmgVector.toStdVector(),dataSets[current_dataSet]->R2effVector.toStdVector(),
                dataSets[current_dataSet]->dyVector.toStdVector(), dataSets[current_dataSet]->a,
                carverrichards,carverrichards};
    fit.ia = { !guess.a[0], !guess.a[1], !guess.a[2], !guess.a[3], false};
    int did_fit_work=fit.fit();

    if (did_fit_work==1 || did_fit_work==2)
        handleFitDidNotWork(did_fit_work, dataSets[current_dataSet]);
    else {
        dataSets[current_dataSet]->addCPMGresults(fit.a, fit.chisq, fit.plot());
        Jackknife jack{dataSets[current_dataSet], fit.ia, mode};
        jack.compute_cpmg();
        ui->plot_area->addGraph();
        graphSets[current_dataSet] = ui->plot_area->graph();
        ui->plot_area->addGraph();
        graphSets[current_dataSet] = ui->plot_area->graph();
        render_graphs();
        write_In_Result();
        ui->menuOptimize_fit->setDisabled(false);
    }
}

void Plis::on_actionRemove_a_fitted_curve_triggered()
{
    if (dataSets[current_dataSet]->has_curve){
        if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software", "Warning: Are you sure that you want to remove the fitted curve?",
                                 QMessageBox::Yes | QMessageBox::No) == QMessageBox::No)
            return;
        QCPGraph * temp=graphSets[current_dataSet];
        graphSets[current_dataSet]=nullptr;
        ui->plot_area->removePlottable(temp);
        ui->menuOptimize_fit->setDisabled(true);
        dataSets[current_dataSet]->removeFittedCurve();
        write_In_Result();
        render_graphs();
    }
}

void Plis::on_actionRemove_all_fitted_curves_triggered()
{
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software", "Warning: Are you sure that you want to remove all fitted curve?",
                             QMessageBox::Yes | QMessageBox::No) == QMessageBox::No)
        return;
    for (int i=0;i<graphSets.size();i++){
        if(graphSets[i]!=nullptr){
            ui->plot_area->removePlottable(graphSets[i]);
            graphSets[i]=nullptr;
            dataSets[i]->removeFittedCurve();
        }
    }
    ui->menuOptimize_fit->setDisabled(true);
    write_In_Result();
    render_graphs();
}

void Plis::on_add_row_button_clicked()
{
    Table table(ui->table);
    table.addRow(dataSets, current_dataSet, mode);
}

void Plis::on_remove_row_button_clicked()
{
    Table table(ui->table);
    table.removeRow(dataSets, current_dataSet, mode);
}

void Plis::on_table_cellChanged(int row, int column)
{
    Table table(ui->table);
    table.cellChanged(dataSets, current_dataSet, mode, row, column);
}

void Plis::on_selected_dataSet_currentIndexChanged(int index)
{
    current_dataSet=index;
    Table table(ui->table);
    table.write_In_Table(dataSets, current_dataSet, mode);
    if (mode=="standard" || mode=="chemshift") {
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
    if(mode=="standard" || mode=="chemshift")
        ui->protein_conc_square->setText(QString::number(dataSets[current_dataSet]->protein_conc));
    else if(mode=="cpmg")
        ui->protein_conc_square->setText(QString::number(dataSets[current_dataSet]->ligand_cpmg));

    ui->menuOptimize_fit->setDisabled(!dataSets[current_dataSet]->has_curve);

    write_In_Result();
}

void Plis::on_protein_conc_square_editingFinished()
{
    if(current_dataSet!=(-1)){
        QString s=ui->protein_conc_square->text();
        s.replace(',','.');
        if(mode=="standard" || mode=="chemshift") {
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
        Table table(ui->table);
        table.write_In_Table(dataSets, current_dataSet, mode);
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
    if (ai!=nullptr) {
        new_name=QInputDialog::getText(this, "PLIS",
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
                QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
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

    if (ui->plot_area->legend->selectTest(pos, false) >= 0) { // context menu on legend requested
        menu->addAction("Move to top left", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignTop|Qt::AlignLeft));
        menu->addAction("Move to top center", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignTop|Qt::AlignHCenter));
        menu->addAction("Move to top right", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignTop|Qt::AlignRight));
        menu->addAction("Move to bottom right", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignBottom|Qt::AlignRight));
        menu->addAction("Move to bottom left", this, SLOT(moveLegend()))->setData(static_cast<int>(Qt::AlignBottom|Qt::AlignLeft));
    } else { // general context menu on graphs requested
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
  if (QAction* contextAction = qobject_cast<QAction*>(sender())) {
    bool ok;
    int dataInt = contextAction->data().toInt(&ok);
    if (ok) {
      ui->plot_area->axisRect()->insetLayout()->setInsetAlignment(0,static_cast<Qt::Alignment>(dataInt));
      ui->plot_area->replot();
    }
  }
}

void Plis::change_title(QMouseEvent *){
    if (QCPTextElement *t = qobject_cast<QCPTextElement*>(sender())) {
        Modify_text new_title{this, title};
        new_title.exec();
        ui->plot_area->replot();
    }
}

void Plis::render_graphs(){
    ui->plot_area->clearGraphs();
    int temp_current_dataSet{current_dataSet};
    for (int i=0; i<dataSets.size();i++) {
        current_dataSet=i;
        if (dataSets[i]->ligand_unit!=dataSets[0]->ligand_unit && (mode=="standard" || mode=="chemshift")){
            QVector<double> temp{dataSets[i]->concVector};
            dataSets[i]->change_unit_of_data(dataSets[0]->ligand_unit, "conc");
            setGraphStyle(ui->plot_area, dataSets[current_dataSet], mode);
            dataSets[i]->concVector=temp;
        }
        else
            setGraphStyle(ui->plot_area, dataSets[current_dataSet], mode);
        if (graphSets[i]!=nullptr)
            draw_curve();
    }

    int curve{0};
    for (int i=0;i<dataSets.size();i++) {
        if (dataSets[i]->data_visible==false) {
            ui->plot_area->legend->item(i+curve)->setTextColor("Grey");
            ui->plot_area->graph(i+curve)->setVisible(false);
        }
        if(dataSets[i]->has_curve)
            curve++;
    }
    curve=0;
    for (int i=0; i<dataSets.size();i++) {
        if (dataSets[i]->has_curve)
            curve++;
        if (dataSets[i]->curve_visible==false && graphSets[i]!=nullptr) {
            ui->plot_area->legend->item(i+curve)->setTextColor("Grey");
            ui->plot_area->graph(i+curve)->setVisible(false);
        }
    }
    current_dataSet=temp_current_dataSet;
    updatePlotArea(ui->plot_area, xlabel, ylabel, xaxis_font, yaxis_font);
}

void Plis::change_penColor(int i)
{
    dataSets[current_dataSet]->setPenColor(i);
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

void Plis::draw_curve(){
    ui->plot_area->addGraph();
    graphSets[current_dataSet]=ui->plot_area->graph();
    graphSets[current_dataSet]->setName("Curve - " + dataSets[current_dataSet]->name);
    graphSets[current_dataSet]->setPen(dataSets[current_dataSet]->pen_curve);
    graphSets[current_dataSet]->setLineStyle(QCPGraph::lsLine);
    graphSets[current_dataSet]->setData(dataSets[current_dataSet]->calc_data[0],dataSets[current_dataSet]->calc_data[1]);
}

void Plis::write_In_Result() {
    Result result(ui->kd1_box, ui->kd1_error, ui->Kd_label, ui->plus_minus,
                  ui->kd2_box,  ui->kd2_error, ui->kd2_label, ui->plus_minus2,
                  ui->kd3_box, ui->kd3_error, ui->Kd3_label, ui->plus_minus3,
                  ui->kd4_box, ui->kd4_error, ui->Kd4_label, ui->plus_minus4,
                  ui->name, ui->model, ui->num_data, ui->chi2_box);
    if (current_dataSet!=-1) {
        result.writeNameModelNdataChi2(mode, dataSets[current_dataSet]);
        if (dataSets[current_dataSet]->num_bind_site==1)
            result.writeKd1(dataSets[current_dataSet]);
        else if (dataSets[current_dataSet]->num_bind_site==2)
            result.writeKd2(dataSets[current_dataSet]);
        else if (dataSets[current_dataSet]->num_bind_site==3)
            result.writeKd3(dataSets[current_dataSet]);
        else if (dataSets[current_dataSet]->num_bind_site==5)
            result.writeKd4(dataSets[current_dataSet]);
        else if (dataSets[current_dataSet]->num_bind_site==4)
            result.writeKd1_CPMG(dataSets[current_dataSet]);
        else if (dataSets[current_dataSet]->num_bind_site==0)
            result.clearKd1();
    }
    else
        result.clearAll();
}

void Plis::delete_dataSets(int i) {
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
        Table table(ui->table);
        table.write_In_Table(dataSets, current_dataSet, mode);
    }
    write_In_Result();
}

QString Plis::get_name(QString name){
    if (ui->selected_dataSet->findText(name) == -1)
        return name;
    for (int i=1;; ++i) // correct to start at i=1
        if (ui->selected_dataSet->findText(name+" (" + QString::number(i)+")") == -1)
            return (name+" ("+QString::number(i)+")");
}


int Plis::make_fitmrq(Fitmrq &curve_fit){
    int did_fit_work = curve_fit.fit();
    if (did_fit_work==0) {
        dataSets[current_dataSet]->setResult(curve_fit.a, curve_fit.chisq, curve_fit.plot()[0], curve_fit.plot()[1]);
        ui->plot_area->addGraph();
        graphSets[current_dataSet]=ui->plot_area->graph();
        render_graphs();
        ui->menuOptimize_fit->setDisabled(false);
        return 0;
    }
    handleFitDidNotWork(did_fit_work, dataSets[current_dataSet]);
    return 1;
}

int Plis::make_fitmrq2(Fitmrq2 &curve_fit) {
    int did_fit_work = curve_fit.fit();
    if (did_fit_work==0) {
        dataSets[current_dataSet]->setResult(curve_fit.a, curve_fit.chisq, curve_fit.plot()[0], curve_fit.plot()[1]);
        ui->plot_area->addGraph();
        graphSets[current_dataSet]=ui->plot_area->graph();
        render_graphs();
        ui->menuOptimize_fit->setDisabled(false);
        return 0;
    }
    handleFitDidNotWork(did_fit_work, dataSets[current_dataSet]);
    return 1;
}

void Plis::set_scatterstyle() {
    if (shapes.at(shape) != QCPScatterStyle::ssCustom){
        dataSets[current_dataSet]->style=shapes[shape];
        dataSets[current_dataSet]->style.setSize(10);
        change_penColor(shape);
    }
    else {
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

void Plis::graph_clicked(QCPAbstractPlottable* a, int) {
    for (int i=0; i<graphSets.size();i++)
        if (graphSets[i]!=nullptr)
            if (a->name()==graphSets[i]->name()) {
                ui->selected_dataSet->setCurrentIndex(i);
                current_dataSet=i;
                break;
            }
    for (int j=0; j<dataSets.size();j++)
        if (a->name()==dataSets[j]->name) {
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
    Optimize opt{this, dataSets[current_dataSet], mode};
    opt.exec();
    write_In_Result();
    render_graphs();
}

void Plis::on_actionStandard_triggered()
{
    if (mode=="standard") return;
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                            "Do you want to change state?\nYour current project will be discarded!",
                            QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
        mode="standard";
        checkMode(true, false, false);
        QStringList names{"Response", "[Ligand] (µM)","Volume"};
        ui->table->setHorizontalHeaderLabels(names);
        setAxisLabels("Concentration of Ligand", "Response (a.u.)");
        ui->menuFit->clear();
        ui->menuFit->addActions(s_action);
        ui->menuFit->setTitle("Fit");
        ui->label_2->setText("Protein conc:");
        clearAndReplot(ui->plot_area);
        int temp{dataSets.size()};
        for (int i=0;i<temp ;i++)
            delete_dataSets(0);
        render_graphs();
    }
    else
        ui->actionStandard->setChecked(false);

}

void Plis::on_actionTitration_NMR_triggered()
{
    if (mode=="chemshift") return;
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                            "Do you want to change state?\nYour current project will be discarded!",
                            QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
        mode="chemshift";
        checkMode(false, true, false);
        QStringList names{"Response", "[Ligand] (µM)","Volume"};
        ui->table->setHorizontalHeaderLabels(names);
        setAxisLabels("Concentration of Ligand", "Response (a.u.)");
        ui->menuFit->clear();
        ui->menuFit->addActions(s_action);
        ui->menuFit->setTitle("Fit");
        ui->label_2->setText("Protein conc:");
        clearAndReplot(ui->plot_area);
        int temp{dataSets.size()};
        for (int i=0;i<temp ;i++)
            delete_dataSets(0);
        render_graphs();
    }
    else
        ui->actionTitration_NMR->setChecked(false);
}

void Plis::on_switch_to_CPGM_triggered()
{
    if (mode=="cpmg") return;
    if (QMessageBox::warning(this, "PLIS - Protein Ligand Interaction Software",
                            "Do you want to change state?\nYour current project will be discarded!",
                            QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
        mode="cpmg";
        checkMode(false, false, true);
        QStringList names{"nu_CPMG (Hz)", "R2,eff (/s)", "Error in R2,eff"};
        ui->table->setHorizontalHeaderLabels(names);
        setAxisLabels("nu_CPMG (Hz)", "R2,eff (/s)");
        ui->menuFit->setTitle("Fit");
        ui->menuFit->clear();
        ui->menuFit->addAction("Fit",this,SLOT(CPMG()));
        ui->label_2->setText("Ligand conc:");
        ui->label_3->hide();
        ui->ligan_conc_unit->hide();
        clearAndReplot(ui->plot_area);
        int temp{dataSets.size()};
        for (int i=0;i<temp ;i++)
            delete_dataSets(0);
        render_graphs();
    }
    else
        ui->switch_to_CPGM->setChecked(false);
}

void Plis::setAxisLabels(QString xlab, QString ylab) {
        xlabel = xlab;
        ylabel = ylab;
}

void Plis::change_mode(QString new_mode) {
    if (new_mode=="standard") {
        mode="standard";
        checkMode(true, false, false);
        QStringList names{"Response", "[Ligand]","Volume"};
        ui->table->setHorizontalHeaderLabels(names);
        setAxisLabels("Concentration of Ligand", "Response");
        ui->menuFit->clear();
        ui->menuFit->addActions(s_action);
        ui->menuFit->setTitle("Fit");
        ui->label_2->setText("Protein\nConcentration");
    }
    else if (new_mode=="chemshift") {
        mode="chemshift";
        checkMode(false, true, false);
        QStringList names{"Response", "[Ligand]","Volume"};
        ui->table->setHorizontalHeaderLabels(names);
        setAxisLabels("Concentration of Ligand", "Response");
        ui->menuFit->clear();
        ui->menuFit->addActions(s_action);
        ui->menuFit->setTitle("Fit");
        ui->label_2->setText("Protein\nConcentration");
    }
    else if (new_mode=="cpmg") {
        mode="cpmg";
        checkMode(false, false, true);
        QStringList names{"nu_CPMG (Hz)", "R2,eff (/s)", "Error in R2,eff"};
        ui->table->setHorizontalHeaderLabels(names);
        setAxisLabels("nu_CPMG (Hz)", "R2,eff (/s)");
        ui->menuFit->setTitle("Fit");
        ui->menuFit->clear();
        ui->menuFit->addAction("Fit",this,SLOT(CPMG()));
        ui->label_2->setText("Ligand conc:");
    }
}

void Plis::checkMode(bool standard, bool chemshift, bool cpmg)
{
    ui->actionStandard->setChecked(standard);
    ui->actionTitration_NMR->setChecked(chemshift);
    ui->switch_to_CPGM->setChecked(cpmg);
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

void Plis::keyReleaseEvent(QKeyEvent *event) {
    if (event->key() != Qt::Key_Delete) return;
    if (ui->plot_area->selectedPlottables().size()>0) {
        for (int i=0;i<dataSets.size();i++)
            if (dataSets[i]->name==ui->plot_area->selectedPlottables().at(0)->name()){
                on_actionRemove_dataset_triggered();
                return;
            }
        for (int i=0;i<graphSets.size();i++)
            if (graphSets[i]!=nullptr) {
                if(graphSets[i]->name()==ui->plot_area->selectedPlottables().at(0)->name()) {
                    on_actionRemove_a_fitted_curve_triggered();
                    return;
                }
            }
    }
}
