#ifndef PLIS_H
#define PLIS_H

#include <QMainWindow>
#include <QFileDialog>
#include "qcustomplot.h"
#include "data.h"
#include "fitmrq.h"
#include "simulatedialog.h"
#include "jackknife.h"
#include "modify_data.h"
#include "modify_text.h"
#include "show_result.h"
#include "guess_parameters.h"
#include "help.h"
#include "optimize.h"
#include <iostream>
#include <QFile>
#include <QTextStream>
#include <QKeyEvent>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace Ui {
class Plis;
}

class Plis : public QMainWindow
{
    Q_OBJECT

public:
    explicit Plis(QWidget *parent = nullptr);
    ~Plis();

    void keyReleaseEvent(QKeyEvent *event);

private slots:
    
    void on_actionImport_data_triggered();

    void on_actionExit_triggered();

    void closeEvent(QCloseEvent *event);

    void on_actionAdd_dataset_triggered();

    void on_actionSimulate_data_triggered();

    void on_actionRemove_dataset_triggered();

    void on_add_row_button_clicked();

    void on_remove_row_button_clicked();

    void on_actionOpen_project_triggered();

    void on_actionSave_image_triggered();

    void on_actionSave_project_triggered();

    void on_actionRemove_all_datasets_triggered();

    void on_actionPaste_data_to_table_triggered();

    void on_table_cellChanged(int row, int column);

    void on_selected_dataSet_currentIndexChanged(int index);

    void on_protein_conc_square_editingFinished();

    void on_unit_currentIndexChanged(const QString &arg1);

    void on_update_graph_button_clicked();

    void selection_changed(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*);

    void legend_dubbel_clicked(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*);

    void on_ligan_conc_unit_currentIndexChanged(const QString &arg1);

    void contextMenuRequest(QPoint pos);

    void moveLegend();

    void change_title(QMouseEvent *event);

    void on_actionOne_bindning_site_triggered();

    void on_actionTwo_binding_sites_triggered();

    void on_actionFour_binding_sites_triggered();

    void on_actionCompetitive_binding_triggered();

    void on_actionRemove_a_fitted_curve_triggered();

    void on_actionRemove_all_fitted_curves_triggered();

    void graph_double_clicked(QCPAbstractPlottable*, int);

    void graph_clicked(QCPAbstractPlottable*, int);

    void axis_double_clicked(QCPAxis *, QCPAxis::SelectablePart, QMouseEvent*);

    void on_pop_result_button_clicked();

    void on_actionOptimize_fit_triggered();

    void on_actionStandard_triggered();

    void on_actionTitration_NMR_triggered();

    void on_switch_to_CPGM_triggered();

    void CPMG();

    void on_actionCopy_data_from_table_triggered();

    void on_actionHelp_triggered();

    void on_actionInformation_about_PLIS_triggered();

private:
    Ui::Plis *ui;
    QDoubleValidator *v=new QDoubleValidator(0.1,99999,2,this);
    SimulateDialog* sim=new SimulateDialog{};
    SimulateDialog* sim_cpmg=new SimulateDialog{this,"cpmg"};
    QList<QAction*> s_action;
    Help *help=new Help{};

    QVector<Data*> dataSets{}; //Data-pointer to the different dataSets
    QVector<QCPGraph*> graphSets{};
    int current_dataSet{};
    QVector<QCPScatterStyle::ScatterShape> shapes{QCPScatterStyle::ssDisc,QCPScatterStyle::ssSquare, QCPScatterStyle::ssDiamond, QCPScatterStyle::ssStar,
                QCPScatterStyle::ssCross, QCPScatterStyle::ssPlus, QCPScatterStyle::ssCircle, QCPScatterStyle::ssTriangle, QCPScatterStyle::ssTriangleInverted,
                QCPScatterStyle::ssCrossSquare, QCPScatterStyle::ssPlusSquare, QCPScatterStyle::ssCrossCircle, QCPScatterStyle::ssPlusCircle, QCPScatterStyle::ssPeace,
                QCPScatterStyle::ssCustom};

    int shape{};
    bool is_disconnected{false};
    void change_penColor(int);
    QVector<double> kd_testVector{0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000};
                                //nM                uM           mM                     M
    QString xlabel{};
    QString ylabel{};
    QCPTextElement *title{};
    QFont xaxis_font{QFont("Helvetica",9)}, yaxis_font{QFont("Helvetica",9)}, title_font{QFont("Helvetica",9)}, legend_font{QFont("Helvetica",9)};
    QString mode{"standard"}; //"standard": normal titration, "chemshift": tiitration where concentration does not affect relevant signal, "CPMG": relaxation dispersion

    //private funtions
    void render_graphs();
    void disconnect_when_no_dataset(bool);
    void write_In_Table();
    void draw_curve();
    void write_In_Result();
    void delete_dataSets(int i=(-1));
    void check_guess(double &y_calc, double &chi2, double &chi2_temp, int &best_index_kd1, int &best_index_kd2);
    QString get_name(QString);
    int make_fitmrq(Fitmrq &curve_fit);
    int make_fitmrq2(Fitmrq2 &curve_fit);
    void set_scatterstyle();
    void change_mode(QString new_mode);
    void setAxisLabels(QString xlab, QString ylab);
    void handleFitDidNotWork(int did_fit_work, Data *dataSet);
    void checkMode(bool standard, bool chemshift, bool cpmg);
};

#endif // PLIS_H
