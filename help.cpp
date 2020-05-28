#include "help.h"
#include "ui_help.h"

Help::Help(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Help)
{
    ui->setupUi(this);
    ui->listWidget->setCurrentRow(0);
    ui->text->setReadOnly(true);
}

Help::~Help()
{
    delete ui;
}

void Help::setup(int i){
    if(i==1)
        ui->listWidget->setCurrentRow(ui->listWidget->count()-1);
}

void Help::on_listWidget_currentTextChanged(const QString &currentText)
{
    if(currentText=="Main window")
        ui->text->setHtml("<p><strong><u>Main window</u></strong></p>"
                          "<hr>"
                          "<p><em><u>Table</u><u>:</u>&nbsp;</em>To the left of the main window the data could be edited with the table view. If more/less rows is needed use add/remove row. And remember to update graph when finished editing data.&nbsp;</p>"
                          "<p>Make sure to fill in the protein concentration (&gt;0).<br><u>Note:</u> Make sure you have as many rows as datapoints.</p>"
                          "<p><em><u>Manage the Plot area:</u>&nbsp;</em>"
                          "<br>- Axis label: Change the font and a label text by dubbel clicking the axis of choice."
                          "<br>- Zoom: Use the scroll function on the computer mouse or touchpad to zoom in/out in the plot area. One can also choose to zoom x-axis or y-axis simply by scrolling when hovering the mouse over the axis of choice. <"
                          "br>- Legend: Right click on the legend to move it around on the plot area.<br>- "
                          "Rename dataset: Dubbel click on the dataset that you wish to rename (click in the legend)."
                          "<br>- Design the plotted curve: Do this by dubbel clicking the plotted curve on the plot area.</p>");
    else if(currentText=="File")
        ui->text->setHtml("<p><u><strong>File</strong></u>:</p>"
                          "<hr>"
                          "<p><u>Import data</u>: <br>Simply choose the file you want to import to the table view, make sure that the file is in txt format and your datapoints is separated by tab or space.</p>"
                          "<p><u>Open project</u>:<br>Choose the correct PLIS project file (.proj) that you want to open.</p>"
                          "<p><u>Save project</u>:<br>Choose the location where you want to save your current project.</p>"
                          "<p><u>Export image:</u><br>Choose the correct folder to save your plot area image to. Choose between different formats to make sure you get the quality of image that you want.</p>");
    else if(currentText=="Edit")
        ui->text->setHtml("<p><strong><u>Edit</u>:</strong></p>"
                          "<hr>"
                          "<p><u>Add dataset</u>:<br>Adds a new dataset to the software, you will be able to rename it.&nbsp;</p>"
                          "<p><u>Simulate data</u>:<br>Simulate data by: <br>1. Choosing the model which the data should follow.<br>2. Fill in the different parameters correlated to that particular model.<br>3. Edit the parameters to the right which is the same for all models. (Ligand range, noise and so on)<br>4. Click simulate</p>"
                          "<p><u>Copy/paste data from/to table</u>:<br>Simply copy/paste data as a normal table.</p>"
                          "<p><u>Remove</u>:<br>- Remove dataset: Removes the current selected dataset. (Which currently can be seen in the table)<br>- Remove all dataset: Removes all the datasets and fitted curves.<br>- Remove selected fitted curve: Removes the currently selected fitted curve. <br>- Remove all fitted curves: Removes all the fitted curves on the plot area.</p>");
    else if(currentText=="Fit")
        ui->text->setHtml("<p><strong>Fit</strong><strong>:</strong></p>"
                          "<hr>"
                          "<p>1. Choose the model which you believe the data follows.<br>2. Either choose to guess the starting parameters or let the program guess them for you.</p>"
                          "<p><u>If</u> you chose to guess the parameters continue here, else you are done.<br>3. Choose if you want du estimate the exact values or choose to guess within a range. </p>"
                          "<p><u>If</u> you chose range: <br>4. Choose the range and the number of datapoints within that range that should be analyzed. <br>5. Click Ok, when done.&nbsp;</p>"
                          "<p><u>Note</u>: If a value should be fixed to a specific value, simply check the checkbox corresponding to that parameter.</p>");
    else if(currentText=="Optimize fit")
        ui->text->setHtml("<p><strong><u>Optimize fit</u>:</strong></p>"
                          "<hr>"
                          "<p>1. Choose the fit which you want to optimize. <br>"
                          "2. Click on optimize fit. <br>"
                          "3. Modify the parameters<br>"
                          "4. Click on optimize and the software calculates the best fit with those parameters as starting guesses. <br>"
                          "5. Redo the optimization until you found the best fit. <br>"
                          "6. Click apply.</p>");
    else if(currentText=="Mode")
        ui->text->setHtml("<p><u><strong>Switch mode</strong></u><strong>:</strong></p>"
                          "<hr>"
                          "<p>Choose the mode which you want to work with, either titration or CPMG.&nbsp;</p>"
                          "<p>The current dataset will all be deleted when switching mode.</p>");
    else if(currentText=="About PLIS")
        ui->text->setHtml("<p><strong>Information about PLIS - Protein Ligand Interaction Software</strong><strong>:</strong></p>"
                          "<hr>"
                          "<p>PLIS was created as a Master Thesis project by Emil Bj&ouml;rklund in 2020 to help the researchers and students of Link&ouml;pings University with the analysis of mainly fluorescence spectroscopy data but also CPMG measurements from NMR.</p>"
                          "<p>The credits for both the id&eacute; and most of the features belong to the examiner of the project, Patrik Lundstr&ouml;m.&nbsp;</p>"
                          "<hr>"
                          "<p><strong>The software:</strong><br>PLIS is written in C++ ISO2009. With a graphical library by Qt and Qt creator will be used to build the program and write the code. The numerical routines was inspired by Numerical Recipes 3rd Edition (Press et al.).</p>"
                          "<p>As an error estimator of the K<sub>d</sub> - values Jackknife will be used.&nbsp;</p>");
}

void Help::on_pushButton_2_clicked()
{
    hide();
}
