#ifndef SAVE_OPEN_H
#define SAVE_OPEN_H

#include "plis.h"

void save_project(QMainWindow *this1, const QString &x, const QString &y,
                  const QVector<Data*> &dataSets, const int &current_dataSet, const bool &is_disconnected, int &shape,
                  QCPTextElement *&title, QFont &xaxis_font, QFont &yaxis_font, QString &mode);


int open_project(QMainWindow *this1, QString &x, QString &y, QVector<Data*> &dataSets,
                  int &current_dataSet, bool &is_disconnected, QVector<QCPGraph *> &graphSets, int &shape,
                  QCPTextElement *&title, QFont &xaxis_font, QFont &yaxis_font, QString &mode);


void line2words(const std::string& line, std::vector<std::string> &word);

void toDouble(const std::string& line, double &temp);

void toDouble_cpmg(std::string line, Data* &data);

void toDouble(std::string line, Data* &dataSets, bool curve=false);

int toInt(const std::string& word);

#endif // SAVE_OPEN_H
