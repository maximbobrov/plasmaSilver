#include "mainwindow.h"
#include "qcustomplot.h"
#include <QGridLayout>
#include <QPushButton>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <stdio.h>

#define NZ 100

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    m_widget=new QWidget();
    m_widget->setFixedHeight(700);
    m_widget->setFixedWidth(900);

    m_widget->setObjectName("central");
    m_customPlot = new QCustomPlot(this);
    m_grid = new QGridLayout(this);
    m_vLayoutCheckBoxes = new QVBoxLayout(this);
    m_hLayout = new QHBoxLayout(this);
    m_textStartTime = new QLineEdit(this);
    m_textDeltaTime = new QLineEdit(this);
    m_textEndTime =   new QLineEdit(this);
    m_scrollBar = new QScrollBar(Qt::Horizontal,this);
    m_simulateButton = new QPushButton("Simulate");
    m_debugButton = new QPushButton("plotDebug");

    m_grid->addWidget(m_customPlot,0,0,3,3);
    m_grid->addLayout(m_vLayoutCheckBoxes,0,3,1,1);
    m_grid->addWidget(m_simulateButton,5,0);
    m_grid->addWidget(m_debugButton,5,1);

    m_hLayout->addWidget(new QLabel("Start time:"));
    m_hLayout->addWidget(m_textStartTime);
    m_hLayout->addWidget(new QLabel("dt:"));
    m_hLayout->addWidget(m_textDeltaTime);
    m_hLayout->addWidget(new QLabel("End time:"));
    m_hLayout->addWidget(m_textEndTime);
    m_grid->addLayout(m_hLayout,4,0,1,2);
    m_grid->addWidget(m_scrollBar,4,2,1,1);


    m_data = new simulationData(NZ);
    m_data->setDt();
    m_data->setDz();
    m_sNe = new solverNe(m_data);
    m_sEn = new solverEnergy(m_data);
    m_sPhi = new solverPhi(m_data);

    m_textStartTime->setText(QString().number(0.0));
    m_textDeltaTime->setText(QString().number(m_data->getDt()));
    m_textEndTime->setText(QString().number(0.0001));


    connect (m_simulateButton, SIGNAL(clicked(bool)), this, SLOT(simulateData(bool)));
    connect (m_debugButton, SIGNAL(clicked(bool)), this, SLOT(drawDebug(bool)));
    connect(m_scrollBar, SIGNAL(valueChanged(int)), this, SLOT(replotGraph(int)));

    m_simulateButton->setCheckable(true);

    m_widget->setLayout(m_grid);
    setCentralWidget(m_widget);
    setWindowTitle("PlasmaSolver");
    m_animStopped=true;
}

MainWindow::~MainWindow()
{


}

void MainWindow::replotGraph(int number)
{
    m_maxY=-1e30;
    m_minY=+1e30;
    m_plots.clear();
    m_plots = m_storage[m_scrollBar->value()].plots;
    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel("x");
    m_customPlot->yAxis->setLabel("y");
    QColor colors[6] = {QColor(255,0,0),QColor(0,255,0),QColor(0,0,255),QColor(255,255,0),QColor(0,255,255),QColor(255,0,255)};
    for (int j = 0; j < m_plots.size(); ++j)
    {
        m_customPlot->addGraph();
        m_customPlot->graph(j)->setName(QString(m_plots[j].name));
        m_customPlot->graph(j)->addToLegend();
        m_customPlot->graph(j)->setData(m_plots[j].x,m_plots[j]. y);
        m_customPlot->graph(j)->setLineStyle((QCPGraph::LineStyle)(1));
        QPen graphPen;
        graphPen.setColor(colors[j]);
        graphPen.setWidthF(2);
        m_customPlot->graph(j)->setPen(graphPen);
        if(!m_checkBoxes[j]->isChecked())
        {
            m_customPlot->graph(j)->setVisible(false);
        }
        else {
            for (int i=0; i < m_plots[j].size; ++i)
            {
                if(m_plots[j].y[i] > m_maxY)
                    m_maxY = m_plots[j].y[i];
                if(m_plots[j].y[i] < m_minY)
                    m_minY = m_plots[j].y[i];
            }

        }
    }
    m_customPlot->xAxis->setRange(-1, 1);
    m_customPlot->yAxis->setRange(m_minY, m_maxY);
    m_customPlot->legend->setVisible(true);
    m_customPlot->replot();
}

void MainWindow::saveInStorage()
{
    for (int j = 0; j < m_plots.size(); ++j)
    {
        for (int i=0; i < m_plots[j].size; ++i)
        {
            m_plots[j].y[i] = m_plots[j].scale * m_plots[j].arr[i];
        }
    }
    storeStruct storeItem;
    storeItem.plots = m_plots;
    storeItem.time = m_time;
    m_storage.push_back(storeItem);

}

void MainWindow::addPlot(double *arr, char *name, int size, double scale)
{
    plotStruct plot;
    plot.x.resize(size);
    plot.y.resize(size);
    plot.arr = arr;
    for (int i = 0; i < size; ++i)
    {
        plot.x[i] = i / (size * 1.0 / 2)-1;
        plot.y[i] = plot.arr[i];
    }
    plot.size = size;
    plot.name = name;
    plot.scale = scale;
    m_plots.push_back(plot);

    QCheckBox* checkBox = new QCheckBox();
    checkBox->setText(QString(name));
    checkBox->setCheckState(Qt::Checked);
    m_checkBoxes.push_back(checkBox);
    m_vLayoutCheckBoxes->addWidget(checkBox);
    //connect(checkBox, SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
}

void MainWindow::addPlotXY(double *arr,double*xx, char *name, int size, double scale)
{
    plotStruct plot;
    plot.x.resize(size);
    plot.y.resize(size);
    plot.arr = arr;
    for (int i = 0; i < size; ++i)
    {
        plot.x[i] = xx[i];//i / (size * 1.0 / 2)-1;
        plot.y[i] = plot.arr[i];
    }
    plot.size = size;
    plot.name = name;
    plot.scale = scale;
    m_plots.push_back(plot);

    QCheckBox* checkBox = new QCheckBox();
    checkBox->setText(QString(name));
    checkBox->setCheckState(Qt::Checked);
    m_checkBoxes.push_back(checkBox);
    m_vLayoutCheckBoxes->addWidget(checkBox);
    //connect(checkBox, SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
}

void MainWindow::initData()
{

    m_fNe = m_data->getFieldNe();
    m_fEnergy = m_data->getFieldEnergy();
    m_fPhi = m_data->getFieldPhi();
    m_numberHeavySpicies = m_data->getNumberHeavySpicies();
    m_time = 0.0;
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_fHeavy.push_back(m_data->getFieldHeavySpicies(j));
        m_sHeavy.push_back(new solverHeavySpicies(m_data, j));
    }

    for (int i = 0; i < NZ; ++i) {
        m_fNe->arr[i] = simulationTools::gauss(i-2*NZ/3, 10);
        m_fNe->arrPrev[i] = simulationTools::gauss(i-2*NZ/3, 10);
        m_fEnergy->arr[i] = simulationTools::gauss(i-NZ/2, 10);
        m_fEnergy->arrPrev[i] = simulationTools::gauss(i-NZ/2, 10);
        m_fPhi->arr[i] = 0.0;
        m_fPhi->arrPrev[i] = 0.0;
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_fHeavy[j]->arr[i] = simulationTools::gauss(i-NZ/3, 2) / m_numberHeavySpicies;
            m_fHeavy[j]->arrPrev[i] = simulationTools::gauss(i-NZ/3, 2) / m_numberHeavySpicies;
        }
    }
    m_plots.clear();
    while (QLayoutItem* item = m_vLayoutCheckBoxes->takeAt(0)) {
        delete item->widget();
        delete item;
    }
    /*addPlot(m_fNe->arr, m_fNe->name ,m_fNe->cellsNumber);
    addPlot(m_fEnergy->arr, m_fEnergy->name, m_fEnergy->cellsNumber);
    addPlot(m_fPhi->arr, m_fPhi->name, m_fPhi->cellsNumber, 100);
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
      addPlot(m_fHeavy[j]->arr, m_fHeavy[j]->name, m_fHeavy[j]->cellsNumber);
    }*/
    // addPlot(m_data->getReactionRate(simulationData::ReactionName::eAr_eAr), "eAr_eAr" ,m_fNe->cellsNumber, 0.5e12);
    // addPlot(m_data->getReactionRate(simulationData::ReactionName::eAr_eArs), "eAr_eArs" ,m_fNe->cellsNumber,0.5e12);

    m_data->calcReaction(simulationData::ReactionName::eAr_2eArp);
    addPlot(m_data->getReactionRate(simulationData::ReactionName::eAr_2eArp), "eAr_2eArp" ,m_fNe->cellsNumber,0.5e12);

    // replotGraph(0);

    // addPlot(m_data->getReactionRate(simulationData::ReactionName::eArs_2eArp), "eArs_2eArp" ,m_fNe->cellsNumber,0.5e12);
    for (int i = 0; i < m_checkBoxes.size(); ++i) {
        connect(m_checkBoxes[i], SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
    }
    m_scrollBar->setRange(0,1);
    m_storage.clear();
    saveInStorage();
}

void MainWindow::updateData()
{
    for (int i = 0; i < 10; ++i)
    {
        m_data->updateParams();
        m_sNe->solve(5);
        m_sEn->solve(5);
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_sHeavy[j]->solve(5);
        }
        m_sPhi->solve(500);
    }
    m_sNe->getStepEuler();
    m_sEn->getStepEuler();
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_sHeavy[j]->getStepEuler();
    }
    saveInStorage();
    replotGraph(m_storage.size()-1);
}

void MainWindow::simulateData(bool status)
{
    if (status==true)
    {
        initData();
        m_animStopped=false;
        while(!m_animStopped &&  m_time <= (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()))
        {
            m_time += m_textDeltaTime->text().toDouble();
            updateData();
            QCoreApplication::processEvents();
            m_scrollBar->setRange(0, m_storage.size() - 1);
            m_scrollBar->setValue(m_storage.size() - 1);
        }
        m_simulateButton->click();
    }
    else
    {
        m_animStopped=true;
    }
}

void MainWindow::drawDebug(bool)
{



    QVector<double> xf,yf;

    FILE* file=fopen( "c:\\user\\devel\\RFFI_petya\\electro\\arpp.txt","r");
    double aa,bb;
    char buf[1024];
    int i=0;
    while (fgets(buf,1024,file)!=NULL)
    {
        sscanf(buf,"%lf %lf", &aa,&bb);
        xf.push_back(aa);
        yf.push_back(bb*1e-10);
        qDebug()<<"i="<<i<<" "<<xf[i]<<" "<<yf[i];
        i++;
    }


    m_data->calcReaction(simulationData::ReactionName::comsol_eAr_2eArp);
    //addPlot(m_data->getReactionRate(simulationData::ReactionName::eAr_2eArp), "eAr_2eArp" ,m_fNe->cellsNumber,0.5e12);

    QVector<double> xv,yv;

    double *r=m_data->getReactionRate(simulationData::ReactionName::comsol_eAr_2eArp);
    for (int i=0;i<m_data->getCellsNumber();i++)
    {
        yv.push_back(r[i]*1e-10);
        xv.push_back(i*100.0/m_data->getCellsNumber());

        //  qDebug()<<"i="<<i<<" "<<xv[i]<<" "<<yv[i];
    }




    double m_maxY=-1e30;
    double m_minY=+1e30;


    double m_maxX=-1e30;
    double m_minX=+1e30;

    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel("x");
    m_customPlot->yAxis->setLabel("y");
    QColor colors[6] = {QColor(255,0,0),QColor(0,255,0),QColor(0,0,255),QColor(255,255,0),QColor(0,255,255),QColor(255,0,255)};

    m_customPlot->addGraph();
    m_customPlot->graph(0)->setName(QString("ARR"));
    m_customPlot->graph(0)->setLineStyle((QCPGraph::LineStyle)(1));
    m_customPlot->graph(0)->addToLegend();
    m_customPlot->graph(0)->setData(xv,yv);
    QPen graphPen;
    graphPen.setColor(colors[0]);
    graphPen.setWidthF(1.5);
    m_customPlot->graph(0)->setPen(graphPen);


    m_customPlot->addGraph();
    m_customPlot->graph(1)->setName(QString("ARR_coms"));
    m_customPlot->graph(1)->setLineStyle((QCPGraph::LineStyle)(1));
    m_customPlot->graph(1)->addToLegend();
    m_customPlot->graph(1)->setData(xf,yf);

    graphPen.setColor(colors[2]);
    graphPen.setWidthF(1.5);
    m_customPlot->graph(1)->setPen(graphPen);

    for (int i=0; i < xv.size(); ++i)
    {
        if(yv[i] > m_maxY)
            m_maxY = yv[i];
        if(yv[i] < m_minY)
            m_minY = yv[i];

        if(xv[i] > m_maxX)
            m_maxX = xv[i];
        if(xv[i] < m_minX)
            m_minX = xv[i];
    }



    for (int i=0; i < xf.size(); ++i)
    {
        if(yf[i] > m_maxY)
            m_maxY = yf[i];
        if(yf[i] < m_minY)
            m_minY = yf[i];

        if(xf[i] > m_maxX)
            m_maxX = xf[i];
        if(xf[i] < m_minX)
            m_minX = xf[i];
    }



    m_customPlot->xAxis->setRange(m_minX,m_maxX);
    m_customPlot->yAxis->setRange(m_minY, m_maxY);
    m_customPlot->legend->setVisible(true);
    m_customPlot->replot();
}



