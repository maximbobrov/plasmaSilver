#include "mainwindow.h"
#include "qcustomplot.h"
#include <QGridLayout>
#include <QPushButton>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <stdio.h>
#include <reactionsolver.h>



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
    m_stepButton = new QPushButton("solve 1");
    m_advectButton = new QPushButton("advect 1");

    m_grid->addWidget(m_customPlot,0,0,3,3);
    m_grid->addLayout(m_vLayoutCheckBoxes,0,3,1,1);
    m_grid->addWidget(m_simulateButton,5,0);
    m_grid->addWidget(m_debugButton,5,1);
    m_grid->addWidget(m_stepButton,5,2);
    m_grid->addWidget(m_advectButton,5,3);

    m_hLayout->addWidget(new QLabel("Start time:"));
    m_hLayout->addWidget(m_textStartTime);
    m_hLayout->addWidget(new QLabel("dt:"));
    m_hLayout->addWidget(m_textDeltaTime);
    m_hLayout->addWidget(new QLabel("End time:"));
    m_hLayout->addWidget(m_textEndTime);
    m_grid->addLayout(m_hLayout,4,0,1,2);
    m_grid->addWidget(m_scrollBar,4,2,1,1);




    m_data = new simulationData(NZ);

    m_rSolver= new reactionSolver(m_data);

    m_data->setDz(0.368/NZ);
    m_data->setDt(1e-11);

    m_sNe = new solverNe(m_data);
    m_sEn = new solverEnergy(m_data);
    m_sPhi = new solverPhi(m_data);

    m_textStartTime->setText(QString().number(0.0));
    m_textDeltaTime->setText(QString().number(m_data->getDt()));
    m_textEndTime->setText(QString().number(1e-7));


    connect (m_simulateButton, SIGNAL(clicked(bool)), this, SLOT(simulateData(bool)));
    connect (m_debugButton, SIGNAL(clicked(bool)), this, SLOT(drawDebug(bool)));
    connect (m_stepButton, SIGNAL(clicked(bool)), this, SLOT(singleStep(bool)));
    connect (m_advectButton, SIGNAL(clicked(bool)), this, SLOT(singleAdvect(bool)));

    connect(m_scrollBar, SIGNAL(valueChanged(int)), this, SLOT(replotGraph(int)));

    m_simulateButton->setCheckable(true);

    m_widget->setLayout(m_grid);
    setCentralWidget(m_widget);
    setWindowTitle("PlasmaSolver");
    m_animStopped=true;
    initData();
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
    m_customPlot->xAxis->setLabel(QString("x, time= %1").arg(m_time));
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

    // m_data->setDz(2e-4/NZ);
    //  m_data->setDt(1e-14);

    simulationData::simulationParameters* pParams=m_data->getParameters();
    for (int i = 0; i < NZ; ++i) {

        double x_=i*m_data->getDz();
        m_fNe->arr[i] =5e11;//1e5+ 2.5e11*simulationTools::gauss(x_-0.15, 0.05);
        m_fNe->arrPrev[i] =m_fNe->arr[i];
        m_fEnergy->arr[i] = 5.0*m_fNe->arr[i];
        m_fEnergy->arrPrev[i] = m_fEnergy->arr[i];

        //for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_fHeavy[simulationData::SpecieName::Ar_plus]->arr[i] =m_fNe->arr[i];
            m_fHeavy[simulationData::SpecieName::Ar_plus]->arrPrev[i] = m_fHeavy[simulationData::SpecieName::Ar_plus]->arr[i];
            m_fHeavy[simulationData::SpecieName::Ar_star]->arr[i] =3.26395e13;//3.5372e13;
            m_fHeavy[simulationData::SpecieName::Ar_star]->arrPrev[i] = m_fHeavy[simulationData::SpecieName::Ar_star]->arr[i];
        }
        m_fPhi->arr[i] = 0.0;
        m_fPhi->arrPrev[i] = 0.0;
    }
    /*m_fNe->arr[0] =1e5+ 1e11*simulationTools::gauss(x_-1e-4, 2e-5);
    m_fNe->arrPrev[0] =m_fNe->arr[0];
    m_fEnergy->arr[0] = 5.0*m_fNe->arr[0];
    m_fEnergy->arrPrev[0] = m_fEnergy->arr[0];

    /*
    m_fNe->arr[NZ-1] =1e6;//1e5+ 1e11*simulationTools::gauss(x_-1e-4, 2e-5);
    m_fNe->arrPrev[NZ-1] =m_fNe->arr[NZ-1];
    m_fEnergy->arr[NZ-1] = 5.0*m_fNe->arr[NZ-1];
    m_fEnergy->arrPrev[NZ-1] = m_fEnergy->arr[NZ-1];
*/

    m_sPhi->solve(100);
    /* m_sPhi->solve(100000);
     m_sPhi->solve(100000);
       m_sPhi->solve(100000);
         m_sPhi->solve(100000);*/
    m_data->updateParams();

    m_plots.clear();
    while (QLayoutItem* item = m_vLayoutCheckBoxes->takeAt(0)) {
        delete item->widget();
        delete item;
    }
    addPlot(m_fNe->arr, m_fNe->name ,m_fNe->cellsNumber);
    addPlot(/*m_fEnergy->arr*/pParams->arrTe, m_fEnergy->name, m_fEnergy->cellsNumber);
    addPlot(/*pParams->arrE*/m_fPhi->arr, m_fPhi->name, m_fPhi->cellsNumber-1, 1.0);
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        addPlot(m_fHeavy[j]->arr, m_fHeavy[j]->name, m_fHeavy[j]->cellsNumber);
    }

    for (int i = 0; i < m_checkBoxes.size(); ++i) {
        connect(m_checkBoxes[i], SIGNAL(stateChanged(int)), this, SLOT(replotGraph(int)));
    }
    m_scrollBar->setRange(0,1);
    m_storage.clear();
    saveInStorage();
    replotGraph(0);
}

void MainWindow::updateData()
{

   double dt= m_data->getDt();
   // qDebug()<<"dt="<<dt;
 // m_data->setDt(1e+5);
//m_sPhi->solve(100);
   dt= 1.e-12;
   m_data->setDt(dt);
 /*    for (int i = 0; i < 2; ++i)
    {
        m_data->updateParams();
        m_sNe->solve(5);
       m_sEn->solve(5);
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_sHeavy[j]->solve(5);
        }
        m_sPhi->solve(5000);

    }*/
 solveNewton();

  /*  for (int i = 0; i < 10; ++i) {
        if(!solveNewton())
        {dt/=2.0;
        m_data->setDt(dt);}

    }
    /* solveNewton();
    solveNewton();
    solveNewton();*/

   m_sNe->getStepEuler();
    m_sEn->getStepEuler();
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_sHeavy[j]->getStepEuler();
    }

}

bool MainWindow::solveNewton()
{


    double *ne=m_data->getFieldNe()->arr;
    double *ars=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;
    double *arp=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_plus)->arr;

    double *en=m_data->getFieldEnergy()->arr;

    double *ne0=m_data->getFieldNe()->arrPrev;
    double *ars0=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev;
    double *arp0=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_plus)->arrPrev;
    double *en0=m_data->getFieldEnergy()->arrPrev;


    //set leftbc;

/*

    m_rSolver->n_n=m_data->getParameters()->N;

    m_rSolver->dt=m_data->getDt();
    m_rSolver->dz=m_data->getDz();

    m_rSolver->nu_ne=m_data->getParameters()->arrDe[1];
    m_rSolver->nu_narp=m_data->getParameters()->arrDomega[1];
    m_rSolver->nu_neps=m_data->getParameters()->arrDeps[1];;
    m_rSolver->nu_nars=m_data->getParameters()->arrDomega[1];;*/

for (int nn=0;nn<40;nn++)
{
    m_sPhi->solve(1);
    m_data->updateParams();

    for (int j=1; j <m_data->getFieldNe()->cellsNumber - 1; j++)
    {
      /*  if (j==1)
        {
        qDebug()<<"dt="<<m_rSolver->dt<<" ne="<<m_rSolver->ne_o<<" nars="<<m_rSolver->nars_o <<" eps="<<m_rSolver->eps_o;
        qDebug()<<" ne0="<<ne[0];
        }

        m_rSolver->ne_i=ne[j];
        m_rSolver->ne_0=ne0[j];
        m_rSolver->nars_i=ars[j];//3.2983e13;
        m_rSolver->nars_0=ars0[j];//3.2983e13;

        m_rSolver->eps_o=en[j]/ne[j];
        //   m_rSolver->n_n=m_data->getParameters()->N;
        m_rSolver->dt=m_data->getDt();

        //x/dt_x - (x0/dt-nu*(xr+xl)/(dz*dz) + rhs) = x/dt_x -rhs_x
        m_rSolver->rhs_neps=m_sEn->getNewtonRhs(j);
        m_rSolver->rhs_nars=m_sHeavy[simulationData::SpecieName::Ar_star]->getNewtonRhs(j);
        m_rSolver->rhs_ne=m_sNe->getNewtonRhs(j);
        m_rSolver->rhs_narp=m_sHeavy[simulationData::SpecieName::Ar_plus]->getNewtonRhs(j);

        if(!m_rSolver->solve_diffuse(4))
            return false;
            */

       /*ne[j]=m_rSolver->ne_o;
        ars[j]=m_rSolver->nars_o;
        en[j]=m_rSolver->neps_o;
        arp[j]=m_rSolver->narp_o;*/

        ne[j]=ne[j]*0.5+0.5*m_sNe->getNewtonRhs(j);
       // ars[j]=m_sNe->getNewtonRhs(j);

        en[j]=en[j]*0.5+0.5*m_sEn->getNewtonRhs(j);
        ars[j]=ars[j]*0.5+0.5*m_sHeavy[simulationData::SpecieName::Ar_star]->getNewtonRhs(j);//*m_pParam-> arrMaskNe[i][j];;
        arp[j]=arp[j]*0.5+0.5*m_sHeavy[simulationData::SpecieName::Ar_plus]->getNewtonRhs(j);//*m_pParam-> arrMaskNe[i][j];;


    }
}
double netot=0.0;
for (int j=1; j <m_data->getFieldNe()->cellsNumber - 1; j++)
{
    netot+=ne[j];
}

   qDebug()<<"ne="<<netot/4.96e8;
m_sPhi->solve(1);

    /*  for ( j=m_data->getFieldNe()->cellsNumber - 3; j >1; j--)
    {
        m_rSolver->ne_i=ne[j];
        m_rSolver->ne_0=ne0[j];
        m_rSolver->nars_i=ars[j];//3.2983e13;
        m_rSolver->nars_0=ars0[j];//3.2983e13;

        m_rSolver->eps_o=en[j]/ne[j];
        //x/dt_x - (x0/dt-nu*(xr+xl)/(dz*dz) + rhs) = x/dt_x -rhs_x
        m_rSolver->rhs_neps=m_sEn->getNewtonRhs(j);
        m_rSolver->rhs_nars=m_sHeavy[simulationData::SpecieName::Ar_star]->getNewtonRhs(j);
        m_rSolver->rhs_ne=m_sNe->getNewtonRhs(j);
        m_rSolver->rhs_narp=m_sHeavy[simulationData::SpecieName::Ar_plus]->getNewtonRhs(j);


        m_rSolver->solve_diffuse(4);

        ne[j]=m_rSolver->ne_o;
        ars[j]=m_rSolver->nars_o;
        en[j]=m_rSolver->neps_o;
        arp[j]=m_rSolver->narp_o;
    }*/
    //
   // m_sPhi->solve(5500);


    /*  qDebug()<<"N="<<m_data->getParameters()->N;
    for (int i = 0; i < 2; ++i)
    {
        m_data->updateParams();
        m_sNe->solve(2);
        m_sEn->solve(2);
        //for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            //m_sHeavy[j]->solve(2);
        }
          m_sHeavy[simulationData::SpecieName::Ar_star]->solve(2);
        m_sPhi->solve(25);

    }

    double *Ne_rhs=m_sNe->getNewtonRhs();
    double *En_rhs=m_sNe->getNewtonRhs();
    double *Ars_rhs=m_sHeavy[simulationData::SpecieName::Ar_star]->getNewtonRhs();

    double *ne=m_data->getFieldNe()->arr;
    double *ars=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;
    double *arp=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_plus)->arr;

    double *en=m_data->getFieldEnergy()->arr;

    double *ne0=m_data->getFieldNe()->arrPrev;
    double *ars0=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev;
    double *arp0=m_data->getFieldHeavySpicies(simulationData::SpecieName::Ar_plus)->arrPrev;
    double *en0=m_data->getFieldEnergy()->arrPrev;

     for (int j=1; j <m_data->getFieldNe()->cellsNumber - 1; j++)
     {
         m_rSolver->rhs_eps=0.0;//En_rhs[j];
         m_rSolver->rhs_nars=0.0;//Ars_rhs[j];
         m_rSolver->rhs_ne=0.0;//Ne_rhs[j];

         m_rSolver->ne_i=ne[j];//1e13;
         m_rSolver->nars_i=ars[j];//3.5372e13;//3.2983e13;
         m_rSolver->n_n=m_data->getParameters()->N;//3.5372e21;
         m_rSolver->eps_i=fabs(en[j])/(fabs(ne[j]+1));//5.0;//3.3333;
         m_rSolver->dt=m_data->getDt();

         m_rSolver->ne_0=ne0[j];//1e13;
         m_rSolver->nars_0=ars0[j];//3.5372e13;//3.2983e13;
         m_rSolver->eps_0=fabs(en0[j])/(fabs(ne0[j]+1));//5.0;//3.3333;


         m_rSolver->solve(4);

         ne[j]+=m_rSolver->ne_o-m_rSolver->ne_i;
         ars[j]+=m_rSolver->nars_o-m_rSolver->nars_i;
         en[j]+=m_rSolver->eps_o*m_rSolver->ne_o-m_rSolver->eps_i*m_rSolver->ne_i;
         arp[j]+=(m_rSolver->ne_o-m_rSolver->ne_i);
     }*/

    return true;

}

void MainWindow::simulateData(bool status)
{
    if (status==true)
    {

        m_animStopped=false;
        while(!m_animStopped &&  m_time <= (m_textEndTime->text().toDouble() - m_textStartTime->text().toDouble()))
        {
            m_time += 10.0*m_textDeltaTime->text().toDouble();
            for (int i=0;i<10;i++)
            {
                updateData();
            }

            saveInStorage();
            replotGraph(m_storage.size()-1);

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

    //reactionSolver *solver=new reactionSolver(m_data);


    /*
    double e=0;
    QVector<double> rr[7], dd[7] , ee;
    while (e < 40)
    {
        ee.push_back(e);
        for (int i=0;i<7;i++)
        {
            rr[i].push_back(solver->m_reactions[i]->getRate(e));
            dd[i].push_back(solver->m_reactions[i]->getDeriv(e));
        }
        e+=0.01;
    }
    double m_maxY=-1e30;
    double m_minY=+1e30;


    double m_maxX=-1e30;
    double m_minX=+1e30;

    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel("x");
    m_customPlot->yAxis->setLabel("y");
    QColor colors[7] = {QColor(255,0,0),QColor(0,255,0),QColor(0,0,255),QColor(255,255,0),QColor(0,255,255),QColor(255,0,255),QColor(0,0,0)};


    for (int i=0;i<7;i++)
    {

        m_customPlot->addGraph();
        m_customPlot->graph(i)->setName(QString("num = ").number(i));
        m_customPlot->graph(i)->setLineStyle((QCPGraph::LineStyle)(1));
        m_customPlot->graph(i)->addToLegend();
        m_customPlot->graph(i)->setData(ee,rr[i]);
        QPen graphPen;
        graphPen.setColor(colors[i]);
        graphPen.setWidthF(1.5);
        m_customPlot->graph(i)->setPen(graphPen);

        for (int j=0; j < rr[i].size(); ++j)
        {
            if(rr[i][j] > m_maxY)
                m_maxY = rr[i][j];
            if(rr[i][j] < m_minY)
                m_minY = rr[i][j];

            if(ee[j] > m_maxX)
                m_maxX = ee[j];
            if(ee[j] < m_minX)
                m_minX = ee[j];
        }
    }
    m_customPlot->xAxis->setRange(m_minX,m_maxX);
    m_customPlot->yAxis->setRange(m_minY, m_maxY);
    m_customPlot->legend->setVisible(true);
    m_customPlot->replot();*/


    qDebug()<<"aaa!";



    m_rSolver->ne_i=1e13;
    m_rSolver->ne_0=1e13;
    m_rSolver->eps_o=5.0;
    m_rSolver->nars_i=3.5372e13;//3.2983e13;
    m_rSolver->nars_0=3.5372e13;//3.2983e13;
    m_rSolver->n_n=3.5372e21;

    m_rSolver->dt=0.5e-9;
    m_rSolver->dz=1e-4;

    m_rSolver->nu_ne=0.0;
    m_rSolver->nu_narp=0.0;
    m_rSolver->nu_neps=0.0;
    m_rSolver->nu_nars=0.0;


    //x/dt_x - (x0/dt-nu*(xr+xl)/(dz*dz) + rhs) = x/dt_x -rhs_x
    m_rSolver->rhs_neps=5.0*m_rSolver->ne_0/m_rSolver->dt;
    m_rSolver->rhs_nars=m_rSolver->nars_0/m_rSolver->dt;;
    m_rSolver->rhs_ne=m_rSolver->ne_0/m_rSolver->dt;
    m_rSolver->rhs_narp=m_rSolver->ne_0/m_rSolver->dt;

    QVector<double> ne,nars,eps,t;

    //ne.push_back(m_rSolver->ne_i);
    //nars.push_back(m_rSolver->nars_i);
    //eps.push_back(m_rSolver->eps_i);
    t.push_back(0.0);
    double t_cur=0.0;
    for (int i=0;i<300;i++)
    {

        //    m_rSolver->dt=1e-9*pow((i+1),0.71);
        t_cur+=m_rSolver->dt;
        //qDebug()<<"t="<<t_cur <<" eps="<<m_rSolver->eps_o<<" ne="<<m_rSolver->ne_o<<" ars="<<m_rSolver->nars_o;
        // for (int j=0;j<100;j++)
        m_rSolver->solve_diffuse(4);
        m_rSolver->ne_i=m_rSolver->ne_o;
        m_rSolver->nars_i=m_rSolver->nars_o;

        m_rSolver->ne_0=m_rSolver->ne_o;
        m_rSolver->nars_0=m_rSolver->nars_o;
        m_rSolver->dt=1e-9*pow((i+1),0.71);

        //x/dt_x - (x0/dt-nu*(xr+xl)/(dz*dz) + rhs) = x/dt_x -rhs_x
        m_rSolver->rhs_neps=m_rSolver->neps_o/m_rSolver->dt;
        m_rSolver->rhs_nars=m_rSolver->nars_o/m_rSolver->dt;
        m_rSolver->rhs_ne=m_rSolver->ne_o/m_rSolver->dt;
        m_rSolver->rhs_narp=m_rSolver->narp_o/m_rSolver->dt;





        ne.push_back(m_rSolver->ne_o/1e13);
        nars.push_back(m_rSolver->nars_o/3.5372e13);
        eps.push_back(m_rSolver->eps_o/5.0);
        t.push_back(t_cur);


        qDebug()<<"dt="<<m_rSolver->dt<<"t="<<t_cur <<" ne="<<m_rSolver->ne_o/1e13<<" nars="<<m_rSolver->nars_o/3.5372e13 <<" eps="<<m_rSolver->eps_o/5.0;


        /* if (m_rSolver->ne_i<1e5) m_rSolver->ne_i=1e5;
        if (m_rSolver->eps_i<0.01) m_rSolver->eps_i=0.01;
        if (m_rSolver->eps_i>40.0) m_rSolver->eps_i=40.0;

        if (m_rSolver->nars_i<1e5) solver->nars_i=1e5;*/
    }



    double m_maxY=-1e50;
    double m_minY=+1e50;


    double m_maxX=-1e50;
    double m_minX=+1e50;

    m_customPlot->clearGraphs();
    m_customPlot->xAxis->setLabel("x");
    m_customPlot->yAxis->setLabel("y");
    QColor colors[6] = {QColor(255,0,0),QColor(0,255,0),QColor(0,0,255),QColor(255,255,0),QColor(0,255,255),QColor(255,0,255)};

    m_customPlot->addGraph();
    m_customPlot->graph(0)->setName(QString("nars"));
    m_customPlot->graph(0)->setLineStyle((QCPGraph::LineStyle)(1));
    m_customPlot->graph(0)->addToLegend();
    m_customPlot->graph(0)->setData(t,nars);
    QPen graphPen;
    graphPen.setColor(colors[0]);
    graphPen.setWidthF(1.5);
    m_customPlot->graph(0)->setPen(graphPen);

    m_customPlot->addGraph();
    m_customPlot->graph(1)->setName(QString("ne"));
    m_customPlot->graph(1)->setLineStyle((QCPGraph::LineStyle)(1));
    m_customPlot->graph(1)->addToLegend();
    m_customPlot->graph(1)->setData(t,ne);

    graphPen.setColor(colors[2]);
    graphPen.setWidthF(1.5);
    m_customPlot->graph(1)->setPen(graphPen);


    m_customPlot->addGraph();
    m_customPlot->graph(2)->setName(QString("eps"));
    m_customPlot->graph(2)->setLineStyle((QCPGraph::LineStyle)(1));
    m_customPlot->graph(2)->addToLegend();
    m_customPlot->graph(2)->setData(t,eps);

    graphPen.setColor(colors[1]);
    graphPen.setWidth(1.5);
    m_customPlot->graph(2)->setPen(graphPen);

    for (int i=0; i < nars.size()-1; ++i)
    {
        if(nars[i] > m_maxY)
            m_maxY = nars[i];
        if(nars[i] < m_minY)
            m_minY = nars[i];

        if(t[i] > m_maxX)
            m_maxX = t[i];
        if(t[i] < m_minX)
            m_minX = t[i];
    }

    qDebug()<<"min="<<m_minY<<" max="<<m_maxY;

    m_customPlot->xAxis->setRange(m_minX,m_maxX);
    m_customPlot->yAxis->setRange(0,1.06);//m_minY, m_maxY);
    m_customPlot->legend->setVisible(true);
    m_customPlot->replot();
}

void MainWindow::singleStep(bool)
{
    m_time += 10.0*m_textDeltaTime->text().toDouble();
    {
        /*  m_data->updateParams();
        m_sNe->solve(5);
        m_sEn->solve(5);
        for (int j = 0; j < m_numberHeavySpicies; ++j)
        {
            m_sHeavy[j]->solve(5);
        }
        m_sPhi->solve(50);*/
        solveNewton();
    }

    /* saveInStorage();
    replotGraph(m_storage.size()-1);

    QCoreApplication::processEvents();
    m_scrollBar->setRange(0, m_storage.size() - 1);
    m_scrollBar->setValue(m_storage.size() - 1);*/
}

void MainWindow::singleAdvect(bool)
{
    m_sNe->getStepEuler();
    m_sEn->getStepEuler();
    for (int j = 0; j < m_numberHeavySpicies; ++j)
    {
        m_sHeavy[j]->getStepEuler();
    }
}



