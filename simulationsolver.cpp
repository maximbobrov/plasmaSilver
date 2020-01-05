#include "simulationsolver.h"
#include "simulationtools.h"
#include "math.h"
#include <QDebug>


simulationSolver::simulationSolver(simulationData* ipData/*= nullptr*/)
{
    m_pData = ipData;
    m_field = nullptr;
    m_aNu = nullptr;

    m_aRHS = nullptr;//new double[m_field ->cellsNumber];
    m_NewtonRHS = nullptr; //new double[ipData->getFieldNe()->cellsNumber];
}

simulationSolver::~simulationSolver()
{

}

double simulationSolver::solve(int iNumberIteration)
{
    double res = 0.0;
    if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dt = m_pData->getDt();
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();
            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = (dt * m_aNu[j] * (m_field->arr[j+1] + m_field->arr[j-1]) / (dz * dz) +  m_aRHS[j] * dt + m_field->arrPrev[j]) / (1.0 + 2.0 * m_aNu[j] * dt / (dz * dz));
            }
            setBc();

        }
        for (int j=1; j < m_field->cellsNumber - 1; j++)
        {
            res += (m_field->arr[j] - m_field->arrPrev[j]) / dt - (m_aNu[j] * (m_field->arr[j+1] - 2.0 * m_field->arr[j] + m_field->arr[j-1]) / (dz * dz) + m_aRHS[j]);
        }
        return res;
    }
    return -1;
}


double simulationSolver::getRhs()
{
    return -1;
}

double simulationSolver::getRhsAt(int j)
{
    return 0;
}

void simulationSolver::getStepEuler()
{
    for (int i = 0; i < m_field->cellsNumber; ++i)
    {
        m_field->arrPrev[i] = m_field->arr[i];
    }
}

void simulationSolver::setBc()
{
    //periodic bc initally
    double tmp=m_field->arr[1];
    int last=m_field->cellsNumber - 1;
    m_field->arr[0]=m_field->arr[m_field->cellsNumber - 2];
    m_field->arr[last]=tmp;

}

double simulationSolver::getNewtonRhs(int j)
{

    if((j==1) || (j==m_field->cellsNumber - 2))
        setBc();
    //x/dt_x - (x0/dt+nu*(xr+xl)/(dz*dz) + rhs) = x/dt_x -rhs_x
        double rhs=getRhsAt(j);
        double dz = m_pData->getDz();
        double dt =m_pData->getDt();
        return   (m_field->arrPrev[j])/dt + m_aNu[j]*(m_field->arr[j+1] + m_field->arr[j-1])/(dz*dz) + rhs;
       // return m_field->arrPrev[j]/dt;
}

solverNe::~solverNe()
{

}

solverNe::solverNe(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldNe();
    m_aNu = m_pData->getParameters()->arrDe;
    m_aRHS = new double[m_field ->cellsNumber];
    m_NewtonRHS = new double[m_field ->cellsNumber];
}

double solverNe::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();


    double dz = m_pData->getDz();
    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i]=
                pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);
    }
    return 1;
}

double solverNe::getRhsAt(int i)
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();


    double dz = m_pData->getDz();

       return     pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddz(pParams->arrE/*, m_field ->cellsNumber*/, dz, i);

}


void solverNe::setBc()
{
    int last=m_field->cellsNumber - 1;

    double te=m_pData->getParameters()->arrTe[1];
    double dz=m_pData->getDz();
    double E0=(m_pData->getFieldPhi()->arr[1]-m_pData->getFieldPhi()->arr[0])/dz;
    double mu0=m_pData->getParameters()->arrMue[1];
    double D0=m_pData->getParameters()->arrDe[1];
    double gam_gam=0.0;//-1.1e19;
    double nuen=6.69e5*sqrt(te);

    m_field->arr[0]= 1e5;//(-gam_gam + m_field->arr[1]*D0/dz)/(mu0*E0 + D0/dz +nuen);
    m_field->arr[last]=m_field->arr[last-1];//1e5;
    /*
    int last=m_field->cellsNumber - 1;
    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double a = pParams->arrMue[last-1] * (pParams->arrE[last-1] + pParams->arrE[last])/2;
    m_field->arr[last]=m_field->arr[last-1] * (pParams->arrDe[last-1] / m_pData->getDz() - a) / (pParams->arrDe[last-1] / m_pData->getDz() + a);

    a = pParams->arrMue[1] * (pParams->arrE[1] + pParams->arrE[0])/2;
        m_field->arr[0] =m_field->arr[1] * (pParams->arrDe[1] / m_pData->getDz() + a) / (pParams->arrDe[1] / m_pData->getDz() - a);*/
}

solverEnergy::~solverEnergy()
{

}

solverEnergy::solverEnergy(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldEnergy();
    m_aNu = m_pData->getParameters()->arrDeps;
    m_aRHS = new double[m_field ->cellsNumber];
     m_NewtonRHS = new double[m_field ->cellsNumber];
}

double solverEnergy::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField* pNe = m_pData->getFieldNe();
   // double* pAr_star = m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    double dz = m_pData->getDz();

    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*simulationTools::ddzCentral(pNe->arr, m_field ->cellsNumber, dz, i);
        //(pNe->arr[i+1]-pNe->arr[i]);

        m_aRHS[i]=
                pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i)
                - electronFlux * pParams->arrE[i];
    }
    return 1;
}

double solverEnergy::getRhsAt(int i)
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField* pNe = m_pData->getFieldNe();
   // double* pAr_star = m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    double dz = m_pData->getDz();

    double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*simulationTools::ddzCentral(pNe->arr, m_field ->cellsNumber, dz, i);

   return     pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
            + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
            + simulationTools::ddzCentral(m_field->arr , m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i)
            - 1.0*electronFlux * pParams->arrE[i];
}

void solverEnergy::setBc()
{

    int last=m_field->cellsNumber - 1;
    //m_field->arr[0]=m_field->arr[1];//1e5;
    //m_field->arr[last]=m_field->arr[last-1];//1e5;

    double te=m_pData->getParameters()->arrTe[1];
    double dz=m_pData->getDz();
    double E0=(m_pData->getFieldPhi()->arr[1]-m_pData->getFieldPhi()->arr[0])/dz;
    double mu0=m_pData->getParameters()->arrMueps[1];
    double D0=m_pData->getParameters()->arrDeps[1];
    double gam_gam=0.0;//-4e20;
    double nuen=6.69e5*sqrt(te);

    m_field->arr[0]=1e5;// (-gam_gam + m_field->arr[1]*D0/dz)/(mu0*E0 + D0/dz +nuen);
    m_field->arr[last]=m_field->arr[last-1];//1e5;


    /*
    int last=m_field->cellsNumber - 1;
    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double a = pParams->arrMueps[last-1] * (pParams->arrE[last-1] + pParams->arrE[last])/2;
    m_field->arr[last]=m_field->arr[last-1] * (pParams->arrDeps[last-1] / m_pData->getDz() - a) / (pParams->arrDeps[last-1] / m_pData->getDz() + a);

    a = pParams->arrMueps[1] * (pParams->arrE[1] + pParams->arrE[0])/2;
        m_field->arr[0] =m_field->arr[1] * (pParams->arrDeps[1] / m_pData->getDz() + a) / (pParams->arrDeps[1] / m_pData->getDz() - a);*/
}

solverPhi::~solverPhi()
{

}

solverPhi::solverPhi(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldPhi();
    m_aRHS = new double[m_field ->cellsNumber];
     m_NewtonRHS = new double[m_field ->cellsNumber];
}

double solverPhi::getRhs()
{
    double q_e= 1.6*10e-19;
    double eps_0=8.85*10e-12;
    double q_over_eps0=q_e/eps_0;

    double* pNe = m_pData->getFieldNe()->arr;
    double* pArp = m_pData->getFieldHeavySpicies(0)->arr;

    simulationData::simulationParameters* pParams=m_pData->getParameters();
    // double mult= (pParams->p*6.022e23)/(pParams->T*8.314);
    //m_fHeavy[j]->arr[i] =(m_fNe->arr[i]*pParams->T*8.314)/(pParams->p*6.022e23);

    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i] = q_over_eps0*0.5*(-(pNe[i]+pNe[i+1])+(pArp[i]+pArp[i+1]));
    }
    return 1;
}

double solverPhi::solve(int iNumberIteration)
{
    double res = 0.0;
    if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dzdz=dz*dz;
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();

            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = 0.5* ((m_field->arr[j+1] + m_field->arr[j-1]) +  m_aRHS[j]*dzdz);
            }
            setBc();
        }
        return res;
    }
    return -1;
}

void solverPhi::setBc()
{
    int last=m_field->cellsNumber - 1;
    m_field->arr[0] = 0;
    m_field->arr[last]=100;
}


solverHeavySpicies::solverHeavySpicies(simulationData *pData, int num)
{
    m_pData = pData;
    m_field = m_pData->getFieldHeavySpicies(num);
    m_specie=m_field->m_specie;
    m_charge = m_pData->getHeavySpiciesCharge(num);
    m_aNu = m_pData->getParameters()->arrDomega;
    m_aRHS = new double[m_field ->cellsNumber];
     m_NewtonRHS = new double[m_field ->cellsNumber];
}

solverHeavySpicies::~solverHeavySpicies()
{

}

double solverHeavySpicies::getRhs()
{

    simulationData::simulationField* pNe = m_pData->getFieldNe();

    double* pAr_star = m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dz = m_pData->getDz();




    if (m_specie==simulationData::SpecieName::Ar_plus)
    {


        for (int i = 0; i < m_field ->cellsNumber-1; ++i)
        {

            m_aRHS[i]=0.0;/*
                    -  m_charge  * pParams->arrMuomega[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                    -  m_charge  * pParams->arrMuomega[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);*/
        }


    }else if (m_specie==simulationData::SpecieName::Ar_star)
    {



        for (int i = 0; i < m_field ->cellsNumber-1; ++i)
        {

            m_aRHS[i]=0.0;/*
                    -  m_charge  * pParams->arrMuomega[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                    -  m_charge  * pParams->arrMuomega[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);*/
        }
    }

    /*if (m_specie==simulationData::SpecieName::Ar_star)
    {
    //m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_eArs);
    m_pData->calcReaction(simulationData::ReactionName::comsol_eArs_2eArp);
    R1=m_pData->getReactionRate(simulationData::ReactionName::comsol_eArs_2eArp);
    mult1=pParams->p / ( pParams->T * 8.314 * 6.022e23 * pParams->rho);
    }*/



    return 1;
}

double solverHeavySpicies::getRhsAt(int j)
{
    return 0.0;
}

void solverHeavySpicies::setBc()
{
    int last=m_field->cellsNumber - 1;
    m_field->arr[0] =m_field->arr[1];//m_field->arr[1]*0.5;//  m_field->arr[1];
    m_field->arr[last]=m_field->arr[last-1] ;
}


