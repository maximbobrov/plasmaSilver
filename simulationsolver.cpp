#include "simulationsolver.h"
#include "simulationtools.h"
#include "math.h"


simulationSolver::simulationSolver(simulationData* ipData/*= nullptr*/)
{
    m_pData = ipData;
    m_aRHS = nullptr;
    m_field = nullptr;
    m_aNu = nullptr;
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

        }
        for (int j=1; j < m_field->cellsNumber - 1; j++)
        {
            res += (m_field->arr[j] - m_field->arrPrev[j]) / dt - (m_aNu[j] * (m_field->arr[j+1] - 2.0 * m_field->arr[j] + m_field->arr[j-1]) / (dz * dz) + m_aRHS[j]);
        }
        return res;
    }
    return -1;
}

double simulationSolver::updateParams()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField* pEn = m_pData->getFieldEnergy();
    simulationData::simulationField* pPhi = m_pData->getFieldPhi();
    double dz = m_pData->getDz();
    for (int j=0; j<m_field->cellsNumber; j++)
    {
        pParams->arrMue[j] =2.6e3;//4e4;
        pParams->arrMueps[j] =5.0 * pParams->arrMue[j] / 3.0;
        pParams->arrDe[j] = pParams->arrMue[j] * 2.0*pEn->arr[j] / 3.0;
        pParams->arrDeps[j] = pParams->arrMueps[j] * 2.0 * pEn->arr[j] / 3.0;
        pParams->arrE[j] = - simulationTools::ddzCentral(pPhi->arr, pPhi->cellsNumber, dz, j);
        pParams->arrMuomega[j]=8.0e3/760.0;
        pParams->arrDomega[j] = 0.15;
    }
}

double simulationSolver::getRhs()
{

}

double simulationSolver::getStepEuler()
{
    for (int i = 0; i < m_field->cellsNumber; ++i)
    {
         m_field->arrPrev[i] = m_field->arr[i];
    }
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

}

double solverNe::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dz = m_pData->getDz();
    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i]= + pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                   + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDe, m_field ->cellsNumber, dz, i);
    }
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
}

double solverEnergy::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField* pNe = m_pData->getFieldNe();
    double dz = m_pData->getDz();
    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*(pNe->arr[i+1]-pNe->arr[i]);
        m_aRHS[i]= + pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                   + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i)
                   + electronFlux * pParams->arrE[i];
    }
}

solverPhi::~solverPhi()
{

}

solverPhi::solverPhi(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldPhi();
    m_aRHS = new double[m_field ->cellsNumber];
}

double solverPhi::getRhs()
{
    simulationData::simulationField* pNe = m_pData->getFieldNe();
    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i] = ((pNe->arr[i])+(pNe->arr[i+1]))/(2*4*M_PI);
    }
}

double solverPhi::solve(int iNumberIteration)
{
    double res = 0.0;
    if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();
            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = ((m_field->arr[j+1] + m_field->arr[j-1]) / (dz * dz) +  m_aRHS[j]) / (2.0 / (dz * dz));
            }
        }
        return res;
    }
    return -1;
}

solverHeavySpicies::solverHeavySpicies(simulationData *pData, int num)
{
    m_pData = pData;
    m_field = m_pData->getFieldHeavySpicies(num);
    m_charge = m_pData->getHeavySpiciesCharge(num);
    m_aNu = m_pData->getParameters()->arrDomega;
    m_aRHS = new double[m_field ->cellsNumber];
}

solverHeavySpicies::~solverHeavySpicies()
{

}

double solverHeavySpicies::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dz = m_pData->getDz();
    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i]= - 10000 * m_charge  * pParams->arrMuomega[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                   - 10000 * m_charge  * pParams->arrMuomega[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);
    }
}
