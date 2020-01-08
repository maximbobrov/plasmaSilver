#include "simulationdata.h"
#include "simulationtools.h"
#include "string.h"
#include "reaction.h"
#include "math.h"
#include <QDebug>

simulationData::simulationData(int iCellsNumber)
{
    m_defaultCellsNumber = iCellsNumber;
    m_fieldNe = new simulationField(iCellsNumber,"electrons",simulationData::SpecieName::e);
    m_fieldEnergy = new simulationField(iCellsNumber,"energy",simulationData::SpecieName::En);
    m_fieldPhi = new simulationField(iCellsNumber+1,"potential",simulationData::SpecieName::phi);
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar+",simulationData::SpecieName::Ar_plus));
    m_chargeHeavySpecies.push_back(1);
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ars",simulationData::SpecieName::Ar_star));
    m_chargeHeavySpecies.push_back(0);


    m_numberHeavySpicies = m_fieldsHeavySpecies.size();
    m_params = new simulationData::simulationParameters(iCellsNumber);

    m_reactions.push_back(new reactionEAr_EAr_comsol(this));
    m_reactions.push_back(new reactionEAr_EArs_comsol(this));
    m_reactions.push_back(new reactionEArs_EAr_comsol(this));
    m_reactions.push_back(new reactionEAr_2EArp_comsol(this));
    m_reactions.push_back(new reactionEArs_2EArp_comsol(this));
    m_reactions.push_back(new reactionArsArs_EArArp_comsol(this));
    m_reactions.push_back(new reactionArsAr_ArAr_comsol(this));

}

simulationData::~simulationData()
{

}

void simulationData::setDt(double idt)
{
    m_dt = idt;
}

void simulationData::setDz(double idz)
{
    m_dz = idz;
}

void simulationData::setCellsNumber(int iCellsNumber)
{
    m_defaultCellsNumber = iCellsNumber;
}

double simulationData::getDt()
{
    return m_dt;
}

double simulationData::getDz()
{
    return m_dz;
}

simulationData::simulationField* simulationData::getFieldNe()
{
    return m_fieldNe;
}

double *simulationData::getArrTe()
{
   return m_params->arrTe;
}


simulationData::simulationField* simulationData::getFieldHeavySpicies(int num)
{
    return m_fieldsHeavySpecies[num];
}

int simulationData::getHeavySpiciesCharge(int num)
{
    return m_chargeHeavySpecies[num];
}

int simulationData::getNumberHeavySpicies()
{
    return m_numberHeavySpicies;
}

double simulationData::getN()
{
   // pv=nkT
   // n/v=p/kT;
   // double pres=13.21;//101505;
   // double T=293.15;
    double k=1.38e-23;


    return m_params->p/(k*m_params->T);
}

double *simulationData::getReactionRate(simulationData::ReactionName reactName)
{
    return m_reactions[reactName]->getR();
}

double simulationData::getReactionDe(simulationData::ReactionName reactName)
{
     return m_reactions[reactName]->getDe();
}

void simulationData::calcReaction(simulationData::ReactionName reactName)
{
    m_reactions[reactName]->calc();
}

simulationData::simulationField* simulationData::getFieldEnergy()
{
    return m_fieldEnergy;
}

simulationData::simulationField* simulationData::getFieldPhi()
{
    return m_fieldPhi;
}

simulationData::simulationParameters *simulationData::getParameters()
{
    return m_params;
}

int simulationData::getCellsNumber()
{
    return m_defaultCellsNumber;
}

simulationData::simulationField::simulationField(int iCellsNumber, char* iName,SpecieName specie)
{

    init(iCellsNumber, iName,specie);
}

void simulationData::simulationField::init(int iCellsNumber, char* iName, SpecieName specie)
{
  cellsNumber =  iCellsNumber;
  arr = new double[iCellsNumber];
  arrPrev = new double[iCellsNumber];
  int len=strlen(iName);
  name = new char[len+1];
  name[len]=0;
  strcpy(name, iName);
  m_specie=specie;
}

void simulationData::simulationParameters::init(int iCellsNumber)
{
    cellsNumber =  iCellsNumber;
    arrDe = new double[iCellsNumber];
    arrDeps = new double[iCellsNumber];
    arrDomega = new double[iCellsNumber];
    arrMue = new double[iCellsNumber];
    arrMueps = new double[iCellsNumber];
    arrMuomega = new double[iCellsNumber];
    arrE = new double[iCellsNumber];
    arrTe = new double[iCellsNumber];

    T=293.15; //K
    p=13.21;//pa
    mAr=39.948/1000.0;//kg/mol
    rho=p*mAr/(8.314*T);//kg/m^3
    double Na=6.022e23; //1/mol
    N=p*Na/(8.314*T);

}

void simulationData::updateParams()
{
    simulationData::simulationParameters* pParams = m_params;
    simulationData::simulationField* pEn = m_fieldEnergy;
    simulationData::simulationField* pNe = m_fieldNe;
    simulationData::simulationField* pPhi = m_fieldPhi;
    double dz = m_dz;


    for (int j=0; j<pParams->cellsNumber; j++)
    {
        pParams->arrMue[j] = 1e25 / m_params->N; ; //4e4; m^2/(V*s)
        pParams->arrMueps[j] =5.0 * pParams->arrMue[j] / 3.0;

        pParams->arrTe[j]=/*3.0;//*/fmax(fmin((2.0/3.0)*(fabs(pEn->arr[j]))/(fabs(pNe->arr[j])+1e8),40),1.5);


        if (pParams->arrTe[j]>60.0) pParams->arrTe[j]=60.0; //some limiter here


        pParams->arrDe[j] = pParams->arrMue[j] * pParams->arrTe[j];
        pParams->arrDeps[j] = pParams->arrMueps[j] * pParams->arrTe[j];
        pParams->arrE[j] = -(pPhi->arr[j+1]-pPhi->arr[j])/dz;//- simulationTools::ddzCentral(pPhi->arr, pPhi->cellsNumber, dz, j);
        pParams->arrDomega[j] = 0.01;// m^2/s
        pParams->arrMuomega[j]=pParams->arrDomega[j]*q/(k_B_const*pParams->T);
    }
}



simulationData::simulationParameters::simulationParameters(int iCellsNumber)
{
    init(iCellsNumber);
}
