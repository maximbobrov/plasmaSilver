#include "simulationdata.h"
#include "simulationtools.h"
#include "string.h"
#include "reaction.h"
#include <QDebug>

simulationData::simulationData(int iCellsNumber)
{
    m_defaultCellsNumber = iCellsNumber;
    m_fieldNe = new simulationField(iCellsNumber,"electrons");
    m_fieldEnergy = new simulationField(iCellsNumber,"energy");
    m_fieldPhi = new simulationField(iCellsNumber+1,"potential");
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar"));
    m_chargeHeavySpecies.push_back(0);
 //   m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ars"));
  //  m_chargeHeavySpecies.push_back(1);
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar+"));
    m_chargeHeavySpecies.push_back(1);
    m_numberHeavySpicies = m_fieldsHeavySpecies.size();
    m_params = new simulationData::simulationParameters(iCellsNumber);

    m_reactions.push_back(new reactionEAr_EAr(this));
    m_reactions.push_back(new reactionEAr_EArs(this));
    m_reactions.push_back(new reactionEAr_2EArp(this));
    m_reactions.push_back(new reactionEArs_2EArp(this));
    m_reactions.push_back(new reactionEAr_2EArp_comsol(this));
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
    double pres=101505;
    double T=300;
    double k=1.38e-23;


    return pres/(k*T);
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

simulationData::simulationField::simulationField(int iCellsNumber, char* iName)
{

    init(iCellsNumber, iName);
}

void simulationData::simulationField::init(int iCellsNumber, char* iName)
{
  cellsNumber =  iCellsNumber;
  arr = new double[iCellsNumber];
  arrPrev = new double[iCellsNumber];
  int len=strlen(iName);
  name = new char[len+1];
  name[len]=0;
  strcpy(name, iName);
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

    T=400; //K
    p=101325;//pa
    mAr=39.948/1000.0;//kg/mol
    rho=p*mAr/(8.314*T);//kg/m^3
    double Na=6.022e23; //1/mol
    N=p*Na/(8.314*T);
}

void simulationData::updateParams()
{
    simulationData::simulationParameters* pParams = m_params;
    simulationData::simulationField* pEn = m_fieldEnergy;
    simulationData::simulationField* pPhi = m_fieldPhi;
    double dz = m_dz;
    for (int j=0; j<pParams->cellsNumber; j++)
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


simulationData::simulationParameters::simulationParameters(int iCellsNumber)
{
    init(iCellsNumber);
}
