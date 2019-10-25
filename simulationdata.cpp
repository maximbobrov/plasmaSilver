#include "simulationdata.h"
#include "string.h"


simulationData::simulationData(int iCellsNumber)
{
    m_defaultCellsNumber = iCellsNumber;
    m_fieldNe = new simulationField(iCellsNumber,"electrons");
    m_fieldEnergy = new simulationField(iCellsNumber,"energy");
    m_fieldPhi = new simulationField(iCellsNumber+1,"potential");
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar"));
    m_chargeHeavySpecies.push_back(0);
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ars"));
    m_chargeHeavySpecies.push_back(1);
    m_fieldsHeavySpecies.push_back(new simulationField(iCellsNumber,"Ar+"));
    m_chargeHeavySpecies.push_back(1);
    m_numberHeavySpicies = m_fieldsHeavySpecies.size();
    m_params = new simulationData::simulationParameters(iCellsNumber);
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
}

simulationData::simulationParameters::simulationParameters(int iCellsNumber)
{
    init(iCellsNumber);
}
