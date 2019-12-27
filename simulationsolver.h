#ifndef SIMULATIONSOLVER_H
#define SIMULATIONSOLVER_H
#include "simulationdata.h"


class simulationSolver
{
public:

    simulationSolver(simulationData* pData = nullptr);


    virtual double solve(int numberIteration);
    virtual double getRhs();
    virtual void getStepEuler();
    virtual void setBc();

    virtual ~simulationSolver();

protected:
    simulationData* m_pData;
    simulationData::simulationField* m_field;
    double* m_aRHS;
    double* m_aNu;
};

class solverNe : public simulationSolver
{
public:
    solverNe(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual void setBc();
    virtual ~solverNe();
};


class solverEnergy : public simulationSolver
{
public:
    solverEnergy(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual void setBc();
    virtual ~solverEnergy();
};

class solverPhi : public simulationSolver
{
public:
    solverPhi(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double solve(int numberIteration);
    virtual void setBc();
    virtual ~solverPhi();
};

class solverHeavySpicies : public simulationSolver
{
public:
    solverHeavySpicies(simulationData* pData = nullptr, int num = 0);
    virtual double getRhs();

    virtual void setBc();
    virtual ~solverHeavySpicies();
private:
    int m_charge;
    simulationData::SpecieName m_specie;
};

#endif // SIMULATIONSOLVER_H
