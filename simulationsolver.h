#ifndef SIMULATIONSOLVER_H
#define SIMULATIONSOLVER_H
#include "simulationdata.h"


class simulationSolver
{
public:

    simulationSolver(simulationData* pData = nullptr);


    virtual double solve(int numberIteration);
    virtual double getRhs();
     virtual double getRhsAt(int j);
    virtual void getStepEuler();
    virtual void setBc();

    virtual double getNewtonRhs(int i);

    virtual ~simulationSolver();

protected:
    simulationData* m_pData;
    simulationData::simulationField* m_field;
    double* m_aRHS;
    double* m_NewtonRHS;
    double* m_aNu;
     double* m_aMu;

    int m_charge;
};

class solverNe : public simulationSolver
{
public:
    solverNe(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double getRhsAt(int j);
    virtual void setBc();
    virtual ~solverNe();
};


class solverEnergy : public simulationSolver
{
public:
    solverEnergy(simulationData* pData = nullptr);
    virtual double getRhs();
    virtual double getRhsAt(int j);
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
    virtual double getRhsAt(int j);
    virtual void setBc();
    virtual ~solverHeavySpicies();

    simulationData::SpecieName m_specie;
};

#endif // SIMULATIONSOLVER_H
