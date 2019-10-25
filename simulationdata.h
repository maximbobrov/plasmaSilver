#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H
#include <iostream>
#include <vector>


class simulationData
{
public:

    simulationData(int cellsNumber);
    ~simulationData();

    struct simulationField
    {
        char* name;
        double* arr;
        double* arrPrev;
        int cellsNumber;
        simulationField(int cellsNumber, char* name);
        private: void init(int cellsNumber, char* name);
    };

    struct simulationParameters
    {
        double* arrDe;
        double* arrDeps;
        double* arrDomega;
        double* arrMue;
        double* arrMueps;
        double* arrMuomega;
        double* arrE;
        int cellsNumber;
        simulationParameters(int cellsNumber);
        private: void init(int cellsNumber);
    };

    void setDt(double dt = 0.000003);
    void setDz(double dz = 0.01);
    void setCellsNumber(int cellsNumber);

    int getCellsNumber();
    double getDt();
    double getDz();
    simulationField* getFieldNe();
    simulationField* getFieldEnergy();
    simulationField* getFieldPhi();
    simulationField* getFieldHeavySpicies(int num);
    int getHeavySpiciesCharge(int num);
    int getNumberHeavySpicies();
    simulationParameters *getParameters();

private:
    double m_dt;
    double m_dz;
    int m_defaultCellsNumber;
    int m_numberHeavySpicies;
    simulationField* m_fieldNe;
    simulationField* m_fieldEnergy;
    simulationField* m_fieldPhi;
    std::vector<simulationField* > m_fieldsHeavySpecies;
    std::vector<int> m_chargeHeavySpecies;
    simulationParameters* m_params;
};

#endif // SIMULATIONDATA_H
