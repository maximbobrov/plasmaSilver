#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H
#include <iostream>
#include <vector>
#include "reaction.h"

class simulationData
{
public:

    simulationData(int cellsNumber);
    ~simulationData();


    enum SpecieName
    {
       Ar,
       Ar_star,
       Ar_plus
    };

    enum ReactionName
    {
       eAr_eAr,
       eAr_eArs,
       eAr_2eArp,
       eArs_2eArp
    };

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
    double getN(); //total number of particles in one m^3
    double* getReactionRate(int num);
    void calcReaction(int num);
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
    std::vector<reaction *> m_reactions;
};

#endif // SIMULATIONDATA_H
