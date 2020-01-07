#ifndef SIMULATIONDATA_H
#define SIMULATIONDATA_H
#include <iostream>
#include <vector>
#include "reaction.h"


#define NZ 2048

class simulationData
{
public:

    simulationData(int cellsNumber);
    ~simulationData();


    enum SpecieName
    {

        Ar_plus,
        Ar_star,

        e,
        En,
        phi,
        Ar
    };

    enum ReactionName
    {
        comsol_eAr_eAr,
        comsol_eAr_eArs,
        comsol_eArs_eAr,
        comsol_eAr_2eArp,
        comsol_eArs_2eArp,
        comsol_ArsArs_EArArp,
        comsol_ArsAr_ArAr
        //comsol_eAr_2eArp
    };

    struct simulationField
    {
        char* name;
        double* arr;
        double* arrPrev;
        int cellsNumber;
        simulationData::SpecieName m_specie;
        simulationField(int cellsNumber, char* name,simulationData::SpecieName specie);
    private: void init(int cellsNumber, char* name,simulationData::SpecieName specie);
    };

    struct simulationParameters
    {
        double* arrDe;
        double* arrDeps;
        double* arrDomega;
        double* arrMue;
        double* arrMueps;
        double* arrMuomega;
        double* arrTe;
        double* arrE;
        int cellsNumber;
        double rho; //mixture density
        double p; //pressure
        double T; //temperature
        double mAr; //argon molar mass
        double N; //neutral number denisty
        simulationParameters(int cellsNumber);
    private: void init(int cellsNumber);
    };

    const double q=1.6022e-19; //coulumbs elementary charge
    const double k_B_const=1.3806e-23;// boltzmann constant

    void setDt(double dt = 0.000003);
    void setDz(double dz = 0.01);
    void setCellsNumber(int cellsNumber);

    int getCellsNumber();
    double getDt();
    double getDz();
    void updateParams();
    double* getArrTe();
    simulationField* getFieldNe();
    simulationField* getFieldEnergy();
    simulationField* getFieldPhi();
    simulationField* getFieldHeavySpicies(int num);
    int getHeavySpiciesCharge(int num);
    int getNumberHeavySpicies();
    double getN(); //total number of particles in one m^3
    double* getReactionRate(simulationData::ReactionName reactName);
    double getReactionDe(simulationData::ReactionName reactName);

    void calcReaction(ReactionName reactName);
    simulationParameters *getParameters();

private:
    friend class reactionSolver;
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
