#include "reaction.h"


static const int g_EAr_EAr_num =67;

static  double g_EAr_EAr[g_EAr_EAr_num][2] = {
    #include "reactions/e+Ar_e+Ar.txt"
};

static const int g_EAr_EArs_num =34;

static double g_EAr_EArs[g_EAr_EArs_num][2] = {
    #include "reactions/e+Ar_e+Ars.txt"
};

static const int g_EAr_2EArp_num =31;

static double g_EAr_2EArp[g_EAr_2EArp_num][2] = {
    #include "reactions/e+Ar_2e+Ar+.txt"
};

static const int g_EArs_2EArp_num =20;

static double g_EArs_2EArp[g_EArs_2EArp_num][2] = {
    #include "reactions/e+Ars_2e+Ar+.txt"
};





reaction::reaction(simulationData* data)
{

    if (data != nullptr)
    {
        int cn=data->getCellsNumber();
        m_R=new double[cn];
        for (int i = 0 ; i < cn ; i++)
        {
            m_R[i] = 0;
        }
    }
}

void reaction::calc()
{
//just a placeholder here
}

reactionEAr_EAr::reactionEAr_EAr(simulationData *data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(30,g_EAr_EAr_num,30);
        m_cs->fillSigmas2(g_EAr_EAr,g_EAr_EAr_num); //fill from static arrays (safer then messing up with external txt files)
    }
}

void reactionEAr_EAr::calc()
{

}

reactionEAr_EArs::reactionEAr_EArs(simulationData *data)
{
    if (data != nullptr)
    {
       m_cs=new crossSection(30,g_EAr_EArs_num,30);
        m_cs->fillSigmas2(g_EAr_EArs,g_EAr_EArs_num);
    }
}

void reactionEAr_EArs::calc()
{

}

reactionEAr_2EArp::reactionEAr_2EArp(simulationData *data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(30,g_EAr_2EArp_num,30);
       m_cs->fillSigmas2(g_EAr_2EArp,g_EAr_2EArp_num);
    }
}

void reactionEAr_2EArp::calc()
{

}

reactionEArs_2EArp::reactionEArs_2EArp(simulationData *data)
{
    if (data != nullptr)
    {
        m_cs=new crossSection(30,g_EArs_2EArp_num,30);
        m_cs->fillSigmas2(g_EArs_2EArp,g_EArs_2EArp_num);
    }
}

void reactionEArs_2EArp::calc()
{

}
