#ifndef CROSSSECTION_H
#define CROSSSECTION_H


class crossSection
{

public:
    crossSection();
    crossSection(int kN,int sigN,int splN);

    void init(int kN,int sigN,int splN);
    void fillSigmas(double* sNew, double*eNew, int n); //copy from 2 1D arrays
    void fillSigmas2(double esNew[][2], int n); //copy from 2D array
    double getSpline(double xp);
    void readFromFile(char* fileName);

    ~crossSection();

private:
    void fillK();
    void buildSpline();
    double integrate(double sm,double sp,double em,double ep, double phi);
    double getKICoarse(double phi);
    double getKIFine(double phi);


private:
    int m_kNum;
    int m_sigNum;
    int m_splNum;
    double* m_sigmas;
    double* m_energs;
    double* m_K;
    double* m_phis;
    double* m_splX, * m_splA, * m_splB, * m_splC, * m_splD; //spl_num+1

};

#endif // CROSSSECTION_H
