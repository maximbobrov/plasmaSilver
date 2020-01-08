#include "simulationsolver.h"
#include "simulationtools.h"
#include "math.h"
#include <QDebug>


simulationSolver::simulationSolver(simulationData* ipData/*= nullptr*/)
{
    m_pData = ipData;
    m_field = nullptr;
    m_aNu = nullptr;
    m_aMu = nullptr;
    m_aRHS = nullptr;//new double[m_field ->cellsNumber];
    m_NewtonRHS = nullptr; //new double[ipData->getFieldNe()->cellsNumber];
    m_charge=-1;
}

simulationSolver::~simulationSolver()
{

}

double simulationSolver::calcE(int j)
{
    double * phi  = m_pData->getFieldPhi()->arr;
    //   double * e  = m_pData->m_params->arrE;


    //dt=1e-20;
    double eu=phi[j+1] - phi[j];//phi[j+2] - phi[j+1];               //j+1
    //double ed=phi[j] - phi[j-1];
}

double simulationSolver::solve(int iNumberIteration)
{
    double res = 0.0;
    if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dt = m_pData->getDt();
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();
            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = (dt * m_aNu[j] * (m_field->arr[j+1] + m_field->arr[j-1]) / (dz * dz) +  m_aRHS[j] * dt + m_field->arrPrev[j]) / (1.0 + 2.0 * m_aNu[j] * dt / (dz * dz));
            }
            setBc();

        }
        for (int j=1; j < m_field->cellsNumber - 1; j++)
        {
            res += (m_field->arr[j] - m_field->arrPrev[j]) / dt - (m_aNu[j] * (m_field->arr[j+1] - 2.0 * m_field->arr[j] + m_field->arr[j-1]) / (dz * dz) + m_aRHS[j]);
        }
        return res;
    }
    return -1;
}


double simulationSolver::getRhs()
{
    return -1;
}

double simulationSolver::getRhsAt(int j)
{
    return 0;
}

void simulationSolver::getStepEuler()
{
    for (int i = 0; i < m_field->cellsNumber; ++i)
    {
        m_field->arrPrev[i] = m_field->arr[i];
    }
}

void simulationSolver::setBc()
{
    //periodic bc initally
    double tmp=m_field->arr[1];
    int last=m_field->cellsNumber - 1;
    m_field->arr[0]=m_field->arr[m_field->cellsNumber - 2];
    m_field->arr[last]=tmp;

}

double simulationSolver::getNewtonRhs(int j)
{

    /*
    if((j==1) || (j==m_field->cellsNumber - 2))
        setBc();
    //x/dt_x - (x0/dt+nu*(xr+xl)/(dz*dz) + rhs) = x/dt_x -rhs_x
        double rhs=getRhsAt(j);
        double dz = m_pData->getDz();
        double dt =m_pData->getDt();
        return   (m_field->arrPrev[j])/dt + m_aNu[j]*(m_field->arr[j+1] + m_field->arr[j-1])/(dz*dz) + rhs;
      */

    double dz = m_pData->getDz();
    double dt = m_pData->getDt();
    double * phi  = m_pData->getFieldPhi()->arr;
    //   double * e  = m_pData->m_params->arrE;

    int last =m_field->cellsNumber - 2;
    /*    double Grij = m_mask[i+1][j]*(-m_aNu[i][j] - m_charge * (phi[i+1][j] - phi[i][j]) * m_aMu[i][j]) / dx;
    double Glij = m_mask[i-1][j]*( m_aNu[i][j] - m_charge * (phi[i][j] - phi[i-1][j]) * m_aMu[i][j]) / dx;
    double Guij = m_mask[i][j+1]*(-m_aNu[i][j] - m_charge * (phi[i][j+1] - phi[i][j]) * m_aMu[i][j]) / dy;
    double Gdij = m_mask[i][j-1]*( m_aNu[i][j] - m_charge * (phi[i][j] - phi[i][j-1]) * m_aMu[i][j]) / dy;
    double Gr = m_mask[i+1][j]*(m_field->arr[i+1][j] * m_aNu[i][j]) / dx;
    double Gl = - m_mask[i-1][j]*(m_field->arr[i-1][j] * m_aNu[i][j]) / dx;
    double Gu = m_mask[i][j+1]*(m_field->arr[i][j+1] * m_aNu[i][j]) / dy;
    double Gd = - m_mask[i][j-1]*(m_field->arr[i][j-1] * m_aNu[i][j]) / dy;
*/

    //dt=1e-20;
    double eu=/*calcE(j+1);*/phi[j+2] - phi[j+1];//phi[j+2] - phi[j+1];               //j+1
    double ed=/*calcE(j-1);*/phi[j+1] - phi[j];
    double au= m_charge * eu *m_aMu[j] ;//(phi[j+1] - phi[j]) * m_aMu[j];
    double ad= m_charge * ed *m_aMu[j];//(phi[j] - phi[j-1]) * m_aMu[j];

    double Ga_u = (j<last)*(-m_aNu[j] /*+ (au>0)*au*/) / dz;
    double Ga_d = (j>1)*( m_aNu[j] /*+0.5*au*/  /*+ad<0)*ad*/) / dz;
    double Gu = (j<last)*(m_field->arr[j+1] * m_aNu[j]   +/*(au<0)*/2.0*(m_field->arrPrev[j+1]  )* au) / dz;
    double Gd = (j>1)*(-m_field->arr[j-1] * m_aNu[j]  +/*(ad>0)*/2.0*(m_field->arrPrev[j]) * ad) / dz;

    //

    double rhs=getRhsAt(j);

    double a = 1.0/dt -  (Ga_u - Ga_d) / dz;

    // qDebug()<<"au"<<m_aNu[j]<<" ad"<<( (Gu - Gd) / dz)/a;
    //return ((rhs-(bp*m_field->arr[i+1][j]+bm*m_field->arr[i-1][j]+cp*m_field->arr[i][j+1]+cm*m_field->arr[i][j-1]))
    //       /a)*m_mask[i][j]+m_maskValue[i][j];

    //  double res0=(m_field->arrPrev[j]/dt+m_aNu[j]*(m_field->arr[j+1] + m_field->arr[j-1])/(dz*dz))/(1.0/dt+ 2.0*m_aNu[j]/(dz*dz));//-m_field->arrPrev[j];
    double res=((rhs + m_field->arrPrev[j] / dt + ( (Gu - Gd) / dz)) / a);
    //qDebug()<<"res1"<<res/m_field->arrPrev[j];


    return res;//((rhs + m_field->arrPrev[j] / dt + ( (Gu - Gd) / dz)) / a);



    // return m_field->arrPrev[j]/dt;
}

solverNe::~solverNe()
{

}

solverNe::solverNe(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldNe();
    m_aNu = m_pData->getParameters()->arrDe;
    m_aRHS = new double[m_field ->cellsNumber];
    m_NewtonRHS = new double[m_field ->cellsNumber];
    m_aMu = m_pData->getParameters()->arrMue;
    m_charge = - 1;
}

double solverNe::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();


    double dz = m_pData->getDz();
    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i]=
                pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDe, m_field ->cellsNumber, dz, i);
    }
    return 1;
}

double solverNe::getRhsAt(int i)
{
    // simulationData::simulationParameters* pParams = m_pData->getParameters();


    //  double dz = m_pData->getDz();

    //return 0.0;
    //pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
    // + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddz(pParams->arrE/*, m_field ->cellsNumber*/, dz, i);




    double ne=m_pData->getFieldNe()->arrPrev[i];
    double nars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev[i];
    double Te=m_pData->getParameters()->arrTe[i];
    static double n_n=m_pData->getParameters()->N;

    static double k[7];


    //for (int i=0;i<7;i++)
    {
        k[3]=m_pData->m_reactions[3]->getRate(fmax(fmin(Te,40),0.0001));
        k[4]=m_pData->m_reactions[4]->getRate(fmax(fmin(Te,40),0.0001));
        k[5]=m_pData->m_reactions[5]->getRate(fmax(fmin(Te,40),0.0001));

        //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
    }

    return ne*(k[3]*n_n + k[4]*nars) + k[5]*nars*nars;


}


void solverNe::setBc()
{
    int last=m_field->cellsNumber - 1;

    double te=m_pData->getParameters()->arrTe[1];
    double dz=m_pData->getDz();
    double E0=(m_pData->getFieldPhi()->arr[1]-m_pData->getFieldPhi()->arr[0])/dz;
    double mu0=m_pData->getParameters()->arrMue[1];
    double D0=m_pData->getParameters()->arrDe[1];
    double gam_gam=0.0;//-1.1e19;
    double nuen=6.69e5*sqrt(te);

    m_field->arr[0]= 1e5;//(-gam_gam + m_field->arr[1]*D0/dz)/(mu0*E0 + D0/dz +nuen);
    m_field->arr[last]=m_field->arr[last-1];//1e5;
    /*
    int last=m_field->cellsNumber - 1;
    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double a = pParams->arrMue[last-1] * (pParams->arrE[last-1] + pParams->arrE[last])/2;
    m_field->arr[last]=m_field->arr[last-1] * (pParams->arrDe[last-1] / m_pData->getDz() - a) / (pParams->arrDe[last-1] / m_pData->getDz() + a);

    a = pParams->arrMue[1] * (pParams->arrE[1] + pParams->arrE[0])/2;
        m_field->arr[0] =m_field->arr[1] * (pParams->arrDe[1] / m_pData->getDz() + a) / (pParams->arrDe[1] / m_pData->getDz() - a);*/
}

solverEnergy::~solverEnergy()
{

}

solverEnergy::solverEnergy(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldEnergy();
    m_aNu = m_pData->getParameters()->arrDeps;
    m_aRHS = new double[m_field ->cellsNumber];
    m_NewtonRHS = new double[m_field ->cellsNumber];
    m_aMu = m_pData->getParameters()->arrMueps;
    m_charge = - 1;
}

double solverEnergy::getRhs()
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField* pNe = m_pData->getFieldNe();
    // double* pAr_star = m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    double dz = m_pData->getDz();

    for (int i = 0; i < m_field ->cellsNumber-1; ++i)
    {
        double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*simulationTools::ddzCentral(pNe->arr, m_field ->cellsNumber, dz, i);
        //(pNe->arr[i+1]-pNe->arr[i]);

        m_aRHS[i]=
                pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i)
                + simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i);
        - electronFlux * pParams->arrE[i];


        /*  m_aRHS[i]=
                    pParams->arrMue[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                    + pParams->arrMue[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);*/

    }
    return 1;
}

double solverEnergy::getRhsAt(int i)
{
    simulationData::simulationParameters* pParams = m_pData->getParameters();
    simulationData::simulationField* pNe = m_pData->getFieldNe();
    // double* pAr_star = m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    double dz = m_pData->getDz();

    double electronFlux =  - pParams->arrMue[i] * pParams->arrE[i] * pNe->arr[i] - pParams->arrDe[i]*simulationTools::ddzCentral(pNe->arr, m_field ->cellsNumber, dz, i);

    // return  0.0 // //  pParams->arrMueps[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
    // + pParams->arrMueps[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);
    //+ simulationTools::ddzCentral(m_field->arr , m_field ->cellsNumber, dz, i) * simulationTools::ddzCentral(pParams->arrDeps, m_field ->cellsNumber, dz, i);
    //+ 1.0*electronFlux * pParams->arrE[i];




    double ne=m_pData->getFieldNe()->arrPrev[i];
    double nars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev[i];
    double Te=m_pData->getParameters()->arrTe[i];
    static double n_n=m_pData->getParameters()->N;

    static double k[7],dE[7];


    for (int i=1;i<6;i++)
    {
        k[i]=m_pData->m_reactions[i]->getRate(fmax(fmin(Te,40),0.0001));
        dE[i]=m_pData->m_reactions[i]->getDe();
    }

    return  1.0*electronFlux * pParams->arrE[i] +
            ne*(n_n*(dE[1]*k[1] + dE[3]*k[3])
            +nars*(dE[2]*k[2] + dE[4]*k[4])) + dE[5]*k[5]*nars*nars;
    /*  x*x0*(
       n_n *( dE[1]*k[1] +  dE[3]*k[3]) +
       y*y0*( dE[2]*k[2] +  dE[4]*k[4])
       ) +
        (dE[5]) * k[5] * y*y0 * nars_i  ;*/
}

void solverEnergy::setBc()
{

    int last=m_field->cellsNumber - 1;
    //m_field->arr[0]=m_field->arr[1];//1e5;
    //m_field->arr[last]=m_field->arr[last-1];//1e5;

    double te=m_pData->getParameters()->arrTe[1];
    double dz=m_pData->getDz();
    double E0=(m_pData->getFieldPhi()->arr[1]-m_pData->getFieldPhi()->arr[0])/dz;
    double mu0=m_pData->getParameters()->arrMueps[1];
    double D0=m_pData->getParameters()->arrDeps[1];
    double gam_gam=0.0;//-4e20;
    double nuen=6.69e5*sqrt(te);

    m_field->arr[0]=5e5;// (-gam_gam + m_field->arr[1]*D0/dz)/(mu0*E0 + D0/dz +nuen);
    m_field->arr[last]=m_field->arr[last-1];//1e5;


    /*
    int last=m_field->cellsNumber - 1;
    simulationData::simulationParameters* pParams = m_pData->getParameters();

    double a = pParams->arrMueps[last-1] * (pParams->arrE[last-1] + pParams->arrE[last])/2;
    m_field->arr[last]=m_field->arr[last-1] * (pParams->arrDeps[last-1] / m_pData->getDz() - a) / (pParams->arrDeps[last-1] / m_pData->getDz() + a);

    a = pParams->arrMueps[1] * (pParams->arrE[1] + pParams->arrE[0])/2;
        m_field->arr[0] =m_field->arr[1] * (pParams->arrDeps[1] / m_pData->getDz() + a) / (pParams->arrDeps[1] / m_pData->getDz() - a);*/
}

solverPhi::~solverPhi()
{

}

solverPhi::solverPhi(simulationData* pData)
{
    m_pData = pData;
    m_field = m_pData->getFieldPhi();
    m_aRHS = new double[m_field ->cellsNumber];
    m_NewtonRHS = new double[m_field ->cellsNumber];
}

double solverPhi::getRhs()
{
    double q_e= 1.6*10e-19;
    double eps_0=8.85*10e-12;
    double q_over_eps0=q_e/eps_0;

    double* pNe = m_pData->getFieldNe()->arr;
    double* pArp = m_pData->getFieldHeavySpicies(0)->arr;

    simulationData::simulationParameters* pParams=m_pData->getParameters();
    // double mult= (pParams->p*6.022e23)/(pParams->T*8.314);
    //m_fHeavy[j]->arr[i] =(m_fNe->arr[i]*pParams->T*8.314)/(pParams->p*6.022e23);

    for (int i = 1; i < m_field ->cellsNumber-1; ++i)
    {
        m_aRHS[i] = -q_over_eps0*0.5*((pNe[i]+pNe[i-1])-(pArp[i]+pArp[i-1]));//q_over_eps0*0.5*(-(pNe[i]+pNe[i+1])+(pArp[i]+pArp[i+1]));
    }
    return 1;
}
/*
double solverPhi::solve(int iNumberIteration)
{
    double res = 0.0;
    if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dzdz=dz*dz;
        for (int i = 0; i < iNumberIteration; ++i)
        {
            getRhs();

            for (int j = 1; j <  m_field->cellsNumber - 1; ++j)
            {
                m_field->arr[j] = 0.5* ((m_field->arr[j+1] + m_field->arr[j-1]) +  m_aRHS[j]*dzdz);
            }
            setBc();
        }
        return res;
    }
    return -1;
}
*/

double solverPhi::solve(int iNumberIteration)
{
    double res = 0.0;

    static double rhs_[20][NZ+1];
    static double field_[20][NZ+1];



    if(m_field != nullptr)
    {
        double dz = m_pData->getDz();
        double dzdz=dz*dz;

        double a=2.0/dzdz;
        double b=-1.0/dzdz;
        getRhs();

        setBc();


        int nz0=NZ+1;
        for (int j=0;j<nz0;j++)
        {
            field_[0][j]=m_field->arr[j] ;
            rhs_[0][j]=m_aRHS[j];

            //  qDebug()<<j<<" "<<field_[0][j]<<m_field->cellsNumber - 1;
        }

        int nz=nz0;
        int N=6;
        int N_intern =2;
        int nn;
        for (int i = 0; i < iNumberIteration; ++i)
        {
            for (nn=0;nn<N;nn++)
            {
                nz=(nz0-1)/pow(2,nn);
                for (int aa=0;aa<N_intern;aa++)
                {
                    for (int j = 1; j < nz; ++j)
                    {
                        field_[nn][j] = (-b*(field_[nn][j+1] + field_[nn][j-1]) +  rhs_[nn][j])/a;
                    }
                }
                for (int j=1;j<(nz0-1)/pow(2,nn+1);j++)
                {
                    int l=j*2;
                    rhs_[nn+1][j]=(rhs_[nn][l] - (a*field_[nn][l]+b*(field_[nn][l+1]+field_[nn][l-1])))*4.0;
                    field_[nn+1][j]=0.0;
                }
            }
            for (int aa=0;aa<N_intern*3;aa++)
            {
                nz=(nz0-1)/pow(2,N);
                for (int j = 1; j < nz; ++j)
                {
                    field_[N][j] = (-b*(field_[N][j+1] + field_[N][j-1]) +  rhs_[N][j])/a;
                }
            }
            for (nn=N-1;nn>=0;nn--)
            {
                nz=(nz0-1)/pow(2,nn);
                int j,l;
                for (l=1; l<(nz)/2; l++)
                {
                    j=l*2;
                    field_[nn][j]+=field_[nn+1][l];
                    field_[nn][j-1]+=0.5*(field_[nn+1][l]+field_[nn+1][l-1]);
                }
                l=(nz)/2; j=l*2;
                field_[nn][j-1]+=0.5*(field_[nn+1][l]+field_[nn+1][l-1]);

                for (int aa=0;aa<N_intern;aa++)
                {
                    for (int j = 1; j <  nz; ++j)
                    {
                        field_[nn][j] = (-b*(field_[nn][j+1] + field_[nn][j-1]) +  rhs_[nn][j])/a;
                    }
                }
            }
        }

        for (int j=1;j<nz0-1;j++)
        {
            m_field->arr[j]=field_[0][j] ;

        }

    }
    return -1;
}

/*
double multigrid_N(INPUT_PARAM par, double** field, double** rhs, double** mask, double** maskValue, int itn,int N)
{
    int i,j,l,k,nn;
    double a,b_p,b_m,c_p,c_m;
    double res=0;
    a=par.a;//((2.0)/(dx*dx)+2.0/(dy*dy));
    c_p=par.cp;//-1.0/(dy*dy);
    c_m=par.cm;//-1.0/(dy*dy);
    static double rhs_[10][N_Z];
    static double field_[10][N_Z];

    for (j=0;j<N_Z;j++)
    {
        field_[0][j]=field[j];
        rhs_[0][j]=rhs[j];
    }


    for (nn=0;nn<N;nn++)
    {
        jacobi_N(par,field_[nn],rhs_[nn],mask_[nn], maskValue_[nn],itn,(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));

        int nz=(N_Z-1)/pow(2,nn);


        for (j=1;j<(N_Z-1)/pow(2,nn+1);j++)
        {
            l=j*2;
            rhs_[nn+1][j]=rhs_[nn][l] - (a*field_[nn][l]+c_p*field_[nn][l+1]+c_m*field_[nn][l-1]);
            field_[nn+1][j]=0.0;
        }
    }

    jacobi_N(par,field_[N],rhs_[N],mask_[N], maskValue_[N],itn,(N_X-1)/pow(2,N),(N_Y-1)/pow(2,N));

    for (nn=N-1;nn>=0;nn--)
    {
        interp_up(field_[nn],field_[nn+1],(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));

        int j,l;
        for (l=1; l<(ny)/2; l++)
        {
            j=l*2;
            field[j]+=field2[l];
            field[j-1]+=0.5*(field2[l]+field2[l-1]);
        }
        l=(ny)/2; j=l*2;
        field[j-1]+=0.5*(field2[l]+field2[l-1]);

        jacobi_N(par,field_[nn],rhs_[nn],mask_[nn], maskValue_[nn],itn,(N_X-1)/pow(2,nn),(N_Y-1)/pow(2,nn));
    }

    for (j=0;j<N_Y;j++)
    {
        field[j]=field_[0][j];
    }


}

double jacobi_N(INPUT_PARAM par, double field[N_X][N_Y], double rhs[N_X][N_Y],double mask[N_X][N_Y],double maskValue[N_X][N_Y], int itn,int nx,int ny)
{
    int i,j,n;
    double eps =1.0;
    double a,c_p,c_m;
    a=par.a;
    c_p=par.cp;
    c_m=par.cm;

    for(n=0;n<itn;n++)
    {
        for (j=1; j<ny; j++)
        {
            if (j==1)
            {
                if (par.s_bc_type==0)//fixed_value
                    field[i][j-1]=par.s_bc_val;
                if (par.s_bc_type==1)//fixed_gradient
                    field[i][j-1]=field[i][1]-par.dy*par.s_bc_val;
                if (par.s_bc_type==2)//cyclic
                    field[i][j-1]=field[i][ny-1];
            }

            if (j==ny-1)
            {
                if (par.n_bc_type==0)//fixed_value
                    field[i][j+1]=par.n_bc_val;
                if (par.n_bc_type==1)//fixed_gradient
                    field[i][j+1]=field[i][j]+par.dy*par.n_bc_val;
                if (par.n_bc_type==2)//cyclic
                    field[i][j+1]=field[i][1];
            }


            field[i][j]=((rhs[i][j]-(b_p*field[i+1][j]+b_m*field[i-1][j]+c_p*field[i][j+1]+c_m*field[i][j-1]))
                    /a);

        }
    }
}
return 0;
}
double interp_up(double field[N_X][N_Y], double field2[N_X][N_Y], int nx,int ny)
{
    int j,l;
    for (l=1; l<(ny)/2; l++)
    {
        j=l*2;
        field[j]+=field2[l];
        field[j-1]+=0.5*(field2[l]+field2[l-1]);
    }
    l=(ny)/2; j=l*2;
    field[j-1]+=0.5*(field2[l]+field2[l-1]);
}
*/



void solverPhi::setBc()
{
    int last=m_field->cellsNumber - 1;
    m_field->arr[0] = 0;
    m_field->arr[NZ]=100;
}


solverHeavySpicies::solverHeavySpicies(simulationData *pData, int num)
{
    m_pData = pData;
    m_field = m_pData->getFieldHeavySpicies(num);
    m_specie=m_field->m_specie;
    m_charge = m_pData->getHeavySpiciesCharge(num);
    m_aNu = m_pData->getParameters()->arrDomega;
    m_aRHS = new double[m_field ->cellsNumber];
    m_NewtonRHS = new double[m_field ->cellsNumber];
    m_aMu = m_pData->getParameters()->arrMuomega;
}

solverHeavySpicies::~solverHeavySpicies()
{

}

double solverHeavySpicies::getRhs()
{

    simulationData::simulationField* pNe = m_pData->getFieldNe();

    double* pAr_star = m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arr;

    simulationData::simulationParameters* pParams = m_pData->getParameters();
    double dz = m_pData->getDz();




    if (m_specie==simulationData::SpecieName::Ar_plus)
    {


        for (int i = 0; i < m_field ->cellsNumber-1; ++i)
        {

            m_aRHS[i]=0.0;/*
                    -  m_charge  * pParams->arrMuomega[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                    -  m_charge  * pParams->arrMuomega[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);*/
        }


    }else if (m_specie==simulationData::SpecieName::Ar_star)
    {



        for (int i = 0; i < m_field ->cellsNumber-1; ++i)
        {

            m_aRHS[i]=0.0;/*
                    -  m_charge  * pParams->arrMuomega[i] * pParams->arrE[i] * simulationTools::ddzCentral(m_field->arr, m_field ->cellsNumber, dz, i)
                    -  m_charge  * pParams->arrMuomega[i] * m_field->arr[i] * simulationTools::ddzCentral(pParams->arrE, m_field ->cellsNumber, dz, i);*/
        }
    }

    /*if (m_specie==simulationData::SpecieName::Ar_star)
    {
    //m_pData->calcReaction(simulationData::ReactionName::comsol_eAr_eArs);
    m_pData->calcReaction(simulationData::ReactionName::comsol_eArs_2eArp);
    R1=m_pData->getReactionRate(simulationData::ReactionName::comsol_eArs_2eArp);
    mult1=pParams->p / ( pParams->T * 8.314 * 6.022e23 * pParams->rho);
    }*/



    return 1;
}

double solverHeavySpicies::getRhsAt(int i)
{
    //return 0.0;



    double ne=m_pData->getFieldNe()->arrPrev[i];
    double nars=m_pData->getFieldHeavySpicies(simulationData::SpecieName::Ar_star)->arrPrev[i];
    double Te=m_pData->getParameters()->arrTe[i];
    static double n_n=m_pData->getParameters()->N;

    static double k[7];

    if (m_specie==simulationData::SpecieName::Ar_plus)
    {
        //for (int i=0;i<7;i++)
        {
            k[3]=m_pData->m_reactions[3]->getRate(fmax(fmin(Te,40),0.0001));
            k[4]=m_pData->m_reactions[4]->getRate(fmax(fmin(Te,40),0.0001));
            k[5]=m_pData->m_reactions[5]->getRate(fmax(fmin(Te,40),0.0001));

            //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }

        return ne*(k[3]*n_n + k[4]*nars) + k[5]*nars*nars;

    }else if (m_specie==simulationData::SpecieName::Ar_star)
    {

        //        ne*(k[1]*n_n -(k[2]+k[4])*nars) -2.0*k[5]*nars*nars;

        //for (int i=0;i<7;i++)
        {
            k[1]=m_pData->m_reactions[1]->getRate(fmax(fmin(Te,40),0.0001));
            k[2]=m_pData->m_reactions[2]->getRate(fmax(fmin(Te,40),0.0001));
            k[4]=m_pData->m_reactions[4]->getRate(fmax(fmin(Te,40),0.0001));
            k[5]=m_pData->m_reactions[5]->getRate(fmax(fmin(Te,40),0.0001));

            //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }

        return ne*(k[1]*n_n -(k[2]+k[4])*nars) -2.0*k[5]*nars*nars;;
    }
    return 0;
}

void solverHeavySpicies::setBc()
{
    int last=m_field->cellsNumber - 1;
    m_field->arr[0] =m_field->arr[1];//m_field->arr[1]*0.5;//  m_field->arr[1];
    m_field->arr[last]=m_field->arr[last-1] ;
}


