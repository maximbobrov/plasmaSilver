#include "reactionsolver.h"
#include "simulationdata.h"
#include <qdebug.h>

reactionSolver::reactionSolver(simulationData* mdata)
{

    for (int i=0;i<mdata->m_reactions.size();i++)
        m_reactions.push_back(mdata->m_reactions[i]);
    qDebug()<<m_reactions.size()<< "AA!!";

}

void reactionSolver::solve(int itn)
{
    oValues[0]=ne_i;
    oValues[1]=eps_i;
    oValues[2]=nars_i;
    double k[7];
    double dk[7];
    double A[7]; //energy gain/loss of reaction

    //for (int nn=0;nn<itn;nn++)
    {
        for (int i=0;i<7;i++)
        {
            /*oValues[0]=ne_i;
            oValues[1]=eps_i;
            oValues[2]=nars_i;*/

           /* if (oValues[0]<1e5) oValues[0]=1e5;
            if (oValues[1]<0.01) oValues[1]=0.01;
            if (oValues[1]>40.0) oValues[1]=40.0;

            if (oValues[2]<1e5) oValues[2]=1e5;*/


            k[i]=m_reactions[i]->getRate(oValues[1]);
            dk[i]=m_reactions[i]->getDeriv(oValues[1]);
            A[i]=m_reactions[i]->getDe();
        }
        //ne:
        /*f[0]=-(oValues[0]-ne_i)/dt+rhs_ne + oValues[0]*(k[3]*n_n +k[4]*oValues[2])+k[5]*oValues[2]*oValues[2];

        df[0][0]=-(1)/dt+ (k[3]*n_n +k[4]*oValues[2]);
        df[0][1]= oValues[0]*(dk[3]*n_n +dk[4]*oValues[2])+dk[5]*oValues[2]*oValues[2];
        df[0][2]= oValues[0]*(k[4])+2.0*k[5]*oValues[2];

        //eps
        f[1]=-(oValues[1]-eps_i)/dt+(rhs_eps-rhs_ne)/oValues[0] +
                n_n*(A[0]*k[0]+A[1]*k[1]+(A[3]-oValues[1])*k[3]) +
                oValues[2]*(A[2]*k[2]+(A[4]-oValues[1])*k[4]) +
                ((A[5]-oValues[1])*k[5]*oValues[2]*oValues[2] + A[6]*k[6]*n_n*oValues[2])/oValues[0];

        df[1][0]=-(rhs_eps-rhs_ne)/(oValues[0]*oValues[0]) -
                ((A[5]-oValues[1])*k[5]*oValues[2]*oValues[2] + A[6]*k[6]*n_n*oValues[2]*oValues[2])/(oValues[0]*oValues[0]);

        df[1][1]=-(1)/dt+
                n_n*(A[0]*dk[0]+A[1]*dk[1]+A[3]*dk[3]-oValues[1]*dk[3]-k[3]) +
                oValues[2]*(A[2]*dk[2]+A[4]*dk[4]-oValues[1]*dk[4]-k[4]) +
                ( (A[5]*k[5]-oValues[1]*dk[5]-k[5])*oValues[2]*oValues[2] + A[6]*dk[6]*n_n*oValues[2])/oValues[0];

        df[1][2]=(A[2]*k[2]+(A[4]-oValues[1])*k[4]) +
                (2.0*(A[5]-oValues[1])*k[5]*oValues[2] + A[6]*k[6]*n_n)/oValues[0];

        //n_ars
        f[2]=-(oValues[2]-nars_i)/dt+rhs_nars + oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2])
                -2.0*k[5]*oValues[2]*oValues[2]-k[6]*n_n*oValues[2];

        df[2][0]= (k[1]*n_n -(k[2]+k[4])*oValues[2]);
        df[2][1]=oValues[0]*(dk[1]*n_n -(dk[2]+dk[4])*oValues[2])
                -2.0*dk[5]*oValues[2]*oValues[2]-dk[6]*n_n*oValues[2];

        df[2][2]=-(1)/dt+ oValues[0]*( -(k[2]+k[4]))
                -4.0*k[5]*oValues[2]-k[6]*n_n;

        solve3x3(df,f,dValues);

        oValues[0]-=dValues[0];
        oValues[1]-=dValues[1];
        oValues[2]-=dValues[2];*/


        f[2]=-(oValues[2]-nars_i)/dt+rhs_nars + oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2])
                -2.0*k[5]*oValues[2]*oValues[2]-k[6]*n_n*oValues[2];


        df[2][2]=-(1)/dt+ oValues[0]*( -(k[2]+k[4]))
                -4.0*k[5]*oValues[2]-k[6]*n_n;

        //solve3x3(df,f,dValues);

        //oValues[0]-=dValues[0];
        //oValues[1]-=dValues[1];
        //oValues[2]-=f[2]/df[2][2];

        oValues[2]=nars_i+dt*(rhs_nars + oValues[0]*(k[1]*n_n) -(k[2]+k[4])*oValues[2]
               /* -2.0*k[5]*oValues[2]*oValues[2]-k[6]*n_n*oValues[2]*/);
    }


    ne_o=oValues[0];
    eps_o=oValues[1];
    nars_o=oValues[2];

}
