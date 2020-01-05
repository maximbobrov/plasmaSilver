#include "reactionsolver.h"
#include "simulationdata.h"
#include <qdebug.h>
#include <stdio.h>
#include <math.h>
reactionSolver::reactionSolver(simulationData* mdata)
{

    for (int i=0;i<mdata->m_reactions.size();i++)
        m_reactions.push_back(mdata->m_reactions[i]);
    qDebug()<<m_reactions.size()<< "AA!!";

}

void reactionSolver::solve(int itn)
{
    /*   oValues[0]=ne_i;
    oValues[1]=eps_i;
    oValues[2]=nars_i;
    double k[7];
    double dk[7];
    double dE[7]; //energy gain/loss of reaction

    double x0=ne_0;
    double y0=nars_0;


    double A,D,C,B,E,F;

    for (int i=0;i<7;i++)
    {
        dE[i]=m_reactions[i]->getDe();
    }


    double dt_cur=dt;
    double x,y,eps;
    x=ne_i/x0;y=nars_i/y0;

    eps=eps_i;
   //  qDebug()<<"x0="<<x0<<" y0="<<y0<<" dt="<<dt<<" eps="<<eps_i<<" nn="<<n_n;

    for (int nn=0;nn<itn;nn++)
    {
        for (int i=0;i<7;i++)
        {
            k[i]=m_reactions[i]->getRate(eps*2.0/3.0);
            dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }
        /*        //n_e
                (oValues[0] - ne_i) / dt = oValues[0]*(k[3]*n_n + k[4]*oValues[2]) + k[5]*oValues[2]*nars_i;
                -(oValues[0]-ne_i)/dt+rhs_ne + oValues[0]*(k[3]*n_n +k[4]*oValues[2])+k[5]*oValues[2]*oValues[2];
                //n_ars
                 (oValues[2]-nars_i)/dt=oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2]) -2.0*k[5]*oValues[2]*nars_i;


     (x - x0) / dt - nu * (xr - 2.0 * x + xl) / (dz*dz) - rhs = ...
     x*(1/dt+2.0*nu/(dz*dz)) -x0/dt-nu*(xr+xl)/(dz*dz) - rhs =...




        a=dt*k[3]*n_n;   b=dt*k[1]*n_n; c=dt*nars_i*k[5];
        d=dt*k[4];  e=dt*(k[2]+k[4]);

        //ne
        x - (ne_i+rhs*dt) = x*(a + d*y) + c*y;
        //n_ars
         y-(nars_i+rhs*dt)=x*(b -e*y) -2.0*c*y;
                */

    /*    A=(dt_cur*k[3]*n_n-1.0);
        D=dt_cur*k[4]*y0;
        C=dt_cur*k[5]*y0*(y0/x0);
        //Ax+Dxy+Cy+(1+rhs*dt/x0)=0  ne eqn;  x is x/x0 y is y/y0

        B=dt_cur*k[1]*n_n*(x0/y0);
        E=-dt_cur*(k[2]+k[4])*x0;
        F=-1.0-2.0*dt_cur*k[5]*y0;
        //Bx+Exy+Fy+(1+rhs*dt/y0)=0  nars eqn;  x is x/x0 y is y/y0


       // qDebug()<<"nn="<<nn<<" A="<<A<<" D="<<D<<" C="<<C;
       // qDebug()<<"B="<<B<<" E="<<E<<" F="<<F;



        double f_ne,f_nars, dfx_ne,dfy_ne, dfx_nars,dfy_nars;

        f_ne=A*x+D*x*y+C*y+(1+rhs_ne*dt/x0);
        dfx_ne=A+D*y;
        dfy_ne=C+D*x;

        f_nars=B*x+E*x*y+F*y+(1+rhs_nars*dt/y0);
        dfx_nars=B+E*y;
        dfy_nars=F+E*x;

        double det = dfx_ne*dfy_nars - dfy_ne*dfx_nars;

        double d_x = -f_ne*dfy_nars + f_nars*dfy_ne;
        double d_y = f_ne*dfx_nars - f_nars*dfx_ne;

        x+=d_x;
        y+=d_y;

        //here we've got ne and ars now for two other species (en and arp)

      //  qDebug()<<"f_ne"<<f_ne<<" f_nars="<<f_nars<<" x="<<x<<" y="<<y;


/*
          eps = eps_i +
                        dt*n_n * dE[1]*k[1] + dt*n_n*k3*dE[3] - dt*n_n*k3*eps  +
                        dt*y*y0 *dE[2]*k[2] + dt*y*y0 *k[4]*dE[4] - dt*y*y0 *k[4]*eps   +
                       (dt/ (x*x0))* k[5] * y*y0 * nars_i*dE[5] - (dt/ (x*x0))* k[5] * y*y0 * nars_i*eps  ;
 */

    /*      eps  = (eps_0 + dt*(rhs_eps-rhs_ne)/(x*x0) +
                      dt*n_n * dE[1]*k[1] + dt*n_n*k[3]*dE[3] +
                      dt*y*y0 *dE[2]*k[2] + dt*y*y0 *k[4]*dE[4]    +
                     (dt/ (x*x0))* k[5] * y*y0 * nars_0*dE[5])/(1.0  + dt*(n_n*k[3] + y*y0 *k[4] +  ( k[5] * y*y0 * nars_0/ (x*x0)) )) ;


    }


    ne_o=x*x0;//oValues[0];
    eps_o=eps;//oValues[1];
    nars_o=y*y0;//oValues[2];*/

}



bool reactionSolver::solve_diffuse(int itn)
{

    double k[7];
    //  double dk[7];
    double dE[7]; //energy gain/loss of reaction

    double x0=ne_0;
    double y0=nars_0;


    double A,D,C,B,E,F;

    for (int i=0;i<7;i++)
    {
        dE[i]=m_reactions[i]->getDe();
    }



    double x,y,z,eps;


    double dt_ne,dt_nars,dt_neps,dt_narp; //ne ar

    dt_ne=1.0/(1/dt+2.0*nu_ne/(dz*dz));
    dt_nars=1.0/(1/dt+2.0*nu_nars/(dz*dz));
    dt_neps=1.0/(1/dt+2.0*nu_neps/(dz*dz));
    dt_narp=1.0/(1/dt+2.0*nu_narp/(dz*dz));

    x=ne_i/x0;y=nars_i/y0;

    eps=eps_o;
    for (int nn=0;nn<itn;nn++)
    {

        // lets get neps (z=neps/x0)
        /*
                z*(1/dt+2.0*nu/(dz*dz)) -z0/dt-nu*(zr+zl)/(dz*dz) - rhs =
                                x*x0*(
                                n_n *( dE[1]*k[1] +  dE[3]*k[3]) +
                                y*y0*( dE[2]*k[2] +  dE[4]*k[4])
                                ) +
                                 (dE[5]) * k[5] * y*y0 * nars_i  ;
        */


        for (int i=0;i<7;i++)
        {
            k[i]=m_reactions[i]->getRate(fmax(fmin(eps*2.0/3.0,40),0.0001));
            //    dk[i]=m_reactions[i]->getDeriv(eps*2.0/3.0);
        }
        /*        //n_e
                (oValues[0] - ne_i) / dt = oValues[0]*(k[3]*n_n + k[4]*oValues[2]) + k[5]*oValues[2]*nars_i;
                -(oValues[0]-ne_i)/dt+rhs_ne + oValues[0]*(k[3]*n_n +k[4]*oValues[2])+k[5]*oValues[2]*oValues[2];
                //n_ars
                 (oValues[2]-nars_i)/dt=oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2]) -2.0*k[5]*oValues[2]*nars_i;


     (x - x0) / dt - nu * (xr - 2.0 * x + xl) / (dz*dz) - rhs = ...
     x*(1/dt+2.0*nu/(dz*dz)) -x0/dt-nu*(xr+xl)/(dz*dz) - rhs =...

       x/dt_x - (x0/dt-nu*(xr+xl)/(dz*dz) + rhs) =
       x/dt_x -rhs_x =


        a=dt*k[3]*n_n;   b=dt*k[1]*n_n; c=dt*nars_i*k[5];
        d=dt*k[4];  e=dt*(k[2]+k[4]);

        //ne
        x - rhs_x*dt_x = x*(a + d*y) + c*y;
        //n_ars
         y- rhs_y*dt_y=x*(b -e*y) -2.0*c*y;
                */

        A=(dt_ne*k[3]*n_n-1.0);
        D=dt_ne*k[4]*y0;
        C=dt_ne*k[5]*y0*(y0/x0);
        //Ax+Dxy+Cy+(1+rhs*dt/x0)=0  ne eqn;  x is x/x0 y is y/y0

        B=dt_nars*k[1]*n_n*(x0/y0);
        E=-dt_nars*(k[2]+k[4])*x0;
        F=-1.0-2.0*dt_nars*k[5]*y0;
        //Bx+Exy+Fy+(1+rhs*dt/y0)=0  nars eqn;  x is x/x0 y is y/y0

        double f_ne,f_nars, dfx_ne,dfy_ne, dfx_nars,dfy_nars;

        f_ne=A*x+D*x*y+C*y+(rhs_ne*dt_ne/x0);
        dfx_ne=A+D*y;
        dfy_ne=C+D*x;

        f_nars=B*x+E*x*y+F*y+(rhs_nars*dt_nars/y0);
        dfx_nars=B+E*y;
        dfy_nars=F+E*x;

        double det = dfx_ne*dfy_nars - dfy_ne*dfx_nars;

        double d_x = -f_ne*dfy_nars + f_nars*dfy_ne;
        double d_y = f_ne*dfx_nars - f_nars*dfx_ne;

        x+=d_x;
        y+=d_y;

        if ((x<0.95) || (x>1.05) || (y<0.95) || (y>1.05) )
            return false;//x=0.9;




        z=(
                    rhs_neps/x0 + x*(n_n *( dE[1]*k[1] +  dE[3]*k[3]) + y*y0*( dE[2]*k[2] +  dE[4]*k[4])) + dE[5]*k[5]* y*y0/x0*nars_0
                )*dt_neps;

        //    z=0.1;
        eps  = z/x;


    } //n

    //now lets get arp
    //narp_o*(1/dt+2.0*nu/(dz*dz)) -x0/dt-nu*(xr+xl)/(dz*dz) - rhs = (ne_o)*(k[3]*n_n + k[4]*y*y0) + k[5]*nars_o*nars_i;


    ne_o=x*x0;//oValues[0];
    neps_o=z*x0;//oValues[1];
    eps_o=z/x;
    nars_o=y*y0;//oValues[2];
    narp_o=dt_narp*( rhs_narp + (ne_o)*(k[3]*n_n + k[4]*nars_o) + k[5]*nars_o*nars_0);
    return true;
} //eof


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



//        0=   -oValues[2]/dt+nars_i/dt+rhs_nars + oValues[0]*k[1]*n_n -oValues[0]*(k[2]+k[4])*oValues[2])
//                    -2.0*k[5]*oValues[2]*oValues[2]-k[6]*n_n*oValues[2];

//            0=   oValues[2]*(-1/dt-oValues[0]*(k[2]+k[4]) -k[6]*n_n)+nars_i/dt+rhs_nars + oValues[0]*k[1]*n_n
//                                   -2.0*k[5]*oValues[2]*oValues[2];


// oValues[2]*(1/dt+oValues[0]*(k[2]+k[4]) +k[6]*n_n+2.0*k[5]*nars_i)=   nars_i/dt+rhs_nars oValues[0]*k[1]*n_n;

// oValues[2]=   (nars_i/dt+rhs_nars + oValues[0]*k[1]*n_n)/(1/dt+oValues[0]*(k[2]+k[4]) +k[6]*n_n+2.0*k[5]*nars_i);


/*
        m_reactions.push_back(new reactionEAr_EAr_comsol(this));
        m_reactions.push_back(new reactionEAr_EArs_comsol(this));
        m_reactions.push_back(new reactionEArs_EAr_comsol(this));
        m_reactions.push_back(new reactionEAr_2EArp_comsol(this));
        m_reactions.push_back(new reactionEArs_2EArp_comsol(this));
        m_reactions.push_back(new reactionArsArs_EArArp_comsol(this));
        m_reactions.push_back(new reactionArsAr_ArAr_comsol(this));
*/

// 0=-(oValues[2]-nars_i)/dt+oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2]);

//oValues[2]=nars_i+dt*oValues[0]*(k[1]*n_n /*-(k[2]+k[4])*oValues[2]*/); //tested good for k[1]

//oValues[2]=nars_i/(1.0 + dt*(k[2]*oValues[0]));  //tested good for k[2]

//oValues[2]=nars_i/(1.0 + dt*(k[4]*oValues[0]));  //test good for k[4]

//0=(oValues[2]-nars_i)/dt + 2.0*k[5]*oValues[2]*oValues[2] +k[6]*n_n*oValues[2];

//nars_i/dt=oValues[2]/dt + 2.0*k[5]*oValues[2]*oValues[2] ;
//oValues[2]=nars_i/(1.0 + dt*2.0*k[5]*nars_i); //tested good for k[5]
//oValues[2]=nars_i/(1.0 + dt*k[6]*n_n); //tested good for k[6]


//we will neglect k[6] and k[0] lets build a simple system:
//////////derivation start:
/*        //n_e
        (oValues[0] - ne_i) / dt = oValues[0]*(k[3]*n_n + k[4]*oValues[2]) + k[5]*oValues[2]*nars_i;
        //n_ars
         (oValues[2]-nars_i)/dt=oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2]) -2.0*k[5]*oValues[2]*nars_i;
        */

/*    x=oValues[0];  y=oValues[2];
        a=dt*k[3]*n_n;   b=dt*k[1]*n_n; c=dt*nars_i*k[5];
        d=dt*k[4];  e=dt*(k[2]+k[4]);

        //ne
        x - ne_i = x*(a + d*y) + c*y;
        //n_ars
         y-nars_i=x*(b -e*y) -2.0*c*y;

         ne_i=x0;
         nars_i=y0;

         //ne
         x*e*(a-1) + e*c*y + d*e*x*y +e*x0=0
         //nars
         d*b*x - e*d*x*y -d*y*(1+2*c) +d*y0 =0


         //sum
               x*e*(a-1) + e*c*y + e*x0 + d*b*x -d*y*(1+2*c) +d*y0 =0
                 x*(e*(a-1)+d*b) + y*(e*c - d*(1+2*c)) + e*x0 +d*y0 =0

                 A=(e*(a-1)+d*b);
                B=(e*c - d*(1+2*c));
                C= e*x0 +d*y0;
                A*x+B*y+C=0.0;

                x= -C/A-B/A*y;
                -C/A=D;
                -B/A=E;
                x=D-E*y;

                //substitute
                d*b*(D-E*y) - e*d*(D-E*y)*y -d*y*(1+2*c) +d*y0 =0*/



/*        //n_ars
         oValues[0]=((oValues[2]-nars_i)/dt + 2.0*k[5]*oValues[2]*nars_i)/(k[1]*n_n -(k[2]+k[4])*oValues[2]);
         oValues[0]=(1.0/dt)*(oValues[2]*(1.0 + dt*2.0*k[5]*nars_i)  -nars_i)/(k[1]*n_n -(k[2]+k[4])*oValues[2]);
        //n_e
        oValues[0]  =(1.0/dt)*(ne_i + dt*k[5]*oValues[2]*nars_i)/( 1.0/dt -(k[3]*n_n + k[4]*oValues[2]));
         //combin:


        (oValues[2]*(1.0 + dt*2.0*k[5]*nars_i)  -nars_i)/(k[1]*n_n -(k[2]+k[4])*oValues[2]) =
         (ne_i + dt*k[5]*oValues[2]*nars_i)/( 1.0/dt -(k[3]*n_n + k[4]*oValues[2]));


        ( 1.0/dt -k[3]*n_n + k[4]*oValues[2])*(oValues[2]*(1.0 + dt*2.0*k[5]*nars_i)  -nars_i) =
         (ne_i + dt*k[5]*oValues[2]*nars_i)*(k[1]*n_n -(k[2]+k[4])*oValues[2]);


        oValues[2]=x;
        1.0/dt -k[3]*n_n = a;
        1.0 + dt*2.0*k[5]*nars_i =b;
        dt*k[5]*nars_i=c;
        (k[2]+k[4])=d;
        k[4]=e;
        nars_i=f;
        ne_i=g;
        k[1]*n_n=h;

        (a + e*x)*(x*b - f) = (g + c*x)*(h - d*x);

        a*b*x -a*f + e*b*x*x -e*f*x = (g*h-d*g*x + c*h*x- c*d*x*x);       */



/*
               //eps
        (oValues[1] - eps_i) / dt=
                        n_n * ( A[1]*k[1] + ( A[3] - oValues[1] ) * k[3]) +
                        oValues[2] * ( A[2]*k[2] + ( A[4] - oValues[1] )*k[4] ) +
                        ( ( A[5] - oValues[1] ) * k[5] * oValues[2] * nars_i ) / oValues[0];
*/

////////////derivation end



//oValues[2]=   (nars_i/dt+oValues[0]*k[1]*n_n)/(1/dt);

/* a=-2.0*k[5];
    b=-1/dt-oValues[0]*(k[2]+k[4]) -k[6]*n_n;
    c=nars_i/dt+ oValues[0]*k[1]*n_n;

    D=b*b-4*a*c;

    oValues[2]=(-b-sqrt(D))/(2*a);

    qDebug()<<"D="<<D<<" X0="<<(-b+sqrt(D))/(2*a)<<" X1="<<(-b-sqrt(D))/(2*a);*/

/*     f[2]=-(oValues[2]-nars_i)/dt+rhs_nars + oValues[0]*(k[1]*n_n -(k[2]+k[4])*oValues[2])
                -2.0*k[5]*oValues[2]*oValues[2]-k[6]*n_n*oValues[2];


        df[2][2]=-(1)/dt+ oValues[0]*( -(k[2]+k[4]))
                -4.0*k[5]*oValues[2]-k[6]*n_n;

        solve3x3(df,f,dValues);*/

//oValues[0]-=dValues[0];
//oValues[1]-=dValues[1];
//oValues[2]-=f[2]/df[2][2];

//        oValues[2]=nars_i+dt*(rhs_nars + oValues[0]*(k[1]*n_n) -(k[2]+k[4])*oValues[2]
//              -2.0*k[5]*oValues[2]*oValues[2]-k[6]*n_n*oValues[2]);
