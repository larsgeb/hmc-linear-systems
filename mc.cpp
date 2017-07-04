#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mc.hpp"
#include "aux.hpp"

/********************************************************************************//**
 * Base class for Monte Carlo sampling.
************************************************************************************/

/* Constructor. -------------------------------------------------------------------*/
mc::mc(int iterations_in, int nt_in, double dt_in, int Nq, bool verbose_in)
{
    /* Local variables. -----------------------------------------------------------*/
    double *dtdqi; /* vsp derivatives*/
    
    /* Set basic parameters. ------------------------------------------------------*/

    nt=nt_in;
    dt=dt_in;
    verbose=verbose_in;
    iterations=iterations_in;

    /* Read input parameters. */
    q.read_input("INPUT/parameters.txt");
    q_new.read_input("INPUT/parameters.txt");
    
    /* Read observed and compute prior data. */
    dat.read_data("DATA/synthetics.txt");
    syn.make_synthetics(q);
    
    /* Initialise random number generator. ----------------------------------------*/
    srand (time(0));
    
    /* Initialise pre-computed matrices and vectors. ------------------------------*/
    B=new double[Nq];
    A=new double*[Nq];
    for (int i=0; i<Nq; i++) A[i]=new double[Nq];
    
    duxi=new double[dat.nt]; duyi=new double[dat.nt]; duzi=new double[dat.nt];
    duxj=new double[dat.nt]; duyj=new double[dat.nt]; duzj=new double[dat.nt];
    
    /* Compute C. -----------------------------------------------------------------*/
    C=0.0;
    for (int rec=0; rec<dat.nrec; rec++)
    {
        for (int it=0; it<dat.nt; it++)
        {
            C+=pow(syn.ux[rec][it]-dat.ux[rec][it],2.0)*dat.cx[rec][it]+pow(syn.uy[rec][it]-dat.uy[rec][it],2.0)*dat.cy[rec][it]+pow(syn.uz[rec][it]-dat.uz[rec][it],2.0)*dat.cz[rec][it];
        }
    }
    C=C*dat.dt;
    
    if (verbose) printf("C=%.2e\n\n",C);
    
    /* Compute B[i]. --------------------------------------------------------------*/
    if (verbose) printf("B=");
    for (int i=0; i<Nq; i++)
    {
        B[i]=0.0;
        for (int rec=0; rec<dat.nrec; rec++)
        {
            dat.du(q,0,i,rec,duxi);
            dat.du(q,1,i,rec,duyi);
            dat.du(q,2,i,rec,duzi);
            
            for (int it=0; it<dat.nt; it++)
            {
                B[i]+=duxi[it]*(syn.ux[rec][it]-dat.ux[rec][it])*dat.cx[rec][it]+duyi[it]*(syn.uy[rec][it]-dat.uy[rec][it])*dat.cy[rec][it]+duzi[it]*(syn.uz[rec][it]-dat.uz[rec][it])*dat.cz[rec][it];
            }
        }
        B[i]=B[i]*dat.dt;
        if (verbose) printf("%.2e ",B[i]);
    }
    if (verbose) printf("\n\n");

    /* Compute A[i][j]. -----------------------------------------------------------*/
    for (int i=0; i<Nq; i++)
    {
        for (int j=0; j<Nq; j++)
        {
            A[i][j]=0.0;
            
            for (int rec=0; rec<dat.nrec; rec++)
            {
                dat.du(q,0,i,rec,duxi); dat.du(q,1,i,rec,duyi); dat.du(q,2,i,rec,duzi);
                dat.du(q,0,j,rec,duxj); dat.du(q,1,j,rec,duyj); dat.du(q,2,j,rec,duzj);
                
                for (int it=0; it<dat.nt; it++)
                {
                    A[i][j]+=duxi[it]*duxj[it]*dat.cx[rec][it]+duyi[it]*duyj[it]*dat.cy[rec][it]+duzi[it]*duzj[it]*dat.cz[rec][it];
                }
            }
            A[i][j]=A[i][j]*dat.dt;
        }
        
        /* Add prior. */
        A[i][i]+=1.0/(q.sigma_q[i]*q.sigma_q[i]);
    }

    if (verbose)
    {
        printf("A=\n");
        for (int i=0; i<Nq; i++)
        {
            for (int j=0; j<Nq; j++)
            {
                if (A[i][j]>=0)
                {
                    printf("+%.2e ",A[i][j]);
                }
                else
                {
                    printf("%.2e ",A[i][j]);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
    
    /* Initialise mass matrix. ----------------------------------------------------*/
    m=new double[Nq];
    for (int i=0; i<Nq; i++) m[i]=A[i][i];
    
    /* Initialise models and momenta. ---------------------------------------------*/
    p=new double[Nq];
    p_new=new double[Nq];
    
    for (int i=0; i<Nq; i++)
    {
        p_new[i]=randn(0.0,sqrt(m[i]));
        q_new.q[i]=randn(q.mean_q[i],q.sigma_q[i]);
        q.q[i]=q_new.q[i];
    }
    
    /* Clean up. ------------------------------------------------------------------*/
    delete[] duxi; delete[] duyi; delete[] duzi;
    delete[] duxj; delete[] duyj; delete[] duzj;
}

/* Destructor. --------------------------------------------------------------------*/
mc::~mc()
{
    if (p) delete[] p;
    if (p_new) delete[] p_new;
    if (m) delete[] m;
    if (B) delete[] B;
    if (A)
    {
        for (int i=0; i<Nq; i++) delete[] A[i];
        delete[] A;
    }
}


/********************************************************************************//**
 * Proposals.
************************************************************************************/

/* Propose test model based on prior. ---------------------------------------------*/
void mc::propose_metropolis()
{
    for (int i=0; i<Nq; i++) q_new.q[i]=randn(q.mean_q[i],q.sigma_q[i]);
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void mc::propose_hamilton()
{
    /* Draw random prior momenta. */
    for (int i=0; i<Nq; i++) p[i]=randn(0.0,sqrt(m[i]));
    
    /* Integrate Hamilton's equations. */
    leap_frog(verbose);
}


/********************************************************************************//**
 * Misfits.
************************************************************************************/

/* Misfit for Metropolis Hastings. ------------------------------------------------*/
double mc::chi()
{
    /* Local variables. */
    double likelihood;
    double a;
    bool linearised=true; /* Misfit linearised or not. */
    
    /* Compute misfit. */
    if (linearised==true)
    {
        likelihood=0.5*C;
        for (int i=0; i<Nq; i++)
        {
            likelihood+=(q_new.q[i]-q.mean_q[i])*B[i];
            for (int j=0; j<Nq; j++)
            {
                a=A[i][j];
                if (j==i) a-=1.0/(q.sigma_q[i]*q.sigma_q[i]); /* Subtract prior to produce likelihood. */
                likelihood+=0.5*(q_new.q[i]-q.mean_q[i])*(q_new.q[j]-q.mean_q[j])*a;
            }
        }
    }
    else
    {
        likelihood=0.0;
        syn.make_synthetics(q_new);
        for (int n=0; n<dat.nrec; n++)
        {
            for (int i=0; i<dat.nt; i++)
            {
                likelihood+=pow(dat.ux[n][i]-syn.ux[n][i],2.0)*dat.cx[n][i]+pow(dat.uy[n][i]-syn.uy[n][i],2.0)*dat.cy[n][i]+pow(dat.uz[n][i]-syn.uz[n][i],2.0)*dat.cz[n][i];
            }
        }
        likelihood=0.5*likelihood*dat.dt;
    }

    return likelihood;
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double mc::energy()
{
    /* Local variables. */
    double H;
    bool linearised=true; /* Misfit linearised or not. */
    
    /* Compute energy - model part. */
    if (linearised)
    {
        H=0.5*C;
    
        for (int i=0; i<Nq; i++)
        {
            H+=(q_new.q[i]-q.mean_q[i])*B[i];
            for (int j=0; j<Nq; j++)
            {
                H+=0.5*(q_new.q[i]-q.mean_q[i])*(q_new.q[j]-q.mean_q[j])*A[i][j];
            }
        }
    }
    else
    {
        H=0.0;
        syn.make_synthetics(q_new);
        for (int n=0; n<dat.nrec; n++)
        {
            for (int i=0; i<dat.nt; i++)
            {
                H+=pow(dat.ux[n][i]-syn.ux[n][i],2.0)*dat.cx[n][i]+pow(dat.uy[n][i]-syn.uy[n][i],2.0)*dat.cy[n][i]+pow(dat.uz[n][i]-syn.uz[n][i],2.0)*dat.cz[n][i];
            }
        }
        H=0.5*H*dat.dt;
        for (int i=0; i<Nq; i++)
        {
            H+=0.5*pow(q_new.q[i]-q.mean_q[i],2.0)/pow(q.sigma_q[i],2.0);
        }
    }
    
    /* Compute energy - momentum part. */
    for (int i=0; i<Nq; i++) H+=0.5*p_new[i]*p_new[i]/m[i];
    
    return H;
}

/*=================================================================================*/
/* Right-hand sides of Hamilton equations. ----------------------------------------*/
/*=================================================================================*/

void mc::dHdp(double *p_in, double *rhs)
{
    for (int i=0; i<Nq; i++) rhs[i]=p_in[i]/m[i];
}

void mc::dHdq(parameters &q_in, double *rhs)
{
    /* March through components. */
    for (int k=0; k<Nq; k++)
    {
        rhs[k]=B[k];
        for (int i=0; i<Nq; i++) rhs[k]+=(q_in.q[i]-q_in.mean_q[i])*A[i][k];
    }
}

/********************************************************************************//**
 * Miscellaneous.
************************************************************************************/

/* Write a sample to an open file. ------------------------------------------------*/
void mc::write_sample(FILE *pfile, double misfit, int iteration)
{
    if (iteration==0) fprintf(pfile,"%d %d\n",Nq,iterations+1);
    
    for (int i=0; i<Nq; i++) fprintf(pfile,"%lg ",q.q[i]);
    fprintf(pfile,"%lg ",misfit);
    fprintf(pfile,"\n");
}

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
void mc::leap_frog(bool verbose)
{
    /* Local variables and setup. -------------------------------------------------*/
    
    double *p_half, *p_init, *out;
    double angle1, angle2;
    parameters q_init;
    
    out=new double[Nq];
    p_half=new double[Nq];
    p_init=new double[Nq];
    
    FILE *pfile;
    if (verbose) pfile=fopen("OUTPUT/trajectory.txt","w");
    
    /* Set initial values. --------------------------------------------------------*/
    
    q_init=q;
    q_new=q;
    for (int i=0; i<Nq; i++) p_init[i]=p[i];
    
    /* March forward. -------------------------------------------------------------*/
    
    if (verbose) fprintf(pfile,"%d %d\n",2*Nq,nt);
    dHdq(q_init,out);
    
    for (int it=0; it<nt; it++)
    {
        /* Some output. */
        if (verbose)
        {
            for (int i=0; i<Nq; i++) fprintf(pfile,"%lg ",q_init.q[i]);
            for (int i=0; i<Nq; i++) fprintf(pfile,"%lg ",p_init[i]);
            fprintf(pfile,"\n");
        }
        
        /* First half step in momentum. */
        for (int i=0; i<Nq; i++)
        {
            p_half[i]=p_init[i]-0.5*dt*out[i];
        }
        
        /* Full step in position. */
        dHdp(p_half,out);
        for (int i=0; i<Nq; i++)
        {
            q_new.q[i]=q_init.q[i]+dt*out[i];
        }
        
        /* Second half step in momentum. */
        dHdq(q_new,out);
        for (int i=0; i<Nq; i++)
        {
            p_new[i]=p_half[i]-0.5*dt*out[i];
        }
        
        /* Update position and momentum. */
        for (int i=0; i<Nq; i++)
        {
            p_init[i]=p_new[i];
            q_init.q[i]=q_new.q[i];
        }
        
        /* Check no-U-turn criterion. */
        angle1=0.0;
        angle2=0.0;
        for (int i=0; i<Nq; i++)
        {
            angle1+=p_new[i]*(q_new.q[i]-q.q[i]);
            angle2+=p[i]*(q.q[i]-q_new.q[i]);
        }
        
        if (angle1<0.0 && angle2<0.0)
        {
            if (verbose) printf("steps: %d\n",it);
            break;
        }
    }
    
    /* Clean up. ------------------------------------------------------------------*/
    
    if (verbose) fclose(pfile);
    
    delete[] p_half;
    delete[] p_init;
    delete[] out;
}




