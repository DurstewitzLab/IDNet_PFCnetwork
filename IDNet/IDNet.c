/*  IDNet.c
MEX function written in C to simulate arbitrary biological neural networks
To be used with the MATLAB wrapper IDNetSim.m
 
*/

/* -------------------------------------------------------------------------------------------------- */
/* ----------------------- Declarations and definitions --------------------------------------------- */
/* -------------------------------------------------------------------------------------------------- */

/* Include libaries */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include "mex.h"

/* Define global parameters */
#define pi 3.14159265418
#define MaxNumSTperN 20
#define SizeHistOutput 10       /* Sets the number spikes that are considered to modify EPSPs */
#define SizeHistInput 1000000   /* Sets the number of possible spikes in the input neurons */
#define NumCtrPar 5
#define NumVar 2                
#define NumNeuPar 12
#define NumSynTypePar 8
#define NumSynPar 7
#define TRUE 1
#define FALSE 0

/* declare data structures */
struct SynType 
{
    /* parameters */
    int No;
    double gmax;
    double tc_on;
    double tc_off;
    double Erev;
    double Mg_gate;
    double Mg_fac;
    double Mg_slope;
    double Mg_half;
    
    /* variables */
    double Gsyn; 
};

struct Neuron 
{ 
    /* parameters */
    double Cm;
    double gL;
    double EL;
    double sf;
    double Vup;
    double tcw;
    double a;
    double b;
    double Vr;
    double Vth;
    double I_ref;
    double v_dep;
    int NumSynType;

    /* variables */
    double Iinj;
    double v[2];
    double dv[2];
    struct SynType *STList[MaxNumSTperN];
    double gfONsyn[MaxNumSTperN];
    double gfOFFsyn[MaxNumSTperN]; 
    double gfONnoise[MaxNumSTperN];
    double gfOFFnoise[MaxNumSTperN]; 
    double SpikeTimes[SizeHistOutput];
    int NumPreSyn;
    int *PreSynList;
    struct SynDepr *SDf; 
};

struct InpNeuron 
{ 
    double SPtrain[SizeHistInput];
    double SpikeTimes[SizeHistInput];
    int SP_ind;
    int NumSynType;
    int NumPreSyn;
    int *PreSynList;
    struct SynDepr *SDf; 
};


struct Synapse 
{
    /* parameters */
    struct SynType *STPtr;
    double dtax;
    double wgt; 
    double p_fail;
    int PreSynIdx; 
};

struct SynDepr 
{
    double use;
    double tc_rec;
    double tc_fac; 
    double Adepr[SizeHistOutput];
    double uprev[SizeHistOutput];
    double Rprev[SizeHistOutput]; 
};

struct SynList 
{
    int NumSyn;
    struct Synapse *Syn; 
};

/* global declarations */
struct SynList NoiseSyn;
double wV;
double D0;
double gsyn_AN;
double gsyn_G;
double I_tot;
double flag_regime_osc;
int flag_dv;
int stop_flag;



/* -------------------------------------------------------------------------------------------------- */
/* --------------------------------------- NR components -------------------------------------------- */
/* -------------------------------------------------------------------------------------------------- */

#define NR_END 1
#define FREE_ARG char*
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

        
long *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	long *v;

	v=(long *)mxMalloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) mexErrMsgTxt("allocation failure in ivector()");
	return v-nl+NR_END;
}


double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)mxMalloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) mexErrMsgTxt("allocation failure in dvector()");
	return v-nl+NR_END;
}


void free_ivector(long *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	mxFree((FREE_ARG) (v+nl-NR_END));
}


void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	mxFree((FREE_ARG) (v+nl-NR_END));
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) mxMalloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) mexErrMsgTxt("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) mxMalloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) mexErrMsgTxt("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) mxMalloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) mexErrMsgTxt("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) mxMalloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) mexErrMsgTxt("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	mxFree((FREE_ARG) (m[nrl]+ncl-NR_END));
	mxFree((FREE_ARG) (m+nrl-NR_END));
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	mxFree((FREE_ARG) (m[nrl]+ncl-NR_END));
	mxFree((FREE_ARG) (m+nrl-NR_END));
}


/* -------------------------------------------------------------------------------------------------- */
/* ------------------------------  Computational routines ------------------------------------------- */
/* -------------------------------------------------------------------------------------------------- */

double round(double x)
/* simple implementation of the round function */ 
{ 
    return floor(x+0.5); 
}


void IDderiv (struct Neuron np, double *v, double dt, double *dv)
/* ODEs that define the model */
/* Computes left-hand side of ODEs for a single neuron */

/* Input: neuron (np), current variables (v), time step (dt) */
/* Output: derivative of variables (dv) */
{
    double z,sgate,Isyn;
    int j;
    struct SynType syn;
    double I_ex,dD0;

    /* Compute synaptic current Isyn from the three synapse types */
    /* implementing a double-exponential decay with time (for spikes, see below) */
    gsyn_AN=0;
    gsyn_G=0;
    for (Isyn=0.0,j=0;j<np.NumSynType;j++) 
    {
	    syn=*np.STList[j];
	    sgate=1.0;
	    if (syn.Mg_gate>0.0)            /* use Mg gate if flag is on (for NMDA only) */
	        sgate=syn.Mg_gate/(1.0+syn.Mg_fac*exp(syn.Mg_slope*(syn.Mg_half-v[0])));
        
	    Isyn+=sgate*(np.gfOFFsyn[j]*exp(-dt/syn.tc_off)-np.gfONsyn[j]*exp(-dt/syn.tc_on))*(syn.Erev-v[0]);
        if (syn.Erev==0.0)
            gsyn_AN=gsyn_AN+sgate*(np.gfOFFsyn[j]*exp(-dt/syn.tc_off)-np.gfONsyn[j]*exp(-dt/syn.tc_on));
        else 
            gsyn_G=gsyn_G+sgate*(np.gfOFFsyn[j]*exp(-dt/syn.tc_off)-np.gfONsyn[j]*exp(-dt/syn.tc_on));
    }
    
    /* Compute noise input */
    for (j=0;j<NoiseSyn.NumSyn;j++)
    {
        syn = *NoiseSyn.Syn[j].STPtr;
	    sgate=1.0;
	    if (syn.Mg_gate>0.0)            /* use Mg gate if flag is on (for NMDA only) */
	        sgate=syn.Mg_gate/(1.0+syn.Mg_fac*exp(syn.Mg_slope*(syn.Mg_half-v[0])));
        
        Isyn+=sgate*(np.gfOFFnoise[j]*exp(-dt/syn.tc_off)-np.gfONnoise[j]*exp(-dt/syn.tc_on))*(syn.Erev-v[0]);
        if (syn.Erev==0.0)
            gsyn_AN=gsyn_AN+sgate*(np.gfOFFnoise[j]*exp(-dt/syn.tc_off)-np.gfONnoise[j]*exp(-dt/syn.tc_on));
        else
            gsyn_G=gsyn_G+sgate*(np.gfOFFnoise[j]*exp(-dt/syn.tc_off)-np.gfONnoise[j]*exp(-dt/syn.tc_on));
    }  
    
 
    /* Exponential term */   
    I_ex=np.gL*np.sf*exp((v[0]-np.Vth)/np.sf);
    
    /* V-nullcline at v[0] */
    wV = np.Iinj + Isyn - np.gL*(v[0]-np.EL) + I_ex;
    
    /* Calculation of D_0 */
    D0 = (np.Cm/np.gL) * wV;
    
    
    /* Compute membrane potential derivative from all currents */
    if((np.Iinj + Isyn)>=np.I_ref && flag_dv==0) /* use flag_dv for restriction of depolarization block to time x after spike */
    {
        dv[0] = -(np.gL/np.Cm) * (v[0]-np.v_dep);
        flag_regime_osc=0;}
    else
    {
        dv[0] = (np.Iinj - np.gL*(v[0]-np.EL) - v[1] + I_ex + Isyn)/np.Cm;
        flag_regime_osc=1;}  
    
    
    /* derivative of D_0 with respect to V */
    dD0=np.Cm*(exp((v[0]-np.Vth)/np.sf)-1);
    
    /* second differential equation */
    if((v[1]>wV-D0/np.tcw) && (v[1]<wV+D0/np.tcw) && v[0]<=np.Vth && (np.Iinj+Isyn)<np.I_ref)
        dv[1]=-(np.gL*(1-exp((v[0]-np.Vth)/np.sf)) + dD0/np.tcw)*dv[0];
    else
        dv[1]=0;
    
    I_tot=Isyn+np.Iinj;
}


void set_gsyn (struct Neuron np, double dt, double v)
/* double set_Isyn (struct Neuron np, double dt, double v)*/
/* Compute synaptic current Isyn from the three synapse types */
/* implementing a double-exponential decay with time (for spikes, see below) */
{
    double sgate,Isyn;
    int j;
    struct SynType syn;
    
    gsyn_AN=0;
    gsyn_G=0;
    
    for (Isyn=0.0, j=0;j<np.NumSynType;j++) 
    {
        syn = *np.STList[j];
        sgate=1.0;
        if (syn.Mg_gate>0.0)            /* use Mg gate if flag is on (for NMDA only) */
            sgate=syn.Mg_gate/(1.0+syn.Mg_fac*exp(syn.Mg_slope*(syn.Mg_half-v))); 
        
        Isyn+=sgate*(np.gfOFFsyn[j]*exp(-dt/syn.tc_off)-np.gfONsyn[j]*exp(-dt/syn.tc_on))*(syn.Erev-v);
        if (syn.Erev==0.0)
            gsyn_AN=gsyn_AN+sgate*(np.gfOFFsyn[j]*exp(-dt/syn.tc_off)-np.gfONsyn[j]*exp(-dt/syn.tc_on));
        else
            gsyn_G=gsyn_G+sgate*(np.gfOFFsyn[j]*exp(-dt/syn.tc_off)-np.gfONsyn[j]*exp(-dt/syn.tc_on));
    }
        
    /* Compute noise input */ 
    for (j=0;j<NoiseSyn.NumSyn;j++)
    {
        syn = *NoiseSyn.Syn[j].STPtr;
	    sgate=1.0;
	    if (syn.Mg_gate>0.0)            /* use Mg gate if flag is on (for NMDA only) */
	        sgate=syn.Mg_gate/(1.0+syn.Mg_fac*exp(syn.Mg_slope*(syn.Mg_half-v)));
        
        Isyn+=sgate*(np.gfOFFnoise[j]*exp(-dt/syn.tc_off)-np.gfONnoise[j]*exp(-dt/syn.tc_on))*(syn.Erev-v);
        if (syn.Erev==0.0)
            gsyn_AN=gsyn_AN+sgate*(np.gfOFFnoise[j]*exp(-dt/syn.tc_off)-np.gfONnoise[j]*exp(-dt/syn.tc_on));
        else
            gsyn_G=gsyn_G+sgate*(np.gfOFFnoise[j]*exp(-dt/syn.tc_off)-np.gfONnoise[j]*exp(-dt/syn.tc_on));
    }    
    
    I_tot=Isyn+np.Iinj;
}


void update (struct Neuron *np, double dt)
/* Integrates ODEs by an explicit Runge-Kutta method of 2nd order */

/* Input: neuron (np), time step (dt) */
/* Output: updated neuron (np) */
{
    double dv1[2],dv2[2],v[2];
    int i,j,nvar;
    
    nvar=2;
    
    for (i=0;i<nvar;i++) 
        v[i]=(*np).v[i];
    
    IDderiv(*np,v,0.0,dv1);
    
    for (i=0;i<nvar;i++) 
        v[i]+=dt*dv1[i];
    IDderiv(*np,v,dt,dv2); 
    
    for (i=0;i<nvar;i++) 
    {
        (*np).v[i]+=dt/2.0*(dv1[i]+dv2[i]);
        (*np).dv[i]=dt/2.0*(dv1[i]+dv2[i]);
    }
    
    /* Make a jump in w when it approaches the nullcline wV from the right to avoid singularities */
    if(((*np).v[1]>wV-D0/(*np).tcw) && ((*np).v[1]<wV+D0/(*np).tcw) && (*np).v[0]<=(*np).Vth) 
        (*np).v[1]=wV-(D0/(*np).tcw);
}


double syndepr (struct SynDepr *Syn, double ISI, int Nsp)
/* Computes variables for short-term synaptic plasticity */
/* See Tsodyks and Markram 1997 for details */

/* Input: Plastic synapse struct (Syn), inter-spike interval (ISI), # presynaptic spikes (Nsp) */
/* Output: synaptic modification factor (R*u) */
{
  double qu, qR, u, R;

  qu=(*Syn).uprev[Nsp]*exp(-ISI/(*Syn).tc_fac);
  qR=exp(-ISI/(*Syn).tc_rec);
  u=qu+(*Syn).use*(1.0-qu);
  R=(*Syn).Rprev[Nsp]*(1.0-(*Syn).uprev[Nsp])*qR+1.0-qR;       /* Corrected: (1.0-u_(k-1)) instead of (1.0-u_k), see Haeussler and Maass 2007 */
  (*Syn).uprev[(Nsp+1) % SizeHistOutput]=u;
  (*Syn).Rprev[(Nsp+1) % SizeHistOutput]=R;  

  return(R*u);
}


/* -------------------------------------------------------------------------------------------------- */
/* ------------------------------  Mex function ----------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------- */


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
/* Simulates the network dynamics */

/* Input:  see IDNetSim.m */
/* Output: see IDNetSim.m */
{
    /******************* Initialization *********************************/
    
    double *ViewList,*CtrPar,*NeuPar,*STypPar,*SynPar,*EvtMtx,*EvtTimes,*InpSTtrains,*NoiseDistr, *V0;
    double *NPList,*SPList,*NumSynInp,*TnextSyn, *SynExpOn,*SynExpOff, *gsyn1, *gsyn2, *Isyn, *flag_osc, *N_osc;
    double Tstart,Tstop,t0,t1,t11,t0_i,t1_i,vp,wp,w_Vup,NextEvtT,EvtOffT;
    double t_display,dt,dt0,ISI,ISI_inp,Aall,rand_num,NoiseStep,k11,k12;
    double *UniqueNum,TauRise,SynNormalizationFactor;
    int UniquePrint;    
    int NumView,NumOutp,WriteST;
    int N,M,L,i,k,k2,NumSynType,MaxNumSyn,nst,STno,NumEvt;
    int EvtNo,numST;
    int eom_ind, print_flag;
    long *NumSpike,kk,j;
    FILE *fpOut,*fpOut2,*fpOut3, *fpISI;
    char InList,fn[20],uno[12],IDNFilename[20],IDN2Filename[20],IDN3Filename[20],IDN4Filename[20];
    struct SynList **ConMtx0;
    struct Neuron *NPtr0,*NPtr;
    struct InpNeuron *InpNPtr0, *InpNPtr;
    struct SynType *SynTPtr0,*STPtr,**STPtr2;
    struct Synapse *SynLPtr;
    double SynSum;
    double *NeuronGroupsSaveArray;
    int NumViewGroups, NumNeuronsPerGroup;
    
    /* Handle wrong number of input parameters */
    if (nrhs!=14) 
	mexErrMsgTxt("14 inputs required: CtrPar,NeuPar,NPList,STypPar,SynPar,SPMtx,EvtMtx,EvtTimes,ViewList,InpSTtrains,NoiseDistr,V0,UniqueNum,NeuronGroupsSaveArray");

    NumOutp=nlhs;

    /* Create unique filename identifier string (used to avoid conflicts if multiple instances are run simultaneously) */
    UniqueNum = mxGetPr(prhs[12]);
    UniquePrint = (int)*UniqueNum;
    
    /* Read-in control parameters from CtrPar*/
    CtrPar=mxGetPr(prhs[0]);
    if (mxGetN(prhs[0])*mxGetM(prhs[0])<NumCtrPar)
    	mexErrMsgTxt("CtrPar should have NumCtrPar elements");
    Tstart=CtrPar[0];
    Tstop=CtrPar[1];
    dt0=CtrPar[2];
    WriteST=(int)CtrPar[4];
    t_display=0;
    stop_flag = 0;
    
    /* Create arrays and parameters related to neuron group synaptic saving */
    NeuronGroupsSaveArray = mxGetPr(prhs[13]);
    NumViewGroups = mxGetN(prhs[13]);
    NumNeuronsPerGroup = mxGetM(prhs[13]);

    /* Read-in initial conditions */
    V0 = mxGetPr(prhs[11]);
    if (mxGetM(prhs[11])<NumVar)
        mexErrMsgTxt("V0 should have #Var rows");
    
    
    /* Initialize neuron parameters NRtr with NeuPar according to NPList (one parameter set for each neuron type) */
    NeuPar=mxGetPr(prhs[1]);
    NPList=mxGetPr(prhs[2]);
    i=mxGetN(prhs[2]);
    j=mxGetM(prhs[2]);
    N=i*j;
    k=mxGetN(prhs[1]);
    if (mxGetM(prhs[1])!=NumNeuPar) 
        mexErrMsgTxt("NeuPar should have NumNeuPar rows");
    
    if (mxGetN(prhs[11])<N)
        mexErrMsgTxt("V0 should have N columns");
    
    
    /* Create NRtr pointer*/
    NPtr0=(struct Neuron *) mxMalloc(N*sizeof(struct Neuron));
    NumSpike=ivector(0,N-1);
    gsyn1=dvector(0,N-1);
    gsyn2=dvector(0,N-1);
    Isyn=dvector(0,N-1);
    flag_osc=dvector(0,N-1);
    for (i=0,NPtr=NPtr0;i<N;i++,NPtr++) 
    {
	    if ((int)NPList[i]>k || (int)NPList[i]<1) 
            mexErrMsgTxt("Illegal entry in NPList");
        
	    (*NPtr).Cm=NeuPar[0+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).gL=NeuPar[1+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).EL=NeuPar[2+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).sf=NeuPar[3+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).Vup=NeuPar[4+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).tcw=NeuPar[5+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).a=NeuPar[6+((int)NPList[i]-1)*NumNeuPar];
        (*NPtr).b=NeuPar[7+((int)NPList[i]-1)*NumNeuPar];
        (*NPtr).Vr=NeuPar[8+((int)NPList[i]-1)*NumNeuPar];
        (*NPtr).Vth=NeuPar[9+((int)NPList[i]-1)*NumNeuPar]; 
        (*NPtr).I_ref=NeuPar[10+((int)NPList[i]-1)*NumNeuPar];
        (*NPtr).v_dep=NeuPar[11+((int)NPList[i]-1)*NumNeuPar];
	    (*NPtr).Iinj=0.0; 
	    (*NPtr).v[0]=V0[0+i*NumVar];
	    (*NPtr).v[1]=V0[1+i*NumVar];
	    (*NPtr).NumSynType=0;
	    (*NPtr).NumPreSyn=0;
	    for (j=0;j<MaxNumSTperN;j++) 
            (*NPtr).STList[j]=NULL; 
	    (*NPtr).PreSynList=NULL;
	    (*NPtr).SDf=NULL;	   
        NumSpike[i]=0; 
        gsyn1[i]=0; 
        gsyn2[i]=0; 
        Isyn[i]=0; 
        flag_osc[i]=0;
    }

    
    /* Read InpNeuron parameters */
    InpSTtrains=mxGetPr(prhs[9]);
    M = mxGetM(prhs[9]);

    /* Create InpNRtr pointer */
    InpNPtr0=(struct InpNeuron *) mxMalloc(M*sizeof(struct InpNeuron));
    for (i=0,InpNPtr=InpNPtr0;i<M;i++,InpNPtr++) 
    {
        eom_ind = SizeHistInput;
        (*InpNPtr).SP_ind = 0;
        for (j=0;j<eom_ind;j++)     /* Set all valid entries */
        {
            if ((eom_ind == SizeHistInput) && (InpSTtrains[i+(j+1)*M] == -1))
                eom_ind = j+1;

            (*InpNPtr).SPtrain[j] = InpSTtrains[i+j*M];
        }

        for (j=eom_ind;j<SizeHistInput;j++)  /* fill up the remaining entries with -1 */
            (*InpNPtr).SPtrain[j] = -1;

        (*InpNPtr).NumSynType=0;
	    (*InpNPtr).NumPreSyn=0;
	    (*InpNPtr).PreSynList=NULL;
	    (*InpNPtr).SDf=NULL;	   
    }

    /* Create vector containing spike counts for both neurons and input neurons */
    NumSpike=ivector(0,N+M-1);
    for (i=0;i<N+M;i++) 
        NumSpike[i]=0; 
    
    /* Initialize synaptic parameter STPtr with STypPar (one parameter set for each synapse type) */
    STypPar=mxGetPr(prhs[3]);
    NumSynType=mxGetN(prhs[3]);
    if (mxGetM(prhs[3])!=NumSynTypePar) 
        mexErrMsgTxt("STypPar should have NumSynTypePar rows");

    /* Create SynTRtr pointer*/
    SynTPtr0=(struct SynType *) mxMalloc(NumSynType*sizeof(struct SynType));
    for (i=0,STPtr=SynTPtr0;i<NumSynType;i++,STPtr++) 
    {
	    (*STPtr).No=i;
	    (*STPtr).gmax=STypPar[0+i*NumSynTypePar];
	    (*STPtr).tc_on=STypPar[1+i*NumSynTypePar];
	    (*STPtr).tc_off=STypPar[2+i*NumSynTypePar];
	    (*STPtr).Erev=STypPar[3+i*NumSynTypePar];
	    (*STPtr).Mg_gate=STypPar[4+i*NumSynTypePar];
	    (*STPtr).Mg_fac=STypPar[5+i*NumSynTypePar];
	    (*STPtr).Mg_slope=STypPar[6+i*NumSynTypePar];
	    (*STPtr).Mg_half=STypPar[7+i*NumSynTypePar];
	    (*STPtr).Gsyn=(*STPtr).gmax*(*STPtr).tc_on*(*STPtr).tc_off/((*STPtr).tc_off-(*STPtr).tc_on); 
    }
 
    
    /* Initialize synapses/connections (ConMtx, and information within NPtr) with SynPar according to SPList (called SPMtx in wrapper) */
    SynPar=mxGetPr(prhs[4]);
    if (mxGetM(prhs[4])!=NumSynPar) 
        mexErrMsgTxt("SynPar should have NumSynPar rows");

    numST=mxGetN(prhs[4]);
    SPList=mxGetPr(prhs[5]);
    if (mxIsEmpty(prhs[5])) 
        mexErrMsgTxt("SPMtx must not be empty");
    
    k=mxGetN(prhs[5]);
    if (mxGetM(prhs[5])!=N || (k%(N+M))!=0)
    	mexErrMsgTxt("SPMtx should be a matrix of dimensions Nunits x (Nunits+InpNunits) x NumSyn");

    MaxNumSyn=k/(N+M);
    
    /* Create ConMtx pointer */
    ConMtx0=(struct SynList **) mxMalloc(N*sizeof(struct SynList *));
    ConMtx0[0]=(struct SynList *) mxMalloc(N*(N+M)*sizeof(struct SynList));
    for (i=1;i<N;i++) 
        ConMtx0[i]=ConMtx0[i-1]+(N+M);
    
    for (i=0;i<N;i++)               /* loop over target neurons */
    {
    	for (j=0;j<N+M;j++)         /* loop over source neurons */
	    { 
            /* Determine number of synapses for each connection pair */
            ConMtx0[i][j].NumSyn=0;
            while ((int)SPList[i+j*N+ConMtx0[i][j].NumSyn*N*(N+M)]>0)  
            {
                ConMtx0[i][j].NumSyn++;
        	    if (ConMtx0[i][j].NumSyn>=MaxNumSyn) 
                    break; 
            }

            /* Allocate synapse pointer, if any */
            if (ConMtx0[i][j].NumSyn>0)
            {
                ConMtx0[i][j].Syn=(struct Synapse *) 
                mxMalloc(ConMtx0[i][j].NumSyn*sizeof(struct Synapse));
            }
	        else 
                ConMtx0[i][j].Syn=NULL;

            /* loop over synapse types */
	        for (k=0,SynLPtr=ConMtx0[i][j].Syn;k<ConMtx0[i][j].NumSyn;k++,SynLPtr++)  
	        {
                nst=(int)SPList[i+j*N+k*N*(N+M)]-1;        
        	    if (nst+1>numST || nst<0) 
                    mexErrMsgTxt("Illegal entry in SPMtx");

                
                /* Set up synaptic information for neurons and synapse types */
                if(j<N)         /********** Regular neurons **********/
                {
            	    InList=FALSE;
                    for (kk=0;kk<NPtr0[j].NumPreSyn;kk++)
                		if (nst==NPtr0[j].PreSynList[kk]) 
                        {
                    	    InList=TRUE; 
                            break; 
                        }
    	            (*SynLPtr).PreSynIdx=kk;
                    if (InList==FALSE) 
                    {
                		NPtr0[j].NumPreSyn++;
                        NPtr0[j].PreSynList=mxRealloc(NPtr0[j].PreSynList,NPtr0[j].NumPreSyn*sizeof(int));
                        NPtr0[j].PreSynList[kk]=nst;
                        NPtr0[j].SDf=mxRealloc(NPtr0[j].SDf,NPtr0[j].NumPreSyn*sizeof(struct SynDepr));
    		            NPtr0[j].SDf[kk].use=SynPar[1+NumSynPar*nst];
        	            NPtr0[j].SDf[kk].tc_rec=SynPar[2+NumSynPar*nst];
                        NPtr0[j].SDf[kk].tc_fac=SynPar[3+NumSynPar*nst]; 
                        for (k2=0;k2<SizeHistOutput;k2++) 
                            NPtr0[j].SDf[kk].Adepr[k2]=1.0;
                        NPtr0[j].SDf[kk].uprev[0]=SynPar[1+NumSynPar*nst];
                        NPtr0[j].SDf[kk].Rprev[0]=1.0; 
                    }

               	    STno=(int)SynPar[0+NumSynPar*nst]-1;
                    
                    if (STno+1>NumSynType || STno<0) 
                        mexErrMsgTxt("Illegal SType-entry in SynPar");
                    (*SynLPtr).STPtr=SynTPtr0+STno;
                    
                    (*SynLPtr).wgt=SynPar[4+NumSynPar*nst];

                    if (SynPar[5+NumSynPar*nst]==0)
                        mexErrMsgTxt("Error: Synaptic delay is zero.");
                    else
                        (*SynLPtr).dtax=SynPar[5+NumSynPar*nst];
                    
                    (*SynLPtr).p_fail=SynPar[6+NumSynPar*nst];
                        
                    InList=FALSE;
                    for (kk=0,STPtr2=NPtr0[i].STList;*STPtr2 && kk<NPtr0[i].NumSynType;STPtr2++,kk++)
                        if (*STPtr2==(*SynLPtr).STPtr) 
                            InList=TRUE;
                    if (NPtr0[i].NumSynType>=MaxNumSTperN) 
                        mexErrMsgTxt("NumSynType exceeds MaxNumSTperN");

                    if (InList==FALSE)     
                    { 
                        *STPtr2=(*SynLPtr).STPtr;
                        NPtr0[i].NumSynType++;
                        NPtr0[i].gfONsyn[kk]=0.0;
                        NPtr0[i].gfOFFsyn[kk]=0.0;
                    }
                }
                else            /*********** Input neurons ***********/
                {
            	    InList=FALSE;
                    for (kk=0;kk<InpNPtr0[j-N].NumPreSyn;kk++)
                		if (nst==InpNPtr0[j-N].PreSynList[kk]) 
                        {
                    	    InList=TRUE; 
                            break; 
                        }
    	            (*SynLPtr).PreSynIdx=kk;
                    if (InList==FALSE) 
                    {
                		InpNPtr0[j-N].NumPreSyn++;
                        InpNPtr0[j-N].PreSynList=mxRealloc(InpNPtr0[j-N].PreSynList,InpNPtr0[j-N].NumPreSyn*sizeof(int));
                        InpNPtr0[j-N].PreSynList[kk]=nst;
                        InpNPtr0[j-N].SDf=mxRealloc(InpNPtr0[j-N].SDf,InpNPtr0[j-N].NumPreSyn*sizeof(struct SynDepr));
    		            InpNPtr0[j-N].SDf[kk].use=SynPar[1+NumSynPar*nst];
        	            InpNPtr0[j-N].SDf[kk].tc_rec=SynPar[2+NumSynPar*nst];
                        InpNPtr0[j-N].SDf[kk].tc_fac=SynPar[3+NumSynPar*nst]; 
                        for (k2=0;k2<SizeHistOutput;k2++) 
                           InpNPtr0[j-N].SDf[kk].Adepr[k2]=1.0;
                        InpNPtr0[j-N].SDf[kk].uprev[0]=SynPar[1+NumSynPar*nst];
                        InpNPtr0[j-N].SDf[kk].Rprev[0]=1.0; 
                    }
                    
                    STno=(int)SynPar[0+NumSynPar*nst]-1;
                    if (STno+1>NumSynType || STno<0) 
                        mexErrMsgTxt("Illegal SType-entry in SynPar");
                    (*SynLPtr).STPtr=SynTPtr0+STno;
                    (*SynLPtr).wgt=SynPar[4+NumSynPar*nst];
                    (*SynLPtr).dtax=SynPar[5+NumSynPar*nst];
                    (*SynLPtr).p_fail=SynPar[6+NumSynPar*nst];
                    InList=FALSE;
                    for (kk=0,STPtr2=NPtr0[i].STList;*STPtr2 && kk<NPtr0[i].NumSynType;STPtr2++,kk++)
                        if (*STPtr2==(*SynLPtr).STPtr) 
                            InList=TRUE;
                    if (NPtr0[i].NumSynType>=MaxNumSTperN) 
                        mexErrMsgTxt("NumSynType exceeds MaxNumSTperN");

                    if (InList==FALSE)     
                    { 
                        *STPtr2=(*SynLPtr).STPtr;
                        NPtr0[i].NumSynType++;
                        NPtr0[i].gfONsyn[kk]=0.0;
                        NPtr0[i].gfOFFsyn[kk]=0.0;
                    } 
                }
            }
        }
    } 

    
    /* Set NoiseSyn */
    NoiseSyn.NumSyn = NumSynType;
    NoiseSyn.Syn=(struct Synapse *) mxMalloc(NoiseSyn.NumSyn*sizeof(struct Synapse));    

    for (i=0;i<N;i++)               /* loop over target neurons */
    {        
        for (j=0;j<NoiseSyn.NumSyn;j++)           /* loop over noise synapses */
        {            
            STno=(int)SynPar[0+NumSynPar*(numST-NoiseSyn.NumSyn+j)]-1;
            if (STno+1>NumSynType || STno<0)
                mexErrMsgTxt("Illegal SType-entry in SynPar");

            NoiseSyn.Syn[j].STPtr   = SynTPtr0+STno;
            NoiseSyn.Syn[j].wgt     = SynPar[4+NumSynPar*(numST-NoiseSyn.NumSyn+j)];
            NoiseSyn.Syn[j].dtax    = SynPar[5+NumSynPar*(numST-NoiseSyn.NumSyn+j)];
            NoiseSyn.Syn[j].p_fail  = SynPar[6+NumSynPar*(numST-NoiseSyn.NumSyn+j)];
            
            NPtr0[i].gfONnoise[j]=0.0;
            NPtr0[i].gfOFFnoise[j]=0.0;
        }
    }
    
    /* Read in NoiseDistr */
    NoiseDistr=mxGetPr(prhs[10]);
    NoiseStep = 1.0/(mxGetN(prhs[10])-1.0);
    
    
    /* Create SynExp pointers */
    SynExpOn=(double *) mxMalloc(NumSynType*sizeof(double));
    SynExpOff=(double *) mxMalloc(NumSynType*sizeof(double));    
    
        /* Set EvtMtx and EvtTimes */
        EvtMtx=mxGetPr(prhs[6]);
        EvtTimes=mxGetPr(prhs[7]);
        NumEvt=mxGetN(prhs[7]);
        if (NumEvt>0)
        { 
            if (mxGetM(prhs[7])!=2) mexErrMsgTxt("EvtTimes should have 2 rows with ev. times and durations");
            if (mxGetM(prhs[6])!=N || mxGetN(prhs[6])!=NumEvt)
	            mexErrMsgTxt("EvtMtx should be a matrix of dimensions N x NumEvt"); 
        }

        /* Set viewlist */
        ViewList=mxGetPr(prhs[8]);
        NumView=mxGetM(prhs[8])*mxGetN(prhs[8]);
        
        sprintf(IDNFilename,"IDN_%i.dat",UniquePrint);
        sprintf(IDN2Filename, "IDN2_%i.dat", UniquePrint);
        sprintf(IDN3Filename, "IDN3_%i.dat", UniquePrint);
        sprintf(IDN4Filename, "IDN4_%i.dat", UniquePrint);
        if (NumView>0) 
        {
            if ((fpOut=fopen(IDNFilename,"w"))==NULL) 
                mexErrMsgTxt("cannot open FILE IDN_#.dat");
            if ((fpOut2=fopen(IDN2Filename,"w"))==NULL) 
                mexErrMsgTxt("cannot open FILE IDN2_#.dat"); 
            if ((fpOut3=fopen(IDN3Filename,"w"))==NULL) 
                mexErrMsgTxt("cannot open FILE IDN3_#.dat");
	        for (i=0;i<NumView;i++)
	            if (ViewList[i]<1 || ViewList[i]>N) 
                    mexErrMsgTxt("Illegal entry in ViewList"); 
        }

        if (CtrPar[3]>NumVar)
        {
            CtrPar[3]=NumVar;
            mexWarnMsgTxt("Warning: ViewList parameters out of bound. Set to NumVar.");
        }
        
        /* Create output matrix and TnextSyn*/
        if (NumOutp>0)
        { 
            plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
            NumSynInp=mxGetPr(plhs[0]); 
            plhs[1]=mxCreateDoubleMatrix(N,1,mxREAL);
            N_osc=mxGetPr(plhs[1]); 
        }

        TnextSyn=dvector(0,N-1);
    
        
        
    /******************* Actual simulation *********************************/

    /* Time loop (without explicit increment) */
    for (t0=Tstart;t0<Tstop;) 
    {
        /* display actual time */
        if (t0>=t_display)            
        {
            mexPrintf("%f percent\n",t0*100/Tstop);
            mexEvalString("drawnow;");
            t_display=t0+100;
        }
            
        /* Set maximal time t1 until which the ODEs will be integrated */
    	t1=t0+dt0; EvtNo=-999;
    	if (t1>Tstop) 
            t1=Tstop;

        /* Set t1 to the next event time if it falls between t0 and t1 */
    	for (i=0;i<NumEvt;i++)
            if (EvtTimes[i*2]>t0 && EvtTimes[i*2]<=t1) 
            { 
                t1=EvtTimes[i*2]; NextEvtT=t1; EvtNo=i*2; 
            }
            else 
            { 
                EvtOffT=EvtTimes[i*2]+EvtTimes[i*2+1];
                    if (EvtOffT>t0 && EvtOffT<=t1) 
                    { 
                        t1=EvtOffT; NextEvtT=t1; EvtNo=i*2+1; 
                    } 
            }

        /* Set t1 to the next spike in any input neuron if it falls between t0 and t1 */
        t11 = t1;
        for (i=0;i<M;i++)       /* extract next spike in input neurons */
        {
            if(InpNPtr0[i].SPtrain[InpNPtr0[i].SP_ind]>t0 && InpNPtr0[i].SPtrain[InpNPtr0[i].SP_ind]<=t11)
            {
                t11 = InpNPtr0[i].SPtrain[InpNPtr0[i].SP_ind];
                print_flag = 1;
            }
            else
                print_flag = 0;
        }
        t1 = t11;

        /* evaluate spikes in all input neurons that fire at that time */
        for (i=0;i<M;i++)       
        {
            if(InpNPtr0[i].SPtrain[InpNPtr0[i].SP_ind]==t1)
            {
                /* Compute ISI from current and previous spike */
                if (InpNPtr0[i].SP_ind>0)
                    ISI_inp=t1-InpNPtr0[i].SPtrain[InpNPtr0[i].SP_ind-1]; 
                else
                    ISI_inp=10.0e8;

                /* Set new spike times */
                InpNPtr0[i].SpikeTimes[InpNPtr0[i].SP_ind]=InpNPtr0[i].SPtrain[InpNPtr0[i].SP_ind];
                InpNPtr0[i].SP_ind++;
                
                /* Update neuron-internal spike history (bounded by SizeHistOutput), */
                j=NumSpike[i+N] % SizeHistOutput;
                
                /* Loop over presynaptic synapses, update short-term synaptic plasticity variables */
                for (kk=0;kk<InpNPtr0[i].NumPreSyn;kk++)
                    if (InpNPtr0[i].SDf[kk].use>0.0)
                        InpNPtr0[i].SDf[kk].Adepr[j]=syndepr(&InpNPtr0[i].SDf[kk], ISI_inp, j);
                
                /* Write spike times to file (spike times, not ISIs, despite the file name!) */
                if (WriteST>0) 
                {
                    sprintf(uno,"%d",i+N);
                    sprintf(fn,"ISIu%s_%i.dat",uno,UniquePrint);
                    if ((fpISI=fopen(fn,"a"))==NULL)
                        mexErrMsgTxt("cannot open FILE ISIu#.dat");
                    fprintf(fpISI,"%lf\n",t1);
                    fclose(fpISI); 
                 }
                
                NumSpike[i+N]++;
            }   
        }
        
        /* Neuron loop */
	    for (i=0;i<N;i++) 
        {
            /* Define integration start time t0_i for the current neuron */
	        t0_i=t0;
            
            /* Integrate up of t1_i until t0_i >= t1 */
            while (t0_i<t1) 
            {
                /* Define integration stopping time t1_i for the current neuron */
		        t1_i=t1;

                /* Set t1_i to TnextSyn if it falls between t0_i and t1_i */
        		if (TnextSyn[i]>t0_i && TnextSyn[i]<t1_i) 
                    t1_i=TnextSyn[i];
                
                /* Set dynamic time step dt and update ODEs accordingly */
        		vp=NPtr0[i].v[0];
                wp=NPtr0[i].v[1];
		        dt=t1_i-t0_i;  

                /* the following section is necessary if one wants to restrict the depolarization block to time x after a spike */
                if (NumSpike[i]>0)
                {
                    if ((t0_i-NPtr0[i].SpikeTimes[(NumSpike[i]-1) % SizeHistOutput])<5)
                        flag_dv=0;
                    else
                        flag_dv=1;
                }
                else
                {
                    flag_dv=1;} 
                /* the section above is necessary if one wants to restrict the depolarization block to time x after a spike */

                
                /* Update differential equations */
		        update(&NPtr0[i],dt);
                
                if (stop_flag>0)
                {
                    mexPrintf("%f %d %f %f\n", t0_i, i, vp, wp);
                    mexErrMsgTxt("nan error");
                }
                
                for (j=0;j<NPtr0[i].NumSynType;j++) {
                    if (NPtr0[i].gfONsyn[j]<0 | NPtr0[i].gfOFFsyn[j]<0) {
                        mexPrintf("%d %d %f %f %f %f\n", i, j, t0_i, t1_i, NPtr0[i].gfONsyn[j], NPtr0[i].gfOFFsyn[j]);
                        mexErrMsgTxt("g is negative");
                    }
                }
                
                if(t1_i==t1)
                {
                    gsyn1[i]=gsyn_AN;
                    gsyn2[i]=gsyn_G;
                    Isyn[i]=I_tot;
                    if(I_tot<NPtr0[i].I_ref*1.01 && I_tot>NPtr0[i].I_ref*0.99)
                        flag_osc[i]++;
                    else
                        flag_osc[i]=0;
                }
                
                /* check if input oscillates around I_ref */
                if(flag_osc[i]>=200 && NumOutp>0) /*after 10 ms */ 
                    N_osc[i]++;

                
                /* Handle spikes (defined by positive crossings of Vup) */
                if ((NPtr0[i].v[0]>=NPtr0[i].Vup) && (vp<NPtr0[i].Vup)) 
                {
                    /* Interpolate to get exact spike time t1_i */
		            t1_i=t0_i+dt*(NPtr0[i].Vup-vp)/(NPtr0[i].v[0]-vp);

                    /* Compute ISI from current and previous spike */
                    if (NumSpike[i]>0) 
                        ISI=t1_i-NPtr0[i].SpikeTimes[(NumSpike[i]-1) % SizeHistOutput];
                    else 
                        ISI=10.0e8;
                    
                    
                    /* Enforce refractory period of 5 ms */
                    if (ISI > 5)
                    {
                        /* reset V and w */
                        w_Vup=wp + ((NPtr0[i].v[1]-wp)/dt) * (t1_i-t0_i);
                        NPtr0[i].v[0]=NPtr0[i].Vr;
                        NPtr0[i].v[1]= w_Vup + NPtr0[i].b;
                        
                        /* Update neuron-internal spike history (bounded by SizeHistOutput), */
                        j=NumSpike[i] % SizeHistOutput;
                        NPtr0[i].SpikeTimes[j]=t1_i;
                                                
                        /* Loop over presynaptic synapses, update short-term synaptic plasticity variables */
                        for (kk=0;kk<NPtr0[i].NumPreSyn;kk++)
                            if (NPtr0[i].SDf[kk].use>0.0)
                                NPtr0[i].SDf[kk].Adepr[j]=syndepr(&NPtr0[i].SDf[kk],ISI,j);
                        
                        /* Write spike times to file (spike times, not ISIs, despite the file name!) */
                        if (WriteST>0)
                        {
                            sprintf(uno,"%d",i);
                            sprintf(fn, "ISIu%s_%i.dat",uno,UniquePrint);
                            if ((fpISI=fopen(fn,"a"))==NULL)
                                mexErrMsgTxt("cannot open FILE ISIu#.dat");
                            fprintf(fpISI,"%lf\n",t1_i);
                            fclose(fpISI);
                        }
                        
                        /* Increment number of spikes */
                        NumSpike[i]++;
                        
                        /* reset time step */
                        dt=t1_i-t0_i;
                    }
                    else
                    {
                        NPtr0[i].v[0]=vp;
                        NPtr0[i].v[1]=wp;
                        if (ISI > 5)
                            mexErrMsgTxt("Warning: ISI>5");
                        if (ISI < 0)
                        {
                            mexPrintf("%d %f %f %f %f\n",i,ISI,t1_i,NPtr0[i].SpikeTimes[NumSpike[i] % SizeHistOutput],NPtr0[i].SpikeTimes[(NumSpike[i]-1) % SizeHistOutput]);
                            mexErrMsgTxt("Warning: ISI<0");
                        }
                        
                        /* reset t1_i */
                        t1_i=dt+t0_i; 
                    }
                        
                    /* recompute synaptic conductances */
                    if(t1_i==t1)
                    {
                        set_gsyn(NPtr0[i], dt, vp);
                        gsyn1[i]=gsyn_AN;
                        gsyn2[i]=gsyn_G;
                        Isyn[i]=I_tot;
                    }
		        }

                
		        /* Application of synaptic noise (deactivated at the moment) */
                for (j=0;j<NoiseSyn.NumSyn;j++)
                {
                    SynExpOn[j] = exp(-dt/(*NoiseSyn.Syn[j].STPtr).tc_on);
                    SynExpOff[j] = exp(-dt/(*NoiseSyn.Syn[j].STPtr).tc_off);                    
                    rand_num = NoiseDistr[(int) round((double)rand()/(double)RAND_MAX/NoiseStep)];

                    NPtr0[i].gfONnoise[j]  = 0.0;
                    NPtr0[i].gfOFFnoise[j] = 0.0;
                    /*NPtr0[i].gfOFFnoise[j] = NoiseSyn.Syn[j].wgt + (NPtr0[i].gfOFFnoise[j]-NoiseSyn.Syn[j].wgt)*(SynExpOff[j]-SynExpOn[j]) + NoiseSyn.Syn[j].dtax*sqrt(1-(SynExpOff[j]-SynExpOn[j])*(SynExpOff[j]-SynExpOn[j]))*rand_num;*/
                }

                
                /* Decay of synaptic conductances */
		        for (j=0;j<NPtr0[i].NumSynType;j++) 
                {
		            NPtr0[i].gfONsyn[j]*=exp(-dt/(*NPtr0[i].STList[j]).tc_on);
		            NPtr0[i].gfOFFsyn[j]*=exp(-dt/(*NPtr0[i].STList[j]).tc_off); 
                    
                }
                
                
    		    /* Process synaptic events (incoming spikes) */
	    	    TnextSyn[i]=Tstop+100.0;
                
                /* Loop over presynaptic neurons*/
		        for (j=0;j<N;j++)
                {
                    /* Loop over synapses between the current pre/postsynaptic neuron pair */
		            for (k=0,SynLPtr=ConMtx0[i][j].Syn;k<ConMtx0[i][j].NumSyn;k++,SynLPtr++)
                    {
                        /* Loop over spikes in the presynaptic neuron (up to SizeHistOutput) */
			            for (kk=NumSpike[j]-1;kk>=0 && (NumSpike[j]-kk)<=SizeHistOutput;kk--)
                        {
                
                            /* Jump to next spike if current spike arrives at current neuron before t0_i */
			                if (t0_i>=(NPtr0[j].SpikeTimes[kk % SizeHistOutput]+(*SynLPtr).dtax)) 
                                break;
			                else 
                            { 
                                /* Use current spike to update synapse if it arrives at the current neuron before t1_i, and if no synaptic failure occurs */
				                if ((t1_i>=NPtr0[j].SpikeTimes[kk % SizeHistOutput]+(*SynLPtr).dtax) && ((double)rand()/(double)RAND_MAX > (*SynLPtr).p_fail))
                                {
				                    for (k2=0;k2<NPtr0[i].NumSynType;k2++)
					                    if (NPtr0[i].STList[k2]==(*SynLPtr).STPtr) 
                                        {
					                        Aall=NPtr0[j].SDf[(*SynLPtr).PreSynIdx].Adepr[kk % SizeHistOutput]*(*SynLPtr).wgt*(*(*SynLPtr).STPtr).Gsyn;
					                        NPtr0[i].gfONsyn[k2]+=Aall;
					                        NPtr0[i].gfOFFsyn[k2]+=Aall; 
					                        if (NumOutp>0) 
                                                NumSynInp[i]=NumSynInp[i]+1.0; 
                                        } 
                                }
				                else
                                    /* If arrival time of next spike is before TnextSyn[i], set TnextSyn[i] to it */
				                    if (NPtr0[j].SpikeTimes[kk % SizeHistOutput]+(*SynLPtr).dtax<TnextSyn[i])
					                    TnextSyn[i]=NPtr0[j].SpikeTimes[kk % SizeHistOutput]+(*SynLPtr).dtax;            
                            }
                        }
                    }
                }

                
                /* Loop over presynaptic input neurons */
                for(j=N;j<N+M;j++)
                {
                    /* Loop over synapses between the current pre/postsynaptic neuron pair */
		            for (k=0,SynLPtr=ConMtx0[i][j].Syn;k<ConMtx0[i][j].NumSyn;k++,SynLPtr++)
                    {
                        /* Loop over spikes in the presynaptic neuron */
			            for (kk=NumSpike[j]-1;kk>=0;kk--)
                        {
                            /* Jump to next spike if current spike arrives at current neuron before t0_i */
			                if (t0_i>=(InpNPtr0[j-N].SpikeTimes[kk]+(*SynLPtr).dtax)) 
                                break;
			                else 
                            { 
                                /* Use current spike to update synapse if it arrives at the current neuron before t1_i and if no synaptic failure occurs */
				                if ((t1_i>=InpNPtr0[j-N].SpikeTimes[kk]+(*SynLPtr).dtax) && ((double)rand()/(double)RAND_MAX > (*SynLPtr).p_fail))
                                {
				                    for (k2=0;k2<NPtr0[i].NumSynType;k2++)
					                    if (NPtr0[i].STList[k2]==(*SynLPtr).STPtr) 
                                        {
                                            Aall=InpNPtr0[j-N].SDf[(*SynLPtr).PreSynIdx].Adepr[kk % SizeHistOutput]*(*SynLPtr).wgt*(*(*SynLPtr).STPtr).Gsyn;
					                        NPtr0[i].gfONsyn[k2]+=Aall;
					                        NPtr0[i].gfOFFsyn[k2]+=Aall; 
                                            if (NumOutp>0) 
                                                NumSynInp[i]=NumSynInp[i]+1.0; 
                                        } 
                                }
				                else
                                    /* If arrival time of next spike is before TnextSyn[i], set TnextSyn[i] to it */
				                    if (InpNPtr0[j-N].SpikeTimes[kk]+(*SynLPtr).dtax<TnextSyn[i])
					                    TnextSyn[i]=InpNPtr0[j-N].SpikeTimes[kk]+(*SynLPtr).dtax; 
                            }
                        }
                    }
                }
                
                /* Update t0_i */
		        t0_i=t1_i;
	        } 
        }

        /* save first set of variables for neurons in viewlist */
	    for (i=0;i<NumView;i++) 
        {
    	    fprintf(fpOut,"%lf %d",t1,(int)ViewList[i]);
	        for (k=0;k<(int)CtrPar[3];k++)
		        fprintf(fpOut," %lf",NPtr0[(int)ViewList[i]-1].v[k]);
        	    fprintf(fpOut,"\n");
        }
                 
        /* save second set of variables for neurons in viewlist */
	    for (i=0;i<NumView;i++) 
            fprintf(fpOut2, " %f %f %f", gsyn1[(int)ViewList[i]-1], gsyn2[(int)ViewList[i]-1], Isyn[(int)ViewList[i]-1]);
        fprintf(fpOut2, "\n");
        
        
    	/* application of events (current injections) */
	    if (EvtNo>=0 && t1>=NextEvtT)
    	{ 
            if ((EvtNo % 2)==0)
	            for (i=0;i<N;i++) 
                    NPtr0[i].Iinj=EvtMtx[i+(EvtNo/2)*N];
	        else
	            for (i=0;i<N;i++) 
                    NPtr0[i].Iinj=0.0; 
        }

        /* Update t0 */
	    t0=t1;
    }

    /* Close output files and free allocated memory */
    if (NumView>0) 
    {
        fclose(fpOut);
        fclose(fpOut2);
        fclose(fpOut3);
        
    }
    free_dvector(TnextSyn,0,N-1);
    free_ivector(NumSpike,0,N-1);
    free_dvector(gsyn1,0,N-1);
    free_dvector(gsyn2,0,N-1);
    free_dvector(Isyn,0,N-1);

    return;
}


/*
% (c) 2016 J. Hass, L. Hertaeg and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
*/
