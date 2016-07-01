/*                ***********NOTICE**********                              */
/*                                                                         */
/*  The EMAP finite element modeling codes were created at the             */
/*  University of Missouri-Rolla Electromagnetic Compatibility Laboratory. */
/*  They are intended for use in teaching or research.  They may be freely */
/*  copied and distributed PROVIDED THE CODE IS NOT MODIFIED IN ANY MANNER.*/
/*                                                                         */
/*  Suggested modifications or questions about these codes can be          */
/*  directed to Dr. Todd Hubing, Department of Electrical Engineering      */
/*  University of Missouri-Rolla, Rolla, MO  65401.  Principal authors     */
/*  of these codes are Mr. Mohammad Ali and Mr. Girish Bhat of the         */
/*  University of Missouri-Rolla.                                          */
/*                                                                         */
/*  Neither the authors nor the University of Missouri makes any warranty, */
/*  express or implied, or assumes any legal responsibility for the        */
/*  accuracy, completeness or usefulness of these codes or any information */
/*  distributed with these codes.                                          */
/*                                                                         */
/* 5/98                                                                    */
/***************************************************************************
*   PROGRAM     :  Emap4.c (Version 1.0)
*   LAST UPDATE :  May 1998 
*   DESCRIPTION :  A 3-D Vector Finite Element Modeling Code for Analyzing
                   Time Varying Complex Electromagnetic Fields.
*   LATEST DEV. :  Inclusion of Isource, modeling of coax cables.
*   INPUT FILE  :  emap4.in or or 1st argument of emap command
*   OUTPUT FILE :  file(s) specified by input file
*   USER'S GUIDE:  http://www.emclab.umr.edu/emap4
*                  (contains all the relevant references)
*
****************************************************************************/


/**************************** Include Files ********************************/

#include <stdio.h>
#include <signal.h>
#include <string.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>


/***********************complex_function*******************************/
struct complex { float x,y; };
typedef struct complex complex;
complex COMplex_Cmplx();
complex COMplex_Conjg();
complex COMplex_Add();
complex COMplex_Add2();
complex COMplex_Sub();
complex COMplex_Mul();
complex Real_Mul();
complex COMplex_Div();
complex Real_Div();
complex COMplex_Null();
float  Real();
float  Aimag();
float  COMplex_Abs();
complex COMplex_Expon();
 
complex COMplex_Cmplx(float x,float y)
{
    complex z; z.x = x; z.y = y; return z;
}
complex COMplex_Conjg(complex z)
{
    return COMplex_Cmplx(z.x, -z.y);
}
complex COMplex_Add(complex a,complex b)
{
    return COMplex_Cmplx(a.x + b.x,a.y + b.y);
}
complex COMplex_Add2(complex a,complex b,complex c)
{
    return COMplex_Cmplx(a.x + b.x + c.x,a.y + b.y + c.y);
}
complex COMplex_Sub(complex a,complex b)
{
    return COMplex_Cmplx(a.x - b.x,a.y - b.y);
}
complex COMplex_Mul(complex a,complex b)
{
    return COMplex_Cmplx(a.x*b.x - a.y*b.y,a.x*b.y + a.y*b.x);
}
complex Real_Mul(float a,complex z)
{
    return COMplex_Cmplx(a*z.x,a*z.y);
}
complex COMplex_Div(complex a,complex b)
{
    float D = b.x*b.x + b.y*b.y;
    return COMplex_Cmplx((a.x*b.x + a.y*b.y)/D,(a.y*b.x - a.x*b.y)/D);
}
complex Real_Div(complex z,float a)
{
    return COMplex_Cmplx(z.x/a,z.y/a);
}
float Real(complex z)
{
    return z.x;
}
float Aimag(complex z)
{
    return z.y;
}
float COMplex_Abs(complex z)
{
    return sqrt(z.x*z.x + z.y*z.y);
}
complex COMplex_Expon(float a,complex z)
{
    float R = exp((a*z.x));
    return COMplex_Cmplx(R*cos(z.y),(a*R)*sin(z.y));
}
complex COMplex_Null()
{
    complex z; z.x=0.0; z.y=0.0; return z;
}

/*******dynamic memory allocation**********************/
 
complex *CMPLX_Vector(int nh)
/* Allocates a complex vector with range [0...nh-1] */
{
	complex  *v;
	v=(complex *)malloc((unsigned) nh*sizeof(complex));
	if (v==NULL) {printf("Unable to allocate memory");}
	return v;
}
float *FLOAT_Vector(int nh)
/* Allocates a float vector with range [0 ... nh-1] */
{
	float *v;
	v=(float *)malloc((unsigned) nh*sizeof(float));
	return v;
}
int *INT_Vector(int nh)
/* Allocates an int vector with range [0 ... nh-1]  */
{
	int *v;
	v=(int *)malloc((unsigned) nh*sizeof(int));
	return v;
}
float **FLOAT_Matrix(int row,int column)
/*Allocates a float matrix with range [0..rw-1][0..cl-1]*/
{
	int i;float **m;
	m=(float **) malloc((unsigned) row*sizeof(float*));
	for(i=0;i<row;i++)
	{m[i]=(float *) malloc((unsigned) column*sizeof(float));}
	return m;
}
int **INT_Matrix(int row,int column)
/*Allocates an int matrix with range [0..rw-1][0..cl-1]*/
{
	int i;int **m;
	m=(int **) malloc((unsigned) row*sizeof(int));
	for(i=0;i<row;i++)
	{m[i]=(int *) malloc((unsigned) column*sizeof(int));}
	return m;
}
complex **CMPLX_Matrix(int row,int column)
/*Allocates a complex matrix with range [0..rw-1][0..cl-1] */
{
	int i;complex **m;
	m=(complex **) malloc((unsigned) row*sizeof(complex*));
	for(i=0;i<row;i++)
	{m[i]=(complex *) malloc((unsigned) column*sizeof(complex));}
	return m;
}

float ***FLOAT_Matrix2(int row,int column,int width)
{
	int i,j;float ***m;
	m=(float ***) malloc((unsigned) row*sizeof(float**));
	for(i=0;i<row;i++)
	{
		m[i]=(float **) malloc((unsigned) column*sizeof(float*));
		for(j=0;j<column;j++)
		{
			m[i][j]=(float *) malloc((unsigned) width*sizeof(float));
		}
	}
	return m;
}

int ***INT_Matrix2(int row,int column,int width)
{
	int i,j;int ***m;
	m=(int ***) malloc((unsigned) row*sizeof(int**));
	for(i=0;i<row;i++)
	{
		m[i]=(int **) malloc((unsigned) column*sizeof(int*));
		for(j=0;j<column;j++)
		{
			m[i][j]=(int *) malloc((unsigned) width*sizeof(int));
		}
	}
	return m;
}

void free_INT_Vector(int* v,int nh)
{
	free((char *) v);
}

void free_CMPLX_Vector(complex* v,int nh)
{
	free((char *) v);
}

void free_FLOAT_Vector(float* v,int nh)
{
	free((char *) v);
}

void free_INT_Matrix(int** m,int row,int column)
{
	int i;
	for(i=row-1;i>=0;i--)
		free((char*) m[i]);
	free((char*) m);
}

void free_FLOAT_Matrix(float** m,int row,int column)
{
	int i;
	for(i=row-1;i>=0;i--)
		free((char*) m[i]);
	free((char*) m);
}

void free_CMPLX_Matrix(complex** m,int row,int column)
{
	int i;
	for(i=row-1;i>=0;i--)
		free((complex*) m[i]);
	free((complex*) m);
}

/************************** Default Parameters *****************************/

#define  GBL_Matrix_MaxColNos  19  /* The FEM is sparse and any row has a 
				      maximum of 19 columns */
#define  FRCD_Matrix_MaxColNos 11  
#define  MAX_ITERATIONS 25000      /* The maximum number of iterations for 
				      the CBCG solver */
#define  TOLERANCE  0.00001        /* The amount of accuracy required in the 
				      code */
#define M_PI 3.14;
/************************** Function Prototypes *****************************/

void Read_Input_Pass_1();
void ParameterInfo();  
void Read_Input_Pass_2();
void AssignGlobalCoorDinates();
void AssignHexHedronNodeNumbering();
void AssignHexHedronEdgeNumbering();
void HexHedraSubDiv1();
void HexHedraSubDiv2();
void EdgeSubDiv1();
void EdgeSubDiv2();
void FindTetraHedronVolume();
void ComputeTetraCoFactor();
void S_Matrix();
void T_Matrix();
void TetraSubMatrix();
void FindGlobalMatrix();
void CountHalfBandWidth();
void PartitionGlobalMatrix();
void ConjugateSolver(int MatrixSize,int* RowIIndex,int** RowICol,complex** RowIDat,complex* RHSVec,int ITmax,float TOL);
void GlobalEdgeEndsDiv1();
void GlobalEdgeEndsDiv2();
void EfieldatNode(int X,int Y,int Z,FILE* Out_F);
void Produce_Output();

/*************************** Global Variables *******************************/

/* The following are the file variables */

FILE *InF, *OutF_0, *OutF_1, *OutF_2, *OutF_3;

FILE * A_matrix,*b_new_matrix,*Ap_values , *Multiples_Matrix, *col_nos_matrix ,*Parameters_matrix , *B_matrix , *X_matrix;
FILE * A_floating_matrix , * B_floating_matrix , *b_new_floatin_matrix ;

int 

TotBoundNodeNum, /* The number of boundary nodes */
TotBoundEdgeNum, /* The number of unforced boundary edges */
TotForcdEdgeNum, /* The number of forced edges */
TotTrngleNum, /* The total number of triangles on the surface (not used) */
TotInnerEdgeNum, /* The number of unforced or "inner" edges */
start, /* Dummy variable used to define the start in edge and 
	  node numbering schemes */
TotNumHexHedron, /* The total number of hexahedra */
t, /* Dummy variable used to specify the tetrahedron number */
XdiM,YdiM,ZdiM, /* The dimensions of the meshed area */
Max_X=0, /* Start and end coordinates of the meshed area */
Max_Y=0, /* The mesh goes from (Min_X,Min_Y,Min_Z) to (Max_X,Max_Y,Max_Z) */
Max_Z=0,
Min_X=0,
Min_Y=0,
Min_Z=0,
TotEdgeNum=0, /* Total number of edges in the mesh */
TotNodeNum=0, /* Total number of nodes in the mesh */
HexNum, /* Index Variable used to specify the hexahedron number */
TotTetElement; /* Total number of tetrahedral elements in the mesh */

int freq_step = 0; /* The variable to indicate the number of the 
		      frequency step in progress */

int

  IEdge, /* Dummy edge index */
  edge_end[6][2], /* temporary variable while calculating the edge-ends */
**HexNode, /* Node numbers : 1st index is the hexahedron number
	       The second index is the local hexahedron node number */
**HexEdgeNum,/* Edge numbers : 1st index is the hexahedron number
	       The second index is the local hexahedron edge number */
**TetGlobalNodeNum,
/* This is a [5:4] matrix. 
		       Node numbers : 1st index is the tetrahedron number [0:4]
		       The second index is the local tetrahedron node number */
**TetEdgeNum, /* Edge numbers : 1st index is the tetrahedron number
	       The second index is the local tetrahedron edge number */
**TetGlobalEdgeEnds,/* Edge Ends : 1st index is the gloal edge number
	       The second index is 0 (beginning) or 1 (end) */
*ForceStat, /* The flag which indicates whether the edge is forced (1) or
	       an unforced edge (0) */
*InnerEdgeStat, /* Contains TotalInnerEdgeNum entries containing the 
		   global edge number of each inner edge */
*BoundEdgeStat, /* Contains TotalBoundEdgeNum entries containing the 
		   global edge number of each unforced boundary edge */
*ForcdEdgeStat, /* Contains TotalForcdEdgeNum entries containing the 
		   global edge number of each forced edge */
*InnerEdgeStat1,/* Contains TotEdgeNum entries containing the 
		   inner edge number of each global edge if applicable
		   otherwise is zero */
*BoundEdgeStat1,/* Contains TotEdgeNum entries containing the 
		   boundary  edge number of each global edge if applicable
		   otherwise is zero */
*ForcdEdgeStat1,/* Contains TotEdgeNum entries containing the 
		   forced edge number of each global edge if applicable
		   otherwise is zero */
/* The above three variables exist for ease and simplicity of the program */

/* The following variables store the Global and the 
   actual matrix to be solved 
   (see chapter on sparse matrices in the user's guide) */
**GBL_Matrix_ColNos,
**LHS_Matrix_ColNos,
**FORC_Matrix_ColNos,
*GBL_Matrix_Index,
*LHS_Matrix_Index;

complex
**GBL_Matrix_Data,
**LHS_Matrix_Data;

complex

A_mat[5][6][6], /* Temp. variable to store the A_mat for each tetrahedron */
S_mat[6][6],/* Temp. variable to store the S_mat for each tetrahedron */
T_mat1[6][6],/* Temp. variable to store the T_mat for each tetrahedron */
*Resist, /* Variable to store the terms arising due to resistors */
*Isource, /* Variable to store the terms arising due to current sources */
*RHSVector, /* The RHS vector in the final equation which is solved */
*ForcdValue, /* The value of the electric fields along the "forced edges */
*RELPerm, /* The relative permittivity of each hexahedron */
*EfieldData; /* The value of the electric field along every edge in the mesh */




int
**xn,
NUM_OF_STEPS=1; /* The number of frequency steps that this code is run */

 
float

*DivisorX, /* This variable stores the value of cell dimension 
	      (length of each brick in the x-direction */
*DivisorY, /* This variable stores the value of cell dimension 
	      (length of each brick in the y-direction */
*DivisorZ, /* This variable stores the value of cell dimension 
	      (length of each brick in the z-direction */
OperateFreq, /* Frequency of operation 
		(the code can run only for one frequency at a time */
FreeSpaceVel, /* "c" */
AbsPermeable, 
AbsPermitt,
WaveLength, /* the free space wavelength */
WaveNumber, /* the free-space wavenumber */
Omega, /* "w"=2*pi*"operatefreq" */
TetVolume, /* Dummy variable to store the volume of each tetrahedron */
CoFactor[4][4], /* Variable which stores the "co-factor matrix" 
		   (the a's, b;s c's and d's.
		   for each tetrahedron (see reference [9] in the 
		   user's guide on the web "http://www.emclab.umr.edu/emap4*/
**NodeCord, /* Stores the node coordinates for each global node 
	       the first index is the global node number, the second index
	       is 0 (x-coord), 1 (y-coord) or 2 (z_coord) */
*Sigma, /* Stores the value of the conductivity of each brick (hexahedron)
	   if applicable, otherwise stores 0 */
*Epsilon, /* Stores the value of the relative permittivity of each brick 
	     (hexahedron) if applicable, otherwise stores 0 */
*EdgeStatus, /* This variable is a vector of length [TotEdgeNum]. 
		If the edge is an inner edge, it stores -1, 
		if it is an unforced boundary edge, it stores -2 
		and if it is a forced edge, it stores the 
		magnitude specified in the input file */
*EdgeStatus1, /* This variable is the same as EdgeStatus and 
		 stores the phase specified in the input file */
FREQ_INC;

complex   /* The variables used in */ 
*X,       /* the conjugate solver  */
*P,       /* routine (see the sparse matrix chapter in the 
	     user's guide for more 
	     information about the variables */
*PP,      
*R,
*RR,
*A_P,
*AH_PP,
Alpha,
Beta,
temp_var;

complex 
*Old_Solution; /* Used in the conjugate solver to feed the previous solution
		  as the starting vector for the next frequency in a 
		  swept frequency simulation */


complex
a_PML,b_PML,c_PML, /* THe values of the PML variables a,b,c inside the specific
		      hexahedron for which the matrix is being computed*/
*a_PML_vec,*b_PML_vec,*c_PML_vec; /* THe values of the PML variables a,b,c 
				     inside all the hexahedra in the PML*/
int *Is_PML; /* Variable which asks the question 
		"Is this hexahedron part of a PML layer?" */ 

char 
Ifile[20], /* Input filename */ 
Out_FileName0[20], /* Output filenames */
Out_FileName1[20], 
Out_FileName2[20], 
Out_FileName3[20];
/***************************************************************************
*  Utility Functions
****************************************************************************/

/* Useful functions when dealing with matrices and vectors 
   and manupulate them */ 

 void VTXsub1(int Count_i,int Count_k,float Buff[4][3],float Buff1[3])

/* subtracts two vectors of the same matrix from each other */
 {
	 int Count_j;
	 for(Count_j=0;Count_j<=2;Count_j++)
	 Buff1[Count_j] = Buff[Count_i][Count_j]-Buff[Count_k][Count_j];
 }
 void VTXcross(float Buff1[3],float Buff2[3],float Buff[3]) 

/* Vector cross product */ 
 {
	 Buff[0] = Buff1[1]*Buff2[2] - Buff2[1]*Buff1[2];
	 Buff[1] = Buff1[2]*Buff2[0] - Buff2[2]*Buff1[0];
	 Buff[2] = Buff1[0]*Buff2[1] - Buff2[0]*Buff1[1];
 }

 float Sign(int Value) /*sign */
 {
	 float Value1,Value2;
	 Value1 = (float)Value ; Value2 = Value1 / fabs(Value1) ;
	 return Value2;
 }

 float VTXmag(float Buff1[3],float Buff2[3]) /* length of vector */
 {
	 float Value,ValueX,ValueY,ValueZ;
	 ValueX = Buff1[0] - Buff2[0];
	 ValueY = Buff1[1] - Buff2[1];
	 ValueZ = Buff1[2] - Buff2[2];
	 Value  = sqrt(ValueX*ValueX + ValueY*ValueY + ValueZ*ValueZ);
	 return Value;
 }

/* calculates the "F" terms used in the FEM matrix 
   (see appendices A and B in reference [9] or see reference [1])
   These terms are slightly different when PMLs are used */
float ff(int i,int j)
{
   float Prod;
   Prod = CoFactor[1][i-1] * CoFactor[1][j-1] 
	+ CoFactor[2][i-1] * CoFactor[2][j-1]
        + CoFactor[3][i-1] * CoFactor[3][j-1]  ;
              return Prod;
}

/* ff1 is used when PML matrix terms are calculated. 
   The PMLs multiply an anisotropic dielectric to the matrix dot-product. 
   So you can see that for each term, the multiplying factor is different */ 

   complex ff1(int i,int j)
   {
	   complex Prod;
	   Prod = COMplex_Add2( Real_Mul(CoFactor[1][i-1]*CoFactor[1][j-1],a_PML),
		  		Real_Mul(CoFactor[2][i-1]*CoFactor[2][j-1],b_PML),
				Real_Mul(CoFactor[3][i-1]*CoFactor[3][j-1],c_PML))  ;
	              return Prod;
   }

/* Searches for non-zero entries in the row and column specified */

int Srch_NZ_Element_Column(int RowNum,int ColNum)
{
   int Count_i;
   {
	   for(Count_i=0;Count_i<=GBL_Matrix_Index[RowNum]-1;Count_i++)
	   {
		   if(ColNum==(GBL_Matrix_ColNos[RowNum][Count_i]))
		   {
			   return Count_i;
		   }
	   } return -1;
   }
 }

int Locate_NZ_Element_Column(int RowNum,int ColNum)
{
	int Count_i,test_i,test_j;
	{
		for(Count_i=0;Count_i<=FRCD_Matrix_MaxColNos-1;Count_i++)
		{
			if(ColNum==(FORC_Matrix_ColNos[RowNum][Count_i]))
			{
				return Count_i;
			}
		} return -1;
	}
}

/**************************** MAIN PROGRAM *********************************/

int main(int argc,char* argv[])

/* MAIN PROGRAM */

{
int con,i,j,k;
float FreeSpaceVel,AbsPermeable,AbsPermitt,WaveLength;

/*Print error messages if the command line is not correct */
        if (argc > 2)
		{
			fprintf(stderr,"Usage: Emap4 [input_file] \n");
			getche();
			exit(1);
		}
        if (argc < 2)
		{
			argv[1] = "emap4.in";
			strcpy(Ifile, argv[1]);
		}
        else
			strcpy(Ifile, argv[1]);

        fprintf(stdout, "\nEMAP4\n\n");
        fprintf(stdout, "The input will be read from the file:  %s\n", Ifile);
        InF = fopen(Ifile, "r");
        if (InF == NULL)
		{
			fprintf(stderr, " Error: Can't open input file.\n");
			getche();
			exit(1);
		}
        
/* Pre-processing begins */

	Read_Input_Pass_1(); /* check for errors and inconsistencies 
				in the input file */
	ParameterInfo(); /* Take in the values of many global parameters */

/* A little bit of initialization */
	//HexNode = INT_Matrix(TotNumHexHedron,8);
	//HexEdgeNum = INT_Matrix(TotNumHexHedron,19);
	RELPerm     = CMPLX_Vector(TotNumHexHedron);
	NodeCord = FLOAT_Matrix(TotNodeNum,3);
	TetGlobalNodeNum  = INT_Matrix(5,4);
	TetEdgeNum  = INT_Matrix(5,6);
	//TetGlobalEdgeEnds = INT_Matrix(TotEdgeNum,2);

	HexNode = new int*[TotNumHexHedron];
	for(int kkk = 0; kkk < TotNumHexHedron; ++kkk)
		HexNode[kkk] = new int[8];

	HexEdgeNum = new int*[TotNumHexHedron];
	for(int kkk = 0; kkk < TotNumHexHedron; ++kkk)
		HexEdgeNum[kkk] = new int[19];

	TetGlobalEdgeEnds = new int*[TotEdgeNum];
	for(int kkk = 0; kkk < TotEdgeNum; ++kkk)
		TetGlobalEdgeEnds[kkk] = new int[2];


	for(int i=0;i<TotNumHexHedron;i++)
	{ 
		for(int j=0;j<8;j++)
		{
			HexNode[i][j] = 0;
		}
		for(int j=0;j<19;j++)
		{
			HexEdgeNum[i][j] = 0;
		}
	RELPerm[i] = COMplex_Null();
	}
	
	
	for(i=0;i<TotNodeNum;i++)  { 
	for(j=0;j<3;j++) {NodeCord[i][j] = 0.0;}
	                           }
	for(i=0;i<TotEdgeNum;i++)  { 
	for(j=0;j<2;j++) {TetGlobalEdgeEnds[i][j] = 0;}
				}
	for(i=0;i<5;i++)  { 
	for(j=0;j<4;j++) {TetGlobalNodeNum[i][j] = 0;}
	for(j=0;j<6;j++) {TetEdgeNum[i][j] = 0;}
				}
/* Mesh Generation begins*/	

	Read_Input_Pass_2(); /* Complete the reading of the input file and 
				complete the initialization */
	AssignGlobalCoorDinates(); /* Assign global coordinates 
				      based on the input file */
	AssignHexHedronNodeNumbering(); /* Global node and edge numbering */
	AssignHexHedronEdgeNumbering();
	
/* More initialization */

	GBL_Matrix_Data=CMPLX_Matrix(TotEdgeNum,GBL_Matrix_MaxColNos);
	//GBL_Matrix_ColNos=INT_Matrix(TotEdgeNum,GBL_Matrix_MaxColNos);
	GBL_Matrix_Index=INT_Vector(TotEdgeNum);
	LHS_Matrix_Data=CMPLX_Matrix(TotInnerEdgeNum,GBL_Matrix_MaxColNos);
	/*LHS_Matrix_ColNos=INT_Matrix(TotInnerEdgeNum,GBL_Matrix_MaxColNos);*/
	LHS_Matrix_Index=INT_Vector(TotInnerEdgeNum);

	GBL_Matrix_ColNos = new int*[TotEdgeNum];
	for(int kkk = 0; kkk < TotEdgeNum; ++kkk)
		GBL_Matrix_ColNos[kkk] = new int[GBL_Matrix_MaxColNos];

	LHS_Matrix_ColNos = new int*[TotInnerEdgeNum];
	for(int kkk = 0; kkk < TotInnerEdgeNum; ++kkk)
		LHS_Matrix_ColNos[kkk] = new int[GBL_Matrix_MaxColNos];


for(freq_step=0;freq_step<NUM_OF_STEPS;freq_step++){ 
/* For Simulations over a range of frequencies 
   (All the routines above are done only the first step) */

for(i=0;i<TotEdgeNum;i++)  
   { 
   GBL_Matrix_Index[i] = 0;
   for(j=0;j<GBL_Matrix_MaxColNos;j++)
      {
      GBL_Matrix_Data[i][j] = COMplex_Null(); 
      /* Recreate matrix for every iteration */ 
      GBL_Matrix_ColNos[i][j] = 0;
      }
   }

for(i=0;i<TotInnerEdgeNum;i++)  
   {
   LHS_Matrix_Index[i] = 0;
   for(j=0;j<GBL_Matrix_MaxColNos;j++)
      {
      LHS_Matrix_Data[i][j] = COMplex_Null();
      LHS_Matrix_ColNos[i][j] = 0;
      }
   }


/* Step for every hexahedron 
   (XdiM, YdiM and ZdiM are the total number of bricks in each direction)*/
for(k=0;k<=ZdiM-1;k++)
for(j=0;j<=YdiM-1;j++)
for(i=0;i<=XdiM-1;i++)
   {
   HexNum = k*YdiM*XdiM + j*XdiM + i; /*Hexahedron number */
   con = pow((float)-1.0,(float)i)*pow((float)-1.0,(float)j)*pow((float)-1.0,(float)k);

   /* Divide hexahedron into 5 tetrahedra in two different ways */
   /* Adjacent hexahedra are divided in different ways (using "con")*/
   if(con==1) {HexHedraSubDiv1(); EdgeSubDiv1(); GlobalEdgeEndsDiv1();}
   else       {HexHedraSubDiv2(); EdgeSubDiv2(); GlobalEdgeEndsDiv2();}

   /* Compute the tetrahedron matrices and 
      assemble them into the final matrix */

   for(t=0;t<=4;t++)
      {
      FindTetraHedronVolume(); /* Finds volume (see user's guide) */
      ComputeTetraCoFactor(); /* Computes a's, b's c's d's */
      S_Matrix(); /* Computes the Eij part */
      T_Matrix(); /* COmputes the Fij part */
      TetraSubMatrix(); /* Complete Tetrahedron matrix */
      FindGlobalMatrix(); 
      /* Assemble tetrahedron matrices into global matrix */
      }
   }

if(freq_step==0) /* This comes frequently : 
		    Call Malloc only when it is the first frequency step */
   {
   RHSVector = CMPLX_Vector(TotInnerEdgeNum);
   Old_Solution = CMPLX_Vector(TotInnerEdgeNum);
   for(i=0;i<TotInnerEdgeNum;i++)
      {
      Old_Solution[i]=COMplex_Null();
      }
   }

for(i=0;i<TotInnerEdgeNum;i++)
   {
   RHSVector[i]=COMplex_Null();
   }
/*  printf("Reached here \n");  */
/*  fflush(stdout);  */

/* Partition the global matrix into the "inner" and "forced" edges */
PartitionGlobalMatrix();
/*  printf("Reached here \n");  */
/* fflush(stdout);  */


/* Solve the final matrix equation */
ConjugateSolver(TotInnerEdgeNum,LHS_Matrix_Index,LHS_Matrix_ColNos,LHS_Matrix_Data,RHSVector,MAX_ITERATIONS,TOLERANCE);
/*  printf("Reached here \n");  */
/*  fflush(stdout);  */

/* The following routine takes the "inner edge" vector along with the 
   "forced edge" vector and calculates and 
   writes out to all the required files */
Produce_Output();
OperateFreq +=FREQ_INC; /* Goto next frequency in the frequency range */
/* Update all the variables used */

AbsPermeable  = 1.25663706144E-06;
AbsPermitt    = 8.8542E-12;
FreeSpaceVel  = 1.0/sqrt(AbsPermeable*AbsPermitt);
WaveLength = FreeSpaceVel/(OperateFreq*1.0E+06);
WaveNumber = 2.0*3.14/WaveLength;
printf("OperateFrequency=%g MHz\n",OperateFreq);
printf("WaveLength=%f m\n",WaveLength );
printf("FreeSpaceVelocity=%E m/s\n",FreeSpaceVel );
printf("WaveNumber=%g\n",WaveNumber);}

return(0);
}

/***************************************************************************/

void ParameterInfo() /* This is done only for the first frequency step */
{

 AbsPermeable  = 1.25663706144E-06;
 AbsPermitt    = 8.8542E-12;
 FreeSpaceVel  = 1.0/sqrt(AbsPermeable*AbsPermitt);
 WaveLength = FreeSpaceVel/(OperateFreq*1.0E+06);
 WaveNumber = 2.0*3.14/WaveLength;
 Omega = WaveNumber*FreeSpaceVel;
 printf("OperateFrequency=%g MHz\n",OperateFreq);
 printf("WaveLength=%f m\n",WaveLength );
 printf("FreeSpaceVelocity=%E m/s\n",FreeSpaceVel );
 printf("WaveNumber=%g\n",WaveNumber);

/* Calculate the necessary parameters */

 TotEdgeNum      = 6*XdiM*YdiM*ZdiM + (XdiM+YdiM+ZdiM)
                 + 3*(XdiM*YdiM+YdiM*ZdiM+ZdiM*XdiM);
 TotNodeNum      = (XdiM+1)*(YdiM+1)*(ZdiM+1);
 TotNumHexHedron = XdiM*YdiM*ZdiM;
 TotTetElement = 5 * TotNumHexHedron;
 TotTrngleNum    = 4*(XdiM*YdiM + YdiM*ZdiM + ZdiM*XdiM);
 TotBoundEdgeNum = 6*(XdiM*YdiM + YdiM*ZdiM + ZdiM*XdiM);
 TotBoundNodeNum = 2*(XdiM*YdiM + YdiM*ZdiM + ZdiM*XdiM+1);

 printf("Total Edge Number = %d \n",TotEdgeNum);
 printf("Total Node Number = %d \n",TotNodeNum);
 printf("Total Cube Number = %d \n",TotNumHexHedron);
 printf("Total Tetrahedron Number = %d \n",TotTetElement);
 printf("Total Boundary Edge Number = %d \n",TotBoundEdgeNum);
 printf("Total Boundary Node Number = %d \n",TotBoundNodeNum);

  /* Local edge orientation of tetrahedron edges */
  edge_end[0][0] = 0; edge_end[0][1] = 1;
  edge_end[1][0] = 0; edge_end[1][1] = 2;
  edge_end[2][0] = 0; edge_end[2][1] = 3;
  edge_end[3][0] = 1; edge_end[3][1] = 2;
  edge_end[4][0] = 3; edge_end[4][1] = 1;
  edge_end[5][0] = 2; edge_end[5][1] = 3;

 }

/***************************************************************************
*   FUNCTION        : AssignGlobalCoorDinates()
****************************************************************************/

void AssignGlobalCoorDinates() 
/* Calculate the node coordinates for each global node */
 {
  int i,j,k,NodeNum; float Xtemp,Ytemp,Ztemp;
   {
     NodeNum=0; Xtemp=0.0; Ytemp=0.0; Ztemp=0.0;
     for(k=0;k<=ZdiM;k++) {
        for(j=0;j<=YdiM;j++) {
           for(i=0;i<=XdiM;i++) {
             NodeCord[NodeNum][0] = Xtemp;
             NodeCord[NodeNum][1] = Ytemp;
             NodeCord[NodeNum][2] = Ztemp;
/*           printf(" \n %f %f %f",NodeCord[NodeNum][0],NodeCord[NodeNum][1],NodeCord[NodeNum][2]);*/
             NodeNum++;   Xtemp += DivisorX[i+1]; 

				}
             Xtemp = 0.0; Ytemp = Ytemp + DivisorY[j+1]; 
			     }
             Ytemp = 0.0; Ztemp = Ztemp + DivisorZ[k+1]; 
			  }
    }
 }

/***************************************************************************
*   FUNCTION        : AssignHexHedronNodeNumbering()
****************************************************************************/

void AssignHexHedronNodeNumbering()
/* Assign global node numbering 
   (see "Mesh generation" chapter in the user's guide 
   to verify these formulas )*/
{
 int i;
  {
   start=1;
   for(i=0;i<=TotNumHexHedron-1;i++)
     {
      HexNode[i][0]=start;
      HexNode[i][1]=HexNode[i][0]+1;
      HexNode[i][2]=HexNode[i][0]+XdiM+1;
      HexNode[i][3]=HexNode[i][2]+1;
      HexNode[i][4]=HexNode[i][0]+(XdiM+1)*(YdiM+1);
      HexNode[i][5]=HexNode[i][4]+1;
      HexNode[i][6]=HexNode[i][4]+(XdiM+1);
      HexNode[i][7]=HexNode[i][6]+1;
      if(((i+1)%((YdiM*XdiM))==0)) {start=start+(XdiM+3);}
      else                     {start=start+1;}
      if((start%(XdiM+1))==0)  { start++; }
     }
  }
}

/***************************************************************************
*   FUNCTION        : AssignHexHedronEdgeNumbering()
****************************************************************************/

void AssignHexHedronEdgeNumbering()
/* Assign global edge numbering 
   (see "Mesh generation" chapter in the user's guide 
   to verify these formulas and try to work it out by hand) */
  {
  int EdgeNum,xx1,xx2,xx3,kk,dy1,dy2,divz,divz1,i,j,k;
  {
        for(k=0;k<=ZdiM-1;k++) 
        for(j=0;j<=YdiM-1;j++)
        for(i=0;i<=XdiM-1;i++) {  
        EdgeNum = k*YdiM*XdiM + j*XdiM + i;
        for(kk=0;kk<=17;kk++)   {HexEdgeNum[EdgeNum][kk] = 0;}
        divz=(6*XdiM*YdiM+3*XdiM+3*YdiM+1)*k;
        divz1=  6*XdiM*YdiM + 3*XdiM + 3*YdiM + 1 ;
        dy1 = (3*XdiM+1)*j; dy2 = (3*XdiM+2)*j;
        xx1 = 1 + (YdiM+1)*XdiM + (2*XdiM+1)*YdiM;
        xx2 = xx1 + (2*XdiM+1);
        xx3 = xx2 + (XdiM + 1);
        HexEdgeNum[EdgeNum][0] = 1 + i + dy1 + divz;
        HexEdgeNum[EdgeNum][1] = 1 + XdiM + 2*i + dy1 + divz;
        HexEdgeNum[EdgeNum][2] = HexEdgeNum[EdgeNum][1] + 1;  
        HexEdgeNum[EdgeNum][3] = HexEdgeNum[EdgeNum][2] + 1;
        HexEdgeNum[EdgeNum][4] = 1 + 3*XdiM + 1 + i + dy1 + divz;
        HexEdgeNum[EdgeNum][5] = xx1 + 2*i + dy2 + divz;
        HexEdgeNum[EdgeNum][6] = HexEdgeNum[EdgeNum][5] + 1; 
        HexEdgeNum[EdgeNum][7] = HexEdgeNum[EdgeNum][6] + 1;
        HexEdgeNum[EdgeNum][8] = xx2 + i + dy2 + divz;
        HexEdgeNum[EdgeNum][9] = HexEdgeNum[EdgeNum][8] + 1;
        HexEdgeNum[EdgeNum][10]= xx3 + 2*i + dy2 + divz;
        HexEdgeNum[EdgeNum][11]= HexEdgeNum[EdgeNum][10] + 1;
        HexEdgeNum[EdgeNum][12]= HexEdgeNum[EdgeNum][11] + 1; 
        HexEdgeNum[EdgeNum][13]= HexEdgeNum[EdgeNum][0] + divz1;
        HexEdgeNum[EdgeNum][14]= HexEdgeNum[EdgeNum][1] + divz1; 
        HexEdgeNum[EdgeNum][15]= HexEdgeNum[EdgeNum][2] + divz1;
        HexEdgeNum[EdgeNum][16]= HexEdgeNum[EdgeNum][3] + divz1; 
        HexEdgeNum[EdgeNum][17]= HexEdgeNum[EdgeNum][4] + divz1;
                               }
        
 }
}

/***************************************************************************
*   FUNCTION        : HexHedraSubDiv1()
****************************************************************************/

void HexHedraSubDiv1()
/* Divides the 8 nodes in the hexahedron into five groups of 4 each */
{
TetGlobalNodeNum[0][0]=HexNode[HexNum][0];
TetGlobalNodeNum[1][0]=HexNode[HexNum][3];
TetGlobalNodeNum[0][1]=HexNode[HexNum][3];
TetGlobalNodeNum[1][1]=HexNode[HexNum][5];
TetGlobalNodeNum[0][2]=HexNode[HexNum][2];
TetGlobalNodeNum[1][2]=HexNode[HexNum][7];
TetGlobalNodeNum[0][3]=HexNode[HexNum][6];
TetGlobalNodeNum[1][3]=HexNode[HexNum][6];
 
TetGlobalNodeNum[2][0]=HexNode[HexNum][0];
TetGlobalNodeNum[3][0]=HexNode[HexNum][0];
TetGlobalNodeNum[2][1]=HexNode[HexNum][1];
TetGlobalNodeNum[3][1]=HexNode[HexNum][6];
TetGlobalNodeNum[2][2]=HexNode[HexNum][3];
TetGlobalNodeNum[3][2]=HexNode[HexNum][4];
TetGlobalNodeNum[2][3]=HexNode[HexNum][5];
TetGlobalNodeNum[3][3]=HexNode[HexNum][5];
 
TetGlobalNodeNum[4][0]=HexNode[HexNum][0];
TetGlobalNodeNum[4][2]=HexNode[HexNum][3];
TetGlobalNodeNum[4][1]=HexNode[HexNum][5];
TetGlobalNodeNum[4][3]=HexNode[HexNum][6];
}

/***************************************************************************
*   FUNCTION        : HexHedraSubDiv2()
****************************************************************************/

void HexHedraSubDiv2()
/* Divides the 8 nodes in the hexahedron into five groups of 4 each 
   in another way*/
{
TetGlobalNodeNum[0][0]=HexNode[HexNum][0];
TetGlobalNodeNum[1][0]=HexNode[HexNum][1];
TetGlobalNodeNum[0][1]=HexNode[HexNum][1];
TetGlobalNodeNum[1][1]=HexNode[HexNum][7];
TetGlobalNodeNum[0][2]=HexNode[HexNum][2];
TetGlobalNodeNum[1][2]=HexNode[HexNum][4];
TetGlobalNodeNum[0][3]=HexNode[HexNum][4];
TetGlobalNodeNum[1][3]=HexNode[HexNum][5];
 
TetGlobalNodeNum[2][0]=HexNode[HexNum][1];
TetGlobalNodeNum[3][0]=HexNode[HexNum][2];
TetGlobalNodeNum[2][1]=HexNode[HexNum][7];
TetGlobalNodeNum[3][1]=HexNode[HexNum][6];
TetGlobalNodeNum[2][2]=HexNode[HexNum][3];
TetGlobalNodeNum[3][2]=HexNode[HexNum][4];
TetGlobalNodeNum[2][3]=HexNode[HexNum][2];
TetGlobalNodeNum[3][3]=HexNode[HexNum][7];

TetGlobalNodeNum[4][0]=HexNode[HexNum][1];
TetGlobalNodeNum[4][2]=HexNode[HexNum][4];
TetGlobalNodeNum[4][1]=HexNode[HexNum][2];
TetGlobalNodeNum[4][3]=HexNode[HexNum][7];
}

/***************************************************************************
*   FUNCTION        : EdgeSubDiv1()
****************************************************************************/

void EdgeSubDiv1()
/* Divides the 18 edges in the hexahedron into five groups of 6 each */
{
TetEdgeNum[0][0] = HexEdgeNum[HexNum][2];  
TetEdgeNum[0][1] = HexEdgeNum[HexNum][1]; 
TetEdgeNum[0][2] = HexEdgeNum[HexNum][8];    
TetEdgeNum[0][3] =-HexEdgeNum[HexNum][4];  
TetEdgeNum[0][4] =-HexEdgeNum[HexNum][11];  
TetEdgeNum[0][5] = HexEdgeNum[HexNum][10];   

TetEdgeNum[1][0] = HexEdgeNum[HexNum][9];   
TetEdgeNum[1][1] = HexEdgeNum[HexNum][12];  
TetEdgeNum[1][2] = HexEdgeNum[HexNum][11]; 
TetEdgeNum[1][3] = HexEdgeNum[HexNum][16];  
TetEdgeNum[1][4] =-HexEdgeNum[HexNum][15];  
TetEdgeNum[1][5] =-HexEdgeNum[HexNum][17];   
 
TetEdgeNum[2][0] = HexEdgeNum[HexNum][0];   
TetEdgeNum[2][1] = HexEdgeNum[HexNum][2];   
TetEdgeNum[2][2] = HexEdgeNum[HexNum][6]; 
TetEdgeNum[2][3] = HexEdgeNum[HexNum][3];  
TetEdgeNum[2][4] =-HexEdgeNum[HexNum][7];  
TetEdgeNum[2][5] = HexEdgeNum[HexNum][9]; 
 
TetEdgeNum[3][0] = HexEdgeNum[HexNum][8];   
TetEdgeNum[3][1] = HexEdgeNum[HexNum][5];  
TetEdgeNum[3][2] = HexEdgeNum[HexNum][6]; 
TetEdgeNum[3][3] =-HexEdgeNum[HexNum][14];  
TetEdgeNum[3][4] = HexEdgeNum[HexNum][15];  
TetEdgeNum[3][5] = HexEdgeNum[HexNum][13]; 
 
TetEdgeNum[4][0] = HexEdgeNum[HexNum][6];   
TetEdgeNum[4][1] = HexEdgeNum[HexNum][2];   
TetEdgeNum[4][2] = HexEdgeNum[HexNum][8];  
TetEdgeNum[4][3] =-HexEdgeNum[HexNum][9]; 
TetEdgeNum[4][4] =-HexEdgeNum[HexNum][15];  
TetEdgeNum[4][5] = HexEdgeNum[HexNum][11];  
}

/***************************************************************************
*   FUNCTION        : EdgeSubDiv2()
****************************************************************************/

void EdgeSubDiv2()
/* Divides the 18 edges in the hexahedron into five groups of 6 each 
 in another way */
{
TetEdgeNum[0][0] = HexEdgeNum[HexNum][0];
TetEdgeNum[0][1] = HexEdgeNum[HexNum][1];
TetEdgeNum[0][2] = HexEdgeNum[HexNum][5];
TetEdgeNum[0][3] = HexEdgeNum[HexNum][2];
TetEdgeNum[0][4] =-HexEdgeNum[HexNum][6];
TetEdgeNum[0][5] = HexEdgeNum[HexNum][8];
 
TetEdgeNum[1][0] = HexEdgeNum[HexNum][9];
TetEdgeNum[1][1] = HexEdgeNum[HexNum][6];
TetEdgeNum[1][2] = HexEdgeNum[HexNum][7];
TetEdgeNum[1][3] =-HexEdgeNum[HexNum][15];
TetEdgeNum[1][4] = HexEdgeNum[HexNum][16];    
TetEdgeNum[1][5] = HexEdgeNum[HexNum][13];   
 
TetEdgeNum[2][0] = HexEdgeNum[HexNum][9];   
TetEdgeNum[2][1] = HexEdgeNum[HexNum][3];     
TetEdgeNum[2][2] = HexEdgeNum[HexNum][2];
TetEdgeNum[2][3] =-HexEdgeNum[HexNum][12];     
TetEdgeNum[2][4] = HexEdgeNum[HexNum][11];     
TetEdgeNum[2][5] =-HexEdgeNum[HexNum][4];  
 
TetEdgeNum[3][0] = HexEdgeNum[HexNum][10];     
TetEdgeNum[3][1] = HexEdgeNum[HexNum][8];
TetEdgeNum[3][2] = HexEdgeNum[HexNum][11];  
TetEdgeNum[3][3] =-HexEdgeNum[HexNum][14];   
TetEdgeNum[3][4] =-HexEdgeNum[HexNum][17];  
TetEdgeNum[3][5] = HexEdgeNum[HexNum][15];  
 
TetEdgeNum[4][0] = HexEdgeNum[HexNum][2];   
TetEdgeNum[4][1] = HexEdgeNum[HexNum][6];   
TetEdgeNum[4][2] = HexEdgeNum[HexNum][9];   
TetEdgeNum[4][3] = HexEdgeNum[HexNum][8];  
TetEdgeNum[4][4] =-HexEdgeNum[HexNum][11];    
TetEdgeNum[4][5] = HexEdgeNum[HexNum][15];  
}

/***************************************************************************
*   FUNCTION        : GlobalEdgeEndsDiv1()
****************************************************************************/

void GlobalEdgeEndsDiv1()
/* Specifies the end-nodes for each edge for the first division 
See the "Division of hexahedra into five tetrahedra" figure*/
{
TetGlobalEdgeEnds[HexEdgeNum[HexNum][0]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][0]-1][1]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][1]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][1]-1][1]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][2]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][2]-1][1]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][3]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][3]-1][1]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][4]-1][0]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][4]-1][1]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][5]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][5]-1][1]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][6]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][6]-1][1]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][7]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][7]-1][1]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][8]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][8]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][9]-1][0]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][9]-1][1]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][10]-1][0]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][10]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][11]-1][0]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][11]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][12]-1][0]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][12]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][13]-1][0]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][13]-1][1]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][14]-1][0]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][14]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][15]-1][0]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][15]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][16]-1][0]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][16]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][17]-1][0]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][17]-1][1]=HexNode[HexNum][7];
   }

/***************************************************************************
*   FUNCTION        : GlobalEdgeEndsDiv2()
****************************************************************************/

void GlobalEdgeEndsDiv2()
/* Specifies the end-nodes for each edge for the second division */
{
TetGlobalEdgeEnds[HexEdgeNum[HexNum][0]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][0]-1][1]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][1]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][1]-1][1]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][2]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][2]-1][1]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][3]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][3]-1][1]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][4]-1][0]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][4]-1][1]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][5]-1][0]=HexNode[HexNum][0];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][5]-1][1]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][6]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][6]-1][1]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][7]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][7]-1][1]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][8]-1][0]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][8]-1][1]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][9]-1][0]=HexNode[HexNum][1];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][9]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][10]-1][0]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][10]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][11]-1][0]=HexNode[HexNum][2];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][11]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][12]-1][0]=HexNode[HexNum][3];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][12]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][13]-1][0]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][13]-1][1]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][14]-1][0]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][14]-1][1]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][15]-1][0]=HexNode[HexNum][4];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][15]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][16]-1][0]=HexNode[HexNum][5];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][16]-1][1]=HexNode[HexNum][7];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][17]-1][0]=HexNode[HexNum][6];
TetGlobalEdgeEnds[HexEdgeNum[HexNum][17]-1][1]=HexNode[HexNum][7];
}


/***************************************************************************
*   FUNCTION        : Find_Volume_Determinent()
****************************************************************************/

void FindTetraHedronVolume()
/* Finds the volume of each tetrahedron to be used further in calculations */

{
   int Count_i, Count_j;
   float Buff[4][3],Buff1[3],Buff2[3],Buff3[3],Buff4[3];
   for(Count_i=0;Count_i<=3;Count_i++)
   for(Count_j=0;Count_j<=2;Count_j++)
   Buff[Count_i][Count_j]=NodeCord[TetGlobalNodeNum[t][Count_i]-1][Count_j];
   VTXsub1(1,0,Buff,Buff1);
   VTXsub1(2,0,Buff,Buff2);
   VTXsub1(3,0,Buff,Buff3);
   VTXcross(Buff1,Buff2,Buff4);
   TetVolume=(1.0/6.0)*(Buff3[0]*Buff4[0]+Buff3[1]*Buff4[1]+Buff3[2]*Buff4[2]);
}

/***************************************************************************
*   FUNCTION        : ComputeTetraCoFactor()
****************************************************************************/

void ComputeTetraCoFactor()
/* This calculates the a's, b's, c's and d's in reference [9] */
/* This just involves some basic matrix simulations */
   {
   int i,j,RowNum,ColNum,Row,Col,Count=0;
   float CoF[4][4],Buff[4][4],CoFs[4][4]; 
     {
    for(i=0;i<=3;i++)  for(j=0;j<=3;j++) {
    if(j==3) Buff[i][j]= 1.0;
    else     Buff[i][j]= NodeCord[TetGlobalNodeNum[t][i]-1][j];
                                         }

/* The main matrix is made up of a column of 1's and 
   the other columns are made of the node coordinates */

/* The co-factor matrices are formed from the main matrix 
   by taking part of the matrix */

     for(RowNum=0;RowNum<=3;RowNum++)
     for(ColNum=0;ColNum<=3;ColNum++) {
     Row=0;Col=0;Count=0;
        for(i=0;i<=3;i++)
        for(j=0;j<=3;j++) {
        if((i!=RowNum) && (j!=ColNum)) {
        Count++;
        if(Count<4)              Row=0;
        if((Count>3)&&(Count<7)) Row=1;
        else if(Count>6)         Row=2;
        CoF[Row][Col]=Buff[i][j];
        Col++; Col=Col%3;              }
                           }

/* The next few steps calculate the co-factors 
   by taking the determinant of the sub-matrix */
     Row=0;Col=0;Count=0;
     CoFs[RowNum][ColNum]
     =CoF[0][0]*((CoF[1][1]*CoF[2][2])-(CoF[1][2]*CoF[2][1]))
     -CoF[0][1]*((CoF[1][0]*CoF[2][2])-(CoF[2][0]*CoF[1][2]))
     +CoF[0][2]*((CoF[1][0]*CoF[2][1])-(CoF[2][0]*CoF[1][1]));
     if((RowNum+ColNum)%2!=0)  CoFs[RowNum][ColNum]*=-1.0;
                                      }
     for(RowNum=0;RowNum<=3;RowNum++) {
     if     (RowNum==0) Col = 3;
     else if(RowNum==1) Col = 0;
     else if(RowNum==2) Col = 1;
     else if(RowNum==3) Col = 2;
     for(ColNum=0;ColNum<=3;ColNum++) {
     CoFactor[RowNum][ColNum] = -CoFs[ColNum][Col];
                                      }}
      }
    }

/***************************************************************************
*   FUNCTION        : S_Matrix()
****************************************************************************/

void S_Matrix()
/* This calculates the Eij part of the tetrahedron matrix 
All the formulas used are in reference [9] or see reference [1]*/

  {
  int i,j;  float Length_i,Length_j,Mult1,Volume;
  float S_mat1[6][6];
   {

if (Is_PML[HexNum] == 0) /* Not a PML */


   {
   Volume = TetVolume; Mult1 = (1296.0*Volume*Volume*Volume);
   for(i=0;i<=5;i++) for(j=0;j<=5;j++)  {
   Length_i = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[i][0]]-1],
                     NodeCord[TetGlobalNodeNum[t][edge_end[i][1]]-1]);
   Length_j = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[j][0]]-1],
                     NodeCord[TetGlobalNodeNum[t][edge_end[j][1]]-1]);
   S_mat1[i][j] = (CoFactor[2][edge_end[i][0]] * CoFactor[3][edge_end[i][1]]
                 - CoFactor[3][edge_end[i][0]] * CoFactor[2][edge_end[i][1]] )
               * ( CoFactor[2][edge_end[j][0]] * CoFactor[3][edge_end[j][1]]
                 - CoFactor[3][edge_end[j][0]] * CoFactor[2][edge_end[j][1]] )
               + ( CoFactor[3][edge_end[i][0]] * CoFactor[1][edge_end[i][1]]
                 - CoFactor[1][edge_end[i][0]] * CoFactor[3][edge_end[i][1]] )
               * ( CoFactor[3][edge_end[j][0]] * CoFactor[1][edge_end[j][1]]
                 - CoFactor[1][edge_end[j][0]] * CoFactor[3][edge_end[j][1]] )
               + ( CoFactor[1][edge_end[i][0]] * CoFactor[2][edge_end[i][1]]
                 - CoFactor[2][edge_end[i][0]] * CoFactor[1][edge_end[i][1]] )
               * ( CoFactor[1][edge_end[j][0]] * CoFactor[2][edge_end[j][1]]
                 - CoFactor[2][edge_end[j][0]] * CoFactor[1][edge_end[j][1]] );
   S_mat1[i][j] = (Length_i * Length_j * S_mat1[i][j])/ Mult1 ;
   S_mat[i][j] = COMplex_Cmplx(S_mat1[i][j],0.0); /* Make it complex */
   }/*
   printf(" %f %f\n",S_mat[0][0].x,S_mat[0][0].y);
  */ }
else /* The hexahedron is a PML */
/* For PMLs, the matrix calculation is different because the 
   dielectric matrix is anisotropic */

   { 

     a_PML = a_PML_vec[HexNum];
     b_PML = b_PML_vec[HexNum];
     c_PML = c_PML_vec[HexNum];
   Volume = TetVolume; Mult1 = (1296.0*Volume*Volume*Volume);
   for(i=0;i<=5;i++) for(j=0;j<=5;j++)  {
   Length_i = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[i][0]]-1], 
                     NodeCord[TetGlobalNodeNum[t][edge_end[i][1]]-1]);
   Length_j = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[j][0]]-1], 
                     NodeCord[TetGlobalNodeNum[t][edge_end[j][1]]-1]);

/* PML calculations (a_PML, b_PML and c_PML are the 
   diagonal entries in the dielectric matrix). 
   This formula comes from a matrix-vector dot-product */
   S_mat[i][j] = COMplex_Add2(
		 Real_Mul(
		 ( CoFactor[2][edge_end[i][0]] * CoFactor[3][edge_end[i][1]]  
                 - CoFactor[3][edge_end[i][0]] * CoFactor[2][edge_end[i][1]] )
               * ( CoFactor[2][edge_end[j][0]] * CoFactor[3][edge_end[j][1]]  
                 - CoFactor[3][edge_end[j][0]] * CoFactor[2][edge_end[j][1]] ),
		   COMplex_Div(COMplex_Cmplx(1.0,0.0),a_PML))
               , Real_Mul(
		 ( CoFactor[3][edge_end[i][0]] * CoFactor[1][edge_end[i][1]]  
                 - CoFactor[1][edge_end[i][0]] * CoFactor[3][edge_end[i][1]] )
               * ( CoFactor[3][edge_end[j][0]] * CoFactor[1][edge_end[j][1]]  
                 - CoFactor[1][edge_end[j][0]] * CoFactor[3][edge_end[j][1]] ),
		   COMplex_Div(COMplex_Cmplx(1.0,0.0),b_PML))
               , Real_Mul(
		 ( CoFactor[1][edge_end[i][0]] * CoFactor[2][edge_end[i][1]]   
                 - CoFactor[2][edge_end[i][0]] * CoFactor[1][edge_end[i][1]] ) 
               * ( CoFactor[1][edge_end[j][0]] * CoFactor[2][edge_end[j][1]]   
                 - CoFactor[2][edge_end[j][0]] * CoFactor[1][edge_end[j][1]] ),
		   COMplex_Div(COMplex_Cmplx(1.0,0.0),c_PML)));
   S_mat[i][j] = Real_Mul((Length_i*Length_j/Mult1),S_mat[i][j]) ;
   }/*
    printf(" %f %f\n",S_mat[0][0].x,S_mat[0][0].y);
fflush(stdout);*/  
 }
   }
   }

/***************************************************************************
*   FUNCTION        : T_Matrix()
****************************************************************************/

void T_Matrix()
/* This calculate the Fij part of the tetrahedron matrix 
See reference [9]*/

  {
  int i,j;
  float Mult1,Length_i,Length_j,T_mat[6][6];
  Mult1 = 1.0/(720.0 * TetVolume);
if (Is_PML[HexNum] == 0) /* Not a PML layer */
   {

/* (see reference [9] for explanation of the following terms */
   T_mat[0][0] =               ff(1,1) + ff(2,2) - ff(1,2);
   T_mat[0][1] = 2 * ff(2,3) - ff(2,1) - ff(1,3) + ff(1,1);
   T_mat[0][2] = 2 * ff(2,4) - ff(2,1) - ff(1,4) + ff(1,1);
   T_mat[0][3] =-2 * ff(1,3) - ff(2,2) + ff(2,3) + ff(1,2);
   T_mat[0][4] = 2 * ff(1,4) - ff(2,4) - ff(1,2) + ff(2,2);
   T_mat[0][5] =     ff(2,4) - ff(2,3) - ff(1,4) + ff(1,3);
   T_mat[1][1] =               ff(1,1) + ff(3,3) - ff(1,3);
   T_mat[1][2] = 2 * ff(3,4) - ff(1,3) - ff(1,4) + ff(1,1);
   T_mat[1][3] = 2 * ff(1,2) - ff(2,3) - ff(1,3) + ff(3,3);
   T_mat[1][4] =     ff(2,3) - ff(3,4) - ff(1,2) + ff(1,4);
   T_mat[1][5] =-2 * ff(1,4) - ff(3,3) + ff(1,3) + ff(3,4);
   T_mat[2][2] =               ff(1,1) + ff(4,4) - ff(1,4);
   T_mat[2][3] =     ff(3,4) - ff(2,4) - ff(1,3) + ff(1,2);
   T_mat[2][4] =-2 * ff(1,2) - ff(4,4) + ff(1,4) + ff(2,4);
   T_mat[2][5] = 2 * ff(1,3) - ff(3,4) - ff(1,4) + ff(4,4);
   T_mat[3][3] =               ff(3,3) + ff(2,2) - ff(2,3);
   T_mat[3][4] =-2 * ff(3,4) - ff(2,2) + ff(2,3) + ff(2,4);
   T_mat[3][5] =-2 * ff(2,4) - ff(3,3) + ff(2,3) + ff(3,4);
   T_mat[4][4] =               ff(2,2) + ff(4,4) - ff(2,4);
   T_mat[4][5] =-2 * ff(2,3) - ff(4,4) + ff(2,4) + ff(3,4);
   T_mat[5][5] =               ff(3,3) + ff(4,4) - ff(3,4);
   for(i=0;i<=5;i++) for(j=0;j<=i-1;j++)
   {T_mat[i][j] = T_mat[j][i];}
 
   for(i=0;i<=5;i++) for(j=0;j<=5;j++) {
   Length_i = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[i][0]]-1],
                     NodeCord[TetGlobalNodeNum[t][edge_end[i][1]]-1]);
   Length_j = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[j][0]]-1],
                     NodeCord[TetGlobalNodeNum[t][edge_end[j][1]]-1]);
   if(i==j) T_mat[i][j] = 2 * Length_i * Length_j * Mult1 * T_mat[i][j];
   else     T_mat[i][j] =     Length_i * Length_j * Mult1 * T_mat[i][j];
   T_mat[i][j] = WaveNumber*WaveNumber * T_mat[i][j] ;
   T_mat1[i][j] = Real_Mul(T_mat[i][j],RELPerm[HexNum]);
                                       }
   }    
else /* Part of PML layer */
/* For PMLs, the matrix calculation is different because the 
   dielectric matrix is anisotropic. The anisotropicity is used in the 
   ff1 calculations. See function above the main routine.
   The calculation is same as above but uses complex numbers 
   and includes the PML dielectric matrix */

   {

     a_PML = a_PML_vec[HexNum];
     b_PML = b_PML_vec[HexNum];
     c_PML = c_PML_vec[HexNum];
   T_mat1[0][0] = COMplex_Add(ff1(1,1),COMplex_Sub(ff1(2,2),ff1(1,2)));
   T_mat1[1][1] = COMplex_Add(ff1(1,1),COMplex_Sub(ff1(3,3),ff1(1,3)));
   T_mat1[2][2] = COMplex_Add(ff1(1,1),COMplex_Sub(ff1(4,4),ff1(1,4)));
   T_mat1[3][3] = COMplex_Add(ff1(3,3),COMplex_Sub(ff1(2,2),ff1(2,3)));
   T_mat1[4][4] = COMplex_Add(ff1(2,2),COMplex_Sub(ff1(4,4),ff1(2,4)));
   T_mat1[5][5] = COMplex_Add(ff1(3,3),COMplex_Sub(ff1(4,4),ff1(3,4)));

   T_mat1[0][5] = COMplex_Sub(COMplex_Sub(ff1(2,4),ff1(2,3)),
			      COMplex_Sub(ff1(1,4),ff1(1,3)));
   T_mat1[1][4] = COMplex_Sub(COMplex_Sub(ff1(2,3),ff1(3,4)),
			      COMplex_Sub(ff1(1,2),ff1(1,4)));
   T_mat1[2][3] = COMplex_Sub(COMplex_Sub(ff1(3,4),ff1(2,4)),
			      COMplex_Sub(ff1(1,3),ff1(1,2)));

   T_mat1[0][1] = COMplex_Sub(COMplex_Sub(
			      COMplex_Add(ff1(2,3),ff1(2,3)),ff1(2,1)),
			      COMplex_Sub(ff1(1,3),ff1(1,1)));
   T_mat1[0][2] = COMplex_Sub(COMplex_Sub(
			      COMplex_Add(ff1(2,4),ff1(2,4)),ff1(2,1)),
			      COMplex_Sub(ff1(1,4),ff1(1,1)));
   T_mat1[0][4] = COMplex_Sub(COMplex_Sub(
			      COMplex_Add(ff1(1,4),ff1(1,4)),ff1(2,4)),
			      COMplex_Sub(ff1(1,2),ff1(2,2)));
   T_mat1[1][2] = COMplex_Sub(COMplex_Sub(
			      COMplex_Add(ff1(3,4),ff1(3,4)),ff1(1,3)),
			      COMplex_Sub(ff1(1,4),ff1(1,1)));
   T_mat1[1][3] = COMplex_Sub(COMplex_Sub(
			      COMplex_Add(ff1(1,2),ff1(1,2)),ff1(2,3)),
			      COMplex_Sub(ff1(1,3),ff1(3,3)));
   T_mat1[2][5] = COMplex_Sub(COMplex_Sub(
			      COMplex_Add(ff1(1,3),ff1(1,3)),ff1(3,4)),
			      COMplex_Sub(ff1(1,4),ff1(4,4)));

   T_mat1[0][3] = COMplex_Sub(COMplex_Add(ff1(2,3),ff1(1,2)),
			      COMplex_Add2(ff1(1,3),ff1(1,3),ff1(2,2)));
   T_mat1[1][5] = COMplex_Sub(COMplex_Add(ff1(1,3),ff1(3,4)),
			      COMplex_Add2(ff1(1,4),ff1(1,4),ff1(3,3)));
   T_mat1[2][4] = COMplex_Sub(COMplex_Add(ff1(1,4),ff1(2,4)),
			      COMplex_Add2(ff1(1,2),ff1(1,2),ff1(4,4)));
   T_mat1[3][4] = COMplex_Sub(COMplex_Add(ff1(2,3),ff1(2,4)),
			      COMplex_Add2(ff1(3,4),ff1(3,4),ff1(2,2)));
   T_mat1[3][5] = COMplex_Sub(COMplex_Add(ff1(2,3),ff1(3,4)),
			      COMplex_Add2(ff1(2,4),ff1(2,4),ff1(3,3)));
   T_mat1[4][5] = COMplex_Sub(COMplex_Add(ff1(2,4),ff1(3,4)),
			      COMplex_Add2(ff1(2,3),ff1(2,3),ff1(4,4)));

   for(i=0;i<=5;i++) for(j=0;j<=i-1;j++)
   {T_mat1[i][j] = T_mat1[j][i];}
   
   for(i=0;i<=5;i++) for(j=0;j<=5;j++) {
   Length_i = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[i][0]]-1],
                     NodeCord[TetGlobalNodeNum[t][edge_end[i][1]]-1]);
   Length_j = VTXmag(NodeCord[TetGlobalNodeNum[t][edge_end[j][0]]-1], 
                     NodeCord[TetGlobalNodeNum[t][edge_end[j][1]]-1]);   
   if(i==j) T_mat1[i][j] = Real_Mul(2*Length_i*Length_j*Mult1,T_mat1[i][j]);
   else     T_mat1[i][j] = Real_Mul(Length_i*Length_j*Mult1,T_mat1[i][j]);
   T_mat1[i][j] = Real_Mul(WaveNumber*WaveNumber,T_mat1[i][j]) ;
                                       }
 } 

}

/***************************************************************************
*   FUNCTION        : TetraSubMatrix()
****************************************************************************/

void TetraSubMatrix()
/* Combine the Eij and Fij matrices into the A matrix */

  {
  int i,j;
  for(i=0;i<=5;i++) for(j=0;j<=5;j++)
  {A_mat[t][i][j]=COMplex_Null();}
  for(i=0;i<=5;i++) for(j=0;j<=5;j++)
  {A_mat[t][i][j] = COMplex_Sub(Real_Mul(4.0,S_mat[i][j]),T_mat1[i][j]);}
  for(i=0;i<=5;i++) for(j=0;j<=5;j++)
  {A_mat[t][i][j] = Real_Mul(Sign(TetEdgeNum[t][j]),A_mat[t][i][j]);}
  for(i=0;i<=5;i++) for(j=0;j<=5;j++)
  {A_mat[t][i][j] = Real_Mul(Sign(TetEdgeNum[t][i]),A_mat[t][i][j]);}
  }

/***************************************************************************
*   FUNCTION        : FindGlobalMatrix()
****************************************************************************/

void FindGlobalMatrix()
/* Assembles the tetrahedron matrices into the global matrices
   and stores the global matrix in the row-indexed format */
{
	int m,n,k,l,rval;
	{
		for(m=0;m<=5;m++)
			for(n=0;n<=5;n++){
				k=abs(TetEdgeNum[t][m])-1;
				l=abs(TetEdgeNum[t][n])-1;
				if(COMplex_Abs(A_mat[t][m][n])>1.0E-10)
				{
					rval=Srch_NZ_Element_Column(k,l); /* Search if there is such an element */
					/* Matrix term already exists,  so add the new term to the old one */
					if(rval>=0)
					{
						GBL_Matrix_Data[k][rval]=COMplex_Add(GBL_Matrix_Data[k][rval],A_mat[t][m][n]);
					}
					/* Create new global matrix element */
					else
					{
						GBL_Matrix_Data[k][GBL_Matrix_Index[k]]=A_mat[t][m][n];
						GBL_Matrix_ColNos[k][GBL_Matrix_Index[k]] = l;
						GBL_Matrix_Index[k]++;
					}
				}
			}
	}
}

/***************************************************************************
*   FUNCTION      : PartitionGlobalMatrix()
****************************************************************************/

void PartitionGlobalMatrix()
{
	int loop_i,loop_j,No_cols,k_row=0,II,JJ,flag=0;
	/* Just for debugging
	for (II=0;II<TotEdgeNum;II++)
	   {
	   for(JJ=0;JJ<GBL_Matrix_Index[II];JJ++)
		  {
		  if ( GBL_Matrix_Data[II][JJ].y > 1000)
			 {

			 printf("\n %d \n", II);
			 printf("%f  ",NodeCord[TetGlobalEdgeEnds[II][0]][0]);
			 printf("%f  ",NodeCord[TetGlobalEdgeEnds[II][0]][1]);
			 printf("%f  ",NodeCord[TetGlobalEdgeEnds[II][0]][2]);
			 printf("%f  ",NodeCord[TetGlobalEdgeEnds[II][1]][0]);
			 printf("%f  ",NodeCord[TetGlobalEdgeEnds[II][1]][1]);
			 printf("%f\n",NodeCord[TetGlobalEdgeEnds[II][1]][2]);
		 }
		  }
	   }
	*/
	for (loop_i=0;loop_i<TotEdgeNum;loop_i++) /* Loop over all the edges */
	{
		flag=0;
		if (ForceStat[loop_i] == 0) /* Inner edge */
		{
			No_cols = GBL_Matrix_Index[loop_i];
			/* printf("%d \n",No_cols); */
			for (loop_j=0;loop_j<No_cols;loop_j++) 
				/*Loop over all non-zero elements in row loop_i */
			{
				if (ForceStat[GBL_Matrix_ColNos[loop_i][loop_j]] == 0) 
					/* Inner edge */
				{ 
					/* Construct the LHS matrix from the elements of the global matrix */
					LHS_Matrix_Data[k_row][LHS_Matrix_Index[k_row]] = GBL_Matrix_Data[loop_i][loop_j];
					if ((loop_i == GBL_Matrix_ColNos[loop_i][loop_j]) && (COMplex_Abs(Resist[loop_i])>TOLERANCE))
					{
						/* Add the term for the resistor */
						/* printf("Inside this resistor loop"); */
						/* fflush(stdout); */
						LHS_Matrix_Data[k_row][LHS_Matrix_Index[k_row]] = COMplex_Add(LHS_Matrix_Data[k_row][LHS_Matrix_Index[k_row]],Resist[loop_i]);
					}
					LHS_Matrix_ColNos[k_row][LHS_Matrix_Index[k_row]]=InnerEdgeStat1[GBL_Matrix_ColNos[loop_i][loop_j]];
					LHS_Matrix_Index[k_row]++;
				}
				else
				{
					/* Construct the RHS Vector */
					RHSVector[k_row]=COMplex_Add(RHSVector[k_row],COMplex_Mul(ForcdValue[ForcdEdgeStat1[GBL_Matrix_ColNos[loop_i][loop_j]]],GBL_Matrix_Data[loop_i][loop_j]));
					if ((COMplex_Abs(Isource[loop_i])>TOLERANCE) && (flag==0))
					{
						/* Add "isource" term to the RHS vector */
						flag=1; /* Add the "isource" term only once */
						/*printf("Inside this isource loop %d %d\n",loop_i,No_cols);*/
						RHSVector[k_row]=COMplex_Add(RHSVector[k_row],Isource[loop_i]);
					}
				}/* end of else*/
			}   /* end of for loop_j*/
			fflush(stdout);
			k_row++;
		}      /* end of if  */
	}         /* end of for loop_i*/
	printf("%d %d\n",TotInnerEdgeNum,k_row);
}


/*************************************************************************
*   SUBROUTINE       : Conjugatesolver()
*************************************************************************/
 
/**********************************************
 FUNCTION     :Vnorm()
 Computes the Euclidean norm of a vector: |*|^2
***********************************************/

float Vnorm(complex* Vec,int Size)
{
 int II; float sum=0;
 for(II=0;II<Size;II++)
    {
    sum += COMplex_Abs(Vec[II])*COMplex_Abs(Vec[II]);
    }
 return sum;
}

/**********************************************
 FUNCTION     :Inner_Prod()
 Computes the Inner product of two vectors : <A*,B>
***********************************************/
 
complex Inner_Prod(complex* AVec,complex* BVec,int Size)
{
 int II; complex I_P;
 I_P = COMplex_Null();
 for(II=0;II<Size;II++)
    {
    I_P = COMplex_Add(I_P,COMplex_Mul(COMplex_Conjg(AVec[II]),BVec[II]));
    }
 return I_P;
} 

 
complex Inner_Prod_Without_Conjugating(complex* AVec,complex* BVec,int Size)
{
 int II; complex I_P;
 I_P = COMplex_Null();
 for(II=0;II<Size;II++)
    {
    I_P = COMplex_Add(I_P,COMplex_Mul(AVec[II],BVec[II]));
    }
 return I_P;
}  

/**********************************************
 FUNCTION     :MaxVecProd()
 Multiply a complex Matrix with a vector
***********************************************/

void MatrixVectorProd(char S,int Size,complex* XVec,complex* AVec,int* RIIndex,int** RICol,complex** RIDat)

{
int II,JJ,KK;

for(II=0;II<Size;II++) 
   {
   AVec[II]=COMplex_Null();
   }

if(S==' ') /* Multiply matrix with vector */
   {
   for(II=0;II<Size;II++)
      {
      for(KK=0;KK<RIIndex[II];KK++) 
         {
         JJ = RICol[II][KK];
         AVec[II] = COMplex_Add(AVec[II],COMplex_Mul(RIDat[II][KK],XVec[JJ]));
         }
      }
   }
else /* Multiply hermitian of matrix with vector */
     /* Here the matrix is symmetric, so we don't need to transpose */
     /* Only conjugate */
   {
   for(II=0;II<Size;II++)
      {
      for(KK=0;KK<RIIndex[II];KK++) 
         {
         JJ = RICol[II][KK];
         AVec[II] = COMplex_Add(AVec[II],COMplex_Mul(COMplex_Conjg(RIDat[II][KK]),XVec[JJ]));
         }
      }
   }
}

/****************************Main routine************************************/


int clog2(int address)
{
	float clog2 = (log(address*1.0f)/log(2.0f)) ;
	float clog22 = clog2-(int)clog2;
	if(clog22*clog22 <.0000000000001)
		return clog2;
	else return (int)(clog2+1);

}



void ConjugateSolver(int MatrixSize,int* RowIIndex,int** RowICol,complex** RowIDat,complex* RHSVec,int ITmax,float TOL)

{


int      II,JJ,IterIndex;
float   Residue,Vnrm;
int max_index_value = - 10;

for(int i = 0;i<MatrixSize;i++)
	{
		if(RowIIndex[i]>max_index_value) max_index_value = RowIIndex[i];
	}


A_matrix = fopen("A.txt","w");
A_floating_matrix = fopen("A_floating.txt","w");
b_new_matrix = fopen("B_new.txt", "w");
Ap_values = fopen("AP.txt", "w");
Multiples_Matrix = fopen("multiples_matrix.txt","w");
col_nos_matrix = fopen("col_nos.txt","w");
Parameters_matrix = fopen("Parameters.txt","w");
B_matrix = fopen("b.txt","w");
X_matrix =fopen("memx.txt","w");
B_floating_matrix = fopen("b_floating.txt","w");
b_new_floatin_matrix = fopen("b_new_floating.txt","w");


printf("\n %d ",freq_step);
if(freq_step == 0) /* Only for the first frequency step */
{
	X     = CMPLX_Vector(MatrixSize);
	P     = CMPLX_Vector(MatrixSize);
	PP    = CMPLX_Vector(MatrixSize);
	R     = CMPLX_Vector(MatrixSize);
	RR    = CMPLX_Vector(MatrixSize);
	A_P   = CMPLX_Vector(MatrixSize);
	AH_PP = CMPLX_Vector(MatrixSize);
 }

for(II=0;II<MatrixSize;II++) 
 {
   X[II]     = COMplex_Null();
   P[II]     = COMplex_Null();
   PP[II]    = COMplex_Null();
   R[II]     = COMplex_Null();
   RR[II]    = COMplex_Null();
   A_P[II]   = COMplex_Null();
   AH_PP[II] = COMplex_Null();
 } 


int counter = 0;
float zero = 0;
int unused = -1;
int col_nos_unused = -1;


// for debugging 




printf("\n\n MAX INDEX VALUE :: %d \n",max_index_value);

if(max_index_value>20)
{ printf("ERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\nERROR\n");}

// VERY IMPORTANT PARAMETERS FOR VERILOG 

int fixed_col_nos_value = 20;
int limit = (max_index_value>fixed_col_nos_value)?max_index_value:fixed_col_nos_value;
int number_of_row_by_vector_modules =4;
int no_of_units = number_of_row_by_vector_modules*2;
int no_of_elements_in_p_emap_output = 8 ;
int memory_height = 1000;
int memory_A_height = memory_height * (no_of_units/number_of_row_by_vector_modules) + no_of_units ;
int A_address_width = ((int) log(1.0*memory_height)/log(1.0*2))+4;

int A_file_index = 0;
int b_file_index = 0;

printf(" A_ADDRESS_WIDTH  :: %d \n",A_address_width);

printf(" :: %d \n ", clog2(1000));

int additional_b =((MatrixSize*1)%(no_of_units) !=0)? (no_of_units-((MatrixSize*1)%(no_of_units))):0;
printf("MATRIX SIZE IS %d",MatrixSize);

int additional_A = (((MatrixSize)%(number_of_row_by_vector_modules))!=0) ? (number_of_row_by_vector_modules-((MatrixSize)%(number_of_row_by_vector_modules))):0;

int total_with_additional_A = MatrixSize + additional_A;



float* matb = new float[MatrixSize];
float** mata;

mata = new float*[MatrixSize];

for (int i = 0; i < MatrixSize; i++)
{
	mata[i] = new float[limit];
}


 for (II=0;II<MatrixSize;II++) 
    { 
		int no_of_multiples_ii = ((RowIIndex[II]%no_of_elements_in_p_emap_output) ==0) ? (RowIIndex[II]/no_of_elements_in_p_emap_output) :((RowIIndex[II]/no_of_elements_in_p_emap_output)+1);
		fprintf(Multiples_Matrix,"%08X",no_of_multiples_ii);
		float testNum1 = 500.5 + II;
		fprintf(B_matrix,"%08X",*(unsigned int*)&RHSVec[II].x/*&testNum1*/);

		fprintf(X_matrix,"%08X",*(unsigned int*)&zero/*&testNum1*/);
		fprintf(B_floating_matrix,"%f ",RHSVec[II].x/*&testNum1*/);
		matb[II] = RHSVec[II].x/*&testNum1*/;
		//printf("mat_b is %f", matb[II]);
		
		counter=counter+1;
    for(JJ=0;JJ<limit;JJ++) 
      {
		  float testNum2 = 1000.5 +II+JJ;
		  if(JJ<RowIIndex[II])
		  {
			  fprintf(A_matrix,"%08X",*(unsigned int*)&RowIDat[II][JJ].x/*&testNum2*/);
			  fprintf(A_floating_matrix,"%f ",RowIDat[II][JJ].x/*&testNum2*/);
			  mata[II][JJ] =RowIDat[II][JJ].x/*&testNum2*/;
			  fprintf(col_nos_matrix,"%08X",RowICol[II][JJ]);
//////////////////////////////////////////////////////////////////////////
			  
		  }
		  else
		 {
			  fprintf(A_matrix,"%08X",*(unsigned int*)&zero);
			  fprintf(A_floating_matrix,"%f ",zero);
			  fprintf(col_nos_matrix,"%08X",col_nos_unused);
			  mata[II][JJ] = zero;
		  }

      } 
		if((II+1)%number_of_row_by_vector_modules ==0)
		{
			fprintf(A_matrix,"\n"); 
			fprintf(A_floating_matrix,"\n");
			fprintf(col_nos_matrix,"\n");
			fprintf(Multiples_Matrix,"\n");
			A_file_index ++;
		}
		if(counter%no_of_units==0)
		{
			b_file_index++;
			fprintf(B_matrix,"\n");
			fprintf(B_floating_matrix,"\n");
			fprintf(X_matrix,"\n");
		}
		if(
			((II+1)%number_of_row_by_vector_modules != 0) && (II == (MatrixSize-1))
			)
		{
			printf("\n :: %d :L \n",(number_of_row_by_vector_modules-(II+1)%number_of_row_by_vector_modules));
			for(int q = 0 ; q < (number_of_row_by_vector_modules-(II+1)%number_of_row_by_vector_modules) ; q++)
			{
				fprintf(Multiples_Matrix,"%08X",1);
			}

		for(int q = 0 ; q <limit*(number_of_row_by_vector_modules-(II+1)%number_of_row_by_vector_modules);q++)
			 {
	 		  fprintf(A_matrix,"%08X",*(unsigned int*)&zero);
			  fprintf(A_floating_matrix,"%f ",zero);
			  fprintf(col_nos_matrix,"%08X",col_nos_unused);
			}
		
			  A_file_index++;
			  fprintf(Multiples_Matrix,"%\n");
			  fprintf(A_matrix,"\n");
			  fprintf(A_floating_matrix,"\n");
			  fprintf(col_nos_matrix,"\n");
		
		}
    } 

 int III = II ;

	


 int third_additional = ((MatrixSize/number_of_row_by_vector_modules)%no_of_elements_in_p_emap_output!=0)
	 ?
	 (no_of_elements_in_p_emap_output-(MatrixSize/number_of_row_by_vector_modules)%no_of_elements_in_p_emap_output)
	 :0
	 ;

 printf("\n THIRD ADDTITIONAL IS :: %d \n\n",third_additional);

 fprintf(Parameters_matrix,"%08X\n",MatrixSize+additional_b);
 fprintf(Parameters_matrix,"%08X\n",MatrixSize+additional_b);



 for (int k = 0 ; k< third_additional ; k++)
 {
	 for(int q = 0 ; q <limit*number_of_row_by_vector_modules;q++)
	 {
	 		//  fprintf(A_matrix,"%08X",*(unsigned int*)&zero);
			//  fprintf(A_floating_matrix,"%f ",zero);
			//  fprintf(col_nos_matrix,"%08X",col_nos_unused);
	 }
	 for(int t = 0 ; t<number_of_row_by_vector_modules;t++)
	 {
		//fprintf(Multiples_Matrix,"%08X",1);
	 }
		//	 A_file_index++;
	   	 	// fprintf(A_matrix,"\n");
			  //fprintf(A_floating_matrix,"\n");
			  //fprintf(col_nos_matrix,"\n");
			  //fprintf(Multiples_Matrix,"\n");
 }

 /*
	fprintf(Multiples_Matrix,"\nIII IS %d \n",III);
		fprintf(A_matrix,"\n"); 
		fprintf(A_floating_matrix,"\n");
		fprintf(col_nos_matrix,"\n");
		fprintf(Multiples_Matrix,"\n");
		*/

 for (II = 0; II<MatrixSize; II++)
 {
	// printf("mat_b is %f", matb[II] );
	 //printf("\n");
     counter = counter + 1;
	 for (JJ = 0; JJ<limit; JJ++)
	 {
		 
		 if (JJ<RowIIndex[II])
		 {
			
			 fprintf(b_new_matrix, "%08X",*(unsigned int*)&matb[RowICol[II][JJ]]);
			 fprintf(b_new_floatin_matrix, "%f ", matb[RowICol[II][JJ]]);
			/* printf("start");
			 printf("\n");
			 printf("index is %d", index_b_new);
			 printf("\n");
			 printf("mat_b is %f", matb[index_b_new], "\n");
			 printf("\n");*/
			
			 //////////////////////////////////////////////////////////////////////////

		 }
		 else
		 {
			 fprintf(b_new_matrix, "%08X", *(unsigned int*)&zero);
			 fprintf(b_new_floatin_matrix, "%f ", zero);
		
		 }

	 }
	 if ((II + 1) % number_of_row_by_vector_modules == 0)
	 {
		 fprintf(b_new_matrix, "\n");
		 fprintf(b_new_floatin_matrix, "\n");
		 
	 }
	 
}

 
 float row_by_vector ;

 for (int II = 0; II<MatrixSize; II++)
 {
	 
	 row_by_vector = 0;
	 for (int KK = 0; KK<limit-1; KK++)
	 {
		 JJ = RowICol[II][KK];
		 row_by_vector = row_by_vector + mata[II][KK] * matb[JJ]; //COMplex_Add(AVec[II], COMplex_Mul(RIDat[II][KK], XVec[JJ]));
	 }
	 fprintf(Ap_values, "%08x", *(unsigned int*)& row_by_vector);
	 fprintf(Ap_values, "  ");
	 if ((II + 1) % number_of_row_by_vector_modules == 0)
	 fprintf(Ap_values, "\n");
 }



 printf("additional_b is %d",additional_b);
 printf("additional_A is %d",additional_A);
 if(additional_b!=0)
 {
 for(int i=0;i<additional_b;i++)
	{
	 fprintf(B_matrix,"%08X",*(unsigned int*)&zero);
	 fprintf(X_matrix,"%08X",*(unsigned int*)&zero);
	 fprintf(B_floating_matrix,"%f ",zero);
	}
 }

 printf("\n\n\nADDITIONAL A IS %d\n\n",additional_A);

 printf("\n\n\nA file index ::  %d\n\n",A_file_index);
 printf("\n\n\nB file index ::  %d\n\n",b_file_index);



 while(A_file_index < 2*(1+b_file_index))
	 {
		 for(int k =0;k<limit*number_of_row_by_vector_modules ; k++)
		{
			fprintf(A_matrix,"%08X",*(unsigned int*)&zero);
			fprintf(A_floating_matrix,"%f ",zero);
			fprintf(col_nos_matrix,"%08X",col_nos_unused);
		}
		 for(int k =0;k<number_of_row_by_vector_modules ; k++)
		 {fprintf(Multiples_Matrix,"%08X",1);}

		 A_file_index ++;
	 }
 


    fclose(A_matrix);
	fclose(A_floating_matrix);
	fclose(b_new_matrix);
	fclose(Ap_values);
	fclose(Multiples_Matrix);
	fclose(col_nos_matrix);
	fclose(B_matrix);
	fclose(X_matrix);
	fclose(B_floating_matrix);
	fclose(b_new_floatin_matrix);

	


/* Initialization */

MatrixVectorProd(' ',MatrixSize,X,R,RowIIndex,RowICol,RowIDat);

for(II=0;II<MatrixSize;II++) 
   {
   R[II]  = COMplex_Sub(RHSVec[II],R[II]);
   P[II]  = R[II];
   RR[II] = COMplex_Conjg(R[II]);
   PP[II] = RR[II];
   } 
//printf("**********************************");
//for(II=0;II<19;II++) 
//   {
//   printf(" R = %f %f",R[II].x,R[II].y);
//   printf(" RR = %f %f \n ",RR[II].x,RR[II].y);
//   printf(" P = %f %f",P[II].x,P[II].y);
//   printf(" PP = %f %f \n",PP[II].x,PP[II].y);
//   } 
//
//printf("***********************************");

Vnrm = Vnorm(RHSVec,MatrixSize);

IterIndex = 0;

/* Actual Iteration */
//cuComplex a =  make_cuComplex(1,0);
//cuComplex b = make_cuComplex(0,0);

for(;;)
   {

   MatrixVectorProd(' ',MatrixSize,P,A_P,RowIIndex,RowICol,RowIDat);

   for(II=0;II<MatrixSize;II++)
      {
		/* fprintf(Ap_values, "%08x", *(unsigned int*)&A_P[II].x);
		 fprintf(Ap_values, "\n");*/
      AH_PP[II] = COMplex_Conjg(A_P[II]); 
      }

  Alpha = Inner_Prod(RR,R,MatrixSize);
   temp_var = Inner_Prod(PP,A_P,MatrixSize);

   float magnOfbast = sqrt(Alpha.x*Alpha.x +Alpha.y*Alpha.y);
   printf("\n  bast.x : %08X ,  bast.y : %08X  ,bast : %08X   magnOfbast : %08X  \n ",  *(unsigned int*)&Alpha.x, *(unsigned int*)&Alpha.y, *(unsigned int*)&Alpha,*(unsigned int*)&magnOfbast );
   printf("\n  bast.x : %f ,  bast.y : %f  ,bast : %f   magnOfbast : %f  \n ", Alpha.x,Alpha.y,Alpha, magnOfbast);
	
   float magnOfMaqam = sqrt(temp_var.x*temp_var.x +temp_var.y*temp_var.y);
   printf("\n  maqam. x : %08X ,  maqam.y : %08X  ,maqam : %08X   magnOfmaqam : %08X  \n ", *(unsigned int*)&temp_var.x,*(unsigned int*)&temp_var.y,*(unsigned int*)&temp_var,*(unsigned int*)&magnOfMaqam );
   printf("\n  maqam. x : %f ,  maqam.y : %f  ,maqam : %f   magnOfmaqam : %f  \n ", temp_var.x,temp_var.y,temp_var, magnOfMaqam);
   Alpha = COMplex_Div(Alpha,temp_var);  

   float magnOfAlpha = sqrt(Alpha.x*Alpha.x +Alpha.y*Alpha.y) ;
   printf("\n  Alpha.x : %08X ,  Alpha.y : %08X  ,Alpha : %08X   magnOfAlpha : %08X  \n ", *(unsigned int*)&Alpha.x,*(unsigned int*)&Alpha.y,*(unsigned int*)&Alpha, *(unsigned int*)&magnOfAlpha);
   printf("\n  Alpha.x : %f ,  Alpha.y : %f  ,Alpha : %f   magnOfAlpha : %f  \n ", Alpha.x,Alpha.y,Alpha, magnOfAlpha);
 
   /* Alpha is the step length parameter */
   for(II=0;II<MatrixSize;II++)  
      {
      /* New Solution Estimate */
      X[II]  = COMplex_Add(X[II],COMplex_Mul(Alpha,P[II]));   
      R[II]  = COMplex_Sub(R[II],COMplex_Mul(Alpha,A_P[II])); /* Residual */
      RR[II] = COMplex_Sub(RR[II],COMplex_Mul(COMplex_Conjg(Alpha),AH_PP[II]));
      /* Bi-residual */
      }

   //Beta = Inner_Prod(AH_PP,R,MatrixSize);
   Beta = Inner_Prod_Without_Conjugating(A_P,R,MatrixSize);


   Beta = COMplex_Div(Beta,temp_var);
   Beta = Real_Mul(-1.0,Beta);   /* Beta is the Bi-Conjugacy coefficient  */

   for(II=0;II<MatrixSize;II++)  
      {
      P[II]  = COMplex_Add(R[II],COMplex_Mul(Beta,P[II]));   /* Direction */
      PP[II] = COMplex_Add(RR[II],COMplex_Mul(COMplex_Conjg(Beta),PP[II])); 
      /* Bi-Direction */
      }

   IterIndex = IterIndex + 1; /* Next iteration */
   Residue = sqrt(Vnorm(R,MatrixSize)/Vnrm);
    printf("Iteration = %d Residue = %e TOL = %e\n",IterIndex,Residue,TOL); 
	

	float magnOfBeta = sqrt(Beta.x*Beta.x +Beta.y*Beta.y) ;
   printf("\n  Beta.x : %08X ,  Beta.y : %08X  ,Beta : %08X   magnOfBeta : %08X  \n ", *(unsigned int*)&Beta.x,*(unsigned int*)&Beta.y,*(unsigned int*)&Beta, *(unsigned int*)&magnOfBeta);
   printf("\n  Beta.x : %f ,  Beta.y : %f  ,Beta : %f   magnOfBeta : %f  \n ", Beta.x,Beta.y,Beta, magnOfBeta);


/*    fflush(stdout); */
   if(Residue<=TOL)  /* Test termination condition */
      {
      printf("Convergence achieved\n");
      printf("Iteration = %d Residue = %e TOL = %e\n",IterIndex,Residue,TOL);
      fflush(stdout);
      break;
      }
   else if(IterIndex==ITmax)  /* No. of iterations has exceeded the limit */
      { 
      printf("Iteration Exceeds\n"); 
      fflush(stdout);
      break;
      }
   else 
      { 
      continue ; 
      }
   }/* End of for loop*/
for(II=0;II<MatrixSize;II++) 
   {
   RHSVec[II] = X[II];
   Old_Solution[II] = X[II]; /* Starting point for next iteration */
   }
}
/***************************************************************************
*   FUNCTION        : Read_Input_Pass_1()
****************************************************************************/
void Read_Input_Pass_1()
{
 int    x1, y1, z1, x2, y2, z2, p, p1, p2, OutF_Count=0, output_flag=0;
 char   type[20], att1[20], att2[18], att3[18], att4[18];
 float  tmp1, tmp2;
 float freq = 0;
 char   buffer[80];

 Min_X = 10000; Min_Y=10000; Min_Z=10000;
 Max_X = 0;     Max_Y=0;     Max_Z=0;

   while (fgets(buffer, 80, InF))
   {
        if (buffer[0] == '#') continue; /* Comment statement */
	if(!strncmp(buffer, "freqstep",8)) /* freqstep keyword */
        	   {
		   sscanf(buffer, "%s%d%s",type, &NUM_OF_STEPS, att1);
		   FREQ_INC = atof(att1); /* Convert from string to float */
		   printf("N0. of steps%d inc%f\n",NUM_OF_STEPS,FREQ_INC);
		   continue;
 		   }
        if(!strncmp(buffer, "default_output",14))/*keyword for default output*/
                   {
                   sscanf(buffer, "%s%s", type, Out_FileName0);
                   OutF_0=fopen(Out_FileName0, "w");
                   if (OutF_0 == NULL) 
		     {
		     fprintf(stderr, 
		     " Error: Can't create default output file %s\n",
		     Out_FileName0); 
		     exit(1);
		     }
                   output_flag=1;
                   continue;
                   }
        if(!strncmp(buffer, "celldim", 7)) { /* celldim keyword */
                                              continue;
                                            }


        if(sscanf(buffer, "%s%d%d%d%d%d%d%s%s%s%s", 
		  type, &x1, &y1, &z1, &x2, &y2, &z2, 
		  att1, att2, att3, att4)>2)
           {
              if(!strcmp(type,"boundary")) {
		fprintf(stderr, 
"\n\nERROR: boundary keyword is not allowed by EMAP4 \n"); 
		exit(1);}
              if(!strcmp(type, "dielectric"))/* dielectric block */
              {
               if (x1<Min_X) Min_X=x1; 
	       if (x2<Min_X) Min_X=x2; 
	       if (x1>Max_X) Max_X=x1; 
	       if (x2>Max_X) Max_X=x2;
               if (y1<Min_Y) Min_Y=y1; 
	       if (y2<Min_Y) Min_Y=y2; 
	       if (y1>Max_Y) Max_Y=y1; 
	       if (y2>Max_Y) Max_Y=y2;
               if (z1<Min_Z) Min_Z=z1; 
	       if (z2<Min_Z) Min_Z=z2; 
	       if (z1>Max_Z) Max_Z=z1; 
	       if (z2>Max_Z) Max_Z=z2;
               if((x1==x2)||(y1==y2)||(z1==z2)) 
		 {fprintf(stderr, "ERROR:  diel thickness is not finite./n");
               exit(1); } continue;
              }
              if(!strcmp(type, "PML")) /* Perfectly matched layer */
              {
               if (x1<Min_X) Min_X=x1; 
	       if (x2<Min_X) Min_X=x2; 
	       if (x1>Max_X) Max_X=x1; 
	       if (x2>Max_X) Max_X=x2;
               if (y1<Min_Y) Min_Y=y1; 
	       if (y2<Min_Y) Min_Y=y2; 
	       if (y1>Max_Y) Max_Y=y1; 
	       if (y2>Max_Y) Max_Y=y2;
               if (z1<Min_Z) Min_Z=z1; 
	       if (z2<Min_Z) Min_Z=z2; 
	       if (z1>Max_Z) Max_Z=z1; 
	       if (z2>Max_Z) Max_Z=z2;
               if((x1==x2)||(y1==y2)||(z1==z2)) 
		 {
		 fprintf(stderr, "ERROR:  PML thickness is not finite./n");
                 exit(1); 
	         } continue;
              }
              if(!strcmp(type, "box")) /* rectangular conducting box */
              {
                   if((x1==x2)||(y1==y2)||(z1==z2)){ 
		     fprintf(stderr, "ERROR: box dimension is invalid.");
                   exit(1);}
                   if (x1<Min_X) Min_X=x1; 
		   if (x2<Min_X) Min_X=x2; 
		   if (x1>Max_X) Max_X=x1; 
		   if (x2>Max_X) Max_X=x2;
                   if (y1<Min_Y) Min_Y=y1; 
		   if (y2<Min_Y) Min_Y=y2; 
		   if (y1>Max_Y) Max_Y=y1; 
		   if (y2>Max_Y) Max_Y=y2;
                   if (z1<Min_Z) Min_Z=z1; 
		   if (z2<Min_Z) Min_Z=z2; 
		   if (z1>Max_Z) Max_Z=z1; 
		   if (z2>Max_Z) Max_Z=z2;
                   continue;
              }
              if(!strcmp(type,"conductor")) /* PEC */
              {
                   if (x1<Min_X) Min_X=x1; 
		   if (x2<Min_X) Min_X=x2; 
		   if (x1>Max_X) Max_X=x1; 
		   if (x2>Max_X) Max_X=x2;
                   if (y1<Min_Y) Min_Y=y1; 
		   if (y2<Min_Y) Min_Y=y2; 
		   if (y1>Max_Y) Max_Y=y1; 
		   if (y2>Max_Y) Max_Y=y2;
                   if (z1<Min_Z) Min_Z=z1; 
		   if (z2<Min_Z) Min_Z=z2; 
		   if (z1>Max_Z) Max_Z=z1; 
		   if (z2>Max_Z) Max_Z=z2;
                   continue;
              }
              if(!strcmp(type,"resistor")) /* resistor */
              {
		if (x1<Min_X) Min_X=x1; 
		if (x2<Min_X) Min_X=x2; 
		if (x1>Max_X) Max_X=x1; 
		if (x2>Max_X) Max_X=x2;
                if (y1<Min_Y) Min_Y=y1; 
		if (y2<Min_Y) Min_Y=y2; 
		if (y1>Max_Y) Max_Y=y1; 
		if (y2>Max_Y) Max_Y=y2;
                if (z1<Min_Z) Min_Z=z1; 
		if (z2<Min_Z) Min_Z=z2; 
		if (z1>Max_Z) Max_Z=z1; 
		if (z2>Max_Z) Max_Z=z2;
                   continue;
              }
              if(!strcmp(type, "aperture")) /* aperture */
              {
                   if((x1!=x2) & (y1!=y2) & (z1!=z2))
		     {
		     fprintf(stderr,"ERROR: aperture dimension is invalid.\n");
                     exit(1); 
		     }
                   if(((x1==x2)&&(y1==y2))||((x1==x2)&&(z1==z2))
		                          ||((z1==z2)&&(y1==y2)))
                   {
		   fprintf(stderr, "ERROR: aperture dimension is invalid.\n"); 
		   exit(1); }
                   continue;
              }
              if(!strcmp(type, "seam")) continue; /* Not available */
              if(!strcmp(type,"esource")) /* electric field source */
              {
                   if (x1<Min_X) Min_X=x1; 
		   if (x2<Min_X) Min_X=x2; 
		   if (x1>Max_X) Max_X=x1; 
		   if (x2>Max_X) Max_X=x2;
                   if (y1<Min_Y) Min_Y=y1; 
		   if (y2<Min_Y) Min_Y=y2; 
		   if (y1>Max_Y) Max_Y=y1; 
		   if (y2>Max_Y) Max_Y=y2;
                   if (z1<Min_Z) Min_Z=z1; 
		   if (z2<Min_Z) Min_Z=z2; 
		   if (z1>Max_Z) Max_Z=z1; 
		   if (z2>Max_Z) Max_Z=z2;
                   freq = atof(att1);
                   if((x1==x2)&&(y1==y2)&&(z1==z2)) {
                   fprintf(stderr, 
			   "ERROR: esource cannot be defined at a point.\n"); 
		   exit(1); }
                   if ((att2[0] != 'x') &&(att2[0] != 'y')&& (att2[0] != 'z'))
		     {
                   fprintf(stderr, 
			   "ERROR: unrecognized polarization of esource statement  %c\n",att2[0]);
                   exit(1); }
                   continue;
              }
              if(!strcmp(type,"coax")) /* coaxial cable interface to system  */
/* 				       under construction */
              {
                   if (x1<Min_X) Min_X=x1; 
		   if (x2<Min_X) Min_X=x2; 
		   if (x1>Max_X) Max_X=x1; 
		   if (x2>Max_X) Max_X=x2;
                   if (y1<Min_Y) Min_Y=y1; 
		   if (y2<Min_Y) Min_Y=y2; 
		   if (y1>Max_Y) Max_Y=y1; 
		   if (y2>Max_Y) Max_Y=y2;
                   if (z1<Min_Z) Min_Z=z1; 
		   if (z2<Min_Z) Min_Z=z2; 
		   if (z1>Max_Z) Max_Z=z1; 
		   if (z2>Max_Z) Max_Z=z2;
                   freq = atof(att1);
                   if((x1==x2)&&(y1==y2)&&(z1==z2)) {
                   fprintf(stderr, 
			"ERROR: coax source cannot be defined at a point.\n");
		   exit(1); }
                   continue;
              }
              if(!strcmp(type,"isource")) /* linear current source */
              {
                   if (x1<Min_X) Min_X=x1; 
		   if (x2<Min_X) Min_X=x2; 
		   if (x1>Max_X) Max_X=x1; 
		   if (x2>Max_X) Max_X=x2;
                   if (y1<Min_Y) Min_Y=y1; 
		   if (y2<Min_Y) Min_Y=y2; 
		   if (y1>Max_Y) Max_Y=y1; 
		   if (y2>Max_Y) Max_Y=y2;
                   if (z1<Min_Z) Min_Z=z1; 
		   if (z2<Min_Z) Min_Z=z2; 
		   if (z1>Max_Z) Max_Z=z1; 
		   if (z2>Max_Z) Max_Z=z2;
                   freq = atof(att1);
                   if((x1==x2)&&(y1==y2)&&(z1==z2)) {
                   fprintf(stderr, 
			"ERROR: isource cannot be defined at a point.\n"); 
		   exit(1); }
                   if ((att2[0] != 'x')&&(att2[0] != 'y') && (att2[0] != 'z'))
		     {
                      fprintf(stderr, 
			   "ERROR: unrecognized polarization of isource statement  %c\n",att2[0]);
                   exit(1); }
                   continue;
              }
              if(!strcmp(type, "msource")) /*not available in EMAP4 */
              {
                   if (x1<Min_X) Min_X=x1; 
		   if (x2<Min_X) Min_X=x2; 
		   if (x1>Max_X) Max_X=x1; 
		   if (x2>Max_X) Max_X=x2;
                   if (y1<Min_Y) Min_Y=y1; 
		   if (y2<Min_Y) Min_Y=y2; 
		   if (y1>Max_Y) Max_Y=y1; 
		   if (y2>Max_Y) Max_Y=y2;
                   if (z1<Min_Z) Min_Z=z1; 
		   if (z2<Min_Z) Min_Z=z2; 
		   if (z1>Max_Z) Max_Z=z1; 
		   if (z2>Max_Z) Max_Z=z2;
                   freq = atof(att1);
                   continue;
              }
              if(!strcmp(type, "efield_output")) /* specified output */
              {    
                   switch (OutF_Count)
                        {
                        case 0:  
			{OutF_Count=1; strcpy(Out_FileName1, att1); 
			 OutF_1=fopen(Out_FileName1,"w"); break;}
                        case 1:  
			{if(strcmp(att1, Out_FileName1))
			   {OutF_Count=2; strcpy(Out_FileName2, att1); 
			    OutF_2=fopen(Out_FileName2,"w");}
                                  break;
                                 }
                        case 2:  
			{if(strcmp(att1, Out_FileName1) 
			    && strcmp(att1, Out_FileName2))
                              {OutF_Count=3; strcpy(Out_FileName3, att1); 
			       OutF_3=fopen(Out_FileName3,"w");}
                                  break;
                                 }
                        case 3:  
			{if((strcmp(att1, Out_FileName1) 
			    && strcmp(att1, Out_FileName2)) 
			    && strcmp(att1, Out_FileName3))
                            {fprintf(stderr, "ERROR: too many different files specified by efield_output commands \n"); 
			     exit(1);}
                                  break;
                                 }
                        }
                   output_flag=1;
                   continue;
              }
              fprintf(stderr,"Warning: Unrecognized Keyword %s\n", type);
           }                 /*end of if*/
           OperateFreq = freq;
    }                    /* end of while */
/* The input file has to contain at least one source statement */
    if(freq == 0) fprintf(stderr,"\nWARNING: The input file does not specify a source, code will not produce output!\n\n");
/* The input file has to contain at one output statement */
    if(output_flag == 0) fprintf(stderr,"\nWARNING: The input file does not contain an output statement, code will not produce output!\n\n");
    rewind(InF);

/* The file is read the first time so that we can calculate how many cell 
blocks there are and the second time it initializes all the cell dimensions  */

/* Calculate the dimensions to be meshed */
    XdiM = Max_X - Min_X; 
    YdiM = Max_Y - Min_Y; 
    ZdiM = Max_Z - Min_Z;
/*     printf(" DiM %d %d %d\n",XdiM,YdiM,ZdiM); */
    DivisorX = FLOAT_Vector(XdiM+2);
    DivisorY = FLOAT_Vector(YdiM+2);
    DivisorZ = FLOAT_Vector(ZdiM+2);

/* Initialize cell dimensions */

    for (p=0;p<XdiM+2;p++)
       {
       DivisorX[p]=0.01;
       }
    for (p=0;p<YdiM+2;p++)
       {
       DivisorY[p]=0.01;
       }
    for (p=0;p<ZdiM+2;p++)
       {
       DivisorZ[p]=0.01;
       }

    OperateFreq = freq;

   while (fgets(buffer, 80, InF))
   {
        if (buffer[0] == '#') continue;

        if(!strncmp(buffer, "default_output",14))
                   {
                   continue;
                   }
/* Calculate cell dimensions using information from the input file */
/* tmp1 is the value specified in the input file, tmp2 is the multiplier 
   due to the unit specified . A non-uniform mesh implies that we can 
   assign different values in different directions for each brick 
   (hexahedron)*/

        if(!strncmp(buffer, "celldim", 7)) { 
                                              sscanf(buffer, "%s%s%d%d%g%s", 
					      type, att2,&p1, &p2,&tmp1, att1);
        				if(!strcmp(att2, "x")) /* x-axis */

						{
                                              if (!strcmp(att1, "cm")) 
						      tmp2=0.01;
                                              else if(!strcmp(att1, "m"))
                                                      tmp2=1;
                                              else if(!strcmp(att1, "mm"))
                                                      tmp2=0.001;
                                              else if(!strcmp(att1, "in"))
                                                      tmp2=0.0254;
                                              else {fprintf(stderr,"ERROR: unrecognized cell dimension.\n");
                                                    exit(1);}
					      for ( p=p1;p<p2;p++)
						 {
                                                 DivisorX[p+1]=tmp1*tmp2;
				    /* printf(" %d %f\n", p,DivisorX[p+1]); */}
                                              continue;
                                            }
        				if(!strcmp(att2, "y")) /* y-axis */
						{ 
                                              if (!strcmp(att1, "cm")) 
						      tmp2=0.01;
                                              else if(!strcmp(att1, "m"))
                                                      tmp2=1;
                                              else if(!strcmp(att1, "mm"))
                                                      tmp2=0.001;
                                              else if(!strcmp(att1, "in"))
                                                      tmp2=0.0254;
                                              else {fprintf(stderr,"ERROR: unrecognized cell dimension.\n");
                                                    exit(1);}
					      for ( p=p1;p<p2;p++)
						 {
                                                 DivisorY[p+1]=tmp1*tmp2;
					  /* printf(" %d %f\n", p,DivisorY[p+1]); */}
                                              continue;
                                            }
        				if(!strcmp(att2, "z")) /* z-axis */
						{ 
                                              if (!strcmp(att1, "cm")) 
						      tmp2=0.01;
                                              else if(!strcmp(att1, "m"))
                                                      tmp2=1;
                                              else if(!strcmp(att1, "mm"))
                                                      tmp2=0.001;
                                              else if(!strcmp(att1, "in"))
                                                      tmp2=0.0254;
                                              else {fprintf(stderr,"ERROR: unrecognized cell dimension.\n");
                                                    exit(1);}
					      for ( p=p1;p<p2;p++)
						 {
                                                 DivisorZ[p+1]=tmp1*tmp2;
						 }
                                              continue;
                                            }
	                                    }

        if(sscanf(buffer, "%s%d%d%d%d%d%d%s%s%s%s", 
		type, &x1, &y1, &z1, &x2, &y2, &z2, att1, att2, att3, att4)>2)
           {
              if(!strcmp(type,"boundary")) {
		fprintf(stderr, "\n\nERROR: boundary keyword is not allowed by EMAP3 \n"); 
		exit(1);}
              if(!strcmp(type, "dielectric"))
              {
              continue;
              }
              if(!strcmp(type, "PML"))
              {
              continue;
              }
              if(!strcmp(type, "box"))
              {
                   continue;
              }
              if(!strcmp(type,"conductor"))
              {
                   continue;
              }
              if(!strcmp(type,"resistor"))
              {
                   continue;
              }
              if(!strcmp(type, "aperture"))
              {
                   continue;
              }
              if(!strcmp(type, "seam")) continue;
              if(!strcmp(type,"esource"))
              {
                   continue;
              }
              if(!strcmp(type,"isource"))
              {
                   continue;
              }
              if(!strcmp(type,"coax"))
              {
                   continue;
              }
              if(!strcmp(type, "msource"))
              {
                  continue;
              }
              if(!strcmp(type, "efield_output"))
              {    
                   continue;
              }
           }                 /*end of if*/

    }                    /* end of while */
    rewind(InF);
  for (p=1;p<XdiM+1;p++) 
    { 
    printf("%f\n",DivisorX[p]); 
    } 
 for (p=1;p<YdiM+1;p++) 
    { 
    printf("%f\n",DivisorY[p]); 
    } 
 for (p=1;p<ZdiM+1;p++) 
    { 
    printf("%f\n",DivisorZ[p]); 
 
    } 
}

/***************************************************************************
*   FUNCTION        : Read_Input_Pass_2()

This routine has two main functions: 

1) Identify the type of each edge (inner, boundary, forced) and once identified
   initialize the relevant variables:
   a) ForceStat : 1 for forced and 0 for unforced
   b) EdgeStatus: -1 for inner edges, -2 for boundary edges and the value 
                  of the electric field (real part) for forced edges.
   c) EdgeStatus1: The value of the imaginary part of the electric field.
   d) Resist    : The value of the resistor term (to be added to the final 
                  matrix) if the edge is a resistor.
   e) Isource   : The value of the isource term (to be added to the final 
                  RHS vector) if the edge is an isource.
   f) for coax edges, both Resist and Isource are initialized.

2) Identify the type of hexahedron (dielectric, PML etc) and initialize the 
   corresponding variables (a_PML_vec, b_PML_vec, c_PML_vec, Epsilon, Sigma).

FOR EDGES:
   The routine reads each keyword and reads the coordinates that need to be 
initialized. It then identifies the hexahedron to which the belongs to. It then
identifies the edge in the hexhedron which to initialize. This is tricky (as 
seen from the length of this routine), because the edge can belong to many 
hexahedra. THe default hexahedra is the one with the lowest number. But if 
the edge is on the boundary, then the hexahedron changes. Having 
figured which edge of which hexahedron the coordinates point to, it then 
initializes the corresponding variables.

THINGS DONE FOR EACH KEYWORD :

1) box : Takes the hollow box and divides into six faces and initializes the 
   edges on these faces one face at a time.
2) conductor : Identifies whether the conductor is 1-D, 2-D or 3-D and proceeds
   from there. The 3-D cond. is assumed to be solid (as opposed to the hollow 
   "box").
3) dielectric : Initializes Epsilon and Sigma for the hexahedron 
   (no edges involved).
4) PML : Initializes a_PML_vec, b_PML_vec, c_PML_vec for each hexahedron. 
   For more on the direction attributes, see user's guide. 
5) aperture : The aperture is usually in a conducting plate. The aperture 
   edges are converted to "inner edges". The important thing here is that the 
   aperture boundary remains conducting.
6) resistor: Resist is initialized for each edge. The edge is made an 
   "inner edge". 
7) isource: Isource is initialized for each edge. The edge is made an "inner
   edge".
8) esource: EgdeStatus and EdgeStatus1 are initialized for each edge. These are
   combined in the routine "CountEdgeType". The edge is now a "forced edge". 

POINTS TO NOTE:

1) The resistor and isource keywords are expected to be 1-D. 
2) In the isource and esource, only the edges in the direction specified in 
   the input file are intialized.

****************************************************************************/

void Read_Input_Pass_2()
{
 void swap(int*,int*); /* For these routines see documentation near the routine */
 void edge(int,int,int,int);
 void deter_edge(int,int,int,int);
 void CountEdgeType();
 int    i, j, k, l, x1, y1, z1, x2, y2, z2, nn, tmp;
 char   type[20], att1[20], att2[18], att3[18], att4[18];
 char   buffer[80];
 float length;
 complex temp1_PML,temp2_PML,temp3_PML,temp4_PML,temp5_PML; 

/* Initialization */

 if(freq_step==0){
        /* Dynamically allocate memory for large variables */
        Is_PML= INT_Vector(TotNumHexHedron);
        Epsilon= FLOAT_Vector(TotNumHexHedron);
        a_PML_vec= CMPLX_Vector(TotNumHexHedron);
        b_PML_vec= CMPLX_Vector(TotNumHexHedron);
        c_PML_vec= CMPLX_Vector(TotNumHexHedron);
        Sigma= FLOAT_Vector(TotNumHexHedron);
        EdgeStatus   = FLOAT_Vector(TotEdgeNum);
        EdgeStatus1  = FLOAT_Vector(TotEdgeNum);
        Resist       = CMPLX_Vector(TotEdgeNum);
        Isource      = CMPLX_Vector(TotEdgeNum);
        ForceStat= INT_Vector( TotEdgeNum);
		/*xn = INT_Matrix (TotNumHexHedron,18);*/

		xn = new int*[TotNumHexHedron];
		for(int kkk = 0; kkk < TotNumHexHedron; ++kkk)
			xn[kkk] = new int[18];

        }
        /*Initialize dynamically allocated variables */
        for (i=0;i<=TotEdgeNum;i++) { ForceStat[i]=0;
                                      EdgeStatus[i]=0.0;
                                      EdgeStatus1[i]=0.0;
				      Resist[i]=COMplex_Null();
				      Isource[i]=COMplex_Null();
				    }
        for(i=0;i<XdiM;i++) for(j=0;j<YdiM;j++) for(k=0;k<ZdiM;k++)
		{
			nn=k*YdiM*XdiM+j*XdiM+i; Epsilon[nn]=1.0; Sigma[nn]=0.0;Is_PML[nn]=0;
			edge(nn,i,j,k);  
			deter_edge(nn,i,j,k);
		}
        
        /*Begin reading the input file for the second time */
        while (fgets(buffer, 80, InF))
        {
            if (buffer[0] == '#') continue;
            if(sscanf(buffer, "%s%d%d%d%d%d%d%s%s%s%s", 
            type, &x1, &y1, &z1, &x2, &y2, &z2, att1, att2, att3, att4)>2)
             {
             if(!strcmp(type, "boundary"))
             {
                   continue;
             }
             if(!strcmp(type, "efield_output"))
             {     
                   if((x1<Min_X || x1>Max_X)||(x2<Min_X || x2>Max_X)) {fprintf(stderr,"ERROR: efield_output parameter out of range.\n"); exit(1);}
                   if((y1<Min_Y || y1>Max_Y)||(y2<Min_Y || y2>Max_Y)) {fprintf(stderr,"ERROR: efield_output parameter out of range.\n"); exit(1);}
                   if((z1<Min_Z || z1>Max_Z)||(z2<Min_Z || z2>Max_Z)) {fprintf(stderr,"ERROR: efield_output parameter out of range.\n"); exit(1);}
                   continue;
             }
             if(!strcmp(type, "dielectric"))
             {
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;
                   for(i=x1;i<x2;i++){
                      for(j=y1;j<y2;j++){
                         for(k=z1;k<z2;k++){
                         nn=k*XdiM*YdiM+j*XdiM+i;
                         Epsilon[nn]=atof(att1); /* Relative permittivity */
                         Sigma[nn]=atof(att2); /* conductivity */
                                            }
                                         }
                                      } 
                   continue;
             }
             if(!strcmp(type, "PML"))
             {
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;
		   temp2_PML = COMplex_Cmplx(atof(att2),atof(att3));
		   temp1_PML = COMplex_Div(COMplex_Cmplx(1.0,0.0),temp2_PML);
		   temp3_PML = COMplex_Mul(temp2_PML,temp2_PML);
		   temp4_PML = COMplex_Mul(temp1_PML,temp2_PML);
		   temp5_PML = COMplex_Mul(temp3_PML,temp1_PML);
		   /* printf("%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n\n",
		           a_PML.x,a_PML.y,b_PML.x,b_PML.y,c_PML.x,c_PML.y); */
                   for(i=x1;i<x2;i++){
                      for(j=y1;j<y2;j++){
                         for(k=z1;k<z2;k++){
                         nn=k*XdiM*YdiM+j*XdiM+i;
                         Is_PML[nn]=1;
/* For more information on the direction attribute, see the PML chapter in the  */
/* user's guide */			 
		   if(!strcmp(att1, "x"))
			{
			a_PML_vec[nn] = temp1_PML;
			b_PML_vec[nn] = temp2_PML;
			c_PML_vec[nn] = temp2_PML;
			}
		   if(!strcmp(att1, "y"))
			{
			b_PML_vec[nn] = temp1_PML;
			a_PML_vec[nn] = temp2_PML;
			c_PML_vec[nn] = temp2_PML;
			}
		   if(!strcmp(att1, "z"))
			{
			c_PML_vec[nn] = temp1_PML;
			b_PML_vec[nn] = temp2_PML;
			a_PML_vec[nn] = temp2_PML;
			}
		   if(!strcmp(att1, "xy"))
			{
			c_PML_vec[nn] = temp3_PML;
			b_PML_vec[nn] = temp4_PML;
			a_PML_vec[nn] = temp4_PML;
			}
		   if(!strcmp(att1, "yz"))
			{
			c_PML_vec[nn] = temp4_PML;
			b_PML_vec[nn] = temp4_PML;
			a_PML_vec[nn] = temp3_PML;
			}
		   if(!strcmp(att1, "xz"))
			{
			c_PML_vec[nn] = temp4_PML;
			b_PML_vec[nn] = temp3_PML;
			a_PML_vec[nn] = temp4_PML;
			}
		   if(!strcmp(att1, "xyz"))
			{
			c_PML_vec[nn] = temp5_PML;
			b_PML_vec[nn] = temp5_PML;
			a_PML_vec[nn] = temp5_PML;
			}
                                            }
                                         }
                                      } 
                   continue;
             }
             if(!strcmp(type, "box"))
             {
              if (x2<x1) swap(&x1, &x2);
              if (y2<y1) swap(&y1, &y2);
              if (z2<z1) swap(&z1, &z2);
              x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
              x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;
              /* left plate */
                           for(k=z1;k<z2;k++){
                             for(j=y1;j<y2;j++){
                                                nn=x2*y2*k+j*x2;
                                                ForceStat[xn[nn][1]]=1;
                                                EdgeStatus[xn[nn][1]]=0.0;
                                                ForceStat[xn[nn][5]]=1;
                                                EdgeStatus[xn[nn][5]]=0.0;
                                                ForceStat[xn[nn][8]]=1;
                                                EdgeStatus[xn[nn][8]]=0.0;
                                                ForceStat[xn[nn][10]]=1;
                                                EdgeStatus[xn[nn][10]]=0.0;
                                                ForceStat[xn[nn][14]]=1;
                                                EdgeStatus[xn[nn][14]]=0.0;
                                               }
                                              }
              /* right plate */
                           for(k=z1;k<z2;k++){
                             for(j=y1;j<y2;j++){
                                                nn=x2*y2*k+(x2-1)+j*x2;
                                                ForceStat[xn[nn][3]]=1;
                                                EdgeStatus[xn[nn][3]]=0.0;
                                                ForceStat[xn[nn][7]]=1;
                                                EdgeStatus[xn[nn][7]]=0.0;
                                                ForceStat[xn[nn][9]]=1;
                                                EdgeStatus[xn[nn][9]]=0.0;
                                                ForceStat[xn[nn][12]]=1;
                                                EdgeStatus[xn[nn][12]]=0.0;
                                                ForceStat[xn[nn][16]]=1;
                                                EdgeStatus[xn[nn][16]]=0.0;
                                               }
                                              }
              /* bottom plate */
                           for(k=z1;k<z2;k++){
                             for(i=x1;i<x2;i++){
                                                nn=y2*x2*k+i;
                                                ForceStat[xn[nn][0]]=1;
                                                EdgeStatus[xn[nn][0]]=0.0;
                                                ForceStat[xn[nn][5]]=1;
                                                EdgeStatus[xn[nn][5]]=0.0;
                                                ForceStat[xn[nn][6]]=1;
                                                EdgeStatus[xn[nn][6]]=0.0;
                                                ForceStat[xn[nn][7]]=1;
                                                EdgeStatus[xn[nn][7]]=0.0;
                                                ForceStat[xn[nn][13]]=1;
                                                EdgeStatus[xn[nn][13]]=0.0;
                                               }
                                              }
              /* top plate */
                           for(k=z1;k<z2;k++){
                             for(i=x1;i<x2;i++){
                                                nn=y2*x2*k+x2*(y2-1)+i;
                                                ForceStat[xn[nn][4]]=1;
                                                EdgeStatus[xn[nn][4]]=0.0;
                                                ForceStat[xn[nn][10]]=1;
                                                EdgeStatus[xn[nn][10]]=0.0;
                                                ForceStat[xn[nn][11]]=1;
                                                EdgeStatus[xn[nn][11]]=0.0;
                                                ForceStat[xn[nn][12]]=1;
                                                EdgeStatus[xn[nn][12]]=0.0;
                                                ForceStat[xn[nn][17]]=1;
                                                EdgeStatus[xn[nn][17]]=0.0;
                                               }
                                              }

              /* front plate */
                           for(j=y1;j<y2;j++){
                             for(i=x1;i<x2;i++){
                                                nn=j*x2+i;
                                                ForceStat[xn[nn][0]]=1;
                                                EdgeStatus[xn[nn][0]]=0.0;
                                                ForceStat[xn[nn][1]]=1;
                                                EdgeStatus[xn[nn][1]]=0.0;
                                                ForceStat[xn[nn][2]]=1;
                                                EdgeStatus[xn[nn][2]]=0.0;
                                                ForceStat[xn[nn][3]]=1;
                                                EdgeStatus[xn[nn][3]]=0.0;
                                                ForceStat[xn[nn][4]]=1;
                                                EdgeStatus[xn[nn][4]]=0.0;
                                               }
                                              }
              /* back plate */
                           for(j=y1;j<y2;j++){
                             for(i=x1;i<x2;i++){
                                                nn=(z2-1)*y2*x2+j*x2+i;
                                                ForceStat[xn[nn][13]]=1;
                                                EdgeStatus[xn[nn][13]]=0.0;
                                                ForceStat[xn[nn][14]]=1;
                                                EdgeStatus[xn[nn][14]]=0.0;
                                                ForceStat[xn[nn][15]]=1;
                                                EdgeStatus[xn[nn][15]]=0.0;
                                                ForceStat[xn[nn][16]]=1;
                                                EdgeStatus[xn[nn][16]]=0.0;
                                                ForceStat[xn[nn][17]]=1;
                                                EdgeStatus[xn[nn][17]]=0.0;
                                               }
                                              }
             }
             if(!strcmp(type, "conductor"))
             {
                   /* fprintf(stderr, "conductor \n"); */
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;

              /* three-dimensional conductor */
                   for(i=x1;i<x2;i++){
                     for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                                           nn=k*XdiM*YdiM+j*XdiM+i;
                                           Epsilon[nn]=0;
                                           /* Sigma[nn]=1.0E+09;*/
                                           for(tmp=0;tmp<18;tmp++){
                                                                  ForceStat[xn[nn][tmp]]=1;
                                                                  EdgeStatus[xn[nn][tmp]]=0.0;
                                                                  }
                                          }
                                       }
                                     }
               /* two-dimensional conductor */
                  if(x1==x2){     /*  y-z plane */
                           for(j=y1;j<y2;j++){
                             for(k=z1;k<z2;k++){
                                                if(x1!=XdiM)   {
                                                nn=k*XdiM*YdiM+j*XdiM+x1;
                                                ForceStat[xn[nn][1]]=1;
                                                EdgeStatus[xn[nn][1]]=0.0;
                                                ForceStat[xn[nn][5]]=1;
                                                EdgeStatus[xn[nn][5]]=0.0;
                                                ForceStat[xn[nn][8]]=1;
                                                EdgeStatus[xn[nn][8]]=0.0;
                                                ForceStat[xn[nn][10]]=1;
                                                EdgeStatus[xn[nn][10]]=0.0;
                                                ForceStat[xn[nn][14]]=1;
                                                EdgeStatus[xn[nn][14]]=0.0;
                                                            }
                                                if(x1!=0)   {
                                                nn=k*XdiM*YdiM+j*XdiM+x1-1;
                                                ForceStat[xn[nn][3]]=1;
                                                EdgeStatus[xn[nn][3]]=0.0;
                                                ForceStat[xn[nn][7]]=1;
                                                EdgeStatus[xn[nn][7]]=0.0;
                                                ForceStat[xn[nn][9]]=1;
                                                EdgeStatus[xn[nn][9]]=0.0;
                                                ForceStat[xn[nn][12]]=1;
                                                EdgeStatus[xn[nn][12]]=0.0;
                                                ForceStat[xn[nn][16]]=1;
                                                EdgeStatus[xn[nn][16]]=0.0;
                                                           }
                                               }
                                              }
                           }
                 if(y1==y2){     /*  x-z plane */
                           for(i=x1;i<x2;i++){
                             for(k=z1;k<z2;k++){
                                                if(y1!=YdiM)   {
                                                nn=k*XdiM*YdiM+y1*XdiM+i;
                                                ForceStat[xn[nn][0]]=1;
                                                EdgeStatus[xn[nn][0]]=0.0;
                                                ForceStat[xn[nn][5]]=1;
                                                EdgeStatus[xn[nn][5]]=0.0;
                                                ForceStat[xn[nn][6]]=1;
                                                EdgeStatus[xn[nn][6]]=0.0;
                                                ForceStat[xn[nn][7]]=1;
                                                EdgeStatus[xn[nn][7]]=0.0;
                                                ForceStat[xn[nn][13]]=1;
                                                EdgeStatus[xn[nn][13]]=0.0;
                                                            }
                                                if(y1!=0)   {
                                                nn=k*XdiM*YdiM+(y1-1)*XdiM+i;
                                                ForceStat[xn[nn][4]]=1;
                                                EdgeStatus[xn[nn][4]]=0.0;
                                                ForceStat[xn[nn][10]]=1;
                                                EdgeStatus[xn[nn][10]]=0.0;
                                                ForceStat[xn[nn][11]]=1;
                                                EdgeStatus[xn[nn][11]]=0.0;
                                                ForceStat[xn[nn][12]]=1;
                                                EdgeStatus[xn[nn][12]]=0.0;
                                                ForceStat[xn[nn][17]]=1;
                                                EdgeStatus[xn[nn][17]]=0.0;
                                                           }
                                               }
                                              }
                           }
                 if(z1==z2){     /*   x-y plane */
                           for(i=x1;i<x2;i++){
                             for(j=y1;j<y2;j++){
                                                if(z1!=ZdiM)   {
                                                nn=z1*XdiM*YdiM+j*XdiM+i;
                                                ForceStat[xn[nn][0]]=1;
                                                EdgeStatus[xn[nn][0]]=0.0;
                                                ForceStat[xn[nn][1]]=1;
                                                EdgeStatus[xn[nn][1]]=0.0;
                                                ForceStat[xn[nn][2]]=1;
                                                EdgeStatus[xn[nn][2]]=0.0;
                                                ForceStat[xn[nn][3]]=1;
                                                EdgeStatus[xn[nn][3]]=0.0;
                                                ForceStat[xn[nn][4]]=1;
                                                EdgeStatus[xn[nn][4]]=0.0;
                                                            }
                                                if(z1!=0)   {
                                                nn=(z1-1)*XdiM*YdiM+j*XdiM+i;
                                                ForceStat[xn[nn][13]]=1;
                                                EdgeStatus[xn[nn][13]]=0.0;
                                                ForceStat[xn[nn][14]]=1;
                                                EdgeStatus[xn[nn][14]]=0.0;
                                                ForceStat[xn[nn][15]]=1;
                                                EdgeStatus[xn[nn][15]]=0.0;
                                                ForceStat[xn[nn][16]]=1;
                                                EdgeStatus[xn[nn][16]]=0.0;
                                                ForceStat[xn[nn][17]]=1;
                                                EdgeStatus[xn[nn][17]]=0.0;
                                                           }
                                               }
                                              }
                           }
 
               /* one-dimensional conductor */
 
                  if ((x1==x2)&&(y1==y2)){
                                          for(k=z1;k<z2;k++){
                                                            if((x1!=XdiM) & (y1!=YdiM)) {
                                                             nn=k*XdiM*YdiM+y1*XdiM+x1;
                                                             ForceStat[xn[nn][5]]=1;
                                                             EdgeStatus[xn[nn][5]]=0.0;
                                                                                  }
                                                            if((x1!=0) & (y1!=YdiM)) {
                                                             nn=k*XdiM*YdiM+y1*XdiM+x1-1;
                                                             ForceStat[xn[nn][7]]=1;
                                                             EdgeStatus[xn[nn][7]]=0.0;
                                                                                  }
                                                            if((x1!=XdiM) & (y1!=0)) {
                                                             nn=k*XdiM*YdiM+(y1-1)*XdiM+x1;
                                                             ForceStat[xn[nn][10]]=1;
                                                             EdgeStatus[xn[nn][10]]=0.0;
                                                                                  }
                                                            if((x1!=0) & (y1!=0)) {
                                                             nn=k*XdiM*YdiM+(y1-1)*XdiM+x1-1;
                                                             ForceStat[xn[nn][12]]=1;
                                                             EdgeStatus[xn[nn][12]]=0.0;
                                                                                  }
                                                            }
                                         }
                  if ((x1==x2)&&(z1==z2)){
                                          for(j=y1;j<y2;j++){
                                                            if((x1!=XdiM) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+j*XdiM+x1;
                                                             ForceStat[xn[nn][1]]=1;
                                                             EdgeStatus[xn[nn][1]]=0.0;
                                                                                  }
                                                            if((x1!=0) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+j*XdiM+x1-1;
                                                             ForceStat[xn[nn][3]]=1;
                                                             EdgeStatus[xn[nn][3]]=0.0;
                                                                                  }
                                                            if((x1!=XdiM) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+j*XdiM+x1;
                                                             ForceStat[xn[nn][14]]=1;
                                                             EdgeStatus[xn[nn][14]]=0.0;
                                                                                  }
                                                            if((x1!=0) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+j*XdiM+x1-1;
                                                             ForceStat[xn[nn][16]]=1;
                                                             EdgeStatus[xn[nn][16]]=0.0;
                                                                                  }
                                                            }
                                         }
                  if ((y1==y2)&&(z1==z2)){
                                          for(i=x1;i<x2;i++){
                                                            if((y1!=YdiM) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+y1*XdiM+i;
                                                             ForceStat[xn[nn][0]]=1;
                                                             EdgeStatus[xn[nn][0]]=0.0;
                                                                                  }
                                                            if((y1!=0) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+(y1-1)*XdiM+i;
                                                             ForceStat[xn[nn][4]]=1;
                                                             EdgeStatus[xn[nn][4]]=0.0;
                                                                                  }
                                                            if((y1!=YdiM) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+y1*XdiM+i;
                                                             ForceStat[xn[nn][13]]=1;
                                                             EdgeStatus[xn[nn][13]]=0.0;
                                                                                  }
                                                            if((y1!=0) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+(y1-1)*XdiM+i;
                                                             ForceStat[xn[nn][17]]=1;
                                                             EdgeStatus[xn[nn][17]]=0.0;
                                                                                  }
                                                            }
                                         }
 
                   continue;
             }

             if(!strcmp(type, "resistor"))
             {
                   /* fprintf(stderr, "resistor \n"); */
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;

               /* The resistor is expected to be one-dimensional */
 
                  if ((x1==x2)&&(y1==y2)){
                                          for(k=z1;k<z2;k++){
							    length=DivisorZ[k+1];
                                                            if((x1!=XdiM) & (y1!=YdiM)) {
                                                             nn=k*XdiM*YdiM+y1*XdiM+x1;
                                                             Resist[xn[nn][5]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((x1!=0) & (y1!=YdiM)) {
                                                             nn=k*XdiM*YdiM+y1*XdiM+x1-1;
                                                             Resist[xn[nn][7]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((x1!=XdiM) & (y1!=0)) {
                                                             nn=k*XdiM*YdiM+(y1-1)*XdiM+x1;
                                                             Resist[xn[nn][10]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((x1!=0) & (y1!=0)) {
                                                             nn=k*XdiM*YdiM+(y1-1)*XdiM+x1-1;
                                                             Resist[xn[nn][12]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            }
                                         }
                  if ((x1==x2)&&(z1==z2)){
                                          for(j=y1;j<y2;j++){
							    length=DivisorY[j+1];
                                                            if((x1!=XdiM) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+j*XdiM+x1;
                                                             Resist[xn[nn][1]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((x1!=0) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+j*XdiM+x1-1;
                                                             Resist[xn[nn][3]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((x1!=XdiM) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+j*XdiM+x1;
                                                             Resist[xn[nn][14]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((x1!=0) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+j*XdiM+x1-1;
                                                             Resist[xn[nn][16]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            }
                                         }
                  if ((y1==y2)&&(z1==z2)){
                                          for(i=x1;i<x2;i++){
							    length=DivisorX[i+1];
                                                            if((y1!=YdiM) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+y1*XdiM+i;
                                                             Resist[xn[nn][0]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((y1!=0) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+(y1-1)*XdiM+i;
                                                             Resist[xn[nn][4]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((y1!=YdiM) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+y1*XdiM+i;
                                                             Resist[xn[nn][13]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            if((y1!=0) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+(y1-1)*XdiM+i;
                                                             Resist[xn[nn][17]]=COMplex_Cmplx(0.0,WaveNumber*377.0*length*length/atof(att1));
                                                                                  }
                                                            }
                                         }
 
		   for (l=0;l<TotEdgeNum;l++)
			{
			if (COMplex_Abs(Resist[l])>0)	printf("Resist[%d]=%f %f",l,Resist[l].x,Resist[l].y);
			}
                   continue;
             }

             if(!strcmp(type, "coax"))
             {
                   /* fprintf(stderr, "coax \n"); */
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;

               /* The coax edges are defined */
 
                  if ((x1==x2)&&(y1==y2)){
                                          for(k=z1;k<z2;k++){
							    length=DivisorZ[k+1];
                                                            if((x1!=XdiM) & (y1!=YdiM)) {
                                                             nn=k*XdiM*YdiM+y1*XdiM+x1;
	                                                      ForceStat[xn[nn][5]]=0;
	                                                      EdgeStatus[xn[nn][5]]=-1;
                                                             Resist[xn[nn][5]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][5]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((x1!=0) & (y1!=YdiM)) {
                                                             nn=k*XdiM*YdiM+y1*XdiM+x1-1;
	                                                      ForceStat[xn[nn][7]]=0;
	                                                      EdgeStatus[xn[nn][7]]=-1;
                                                             Resist[xn[nn][7]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][7]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((x1!=XdiM) & (y1!=0)) {
                                                             nn=k*XdiM*YdiM+(y1-1)*XdiM+x1;
	                                                      ForceStat[xn[nn][10]]=0;
	                                                      EdgeStatus[xn[nn][10]]=-1;
                                                             Resist[xn[nn][10]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][10]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((x1!=0) & (y1!=0)) {
                                                             nn=k*XdiM*YdiM+(y1-1)*XdiM+x1-1;
	                                                      ForceStat[xn[nn][12]]=0;
	                                                      EdgeStatus[xn[nn][12]]=-1;
                                                             Resist[xn[nn][12]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][12]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            }
                                         }
                  if ((x1==x2)&&(z1==z2)){
printf("Inside this loop");
                                          for(j=y1;j<y2;j++){
							    length=DivisorY[j+1];
                                                            if((x1!=XdiM) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+j*XdiM+x1;
	                                                      ForceStat[xn[nn][1]]=0;
	                                                      EdgeStatus[xn[nn][1]]=-1;
                                                             Resist[xn[nn][1]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][1]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((x1!=0) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+j*XdiM+x1-1;
	                                                      ForceStat[xn[nn][3]]=0;
	                                                      EdgeStatus[xn[nn][3]]=-1;
                                                             Resist[xn[nn][3]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][3]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((x1!=XdiM) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+j*XdiM+x1;
	                                                      ForceStat[xn[nn][14]]=0;
	                                                      EdgeStatus[xn[nn][14]]=-1;
                                                             Resist[xn[nn][14]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][14]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((x1!=0) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+j*XdiM+x1-1;
	                                                      ForceStat[xn[nn][16]]=0;
	                                                      EdgeStatus[xn[nn][16]]=-1;
                                                             Resist[xn[nn][16]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][16]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            }
                                         }
                  if ((y1==y2)&&(z1==z2)){
printf("Inside this loop");
                                          for(i=x1;i<x2;i++){
							    length=DivisorX[i+1];
                                                            if((y1!=YdiM) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+y1*XdiM+i;
	                                                      ForceStat[xn[nn][0]]=0;
	                                                      EdgeStatus[xn[nn][0]]=-1;
                                                             Resist[xn[nn][0]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][0]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((y1!=0) & (z1!=ZdiM)) {
                                                             nn=z1*XdiM*YdiM+(y1-1)*XdiM+i;
	                                                      ForceStat[xn[nn][4]]=0;
	                                                      EdgeStatus[xn[nn][4]]=-1;
                                                             Resist[xn[nn][4]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][4]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((y1!=YdiM) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+y1*XdiM+i;
	                                                      ForceStat[xn[nn][13]]=0;
	                                                      EdgeStatus[xn[nn][13]]=-1;
                                                             Resist[xn[nn][13]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][13]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            if((y1!=0) & (z1!=0)) {
                                                             nn=(z1-1)*XdiM*YdiM+(y1-1)*XdiM+i;
	                                                      ForceStat[xn[nn][17]]=0;
	                                                      EdgeStatus[xn[nn][17]]=-1;
                                                             Resist[xn[nn][17]]=COMplex_Cmplx(atof(att2),0.0);
                                                             Isource[xn[nn][17]]=COMplex_Cmplx(atof(att3),0.0);
                                                                                  }
                                                            }
                                         }
 
		   for (l=0;l<TotEdgeNum;l++)
			{
			if (COMplex_Abs(Isource[l])>0)	printf("Isource[%d]=%f %f\n",l,Isource[l].x,Isource[l].y);
			if (COMplex_Abs(Resist[l])>0)	printf("Resist[%d]=%f %f\n",l,Resist[l].x,Resist[l].y);
			}
                   continue;
             }

             if(!strcmp(type, "aperture"))
             {
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;
 
                  if(x1==x2){                       /*  y-z plane */
                           for(j=y1;j<y2;j++){
                             for(k=z1;k<z2;k++){
                                                if(x1!=XdiM)   {
                                                nn=k*XdiM*YdiM+j*XdiM+x1;
                                                      if(k!=z1)     {
                                                      ForceStat[xn[nn][1]]=0;
                                                      EdgeStatus[xn[nn][1]]=-1;
                                                                    }
                                                      if(j!=y1)     {
                                                      ForceStat[xn[nn][5]]=0;
                                                      EdgeStatus[xn[nn][5]]=-1;
                                                                    }
                                                      ForceStat[xn[nn][8]]=0;
                                                      EdgeStatus[xn[nn][8]]=-1;
                                                      if(j!=(y2-1)) {
                                                      ForceStat[xn[nn][10]]=0;
                                                      EdgeStatus[xn[nn][10]]=-1;
                                                                    }
                                                      if(k!=(z2-1)) {
                                                      ForceStat[xn[nn][14]]=0;
                                                      EdgeStatus[xn[nn][14]]=-1;
                                                                    }
                                                             }
                                                if(x1!=0)   {
                                                nn=k*XdiM*YdiM+j*XdiM+x1-1;
                                                      if(k!=z1)     {
                                                      ForceStat[xn[nn][3]]=0;
                                                      EdgeStatus[xn[nn][3]]=-1;
                                                                    }
                                                      if(j!=y1)     {
                                                      ForceStat[xn[nn][7]]=0;
                                                      EdgeStatus[xn[nn][7]]=-1;
                                                                    }
                                                      ForceStat[xn[nn][9]]=0;
                                                      EdgeStatus[xn[nn][9]]=-1;
                                                      if(j!=(y2-1)) {
                                                      ForceStat[xn[nn][12]]=0;
                                                      EdgeStatus[xn[nn][12]]=-1;
                                                                    }
                                                      if(k!=(z2-1)) {
                                                      ForceStat[xn[nn][16]]=0;
                                                      EdgeStatus[xn[nn][16]]=-1;
                                                                    }
                                                           }
                                               }
                                              }
                           }
                 if(y1==y2){                             /*  x-z plane */
                           for(i=x1;i<x2;i++){
                             for(k=z1;k<z2;k++){
                                                if(y1!=YdiM)   {
                                                nn=k*XdiM*YdiM+y1*XdiM+i;
                                                      if(k!=z1)     {
                                                      ForceStat[xn[nn][0]]=0;
                                                      EdgeStatus[xn[nn][0]]=-1;
                                                                    }
                                                      if(i!=x1)     {
                                                      ForceStat[xn[nn][5]]=0;
                                                      EdgeStatus[xn[nn][5]]=-1;
                                                                    }
                                                      ForceStat[xn[nn][6]]=0;
                                                      EdgeStatus[xn[nn][6]]=-1;
                                                      if(i!=(x2-1)) {
                                                      ForceStat[xn[nn][7]]=0;
                                                      EdgeStatus[xn[nn][7]]=-1;
                                                                    }
                                                      if(k!=(z2-1)) {
                                                      ForceStat[xn[nn][13]]=0;
                                                      EdgeStatus[xn[nn][13]]=-1;
                                                                    }
                                                            }
                                                if(y1!=0)   {
                                                nn=k*XdiM*YdiM+(y1-1)*XdiM+i;
                                                      if(k!=z1)     {
                                                      ForceStat[xn[nn][4]]=0;
                                                      EdgeStatus[xn[nn][4]]=-1;
                                                                    }
                                                      if(i!=x1)     {
                                                      ForceStat[xn[nn][10]]=0;
                                                      EdgeStatus[xn[nn][10]]=-1;
                                                                    }
                                                      ForceStat[xn[nn][11]]=0;
                                                      EdgeStatus[xn[nn][11]]=-1;
                                                      if(i!=(x2-1)) {
                                                      ForceStat[xn[nn][12]]=0;
                                                      EdgeStatus[xn[nn][12]]=-1;
                                                                    }
                                                      if(k!=(z2-1)) {
                                                      ForceStat[xn[nn][17]]=0;
                                                      EdgeStatus[xn[nn][17]]=-1;
                                                                    }
                                                           }
                                               }
                                              }
                           }
                 if(z1==z2){                           /*   x-y plane */
                           for(i=x1;i<x2;i++){
                             for(j=y1;j<y2;j++){
                                                if(z1!=ZdiM)   {
                                                nn=z1*XdiM*YdiM+j*XdiM+i;
                                                     if(j!=y1)     {
                                                      ForceStat[xn[nn][0]]=0;
                                                      EdgeStatus[xn[nn][0]]=-1;
                                                                    }
                                                      if(i!=x1)     {
                                                      ForceStat[xn[nn][1]]=0;
                                                      EdgeStatus[xn[nn][1]]=-1;
                                                                    }
                                                      ForceStat[xn[nn][2]]=0;
                                                      EdgeStatus[xn[nn][2]]=-1;
                                                      if(i!=(x2-1)) {
                                                      ForceStat[xn[nn][3]]=0;
                                                      EdgeStatus[xn[nn][3]]=-1;
                                                                    }
                                                      if(j!=(y2-1)) {
                                                      ForceStat[xn[nn][4]]=0;
                                                      EdgeStatus[xn[nn][4]]=-1;
                                                                    }
                                                            }
                                                if(z1!=0)   {
                                                nn=(z1-1)*XdiM*YdiM+j*XdiM+i;
                                                     if(j!=y1)     {
                                                      ForceStat[xn[nn][13]]=0;
                                                      EdgeStatus[xn[nn][13]]=-1;
                                                                    }
                                                      if(i!=x1)     {
                                                      ForceStat[xn[nn][14]]=0;
                                                      EdgeStatus[xn[nn][14]]=-1;
                                                                    }
                                                      ForceStat[xn[nn][15]]=0;
                                                      EdgeStatus[xn[nn][15]]=-1;
                                                      if(i!=(x2-1)) {
                                                      ForceStat[xn[nn][16]]=0;
                                                      EdgeStatus[xn[nn][16]]=-1;
                                                                    }
                                                      if(j!=(y2-1)) {
                                                      ForceStat[xn[nn][17]]=0;
                                                      EdgeStatus[xn[nn][17]]=-1;
                                                                    }
                                                           }
                                               }
                                              }
                           }
 
                   continue;
             }
             if(!strcmp(type, "seam"))
             {
                   /* fprintf(stderr, "seam \n"); */
                   /*  need to define seam here  */
                   continue;
             }
             if(!strcmp(type, "esource"))
             {
                   /* fprintf(stderr, "esource \n"); */
                   if(OperateFreq != atof(att1)) { fprintf(stderr, "Error: all sources must have the same frequency\n");
                                            exit(1);
                                          }
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;

                  /* three dimensional source blocks */       
                   if (att2[0] == 'y')
                     {
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                       nn=k*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][1]] = 1;
                       EdgeStatus[xn[nn][1]] = atof(att3);
                       EdgeStatus1[xn[nn][1]] = atof(att4);
                       ForceStat[xn[nn][3]] = 1;
                       EdgeStatus[xn[nn][3]] = atof(att3);
                       EdgeStatus1[xn[nn][3]] = atof(att4);
                       ForceStat[xn[nn][14]] = 1;
                       EdgeStatus[xn[nn][14]] = atof(att3);
                       EdgeStatus1[xn[nn][14]] = atof(att4);
                       ForceStat[xn[nn][16]] = 1;
                       EdgeStatus[xn[nn][16]] = atof(att3);
                       EdgeStatus1[xn[nn][16]] = atof(att4);
                                           }
                                         }
                                         }
                     }

                   if (att2[0] == 'z')
                     { 
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                       nn=k*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][5]] = 1;
                       EdgeStatus[xn[nn][5]] = atof(att3);
                       EdgeStatus1[xn[nn][5]] = atof(att4);
                       ForceStat[xn[nn][7]] = 1;
                       EdgeStatus[xn[nn][7]] = atof(att3);
                       EdgeStatus1[xn[nn][7]] = atof(att4);
                       ForceStat[xn[nn][10]] = 1;
                       EdgeStatus[xn[nn][10]] = atof(att3);
                       EdgeStatus1[xn[nn][10]] = atof(att4);
                       ForceStat[xn[nn][12]] = 1;
                       EdgeStatus[xn[nn][12]] = atof(att3);
                       EdgeStatus1[xn[nn][12]] = atof(att4);
                                          }
                                         }
                                        }
                     }

                  if (att2[0] == 'x')
                     {
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                       nn=k*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][0]] = 1;
                       EdgeStatus[xn[nn][0]] = atof(att3);
                       EdgeStatus1[xn[nn][0]] = atof(att4);
                       ForceStat[xn[nn][4]] = 1;
                       EdgeStatus[xn[nn][4]] = atof(att3);
                       EdgeStatus1[xn[nn][4]] = atof(att4);
                       ForceStat[xn[nn][13]] = 1;
                       EdgeStatus[xn[nn][13]] = atof(att3);
                       EdgeStatus1[xn[nn][13]] = atof(att4);
                       ForceStat[xn[nn][17]] = 1;
                       EdgeStatus[xn[nn][17]] = atof(att3);
                       EdgeStatus1[xn[nn][17]] = atof(att4);
                                          }
                                         }
                                        }
                     }

                  /* two dimensional source blocks */       
                   if ((att2[0] == 'y')&&(z1==z2))
                     {
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                       if(z1!=ZdiM){
                       nn=z1*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][1]] = 1;
                       EdgeStatus[xn[nn][1]] = atof(att3);
                       EdgeStatus1[xn[nn][1]] = atof(att4);
                       ForceStat[xn[nn][3]] = 1;
                       EdgeStatus[xn[nn][3]] = atof(att3);
                       EdgeStatus1[xn[nn][3]] = atof(att4); 
                       ForceStat[xn[nn][2]] = 1;
                       EdgeStatus[xn[nn][2]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][2]] = atof(att4); 
                                   }
                       else       { 
                       nn=(z1-1)*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][14]] = 1; 
                       EdgeStatus[xn[nn][14]] = atof(att3);
                       EdgeStatus1[xn[nn][14]] = atof(att4);
                       ForceStat[xn[nn][16]] = 1;
                       EdgeStatus[xn[nn][16]] = atof(att3);
                       EdgeStatus1[xn[nn][16]] = atof(att4);
                       ForceStat[xn[nn][15]] = 1;
                       EdgeStatus[xn[nn][15]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][15]] = atof(att4); 
                                   }
                                         }
                                         }
                     }
                   if ((att2[0] == 'y')&&(x1==x2))
                     {
                      for(k=z1;k<z2;k++){ 
                       for(j=y1;j<y2;j++){
                       if(x1!=XdiM){
                       nn=k*XdiM*YdiM+j*XdiM+x1;
                       ForceStat[xn[nn][1]] = 1;
                       EdgeStatus[xn[nn][1]] = atof(att3); 
                       EdgeStatus1[xn[nn][1]] = atof(att4);
                       ForceStat[xn[nn][14]] = 1;
                       EdgeStatus[xn[nn][14]] = atof(att3);
                       EdgeStatus1[xn[nn][14]] = atof(att4);
                       ForceStat[xn[nn][8]] = 1;
                       EdgeStatus[xn[nn][8]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][8]] = atof(att4); 
                                   }
                       else        {
                       nn=k*XdiM*YdiM+j*XdiM+x1-1;
                       ForceStat[xn[nn][3]] = 1;
                       EdgeStatus[xn[nn][3]] = atof(att3);  
                       EdgeStatus1[xn[nn][3]] = atof(att4);
                       ForceStat[xn[nn][16]] = 1;
                       EdgeStatus[xn[nn][16]] = atof(att3); 
                       EdgeStatus1[xn[nn][16]] = atof(att4);
                       ForceStat[xn[nn][9]] = 1;
                       EdgeStatus[xn[nn][9]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][9]] = atof(att4); 
                                   }
                                         }
                                         }
                     }
                   if ((att2[0] == 'z')&&(y1==y2))
                     {
                      for(i=x1;i<x2;i++){ 
                        for(k=z1;k<z2;k++){
                       if(y1!=YdiM){
                       nn=k*XdiM*YdiM+y1*XdiM+i;
                       ForceStat[xn[nn][5]] = 1;
                       EdgeStatus[xn[nn][5]] = atof(att3); 
                       EdgeStatus1[xn[nn][5]] = atof(att4);
                       ForceStat[xn[nn][7]] = 1;
                       EdgeStatus[xn[nn][7]] = atof(att3);  
                       EdgeStatus1[xn[nn][7]] = atof(att4);
                       ForceStat[xn[nn][6]] = 1;
                       EdgeStatus[xn[nn][6]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][6]] = atof(att4); 
                                   }
                       else        {
                       nn=k*XdiM*YdiM+(y1-1)*XdiM+i;
                       ForceStat[xn[nn][10]] = 1;
                       EdgeStatus[xn[nn][10]] = atof(att3);
                       EdgeStatus1[xn[nn][10]] = atof(att4);
                       ForceStat[xn[nn][12]] = 1;
                       EdgeStatus[xn[nn][12]] = atof(att3); 
                       EdgeStatus1[xn[nn][12]] = atof(att4);
                       ForceStat[xn[nn][11]] = 1;
                       EdgeStatus[xn[nn][11]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][11]] = atof(att4); 
                                   }
                                          }
                                        }
                     }
                   if ((att2[0] == 'z')&&(x1==x2))
                     {
                      for(j=y1;j<y2;j++){  
                        for(k=z1;k<z2;k++){
                       if(x1!=XdiM){
                       nn=k*XdiM*YdiM+j*XdiM+x1;
                       ForceStat[xn[nn][5]] = 1;
                       EdgeStatus[xn[nn][5]] = atof(att3); 
                       EdgeStatus1[xn[nn][5]] = atof(att4);
                       ForceStat[xn[nn][10]] = 1;
                       EdgeStatus[xn[nn][10]] = atof(att3);
                       EdgeStatus1[xn[nn][10]] = atof(att4);
                       ForceStat[xn[nn][8]] = 1;
                       EdgeStatus[xn[nn][8]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][8]] = atof(att4); 
                                    }
                       else         {
                       nn=k*XdiM*YdiM+j*XdiM+x1-1;
                       ForceStat[xn[nn][7]] = 1;
                       EdgeStatus[xn[nn][7]] = atof(att3);  
                       EdgeStatus1[xn[nn][7]] = atof(att4);
                       ForceStat[xn[nn][12]] = 1;
                       EdgeStatus[xn[nn][12]] = atof(att3); 
                       EdgeStatus1[xn[nn][12]] = atof(att4);
                       ForceStat[xn[nn][9]] = 1;
                       EdgeStatus[xn[nn][9]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][9]] = atof(att4); 
                                    }
                                          }
                                        }
                     }
                  if ((att2[0] == 'x')&&(y1==y2))
                     {
                      for(i=x1;i<x2;i++){  
                        for(k=z1;k<z2;k++){
                       if(y1!=YdiM){
                       nn=k*XdiM*YdiM+y1*XdiM+i;
                       ForceStat[xn[nn][0]] = 1;
                       EdgeStatus[xn[nn][0]] = atof(att3); 
                       EdgeStatus1[xn[nn][0]] = atof(att4);
                       ForceStat[xn[nn][13]] = 1;
                       EdgeStatus[xn[nn][13]] = atof(att3); 
                       EdgeStatus1[xn[nn][13]] = atof(att4);
                       ForceStat[xn[nn][6]] = 1;
                       EdgeStatus[xn[nn][6]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][6]] = atof(att4); 
                                    }
                       else         {
                       nn=k*XdiM*YdiM+(y1-1)*XdiM+i;
                       ForceStat[xn[nn][4]] = 1;
                       EdgeStatus[xn[nn][4]] = atof(att3); 
                       EdgeStatus1[xn[nn][4]] = atof(att4);
                       ForceStat[xn[nn][17]] = 1;
                       EdgeStatus[xn[nn][17]] = atof(att3); 
                       EdgeStatus1[xn[nn][17]] = atof(att4);
                       ForceStat[xn[nn][11]] = 1;
                       EdgeStatus[xn[nn][11]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][11]] = atof(att4); 
                                    }
                                          }
                                        }
                     }
                  if ((att2[0] == 'x')&&(z1==z2))
                     {
                      for(i=x1;i<x2;i++){
                        for(j=y1;j<y2;j++){
                       if(z1!=ZdiM){
                       nn=z1*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][0]] = 1;
                       EdgeStatus[xn[nn][0]] = atof(att3);
                       EdgeStatus1[xn[nn][0]] = atof(att4);
                       ForceStat[xn[nn][4]] = 1;
                       EdgeStatus[xn[nn][4]] = atof(att3);
                       EdgeStatus1[xn[nn][4]] = atof(att4);
                       ForceStat[xn[nn][2]] = 1;
                       EdgeStatus[xn[nn][2]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][2]] = atof(att4); 
                                    }
                       else         {
                       nn=(z1-1)*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][13]] = 1;
                       EdgeStatus[xn[nn][13]] = atof(att3);
                       EdgeStatus1[xn[nn][13]] = atof(att4);
                       ForceStat[xn[nn][17]] = 1;
                       EdgeStatus[xn[nn][17]] = atof(att3);
                       EdgeStatus1[xn[nn][17]] = atof(att4);
                       ForceStat[xn[nn][15]] = 1;
                       EdgeStatus[xn[nn][15]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][15]] = atof(att4); 
                                    }
                                          }
                                        }
 
                     }

                  /* one dimensional source blocks */       
                   if ((att2[0] == 'x')&&(z1==z2)&&(y1==y2))
                     {
                       for(i=x1;i<x2;i++){
                                                  if((y1!=YdiM) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+y1*XdiM+i;
                                                   ForceStat[xn[nn][0]]=1;
                                                   EdgeStatus[xn[nn][0]]=atof(att3);
                                                   EdgeStatus1[xn[nn][0]]=atof(att4);
                                                                              }
                                                  if((y1!=0) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+(y1-1)*XdiM+i;
                                                   ForceStat[xn[nn][4]]=1;
                                                   EdgeStatus[xn[nn][4]]=atof(att3);
                                                   EdgeStatus1[xn[nn][4]]=atof(att4);
                                                                            }
                                                  if((y1!=YdiM) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+y1*XdiM+i;
                                                   ForceStat[xn[nn][13]]=1;
                                                   EdgeStatus[xn[nn][13]]=atof(att3);
                                                   EdgeStatus1[xn[nn][13]]=atof(att4);
                                                                            }
                                                  if((y1!=0) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+(y1-1)*XdiM+i;
                                                   ForceStat[xn[nn][17]]=1;
                                                   EdgeStatus[xn[nn][17]]=atof(att3);
                                                   EdgeStatus1[xn[nn][17]]=atof(att4);
                                                                         }
                                                  
                                         }
                     }
                   if ((att2[0] == 'y')&&(z1==z2)&&(x1==x2))
                     {
                       for(j=y1;j<y2;j++){
                                                  if((x1!=XdiM) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+j*XdiM+x1;
                                                   ForceStat[xn[nn][1]]=1;
                                                   EdgeStatus[xn[nn][1]]=atof(att3);
                                                   EdgeStatus1[xn[nn][1]]=atof(att4);
                                                                               }
                                                  if((x1!=0) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+j*XdiM+x1-1;  
                                                   ForceStat[xn[nn][3]]=1;
                                                   EdgeStatus[xn[nn][3]]=atof(att3);
                                                   EdgeStatus1[xn[nn][3]]=atof(att4);
                                                                           } 
                                                  if((x1!=XdiM) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+j*XdiM+x1;
                                                   ForceStat[xn[nn][14]]=1;
                                                   EdgeStatus[xn[nn][14]]=atof(att3);
                                                   EdgeStatus1[xn[nn][14]]=atof(att4);
                                                                            }
                                                  if((x1!=0) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+j*XdiM+x1-1;
                                                   ForceStat[xn[nn][16]]=1;
                                                   EdgeStatus[xn[nn][16]]=atof(att3);
                                                   EdgeStatus1[xn[nn][16]]=atof(att4);
                                                                         }
                                                  
                                         }
                     }
                   if ((att2[0] == 'z')&&(y1==y2)&&(x1==x2))
                     {
                       for(k=z1;k<z2;k++){
                                                  if((x1!=XdiM) & (y1!=YdiM)) {
                                                   nn=k*XdiM*YdiM+y1*XdiM+x1;    
                                                   ForceStat[xn[nn][5]]=1;
                                                   EdgeStatus[xn[nn][5]]=atof(att3);
                                                   EdgeStatus1[xn[nn][5]]=atof(att4);
                                                                              }
                                                  if((x1!=0) & (y1!=YdiM)) {
                                                   nn=k*XdiM*YdiM+y1*XdiM+x1-1;  
                                                   ForceStat[xn[nn][7]]=1; 
                                                   EdgeStatus[xn[nn][7]]=atof(att3); 
                                                   EdgeStatus1[xn[nn][7]]=atof(att4); 
                                                                            }
                                                  if((x1!=XdiM) & (y1!=0)) {
                                                   nn=k*XdiM*YdiM+(y1-1)*XdiM+x1;  
                                                   ForceStat[xn[nn][10]]=1;
                                                   EdgeStatus[xn[nn][10]]=atof(att3);
                                                   EdgeStatus1[xn[nn][10]]=atof(att4);
                                                                            }
                                                  if((x1!=0) & (y1!=0)) {
                                                   nn=k*XdiM*YdiM+(y1-1)*XdiM+x1-1;
                                                   ForceStat[xn[nn][12]]=1;
                                                   EdgeStatus[xn[nn][12]]=atof(att3);
                                                   EdgeStatus1[xn[nn][12]]=atof(att4);
                                                                         }
                                                  
                                         }
                                                   
                     }

                 continue;
              }
             if(!strcmp(type, "isource"))
             {
                   /* fprintf(stderr, "isource \n"); */
                   if(OperateFreq != atof(att1)) { fprintf(stderr, "Error: all sources must have the same frequency\n");
                                            exit(1);
                                          }
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   x1=x1-Min_X;y1=y1-Min_Y;z1=z1-Min_Z;
                   x2=x2-Min_X;y2=y2-Min_Y;z2=z2-Min_Z;

                  /* three dimensional source blocks        
                   if (att2[0] == 'y')
                     {
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                       nn=k*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][1]] = 0;
                       EdgeStatus[xn[nn][1]] = -1;
                       Isource[xn[nn][1]] = COMplex_Cmplx(atof(att3)*cos(atof(att4),atof(att3)*sin(atof(att4)));
                       ForceStat[xn[nn][3]] = 0;
                       EdgeStatus[xn[nn][3]] = atof(att3);
                       EdgeStatus1[xn[nn][3]] = atof(att4);
                       ForceStat[xn[nn][14]] = 0;
                       EdgeStatus[xn[nn][14]] = atof(att3);
                       EdgeStatus1[xn[nn][14]] = atof(att4);
                       ForceStat[xn[nn][16]] = 0;
                       EdgeStatus[xn[nn][16]] = atof(att3);
                       EdgeStatus1[xn[nn][16]] = atof(att4);
                                           }
                                         }
                                         }
                     }

                   if (att2[0] == 'z')
                     { 
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                       nn=k*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][5]] = 0;
                       EdgeStatus[xn[nn][5]] = atof(att3);
                       EdgeStatus1[xn[nn][5]] = atof(att4);
                       ForceStat[xn[nn][7]] = 0;
                       EdgeStatus[xn[nn][7]] = atof(att3);
                       EdgeStatus1[xn[nn][7]] = atof(att4);
                       ForceStat[xn[nn][10]] = 0;
                       EdgeStatus[xn[nn][10]] = atof(att3);
                       EdgeStatus1[xn[nn][10]] = atof(att4);
                       ForceStat[xn[nn][12]] = 0;
                       EdgeStatus[xn[nn][12]] = atof(att3);
                       EdgeStatus1[xn[nn][12]] = atof(att4);
                                          }
                                         }
                                        }
                     }

                  if (att2[0] == 'x')
                     {
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                        for(k=z1;k<z2;k++){
                       nn=k*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][0]] = 0;
                       EdgeStatus[xn[nn][0]] = atof(att3);
                       EdgeStatus1[xn[nn][0]] = atof(att4);
                       ForceStat[xn[nn][4]] = 0;
                       EdgeStatus[xn[nn][4]] = atof(att3);
                       EdgeStatus1[xn[nn][4]] = atof(att4);
                       ForceStat[xn[nn][13]] = 0;
                       EdgeStatus[xn[nn][13]] = atof(att3);
                       EdgeStatus1[xn[nn][13]] = atof(att4);
                       ForceStat[xn[nn][17]] = 0;
                       EdgeStatus[xn[nn][17]] = atof(att3);
                       EdgeStatus1[xn[nn][17]] = atof(att4);
                                          }
                                         }
                                        }
                     }

                   two dimensional source blocks        
                   if ((att2[0] == 'y')&&(z1==z2))
                     {
                      for(i=x1;i<x2;i++){
                       for(j=y1;j<y2;j++){
                       if(z1!=ZdiM){
                       nn=z1*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][1]] = 0;
                       EdgeStatus[xn[nn][1]] = atof(att3);
                       EdgeStatus1[xn[nn][1]] = atof(att4);
                       ForceStat[xn[nn][3]] = 0;
                       EdgeStatus[xn[nn][3]] = atof(att3);
                       EdgeStatus1[xn[nn][3]] = atof(att4); 
                       ForceStat[xn[nn][2]] = 0;
                       EdgeStatus[xn[nn][2]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][2]] = atof(att4); 
                                   }
                       else       { 
                       nn=(z1-1)*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][14]] = 0; 
                       EdgeStatus[xn[nn][14]] = atof(att3);
                       EdgeStatus1[xn[nn][14]] = atof(att4);
                       ForceStat[xn[nn][16]] = 0;
                       EdgeStatus[xn[nn][16]] = atof(att3);
                       EdgeStatus1[xn[nn][16]] = atof(att4);
                       ForceStat[xn[nn][15]] = 0;
                       EdgeStatus[xn[nn][15]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][15]] = atof(att4); 
                                   }
                                         }
                                         }
                     }
                   if ((att2[0] == 'y')&&(x1==x2))
                     {
                      for(k=z1;k<z2;k++){ 
                       for(j=y1;j<y2;j++){
                       if(x1!=XdiM){
                       nn=k*XdiM*YdiM+j*XdiM+x1;
                       ForceStat[xn[nn][1]] = 0;
                       EdgeStatus[xn[nn][1]] = atof(att3); 
                       EdgeStatus1[xn[nn][1]] = atof(att4);
                       ForceStat[xn[nn][14]] = 0;
                       EdgeStatus[xn[nn][14]] = atof(att3);
                       EdgeStatus1[xn[nn][14]] = atof(att4);
                       ForceStat[xn[nn][8]] = 0;
                       EdgeStatus[xn[nn][8]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][8]] = atof(att4); 
                                   }
                       else        {
                       nn=k*XdiM*YdiM+j*XdiM+x1-1;
                       ForceStat[xn[nn][3]] = 0;
                       EdgeStatus[xn[nn][3]] = atof(att3);  
                       EdgeStatus1[xn[nn][3]] = atof(att4);
                       ForceStat[xn[nn][16]] = 0;
                       EdgeStatus[xn[nn][16]] = atof(att3); 
                       EdgeStatus1[xn[nn][16]] = atof(att4);
                       ForceStat[xn[nn][9]] = 0;
                       EdgeStatus[xn[nn][9]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][9]] = atof(att4); 
                                   }
                                         }
                                         }
                     }
                   if ((att2[0] == 'z')&&(y1==y2))
                     {
                      for(i=x1;i<x2;i++){ 
                        for(k=z1;k<z2;k++){
                       if(y1!=YdiM){
                       nn=k*XdiM*YdiM+y1*XdiM+i;
                       ForceStat[xn[nn][5]] = 0;
                       EdgeStatus[xn[nn][5]] = atof(att3); 
                       EdgeStatus1[xn[nn][5]] = atof(att4);
                       ForceStat[xn[nn][7]] = 0;
                       EdgeStatus[xn[nn][7]] = atof(att3);  
                       EdgeStatus1[xn[nn][7]] = atof(att4);
                       ForceStat[xn[nn][6]] = 0;
                       EdgeStatus[xn[nn][6]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][6]] = atof(att4); 
                                   }
                       else        {
                       nn=k*XdiM*YdiM+(y1-1)*XdiM+i;
                       ForceStat[xn[nn][10]] = 0;
                       EdgeStatus[xn[nn][10]] = atof(att3);
                       EdgeStatus1[xn[nn][10]] = atof(att4);
                       ForceStat[xn[nn][12]] = 0;
                       EdgeStatus[xn[nn][12]] = atof(att3); 
                       EdgeStatus1[xn[nn][12]] = atof(att4);
                       ForceStat[xn[nn][11]] = 0;
                       EdgeStatus[xn[nn][11]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][11]] = atof(att4); 
                                   }
                                          }
                                        }
                     }
                   if ((att2[0] == 'z')&&(x1==x2))
                     {
                      for(j=y1;j<y2;j++){  
                        for(k=z1;k<z2;k++){
                       if(x1!=XdiM){
                       nn=k*XdiM*YdiM+j*XdiM+x1;
                       ForceStat[xn[nn][5]] = 0;
                       EdgeStatus[xn[nn][5]] = atof(att3); 
                       EdgeStatus1[xn[nn][5]] = atof(att4);
                       ForceStat[xn[nn][10]] = 0;
                       EdgeStatus[xn[nn][10]] = atof(att3);
                       EdgeStatus1[xn[nn][10]] = atof(att4);
                       ForceStat[xn[nn][8]] = 0;
                       EdgeStatus[xn[nn][8]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][8]] = atof(att4); 
                                    }
                       else         {
                       nn=k*XdiM*YdiM+j*XdiM+x1-1;
                       ForceStat[xn[nn][7]] = 0;
                       EdgeStatus[xn[nn][7]] = atof(att3);  
                       EdgeStatus1[xn[nn][7]] = atof(att4);
                       ForceStat[xn[nn][12]] = 0;
                       EdgeStatus[xn[nn][12]] = atof(att3); 
                       EdgeStatus1[xn[nn][12]] = atof(att4);
                       ForceStat[xn[nn][9]] = 0;
                       EdgeStatus[xn[nn][9]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][9]] = atof(att4); 
                                    }
                                          }
                                        }
                     }
                  if ((att2[0] == 'x')&&(y1==y2))
                     {
                      for(i=x1;i<x2;i++){  
                        for(k=z1;k<z2;k++){
                       if(y1!=YdiM){
                       nn=k*XdiM*YdiM+y1*XdiM+i;
                       ForceStat[xn[nn][0]] = 0;
                       EdgeStatus[xn[nn][0]] = atof(att3); 
                       EdgeStatus1[xn[nn][0]] = atof(att4);
                       ForceStat[xn[nn][13]] = 0;
                       EdgeStatus[xn[nn][13]] = atof(att3); 
                       EdgeStatus1[xn[nn][13]] = atof(att4);
                       ForceStat[xn[nn][6]] = 0;
                       EdgeStatus[xn[nn][6]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][6]] = atof(att4); 
                                    }
                       else         {
                       nn=k*XdiM*YdiM+(y1-1)*XdiM+i;
                       ForceStat[xn[nn][4]] = 0;
                       EdgeStatus[xn[nn][4]] = atof(att3); 
                       EdgeStatus1[xn[nn][4]] = atof(att4);
                       ForceStat[xn[nn][17]] = 0;
                       EdgeStatus[xn[nn][17]] = atof(att3); 
                       EdgeStatus1[xn[nn][17]] = atof(att4);
                       ForceStat[xn[nn][11]] = 0;
                       EdgeStatus[xn[nn][11]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][11]] = atof(att4); 
                                    }
                                          }
                                        }
                     }
                  if ((att2[0] == 'x')&&(z1==z2))
                     {
                      for(i=x1;i<x2;i++){
                        for(j=y1;j<y2;j++){
                       if(z1!=ZdiM){
                       nn=z1*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][0]] = 0;
                       EdgeStatus[xn[nn][0]] = atof(att3);
                       EdgeStatus1[xn[nn][0]] = atof(att4);
                       ForceStat[xn[nn][4]] = 0;
                       EdgeStatus[xn[nn][4]] = atof(att3);
                       EdgeStatus1[xn[nn][4]] = atof(att4);
                       ForceStat[xn[nn][2]] = 0;
                       EdgeStatus[xn[nn][2]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][2]] = atof(att4); 
                                    }
                       else         {
                       nn=(z1-1)*XdiM*YdiM+j*XdiM+i;
                       ForceStat[xn[nn][13]] = 0;
                       EdgeStatus[xn[nn][13]] = atof(att3);
                       EdgeStatus1[xn[nn][13]] = atof(att4);
                       ForceStat[xn[nn][17]] = 0;
                       EdgeStatus[xn[nn][17]] = atof(att3);
                       EdgeStatus1[xn[nn][17]] = atof(att4);
                       ForceStat[xn[nn][15]] = 0;
                       EdgeStatus[xn[nn][15]] = 0.707 * atof(att3);
                       EdgeStatus1[xn[nn][15]] = atof(att4); 
                                    }
                                          }
                                        }
 
                     }

                  one dimensional source blocks */       
                   if ((att2[0] == 'x')&&(z1==z2)&&(y1==y2))
                     {
                       for(i=x1;i<x2;i++){
  					          length=DivisorX[i+1];
                                                  if((y1!=YdiM) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+y1*XdiM+i;
                                                   ForceStat[xn[nn][0]]=0;
						   Isource[xn[nn][0]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                              }
                                                  if((y1!=0) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+(y1-1)*XdiM+i;
                                                   ForceStat[xn[nn][4]]=0;
						   Isource[xn[nn][4]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                            }
                                                  if((y1!=YdiM) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+y1*XdiM+i;
                                                   ForceStat[xn[nn][13]]=0;
						   Isource[xn[nn][13]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                            }
                                                  if((y1!=0) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+(y1-1)*XdiM+i;
                                                   ForceStat[xn[nn][17]]=0;
						   Isource[xn[nn][17]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                         }
                                                  
                                         }
                     }
                   if ((att2[0] == 'y')&&(z1==z2)&&(x1==x2))
                     {
                       for(j=y1;j<y2;j++){
  					          length=DivisorY[j+1];
                                                  if((x1!=XdiM) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+j*XdiM+x1;
                                                   ForceStat[xn[nn][1]]=0;
						   Isource[xn[nn][1]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                               }
                                                  if((x1!=0) & (z1!=ZdiM)) {
                                                   nn=z1*XdiM*YdiM+j*XdiM+x1-1;  
                                                   ForceStat[xn[nn][3]]=0;
						   Isource[xn[nn][3]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                           } 
                                                  if((x1!=XdiM) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+j*XdiM+x1;
                                                   ForceStat[xn[nn][14]]=0;
						   Isource[xn[nn][14]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                            }
                                                  if((x1!=0) & (z1!=0)) {
                                                   nn=(z1-1)*XdiM*YdiM+j*XdiM+x1-1;
                                                   ForceStat[xn[nn][16]]=0;
						   Isource[xn[nn][16]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                         }
                                                  
                                         }
                     }
                   if ((att2[0] == 'z')&&(y1==y2)&&(x1==x2))
                     {
                       for(k=z1;k<z2;k++){
  					          length=DivisorZ[k+1];
                                                  if((x1!=XdiM) & (y1!=YdiM)) {
                                                   nn=k*XdiM*YdiM+y1*XdiM+x1;    
                                                   ForceStat[xn[nn][5]]=0;
						   Isource[xn[nn][5]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                              }
                                                  if((x1!=0) & (y1!=YdiM)) {
                                                   nn=k*XdiM*YdiM+y1*XdiM+x1-1;  
                                                   ForceStat[xn[nn][7]]=0; 
						   Isource[xn[nn][7]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                            }
                                                  if((x1!=XdiM) & (y1!=0)) {
                                                   nn=k*XdiM*YdiM+(y1-1)*XdiM+x1;  
                                                   ForceStat[xn[nn][10]]=0;
						   Isource[xn[nn][10]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                            }
                                                  if((x1!=0) & (y1!=0)) {
                                                   nn=k*XdiM*YdiM+(y1-1)*XdiM+x1-1;
                                                   ForceStat[xn[nn][12]]=0;
						   Isource[xn[nn][12]]=Real_Mul(WaveNumber*377.0*length,COMplex_Cmplx(atof(att3)*sin(atof(att4)),(-1)*atof(att3)*cos(atof(att4)))); 
                                                                         }
                                                  
                                         }
                                                   
                     }
		   for (l=0;l<TotEdgeNum;l++)
			{
			if (COMplex_Abs(Isource[l])!=0)	printf("Isource[%d]=%f %f",l,Isource[l].x,Isource[l].y);
			}

                 continue;
              }
              if(!strcmp(type, "celldim"))
              {continue;
              }
           }              /*end of if*/ 
        }                 /* end of while */
        rewind(InF);
        remove("source_nodes.tmp");
        CountEdgeType(); 

}

/***************************************************************************
*   FUNCTION        : absol()
****************************************************************************/

int    absol(int x)
{
	if (x<0) return(-x); else  return(x);
}

/***************************************************************************
*
*   FUNCTION    : swap(x,y)
*   DESCRIPTION : swap the values of x and y
*
****************************************************************************/

void swap(int* x,int* y)
{
	int temp; temp = *x; *x=*y; *y=temp;
}

/***************************************************************************
*   FUNCTION        : CountEdgeType()

Once the assignment of all the edges are done, this routine calculates how 
many inner, forced and boundary edges are there. Boundary edges are not used
in EMAP4 and they are treated as though they are electric conductors 
(if not already specified) and the
value of electric field is set to zero. But if you specify that an edge on 
the boundary is part of an aperture, then it treats it as an "inner" edge. 
This is not encouraged as it might give some erroneous results. 
****************************************************************************/

void CountEdgeType()
{

 float temp;

 TotInnerEdgeNum = 0;
 TotBoundEdgeNum = 0;
 TotForcdEdgeNum = 0;   
 
 for(IEdge=0;IEdge<TotEdgeNum;IEdge++)
 {
 if(ForceStat[IEdge] == 0)   {
    if(EdgeStatus[IEdge]==-1)       {TotInnerEdgeNum++;}
    else if(EdgeStatus[IEdge]==-2)  {TotBoundEdgeNum++;}
                         }
 else                    {           TotForcdEdgeNum++;}
  }

if(freq_step==0){
 InnerEdgeStat = INT_Vector(TotInnerEdgeNum);
 BoundEdgeStat = INT_Vector(TotBoundEdgeNum);
 ForcdEdgeStat = INT_Vector(TotForcdEdgeNum);
 InnerEdgeStat1 = INT_Vector(TotEdgeNum);
 ForcdEdgeStat1 = INT_Vector(TotEdgeNum);
 BoundEdgeStat1 = INT_Vector(TotEdgeNum);
 ForcdValue    = CMPLX_Vector(TotForcdEdgeNum);
}
 TotInnerEdgeNum = 0;
 TotBoundEdgeNum = 0;
 TotForcdEdgeNum = 0;

 for(IEdge=0;IEdge<TotEdgeNum;IEdge++)
 {
     if(ForceStat[IEdge] == 0)   {
        if(EdgeStatus[IEdge]==-1)                   {
        InnerEdgeStat[TotInnerEdgeNum] = IEdge;
        InnerEdgeStat1[IEdge] = TotInnerEdgeNum;
        TotInnerEdgeNum++;
                                        }
        else if(EdgeStatus[IEdge]==-2)              {
        BoundEdgeStat[TotBoundEdgeNum]= IEdge;
        BoundEdgeStat1[IEdge] = TotBoundEdgeNum;
        TotBoundEdgeNum++;
                                        }
                       }
     else              {
        ForcdEdgeStat[TotForcdEdgeNum]= IEdge;
        ForcdEdgeStat1[IEdge] = TotForcdEdgeNum;
        ForcdValue[TotForcdEdgeNum] = COMplex_Cmplx(
				    EdgeStatus[IEdge]*cos(EdgeStatus1[IEdge]),
				    EdgeStatus[IEdge]*sin(EdgeStatus1[IEdge]));
        TotForcdEdgeNum++;
                           }
       }
for(IEdge=0;IEdge<TotNumHexHedron;IEdge++)
{
temp=(-1)*(Sigma[IEdge]/(AbsPermitt*Omega));
RELPerm[IEdge] = COMplex_Cmplx(Epsilon[IEdge],temp);
}

printf("TotInnerEdgeNum = %d\n",TotInnerEdgeNum);
printf("TotBoundEdgeNum = %d\n",TotBoundEdgeNum);
printf("TotForcdEdgeNum = %d\n",TotForcdEdgeNum);
}

/***************************************************************************
*   FUNCTION        : edge()
This routine assigns global numbers to the local edges in the hexahedron.
This is the same as AssignHexHedronEdgeNumbering. 

****************************************************************************/

void edge(int nn,int i,int j,int k)
{
int kk,xx,dy1,dy2,divz,divz1;
   
    for (kk=0;kk<18;kk++) 
       xn[nn][kk] = 0;
    divz=(6*XdiM*YdiM+3*XdiM+3*YdiM+1)*k;
    divz1=XdiM*YdiM*6 + XdiM*3 +YdiM*3 + 1 ;
    dy1 = (3*XdiM+1)*j; dy2 = (3*XdiM+2)*j;
    xn[nn][0] = i + dy1 + divz;
    xn[nn][1] = XdiM + 2*i + dy1 + divz;
    xn[nn][2] = xn[nn][1] + 1; 
    xn[nn][3] = xn[nn][2] + 1;
    xn[nn][4] = 3*XdiM + 1 + i + dy1 + divz;
    xx = (YdiM+1)*XdiM + (2*XdiM+1)*YdiM;
    xn[nn][5] = xx + 2*i + dy2 + divz;
    xn[nn][6] = xn[nn][5] + 1; 
    xn[nn][7] = xn[nn][6] + 1;
    xx = xx + (2*XdiM+1); 
    xn[nn][8] = xx + i + dy2 + divz;
    xn[nn][9]= xn[nn][8] + 1;
    xx = xx + XdiM + 1; 
    xn[nn][10]= xx + 2*i + dy2 + divz;
    xn[nn][11]= xn[nn][10] + 1;
    xn[nn][12]= xn[nn][11] + 1; 
    xn[nn][13]= xn[nn][0] + divz1;
    xn[nn][14]= xn[nn][1] + divz1; 
    xn[nn][15]= xn[nn][2] + divz1;
    xn[nn][16]= xn[nn][3] + divz1; 
    xn[nn][17]= xn[nn][4] + divz1;
  return;
}

/***************************************************************************
*   FUNCTION        : deter_edge()
This looks at the position of each edge and decides whether the edge is inside
the region or on the boundary. If on the boundary, it assigns the value of 
EdgeStatus to -2 (boundary edge) or sets it to -1 (inner edge). This is done 
hexahedron by hexahedron.

****************************************************************************/

void deter_edge(int nn,int i,int j,int k)
{
int scratch;

     scratch=xn[nn][0];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 1\n",scratch);  fflush(stdout);
        }
     if ((j==0)||(k==0)) EdgeStatus[scratch] = -2; /* Is it on the boundary? */
     else                EdgeStatus[scratch] = -1;

     scratch=xn[nn][1];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 2\n",scratch); fflush(stdout);
        }
     if ((i==0)||(k==0))   EdgeStatus[scratch] = -2;
     else                  EdgeStatus[scratch] = -1;

     scratch=xn[nn][2];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 3\n",scratch); fflush(stdout);
         fflush(stdout);
       }
     if (k==0)            EdgeStatus[scratch] = -2;
     else                 EdgeStatus[scratch] = -1;            
           
     scratch = xn[nn][3];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 4\n",scratch); fflush(stdout);
       }
     if ((i==(XdiM-1))||(k==0))    EdgeStatus[scratch] = -2; 
     else                          EdgeStatus[scratch] = -1;

     scratch = xn[nn][4];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 5\n",scratch); fflush(stdout);
       }
     if ((j==(YdiM-1))||(k==0))   EdgeStatus[scratch] = -2; 
     else                         EdgeStatus[scratch] = -1; 

     scratch = xn[nn][5];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 6\n",scratch); fflush(stdout);
       }
     if ((j==0)||(i==0))   EdgeStatus[scratch] = -2;
     else                  EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][6];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 7\n",scratch); fflush(stdout);
       }
     if (j==0)   EdgeStatus[scratch] = -2;
     else        EdgeStatus[scratch] = -1;

     scratch = xn[nn][7];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 8\n",scratch); fflush(stdout);
       }
     if ((i==(XdiM-1))||(j==0))  EdgeStatus[scratch] = -2;
     else                        EdgeStatus[scratch] = -1;

     scratch = xn[nn][8];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 9\n",scratch); fflush(stdout);
       }
     if (i==0)    EdgeStatus[scratch] = -2;
     else         EdgeStatus[scratch] = -1;

     scratch = xn[nn][9];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 10\n",scratch); fflush(stdout);
       }
     if (i==(XdiM-1))   EdgeStatus[scratch] = -2;
     else               EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][10];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 11\n",scratch); fflush(stdout);
       }
     if ((j==(YdiM-1))||(i==0))  EdgeStatus[scratch] = -2;
     else                        EdgeStatus[scratch] = -1;

     scratch = xn[nn][11];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 12\n",scratch); fflush(stdout);
       }
     if (j==(YdiM-1))    EdgeStatus[scratch] = -2;
     else                EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][12];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 13\n",scratch); fflush(stdout);
       }
     if ((i==(XdiM-1))||(j==(YdiM-1))) EdgeStatus[scratch] = -2;
     else                              EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][13];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
        {
          printf("scratch=%d out of bounds spot 14\n",scratch); fflush(stdout);
        }
     if ((k==(ZdiM-1))||(j==0))   EdgeStatus[scratch] = -2;
     else                         EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][14];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
        {
          printf("scratch=%d out of bounds spot 15\n",scratch); fflush(stdout);
        }
     if ((k==(ZdiM-1))||(i==0))    EdgeStatus[scratch] = -2;
     else                          EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][15];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
        {
          printf("scratch=%d out of bounds spot 16\n",scratch); fflush(stdout);
        }
     if (k==(ZdiM-1))  EdgeStatus[scratch] = -2;
     else              EdgeStatus[scratch] = -1;
 
     scratch = xn[nn][16];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
       {
         printf("scratch=%d out of bounds spot 17\n",scratch); fflush(stdout);
       }
     if ((i==(XdiM-1))||(k==(ZdiM-1)))  EdgeStatus[scratch] = -2;
     else                               EdgeStatus[scratch] = -1;

     scratch = xn[nn][17];
     if ((scratch >= TotEdgeNum) || (scratch < 0))
        {
          printf("scratch=%d out of bounds spot 18\n",scratch); fflush(stdout);
        }
     if ((j==(YdiM-1))||(k==(ZdiM-1)))  EdgeStatus[scratch] = -2;
     else                               EdgeStatus[scratch] = -1;
  return;
 }

/***************************************************************************
*   FUNCTION        : Produce_Output()
This routine writes out the required electric field values to the output files.
It first writes out the default output file (if required) and then calls the 
EfieldatNode routine to write out electric at nodes as given in the input 
files.  

****************************************************************************/

void Produce_Output()
{
 int    i, j, k, x1, y1, z1, x2, y2, z2, Count_i, dflag=0;
 char   type[20], att1[20], att2[8], att3[8], att4[8];
 char   buffer[80];
 FILE   *OutF;
if(freq_step==0){
   EfieldData = CMPLX_Vector( TotEdgeNum);
          }
/* Open output files for writing */

if(freq_step>0){
InF = fopen(Ifile, "r");
OutF_0=fopen(Out_FileName0,"a");
OutF_1=fopen(Out_FileName1,"a");
OutF_2=fopen(Out_FileName2,"a");
OutF_3=fopen(Out_FileName3,"a");
}

/* Create the complete electric field vector */

   for(i=0;i<TotEdgeNum;i++) {EfieldData[i]=COMplex_Null();}
   for(Count_i=0;Count_i<TotInnerEdgeNum;Count_i++)
   {EfieldData[InnerEdgeStat[Count_i]] = RHSVector[Count_i]; }
   for(Count_i=0;Count_i<TotForcdEdgeNum;Count_i++)
   {EfieldData[ForcdEdgeStat[Count_i]] = ForcdValue[Count_i];}

/* Read input file again to read the filenames and coordinates */

        while (fgets(buffer, 80, InF))
        {
         if (buffer[0] == '#') {
                              if(OutF_0 != NULL) fprintf(OutF_0, "%s", buffer);
                              if(OutF_1 != NULL) fprintf(OutF_1, "%s", buffer);
                              if(OutF_2 != NULL) fprintf(OutF_2, "%s", buffer);
                              if(OutF_3 != NULL) fprintf(OutF_3, "%s", buffer);
                              continue;
                              }

/* default_output keyword */
        if(!strncmp(buffer, "default_output",14)) {dflag=1; continue;}

/* efield_output keyword */
        if(!strncmp(buffer, "efield_output",13))
                   { 
                   sscanf(buffer, "%s%d%d%d%d%d%d%s%s%s%s", type, &x1, &y1, &z1, &x2, &y2, &z2, att1, att2, att3, att4);
                   if (x2<x1) swap(&x1, &x2);
                   if (y2<y1) swap(&y1, &y2);
                   if (z2<z1) swap(&z1, &z2);
                   if (!strcmp(att1, Out_FileName1)) OutF=OutF_1;
                   if (!strcmp(att1, Out_FileName2)) OutF=OutF_2;
                   if (!strcmp(att1, Out_FileName3)) OutF=OutF_3;
                   for(k=z1;k<=z2;k++){
                     for(j=y1;j<=y2;j++){
                        for(i=x1;i<=x2;i++){
                                            EfieldatNode(i, j, k, OutF);
					/* write out E field values at nodes */
                                           }
                                        }
                                      } 
                   continue;
                   }

        }                  /*end of while */
        if(dflag == 1)     /*prints default output data*/ 
              {
               for(Count_i=0;Count_i<TotEdgeNum;Count_i++){
               fprintf(OutF_0,"%d\t%f\t%f\n",Count_i,EfieldData[Count_i].x,EfieldData[Count_i].y);}
              }

       /* fclose(OutF_0);
        fclose(OutF_1);
        fclose(OutF_2);
        fclose(OutF_3);
        fclose(InF);*/
        fprintf(stdout,"\nProgram Successfully Completed: \n\n");
}


/***************************************************************************
*   FUNCTION        : EfieldatNode()
Writes out electric fields at nodes. It first takes the node coordinates and 
converts it to a global number. It then uses this global number to calculate 
the edges which have this node as one of their end-nodes. It then discards the 
diagonal edges. The value at the node is calculated by averaging the values 
along the two edges in each direction (2 for x, y and z directions).
****************************************************************************/

void EfieldatNode(int X,int Y,int Z,FILE* Out_F)
{
int i,j,p,INode,JNode,IEdge,NodeNum[3][2],CountX=0,CountY=0,CountZ=0;
float Xcord=0.0,Ycord=0.0,Zcord=0.0;
complex ValX,ValY,ValZ;
for (p=0;p<X;p++)
   {
   Xcord = Xcord + DivisorX[p+1];
   }
for (p=0;p<Y;p++)
   {
   Ycord = Ycord + DivisorY[p+1];
   }
for (p=0;p<Z;p++)
   {
   Zcord = Zcord + DivisorZ[p+1];
   }

for(i=0;i<3;i++) for(j=0;j<2;j++) {NodeNum[i][j]=0;}

/* Calculate global node number */ 

for(JNode=0;JNode<TotNodeNum;JNode++)
 {
 if((NodeCord[JNode][0]==Xcord) && (NodeCord[JNode][1]==Ycord) && (NodeCord[JNode][2]==Zcord))
 {INode = JNode; break;} else {continue;}
 }
/* Find out which edges are incident on the node  */
   for (IEdge=0;IEdge<TotEdgeNum;IEdge++)
    {
 
    if(((TetGlobalEdgeEnds[IEdge][0])-1==INode)
    || ((TetGlobalEdgeEnds[IEdge][1])-1==INode) )
    {

     if((NodeCord[TetGlobalEdgeEnds[IEdge][0]-1][1]==NodeCord[TetGlobalEdgeEnds[IEdge][1]-1][1]) &&
        (NodeCord[TetGlobalEdgeEnds[IEdge][0]-1][2]==NodeCord[TetGlobalEdgeEnds[IEdge][1]-1][2]))
        {NodeNum[0][CountX] = IEdge; CountX++;}
 
      else if((NodeCord[TetGlobalEdgeEnds[IEdge][0]-1][2]==NodeCord[TetGlobalEdgeEnds[IEdge][1]-1][2]) &&
              (NodeCord[TetGlobalEdgeEnds[IEdge][0]-1][0]==NodeCord[TetGlobalEdgeEnds[IEdge][1]-1][0]))
              {NodeNum[1][CountY] = IEdge; CountY++;}
 
      else if((NodeCord[TetGlobalEdgeEnds[IEdge][0]-1][0]==NodeCord[TetGlobalEdgeEnds[IEdge][1]-1][0]) &&
              (NodeCord[TetGlobalEdgeEnds[IEdge][0]-1][1]==NodeCord[TetGlobalEdgeEnds[IEdge][1]-1][1]))
              {NodeNum[2][CountZ] = IEdge; CountZ++;}
      else
              continue;
    }
   }

/* Average the two edge values to get the value at each node */
 
  ValX = Real_Mul(0.5,COMplex_Add(EfieldData[NodeNum[0][0]],EfieldData[NodeNum[0][1]]));
  ValY = Real_Mul(0.5,COMplex_Add(EfieldData[NodeNum[1][0]],EfieldData[NodeNum[1][1]]));
  ValZ = Real_Mul(0.5,COMplex_Add(EfieldData[NodeNum[2][0]],EfieldData[NodeNum[2][1]]));

/* Print out to file */
  
fprintf(Out_F,"%6.3f\t%6.3f\t%6.3f\t%f\t%f\t%f\t%f\t%f\t%f\n",Xcord,Ycord,Zcord,ValX.x,ValX.y,ValY.x,ValY.y,ValZ.x,ValZ.y);
}

/*************************** END OF PROGRAM ********************************/
