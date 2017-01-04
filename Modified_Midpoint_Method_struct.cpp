#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <memory.h>
#include <iostream>
#include <cstdlib>

const int SIZE = 1;
const float xs	  = 0.0; //valore da cui partire
const float htot  = 3.0;  //H
const int   nstep = 10; //dimensione step
const size_t nvar = 300000; //dimensione vettore
const size_t nvarproc = (nvar+SIZE-1) / SIZE;//ceil(nvar/SIZE); //dimensione problema
const size_t nvarproc_r = nvar % SIZE;
const size_t nvarproc_max = nvarproc + nvarproc_r;

//struct per inviare i dati
struct mul_data
{ 
	int   local_nvarproc;
	float h;
	float ym[nvarproc_max] { 0 };
	float dydx[nvarproc_max] { 0 };
};

float g(float x)
{
	return x*x;
	//return x*x*x;
}

//derivata della funzione g
float devg(float x)
{
	return 2.0*x ;
	//return 3.*x*x;
}

//derivata rispetto a x di g
void derivs(float x, float y[], float dydx[], int size)
{
	for (int i = 0; i < size; i++)
	{
		dydx[i] = devg(x);
	}
}


int main (int argc, char*argv[]) 
{

	int size = 0; 
	int rank = 0;
	MPI_Init(&argc, &argv);
	double start_time = MPI_Wtime ();
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	int i;
	
	MPI_Datatype mul_data_type;
	{
		//dati struct
		MPI_Datatype oldtypes[4];
		int          blockcounts[4]; 
		MPI_Aint     offsets[4];
		//attributi
		
		offsets[0]    = offsetof(mul_data,local_nvarproc);
		oldtypes[0]   = MPI_INT;
		blockcounts[0]= 1;
		
		offsets[1]    = offsetof(mul_data,h);
		oldtypes[1]   = MPI_FLOAT;
		blockcounts[1]= 1;


		offsets[2]    = offsetof(mul_data,ym);
		oldtypes[2]   = MPI_FLOAT;
		blockcounts[2]= nvarproc_max;

		offsets[3]    = offsetof(mul_data,dydx);
		oldtypes[3]   = MPI_FLOAT;
		blockcounts[3]= nvarproc_max;

		//crea tipo mpi
		MPI_Type_create_struct (4, blockcounts, offsets, oldtypes, &mul_data_type);
		MPI_Type_commit (&mul_data_type);
	}


	//---------------------------------------
	if (rank == 0)
	{
		//creo e riempo i vettori 
		float x[nvar];
		for (i = 0; i < sizeof(x)/sizeof(*x); i++) x[i] = i;

		float y[nvar];
		for (i = 0; i < sizeof y / sizeof *y; i++) y[i]= g(x[i]);
			
		float dydx[nvar];
		for (i = 0; i < sizeof dydx / sizeof *dydx; i++) dydx[i]= devg(x[i]);

		float ym [nvar];
		for (i=0; i < nvar; i++) ym [i] = 0;

		float yn [nvar];
		for (i=0; i < nvar; i++) yn [i] = 0;
		
		float h = htot/nstep;
		for (i = 0; i< nvar; i++) ym [i] = y[i];

		float* youtput = new float[nvarproc_max*size];

		memset(youtput,0,sizeof(float)*nvarproc_max*size);

		//riempo la struct
		mul_data data_to_send;
 		data_to_send.h = h;
		//-------------------------------------------------------------------
		for (i = 1; i < size; i++)
		{	
			//per gestire il resto controllo l'ultimo
			if (data_to_send.local_nvarproc = i == size-1) data_to_send.local_nvarproc = nvarproc_max;	
			else data_to_send.local_nvarproc = nvarproc;

			size_t offset = nvarproc*(i-1);
 			memcpy(data_to_send.ym, &ym[offset],sizeof(float)*data_to_send.local_nvarproc);
 			memcpy(data_to_send.dydx, &dydx[offset],sizeof(float)*data_to_send.local_nvarproc);
		    MPI_Send(&data_to_send,1,mul_data_type,i,0,MPI_COMM_WORLD);
		}
		//--------------------------------------------------------------------

		for (i = 1; i < size; i++)
		{	
			//per gestire il resto controllo l'ultimo
			int local_nvarproc;
			if (local_nvarproc = i == size-1) local_nvarproc = nvarproc_max;
			else local_nvarproc = nvarproc;
			
			MPI_Recv(&youtput[nvarproc*(i-1)],local_nvarproc,MPI_FLOAT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}

		//--------------------------- --
	    printf("in tempo[0]: %f \n", MPI_Wtime() - start_time);	
		for (i = 0; i < nvar; ++i)
		{
			printf ("\tyoutput[%i]: %f\n ",i, youtput[i]);
		}

		delete youtput;
	}
	else
	{
		mul_data data_to_receve;

		float vecD  [nvarproc_max];
		float swap  [nvarproc_max];
		float yout 	[nvarproc_max];
		float h2 = 0;
		memset(vecD,0,sizeof(float)*nvarproc_max);
		memset(swap,0,sizeof(float)*nvarproc_max);
		memset(yout,0,sizeof(float)*nvarproc_max);
	    
		MPI_Recv(&data_to_receve,1,mul_data_type,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		for(i=0;i<data_to_receve.local_nvarproc;i++)
		{
			vecD [i] = data_to_receve.ym[i] + data_to_receve.h*data_to_receve.dydx[i];
		}
	                
		float x1;
	    x1=xs+data_to_receve.h;
	    
		derivs(x1, vecD, yout, data_to_receve.local_nvarproc);
		h2 = 2. * data_to_receve.h;
		int n = 0;

	    for (n =2; n <=  nstep; n++)
	    {	
			for(i=0;i<data_to_receve.local_nvarproc;i++)
	        {
	            swap[i]  = data_to_receve.ym[i] + h2*yout[i];
	            data_to_receve.ym[i]  = vecD [i];
	            vecD [i] = swap[i];
	        }
		    x1 += data_to_receve.h;
		    derivs(x1, vecD, yout, data_to_receve.local_nvarproc);
			    
		}
		
		for(i=0;i<data_to_receve.local_nvarproc;i++)
		{
			yout[i] = 0.5*(data_to_receve.ym[i] + vecD[i] + data_to_receve.h*yout[i]);
		}

		MPI_Send(&yout,data_to_receve.local_nvarproc,MPI_FLOAT,0,0,MPI_COMM_WORLD);
		printf("in tempo[%d]: %f \n", rank,MPI_Wtime () - start_time);
	}


		

	MPI_Finalize();
	return 0;

}
