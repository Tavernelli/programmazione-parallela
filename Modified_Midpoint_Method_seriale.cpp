#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*
const size_t nvar = 100;
const float xs = 0.0;
const float htot = 3.0;
const int   nstep = 10; //dimensione step

//funzioni che mi servono di  nrutil.h ---------------

	void nrerror(char error_text[])
	/* Numerical Recipes standard error handler */
	{
		fprintf(stderr, "Numerical Recipes run-time error...\n");
		fprintf(stderr, "%s\n", error_text);
		fprintf(stderr, "...now exiting to system...\n");
		exit(1);
	}


	float *vector(long nl, long nh)
	/* allocate a float vector with subscript range v[nl..nh] */
	{
		float *v;

		v = (float *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(float)));
		if (!v) nrerror("allocation failure in vector()");
		return v - nl + NR_END;
	}


	void free_vector(float *v, long nl, long nh)

	/* free a float vector allocated with vector() */
	{
		free((FREE_ARG)(v + nl - NR_END));
	}

//---------------------------------------

	void mmid
	(
		float y[],
		float dydx[],
		int nvar,
		float xs,
		float htot,
		int nstep,
		float yout[],
		void(*derivs)(float, float[], float[])
	)
	{
		int n, i;
		float x, swap, h2, h, *ym, *yn;

		ym = vector(1, nvar);
		yn = vector(1, nvar);
		h = htot / nstep;




		for (i = 1; i <= nvar; i++)
		{
			ym[i] = y[i];
			yn[i] = y[i] + h*dydx[i];

		}


		x = xs + h;
		(*derivs)(x, yn, yout);
		h2 = 2.0*h;
		for (n = 2; n <= nstep; n++)
		{

			for (i = 1; i <= nvar; i++)
			{
				swap = ym[i] + h2*yout[i];
				ym[i] = yn[i];
				yn[i] = swap;
			}
			x += h;
			(*derivs)(x, yn, yout);


		}

		for (i = 1; i <= nvar; i++)
			yout[i] = 0.5*(ym[i] + yn[i] + h*yout[i]);

		free_vector(yn, 1, nvar);
		free_vector(ym, 1, nvar);
	}


	//funzione g in ingresso
	float g(float x)
	{
		return x*x;
		//return x*x*x;
	}

	//derivata della funzione g
	float devg(float x)
	{
		return 2.0*x;
		//return 3.*x*x;
	}


	//derivata rispetto a x di g
	void derivs(float x, float y[], float dydx[])
	{
		for (int i = 1; i <= nvar; ++i)
		{
			dydx[i] = devg(x);
		}
	}


int main()
{	
	double start_time = MPI_Wtime ();
	float x[nvar + 1];

	for (int i = 0; i < sizeof x / sizeof *x; i++)
	{
		x[i] = i;

	}

	float y[nvar + 1];

	for (int i = 0; i < sizeof y / sizeof *y; i++)
	{
		y[i] = g(x[i]);
	}

	float dydx[nvar + 1];

	for (int i = 0; i < sizeof dydx / sizeof *dydx; i++)
	{
		dydx[i] = devg(x[i]);
	}


	//valore da cui partire
	

	float youtput[nvar + 1] = { 0 };

	mmid(y, dydx, nvar, xs, htot, nstep, youtput, derivs);

	//print output


	for (int i = 1; i < nvar; ++i)
	{
		printf("youtput[%d]=%f\n", i, youtput[i]);
	}

	printf("in tempo: %f \n", MPI_Wtime () - start_time);

	return (0);
}
