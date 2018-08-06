//
// NNLS_Double.cpp
//
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// support functions
//#include "Message.h"
#include "myCallback.h"
#include "NNLS_Double.h"
#include "NNLSDoubleCore.h"

static int nMaxIterations = 12;
static double dEpsilon = 0.0001;

//
// Do the WEIGHTED iterative NNLS regression process 
//
// Takes a nRows X nCols matrix pEnvDataMatrix and a vector pRespVector
// of length nRows and iteratively produces a non-negative
// vector of coefficients representing the best fit of the matrix 
// against the vector using a default weights vector of ones.
//
double *WeightedNNLSRegression( double *pEnvDataMatrix, int nRows, int nCols, 
								double *pRespVector, double *pDeviance, double *pWeights )
{
	int i;

	//
	// allocate some temporary working vectors
	//
	double *pNN = new double [ nRows ];
	double *pUU = new double [ nRows ];
	double *pWW = new double [ nRows ];
	double *pZZ = new double [ nRows ];

    // allocate copy of environment data, as the inner functions modify it
    // replaces old file-based preservation
    size_t nElem = (size_t)nRows * (size_t)nCols;
    double *pCopyOfEnvDataMatrix = new double [nElem];
    memcpy(pCopyOfEnvDataMatrix, pEnvDataMatrix, nElem);


	//
	// initialise pNN prior to regression loop
	//
	
	// Modified by DNL:
	//int nCurrent = 0;
	for ( i=0; i<nRows; i++ )
	{
		pNN[i] = -(log( 1.0 - ( ( ( pRespVector[i] * pWeights[i] ) + 0.5 ) / ( pWeights[i] + 1.0 ) ) ) );
	}


	//
	// do the regression...
	//
	double *pCoeff = NULL;
	double dOldDev = 0.0;
	double dNewDev = 0.0;
	// Modified by DNL:
	//nCurrent = 5;
	for ( int nIter=0; nIter<nMaxIterations; nIter++ )
	{
		for ( i=0; i<nRows; i++ ) 
		{
			// do the initial transformation...
			pUU[i] = 1.0 - exp( -pNN[i] );

			// recalculate the weights
			pWW[i] = sqrt( pWeights[i] * (( 1.0 - pUU[i] ) / pUU[i]) );

			// calculate the new weighted transform to send to nnls
			pZZ[i] = pNN[i] + ( ( pRespVector[i] - pUU[i] ) / ( 1.0 - pUU[i] ) );
		}

		////////////////////////////////////////////////////////////////////////////////
		//
		// check the deviance of the fit, compare to the previous deviance,
		// and if convergence conditions are reached, exit the loop.
		//
		// NOTE: we send a copy of the environmental data matrix as 
		//	     this data gets modified inside the called function
		//
		if (pCoeff) delete[] pCoeff;
		pCoeff = nnlsFITDouble( pEnvDataMatrix, nRows, nCols, pZZ, pWW );
		if ( NULL == pCoeff )
		{
			if ( pNN ) delete[] pNN;
			if ( pUU ) delete[] pUU;
			if ( pWW ) delete[] pWW;
			if ( pZZ ) delete[] pZZ;
            if ( pCopyOfEnvDataMatrix ) delete[] pCopyOfEnvDataMatrix;
			return( NULL );
		}

		// 
		// restore values in the matrix data block
		//
        memcpy(pEnvDataMatrix, pCopyOfEnvDataMatrix, nElem);
		
		//
		// calculate the deviance for this fit 
		//
		dNewDev = CalcGDMDevianceDouble( pRespVector, pUU, pWeights, nRows );
		double dExitVal = ( dNewDev - dOldDev ) / ( dOldDev + dEpsilon );

		// check for convergence...
		if ( dEpsilon > fabs( dExitVal ) )
		{
			// adequate convergence reached
			break;	
		}

		else
		{			
			// setup for another loop iteration...
			dOldDev = dNewDev;

			// apply the new transformation to vector pNN
			for ( i=0; i<nRows; i++ )
			{
				double dVal = 0.0;

				for ( int j=0; j<nCols; j++ )
				{
					// do calculation for this term
					dVal += pCoeff[j] * pEnvDataMatrix[ ( j * nRows ) + i ];
				}

				// set final solution for this row item
				pNN[i] = dVal;				
			}			
		} 	
	} // for ( int nIter=0; nIter<nMaxIterations; nIter++ )


	// clean up
	if ( pNN ) delete[] pNN;
	if ( pUU ) delete[] pUU;
	if ( pWW ) delete[] pWW;
	if ( pZZ ) delete[] pZZ;
    if ( pCopyOfEnvDataMatrix ) delete[] pCopyOfEnvDataMatrix;
	*pDeviance = dNewDev;
	return( pCoeff );
}


//
// Create a NULL model and return the WEIGHTED deviance
//
double GetWeightedNULLDeviance( int nRows, double *pRespVector, double *pWeights )
{
	//
	// The data matrix is just a single vector of ONES
	//
	double *pNullVector = new double [ nRows ];
	for ( int i=0; i<nRows; i++ ) 
	{
		pNullVector[i] = 1.0;
	}

	//
	// allocate some temporary working vectors
	//
	double *pNN = new double [ nRows ];
	double *pUU = new double [ nRows ];
	double *pWW = new double [ nRows ];
	double *pZZ = new double [ nRows ];

	//
	// initialise pNN prior to regression loop
	//
	for ( int i=0; i<nRows; i++ )
	{
		pNN[i] = -(log( 1.0 - ( ( ( pRespVector[i] * pWeights[i] ) + 0.5 ) / ( pWeights[i] + 1.0 ) ) ) );
	}

	//
	// do the regression...
	//
	double *pCoeff = NULL;
	double dOldDev = 0.0;
	double dNewDev = 0.0;
	for ( int nIter=0; nIter<nMaxIterations; nIter++ )
	{
		for ( int i=0; i<nRows; i++ ) 
		{
			// do the initial transformation...
			pUU[i] = 1.0 - exp( -pNN[i] );

			// recalculate the weights
			pWW[i] = sqrt( pWeights[i] * (( 1.0 - pUU[i] ) / pUU[i]) );

			// calculate the new weighted transform to send to nnls
			pZZ[i] = pNN[i] + ( ( pRespVector[i] - pUU[i] ) / ( 1.0 - pUU[i] ) );
		}

		////////////////////////////////////////////////////////////////////////////////
		//
		// check the deviance of the fit, compare to the previous deviance,
		// and if convergence conditions are reached, exit the loop.
		//
		// NOTE: we send a copy of the environmental data matrix as 
		//	     this data gets modified inside the called function
		//
		double *pCopyOfNullVector = CopyEnvMatrixDouble( pNullVector, nRows, 1 );
		if (pCoeff) delete[] pCoeff;
		pCoeff = nnlsFITDouble( pCopyOfNullVector, nRows, 1, pZZ, pWW );
		if ( pCopyOfNullVector ) delete[] pCopyOfNullVector;
		////////////////////////////////////////////////////////////////////////////////

		//
		// calculate the deviance for this fit 
		//
		dNewDev = CalcGDMDevianceDouble( pRespVector, pUU, pWeights, nRows );
		double dExitVal = ( dNewDev - dOldDev ) / ( dOldDev + dEpsilon );

		// check for convergence...
		if ( dEpsilon > fabs( dExitVal ) )
		{
			// adequate convergence reached
			break;	
		}

		else
		{
			// setup for another loop iteration...
			dOldDev = dNewDev;

			// apply the new transformation to vector pNN
			for ( int i=0; i<nRows; i++ )
			{
				// set final solution for this row item
				pNN[i] = pCoeff[0];				
			}
		} 	
	} // for ( int nIter=0; nIter<nMaxIterations; nIter++ )


	// clean up
	if ( pNN ) delete[] pNN;
	if ( pUU ) delete[] pUU;
	if ( pWW ) delete[] pWW;
	if ( pZZ ) delete[] pZZ;
	if ( pCoeff ) delete[] pCoeff;
	if ( pNullVector ) delete[] pNullVector;
	return( dNewDev );
}



//
// Function for calculating a single iteration deviance
//
double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, int nLen )
{
	double dTotal = 0.0;
	double t1, t2;
	for ( int i=0; i<nLen; i++ )
	{
		if ( pU[i] == 0.0 ) 
			t1 = pY[i];

		else if ( pY[i] == 0.0 )
			t1 = 0.0;

		else
			t1 = pY[i] * log( pY[i] / pU[i] );


		if ( pU[i] == 1.0 ) 
			t2 = 1.0 - pY[i];

		else if ( pY[i] == 1.0 )
			t2 = 1.0 - pY[i];

		else
			t2 = ( 1.0 - pY[i] ) * log( ( 1.0 - pY[i] ) / ( 1.0 - pU[i] ) );

		// accumulate the running sum
		dTotal += pW[i] * ( t1 + t2 );
	}
	// final multiply by 2
	return( dTotal * 2.0 );
}


//
// create a copy of the Environmental Data Matrix to send to nnls
// It gets modified inside nnls, so we need to give it a copy at each iteration
//
double *CopyEnvMatrixDouble( double *pMatrix, int nRows, int nCols )
{
	int nTotal = nRows * nCols;
	double *pNew = new double [ nTotal ];
	for ( int i=0; i<nTotal; i++ )
	{
		pNew[i] = pMatrix[i];
	}
	return( pNew );
}



//
// This function calles the JPL code from "Solving Least Squares Problems," by C. Lawson and R. Hanson
//
// The WEIGHTED version of the NNLS Algorithm
//
// Returns the coefficients[nCols] for the best fit of the data to the reponse.
//
double *nnlsFITDouble( double *pEnvDataMatrix, int nRows, int nCols, double *pRespVector, double *pWeights )
{
	//
	// get the data matrix memory block
	//
	double *a = pEnvDataMatrix;
	if ( NULL == a )
	{
		//Message( "nnlsFITDouble() - NULL EnvDataMatrix!", "ERROR" );
		return(NULL);
	}


	//
	// get the response vector memory block
	//
	double *b = pRespVector;
	if ( NULL == b )
	{
		//Message( "nnlsFITDouble() - NULL pRespVector!", "ERROR" );
		return(NULL);
	}


	//
	// setup the parameters for the call to nnls
	// 
	int mda = nRows;
	int m = nRows;
	int n = nCols;

	double rnorm;
	double *x  = new double [ nCols ];
	double *w  = new double [ nCols ];
	double *zz = new  double [ nRows ];

	int *indx  = new int [ nCols + 10 ];
	int mode;

	//////////////////////////////////////////////////////////////////////////
	//  do the matrix regression
	//
	// this is the new version with WEIGHTING
	nnls_Weighted( a, &mda, &m, &n, b, pWeights, x, &rnorm, w, zz, indx, &mode );
	//////////////////////////////////////////////////////////////////////////

	//
	// show results
	//
	double *pCoeffs = NULL;
	//char buff[64];
	if ( mode == 1 )
	{
		//sprintf( buff, "THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY." );
		//DEBUG( buff, "INFO" );

		pCoeffs = new double [ n ];
		for ( int i=0; i<n; i++ )
		{
			//sprintf( buff, "item %d of %d: %lf", i+1, n, x[i] );
			//DEBUG( buff, "INFO" );
			pCoeffs[ i ] = x[ i ];
		}
	}

	else if ( mode == 2 )
	{
		//sprintf( buff, "THE DIMENSIONS OF THE PROBLEM ARE BAD." );
		//Message( buff, "INFO" );
	}

	else if ( mode == 3 )
	{
		//sprintf( buff, "ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS." );
		//Message( buff, "INFO" );
	}

	
	//
	// clean up local allocations...
	//
	if ( x )    delete[] x;
	if ( w )    delete[] w;
	if ( zz )   delete[] zz;
	if ( indx ) delete[] indx;

	return( pCoeffs );
}

