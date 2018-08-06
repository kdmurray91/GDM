//
// NNLS_Double.h
//
#ifndef __NNLS_DOUBLE_H__
#define __NNLS_DOUBLE_H__


double *WeightedNNLSRegression( double *pEnvDataMatrix, int nRows, int nCols, 
								double *pRespVector, double *pDeviance, double *pWeights);

double CalcGDMDevianceDouble( double *pY, double *pU, double *pW, int nLen );

double GetWeightedNULLDeviance( int nRows, double *pRespVector, double *pWeights );

double *CopyEnvMatrixDouble( double *pMatrix, int nRows, int nCols );

double *nnlsFITDouble( double *pEnvDataMatrix, int nRows, int nCols, double *pRespVector, double *pWeights );

#endif  // __NNLS_DOUBLE_H__
