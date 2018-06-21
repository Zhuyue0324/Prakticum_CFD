#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int **Flag,
  int NameQuestion,
  double Th,
  double Tc,
  double **Temp,
  double INUI,
  double INVI,
  double INTI
);

//set boundary for inflow
void spec_boundary_vals(
  int i,
  int j,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **Temp,
  double INUI,
  double INVI,
  double INTI
);

#endif
