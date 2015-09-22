/////////////////////////////////////////////////////////////////////

#ifndef _limioptic_
#define _limioptic_

#include <cstdlib>
#include <cstdio>

#include "climioptic.h"

#ifdef __cplusplus
extern "C" {
#endif

void Clear();
void ClearParticles();
void PrintParticles();
void AddParticle(double,double,double,double,double,double);
void AddGaussBeam(double,double,double,double,double,double,double,double,double,double,double,double,double);
int GetParticleNum();
int GetParticleSize();
void ClearBeamline();
void PrintBeamline();
void AddMatrix(int,double*,double);
void AddDrift(int,double,double);
void AddBeamProfile(int);
void AddWaist();
void AddBeamWaist();
void AddSlit(double,double,double,double,int);
void AddModifyEmittance(double,double);
void ChangeBeamParameters(double, double, double, double, double, double, double, double, double);
void ChangeBeamParameters2(double, double, double, double);
void AddThinLens(int,double,double,double);
void AddQuadrupolRadFoc(int,double,double,double);
void AddAMSQuadrupolRadFoc(int,double,double,double,double);
void AddAMSQuadrupolAxFoc(int,double,double,double,double);
void AddQuadrupolAxFoc(int,double,double,double);
void AddESD(int,double,double,double,double,double,double);
void AddEdgeFocusing(int,double,double,double);
void AddHomDeflectingMagnet(int,double,double,double,double);
void AddInhomDeflectingMagnet(int,double,double,double);
void ClearTrajectories();
void PrintTrajectories();
void CalculateTrajectories();
int GetTrajectoriesSize();
double GetSpotSize();
double GetSigmaX();
double GetSigmaY();
void GetTrajectories(double *);
int GetTrajectorySize();
void GetTrajectory(int,int,double *);
void ApplyDrift(double*,int,double,double);
void ApplyBeamProfile(double*, int);
void ApplyWaist(double*);
void ApplySlit(double*,double,double,double,double,int);
void ApplyModifyEmittance(double*,int,double,double);
void ApplyChangeBeamParameters(double*, int, double, double, double, double, double, double, double, double, double);
void ApplyChangeBeamParameters2(double*, int, double, double, double, double);
void ApplyESD(double*,int,double,double,double,double,double,double);
void ApplyEdgeFocusing(double*,int,double,double,double);
void ApplyHomDeflectingMagnet(double*,int,double,double,double,double);
void ApplyInhomDeflectingMagnet(double*,int,double,double,double);
void ApplyQuadrupolRadFoc(double*,int,double,double,double);
void ApplyAMSQuadrupolRadFoc(double*,int,double,double,double,double);
void ApplyQuadrupolAxFoc(double*,int,double,double,double);
void ApplyAMSQuadrupolAxFoc(double*,int,double,double,double,double);

#ifdef __cplusplus
}
#endif

#endif

