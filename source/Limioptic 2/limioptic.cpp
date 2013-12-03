/////////////////////////////////////////////////////////////////////

#include "limioptic.h"

CLimioptic optic;

void Clear()
{
   optic.Clear();
}

void ClearParticles()
{
   optic.ClearParticles();
}

void PrintParticles()
{
   optic.PrintParticles();
}

void AddParticle(double xdiff, double xangle, double ydiff, double yangle,
   double deltak, double deltam)
{
   optic.AddParticle(xdiff, xangle, ydiff, yangle, deltak, deltam);
}

void AddGaussBeam(double x, double sx, double a, double sa,
   double y, double sy,double b,double sb,double k,double sk,double m,double sm,double number)
{
   optic.AddGaussBeam(x, sx, a, sa, y, sy, b, sb, k, sk, m, sm, number);
}


int GetParticleNum()
{
   return optic.GetParticleNum();
}

int GetParticleSize()
{
   return optic.GetParticleSize();
}

void ClearBeamline()
{
   optic.ClearBeamline();
}

void PrintBeamline()
{
   optic.PrintBeamline();
}

void AddMatrix(int num,double *mat,double length)
{
   optic.AddMatrix(num,mat,length);
}

void AddDrift(int num,double gamma2,double length)
{
   optic.AddDrift(num,gamma2,length);
}

void AddBeamProfile()
{
   optic.AddBeamProfile();
}

void AddWaist()
{
   optic.AddWaist();
}

void AddSlit(double x,double dx,double y,double dy)
{
   optic.AddSlit(x,dx,y,dy);
}

void AddModifyEmittance(double factor1,double factor2)
{
   optic.AddModifyEmittance(factor1, factor2);
}

void ChangeBeamParameters(double dk, double dm, double strag_k, double strag_m, double strag_x, double strag_y, double strag_dx, double strag_dy, double percentage)
{
   optic.ChangeBeamParameters(dk, dm, strag_k, strag_m, strag_x, strag_y, strag_dx, strag_dy, percentage);
}

void ChangeBeamParameters2(double dk, double dm, double strag_k, double strag_m)
{
   optic.ChangeBeamParameters2(dk, dm, strag_k, strag_m);
}

void AddThinLens(int num,double fx,double fy,double length)
{
   optic.AddThinLens(num,fx,fy,length);
}

void AddESD(int num,double gamma2,double alpha,double rho0,double r0,double beta0, double korrektur)
{
   optic.AddESD(num,gamma2,alpha,rho0,r0,beta0,korrektur);
}

void AddEdgeFocusing(int num,double r,double beta,double betaeff)
{
   optic.AddEdgeFocusing(num,r,beta,betaeff);
}

void AddHomDeflectingMagnet(int num,double gamma2,double r,double alpha, double korrektur)
{
   optic.AddHomDeflectingMagnet(num,gamma2,r,alpha, korrektur);
}

void AddInhomDeflectingMagnet(int num,double rho,double phi,double n1)
{
   optic.AddInhomDeflectingMagnet(num,rho,phi,n1);
}

void AddQuadrupolRadFoc(int num,double gamma2,double k,double l)
{
   optic.AddQuadrupolRadFoc(num,gamma2,k,l);
}

void AddAMSQuadrupolRadFoc(int num,double gamma2,double kx,double ky,double l)
{
   optic.AddAMSQuadrupolRadFoc(num,gamma2,kx,ky,l);
}

void AddQuadrupolAxFoc(int num,double gamma2,double k,double l)
{
   optic.AddQuadrupolAxFoc(num,gamma2,k,l);
}

void AddAMSQuadrupolAxFoc(int num,double gamma2,double kx,double ky,double l)
{
   optic.AddAMSQuadrupolAxFoc(num,gamma2,kx,ky,l);
}

void ClearTrajectories()
{
   optic.ClearTrajectories();
}

void PrintTrajectories()
{
   optic.PrintTrajectories();
}

void CalculateTrajectories()
{
   optic.CalculateTrajectories();
}

int GetTrajectoriesSize()
{
   return optic.GetTrajectoriesSize();
}

double GetSpotSize()
{
   return optic.GetSpotSize();
}

double GetSigmaX()
{
   return optic.GetSigmaX();
}

double GetSigmaY()
{
   return optic.GetSigmaY();
}

void GetTrajectories(double *dst)
{
   optic.GetTrajectories(dst);
}

int GetTrajectorySize()
{
   return optic.GetTrajectorySize();
}

void GetTrajectory(int iparticle,int iproperty,double *dst)
{
   optic.GetTrajectory(iparticle,iproperty,dst);
}

void ApplyDrift(double *p,int n,double gamma2,double length)
{
   optic.ApplyDrift(p,n,gamma2,length);
}

void ApplyBeamProfile(double *p)
{
   optic.ApplyBeamProfile(p);
}

void ApplyWaist(double *p)
{
   optic.ApplyWaist(p);
}

void ApplySlit(double *p,double x,double dx,double y,double dy)
{
   optic.ApplySlit(p,x,dx,y,dy);
}

void ApplyModifyEmittance(double *p,int n, double factor1, double factor2)
{
   optic.ApplyModifyEmittance(p,n,factor1,factor2);
}

void ApplyChangeBeamParameters(double *p, int n, double dk, double dm, double strag_k, double strag_m, double strag_x, double strag_y, double strag_dx, double strag_dy, double percentage)
{
   optic.ApplyChangeBeamParameters(p, n, dk, dm, strag_k, strag_m, strag_x, strag_y, strag_dx, strag_dy, percentage);
}

void ApplyChangeBeamParameters2(double *p, int n, double dk, double dm, double strag_k, double strag_m)
{
   optic.ApplyChangeBeamParameters2(p, n, dk, dm, strag_k, strag_m);
}

void ApplyESD(double *p,int n,double gamma2,double alpha,double rho0,
   double r0,double beta0, double korrektur)
{
   optic.ApplyESD(p,n,gamma2,alpha,rho0,r0,beta0, korrektur);
}

void ApplyEdgeFocusing(double *p,int n,double r,double beta,double betaeff)
{
   optic.ApplyEdgeFocusing(p,n,r,beta,betaeff);
}

void ApplyHomDeflectingMagnet(double *p,int n,double gamma2,double r,double alpha, double korrektur)
{
   optic.ApplyHomDeflectingMagnet(p,n,gamma2,r,alpha, korrektur);
}

void ApplyInhomDeflectingMagnet(double *p,int num,double rho,double phi,double n1)
{
   optic.ApplyInhomDeflectingMagnet(p,num,rho,phi,n1);
}

void ApplyQuadrupolRadFoc(double *p,int n,double gamma2,double k,double l)
{
   optic.ApplyQuadrupolRadFoc(p,n,gamma2,k,l);
}

void ApplyAMSQuadrupolRadFoc(double *p,int n,double gamma2,double kx,double ky,double l)
{
   optic.ApplyAMSQuadrupolRadFoc(p,n,gamma2,kx,ky,l);
}

void ApplyQuadrupolAxFoc(double *p,int n,double gamma2,double k,double l)
{
   optic.ApplyQuadrupolAxFoc(p,n,gamma2,k,l);
}

void ApplyAMSQuadrupolAxFoc(double *p,int n,double gamma2,double kx,double ky,double l)
{
   optic.ApplyAMSQuadrupolAxFoc(p,n,gamma2,kx,ky,l);
}

/////////////////////////////////////////////////////////////////////

