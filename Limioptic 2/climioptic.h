/////////////////////////////////////////////////////////////////////

#ifndef _climioptic_
#define _climioptic_

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>

using namespace std;

class CLimioptic {
public:
   CLimioptic();
   CLimioptic(const CLimioptic&); // Copy-Constructor
   virtual ~CLimioptic();

   CLimioptic& operator=(const CLimioptic&); // Zuweisungsoperator

   void Clear();

   void ClearParticles(); // alle Startteilchen loeschen
   void PrintParticles(); // alle Startteilchen ausdrucken
   // Startteilchen hinzufuegen
   void AddParticle(double,double,double,double,double,double);
   int GetParticleNum(); // Anzahl der Teilchen zurueckgeben
   int GetParticleSize(); // Anzahl der Werte, die ein Teilchen definieren, zurueckgeben

   void ClearBeamline(); // Alle ionenoptischen Elemente der Beamline loeschen
   void PrintBeamline();

   void AddMatrix(int,double*,double); // Allgemeine Matrix zur Beamline hinzufuegen
   void AddDrift(int,double,double); // Driftstrecke zur Beamline hinzufuegen
   void AddBeamProfile();
   void AddSlit(double,double,double,double);
   void AddModifyEmittance(double, double);
   void ChangeBeamParameters(double, double, double, double);
   void ChangeBeamParameters2(double, double, double, double);
   void AddThinLens(int,double,double,double); // Duenne Linse zur Beamline hinzufuegen
   // radial fokussierenden Quadrupol hinzufuegen
   void AddQuadrupolRadFoc(int,double,double,double);
   void AddAMSQuadrupolRadFoc(int,double,double,double,double);
   void AddAMSQuadrupolAxFoc(int,double,double,double,double);
   // axial fokussierenden Quadrupol hinzufuegen
   void AddQuadrupolAxFoc(int,double,double,double);
   // Elektrostatischer Deflektor zur Beamline hinzufuegen
   void AddESD(int,double,double,double,double,double,double);
   // Kantenfokussierung eines Dipolmagneten hinzufuegen
   void AddEdgeFocusing(int,double,double,double);
   // homogener Ablenkmagnet hinzufuegen
   void AddHomDeflectingMagnet(int,double,double,double,double);
   void AddInhomDeflectingMagnet(int,double,double,double);

   void ClearTrajectories(); // Teilchen-Trajektorien loeschen
   void PrintTrajectories(); // alle Teilchen-Trajectorien ausgeben
   void CalculateTrajectories(); // Teilchen-Trajektorien berechnen
   int GetTrajectoriesSize(); // Groesse des Arrays mit den Trajektorien zurueckgeben
   void GetTrajectories(double *); // Teilchen-Trajektorien in ein externes Array kopieren
   // Groesse einer Eigenschaft der Trajektorie eines Teilchens zurueckgeben
   int GetTrajectorySize();
   // Eine Eigenschaft der Trajektorie eines Teilchens kopieren
   void GetTrajectory(int,int,double *);

   void ApplyMatrix(double*,int,double*);
   void ApplyDrift(double*,int,double,double);
   void ApplyBeamProfile(double*);
   void ApplySlit(double*,double,double,double,double);
   void ApplyModifyEmittance(double*, int, double, double);
   void ApplyChangeBeamParameters(double*, int, double, double, double, double);
   void ApplyChangeBeamParameters2(double*, int, double, double, double, double);
   void ApplyThinLens(double*,int,double,double,double);
   void ApplyQuadrupolRadFoc(double*,int,double,double,double);
   void ApplyAMSQuadrupolRadFoc(double*,int,double,double,double,double);
   void ApplyAMSQuadrupolAxFoc(double*,int,double,double,double,double);
   void ApplyQuadrupolAxFoc(double*,int,double,double,double);
   void ApplyESD(double*,int,double,double,double,double,double,double);
   void ApplyEdgeFocusing(double*,int,double,double,double);
   void ApplyHomDeflectingMagnet(double*,int,double,double,double,double);
   void ApplyInhomDeflectingMagnet(double*,int,double,double,double);

   // die Startteilchen:
   // A particle is represented by 8 doubles.
   // particle[0] = x-Abstand zu Sollbahn (radiale Ortsabweichung) in [mm]
   // particle[1] = x' bezueglich Sollbahn (radiale Richtungsabweichung) in [mrad]
   // particle[2] = y-Abstand zu Sollbahn (axiale Ortsabweichung) in [mm]
   // particle[3] = y' bezueglich Sollbahn (axiale Richtungsabweichung) in [mrad]
   // particle[4] = longitudinale Ortsabweichung in [mm]
   // particle[5] = relative Impulsabweichung in [promille]
   // particle[6] = Bisher zurueckgelegte Strecke auf der Sollbahn des Teilchens in [m]
   // particle[7] = Index des ionenoptischen Elements, in dem sich das Teilchen befindet (-1 = Start)
   // Die naechsten 8 Eintraege sind dann fuer das naechste Teilchen.
   // x-Achse (auch sog. radiale Richtung) zeigt in Strahlrichtung gesehen nach links
   // y-Achse (auch sog. axiale Richtung) zeigt nach oben
   vector<double> particles;
   static const int particlesize = 8; // enthaelt die Anzahl der Parameter pro Teilchen

   // Die beamline besteht aus einer Anzahl von ionenoptischen Elementen.
   // Jedes ionenoptische Element wird durch einen vector<double> repraesentiert.
   // Diser vector mit Eintraegen vom Typ double enthalt eine ID,
   // welche des Typ bestimmt (z.B. Drift), dann die Anzahl (haeufig braucht man
   // z.B. 100 mal eine identische Transfermatrix), sowie weitere Parameter. Die Anzahl und
   // Bedeutung der weiteren Parameter ist vom Typ abhaengig.
   //  ID  |  Typ
   // -----------
   //  1.0    allgemeine Matrix
   //  2.0    Drift
   //  3.0    Duenne Linse
   //  4.0    radial fokussierender Quadrupol
   //  5.0    axial fokussierender Quadrupol
   //  6.0    ESD   (Elektrostatischer Deflektor)
   //  7.0    Kantenfokussierung eines Dipolmagneten
   //  8.0    homogener Ablenkmagnet
   vector<vector<double> > beamline;

   // Die Teilchen-Trajectorien:
   // Hier werden die 'particles' nach Durchlaufen jedes einzelnen Elements
   // der 'beamline' gespeichert. Zuerst alle Teilchen beim Start, dann alle Teilchen
   // nach der ersten Transfermatrix usw.
   vector<double> trajectories;
};

#endif

/////////////////////////////////////////////////////////////////////
