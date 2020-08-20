/////////////////////////////////////////////////////////////////////

#include "climioptic.h"
#include <iostream>
//#include <fstream>
#include <random>
#include <stdio.h>

#include <sstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()



CLimioptic::CLimioptic()
{
    ClearParticles();
    ClearBeamline();
    ClearTrajectories();
    spotsize = 0.;
    nottransmitted = 0;
    x_verteilung = 0.;
    y_verteilung = 0.;
}

CLimioptic::CLimioptic(const CLimioptic &copyobj)
{
    particles      = copyobj.particles;
    beamline       = copyobj.beamline;
    trajectories   = copyobj.trajectories;
    spotsize       = copyobj.spotsize;
    nottransmitted = copyobj.nottransmitted;
    x_verteilung   = copyobj.x_verteilung;
    y_verteilung   = copyobj.y_verteilung;
}

CLimioptic::~CLimioptic()
{
}

CLimioptic &CLimioptic::operator=(const CLimioptic &assignobj)
{
    if (this != &assignobj) // no self assignment
    {
        particles      = assignobj.particles;
        beamline       = assignobj.beamline;
        trajectories   = assignobj.trajectories;
        spotsize       = assignobj.spotsize;
        nottransmitted = assignobj.nottransmitted;
        x_verteilung   = assignobj.x_verteilung;
        y_verteilung   = assignobj.y_verteilung;
    }

    return *this;
}

void CLimioptic::Clear()
{
    ClearParticles();
    ClearBeamline();
    ClearTrajectories();
    spotsize = 0.;
    nottransmitted = 0;
    x_verteilung = 0.;
    y_verteilung = 0.;
}

void CLimioptic::ClearParticles()
{
    // Alle Startteilchen loeschen.
    particles.clear();
}

void CLimioptic::PrintParticles()
{
    unsigned int i;
    for (i = 0; i < particles.size(); i = i + 7)
    {
        printf("%f %f %f %f %f %f %f\n", particles[i], particles[i + 1], particles[i + 2],
               particles[i + 3], particles[i + 4], particles[i + 5], particles[i + 6]);
    }
}

void CLimioptic::AddParticle(double xdiff, double xangle, double ydiff, double yangle,
                             double deltaz, double deltam)
{
    double z    =  0.0;
    double iele = -1.0;

    particles.push_back(xdiff);
    particles.push_back(xangle);
    particles.push_back(ydiff);
    particles.push_back(yangle);
    particles.push_back(deltaz);
    particles.push_back(deltam);
    particles.push_back(z);
    particles.push_back(iele);
}

void CLimioptic::AddGaussBeam(double x, double sx, double a, double sa,
                              double y, double sy, double b, double sb, double k, double sk, double m, double sm, double number)
{
    double z    =  0.0;
    double iele = -1.0;

    int num = (int)number;
    //std::cout << num << endl;

    std::default_random_engine generator;
    std::normal_distribution<double> x_distribution(x, sx);
    std::normal_distribution<double> a_distribution(a, sa);
    std::normal_distribution<double> y_distribution(y, sy);
    std::normal_distribution<double> b_distribution(b, sb);
    std::normal_distribution<double> k_distribution(k, sk);
    std::normal_distribution<double> m_distribution(m, sm);

    int i;
    for (i = 0; i < num; i++)
    {
        particles.push_back(x_distribution(generator));
        particles.push_back(a_distribution(generator));
        particles.push_back(y_distribution(generator));
        particles.push_back(b_distribution(generator));
        particles.push_back(k_distribution(generator));
        particles.push_back(m_distribution(generator));
        particles.push_back(z);
        particles.push_back(iele);
    }
}

int CLimioptic::GetParticleNum()
{
    return (int)(particles.size() / particlesize);
}

int CLimioptic::GetParticleSize()
{
    return particlesize;
}

void CLimioptic::ClearBeamline()
{
    // Alle ionenoptischen Elemente der Beamline loeschen.
    beamline.clear();
}

void CLimioptic::PrintBeamline()
{
    unsigned int i, j;
    for (i = 0; i < beamline.size(); i++)
    {
        for (j = 0 ; j < beamline[i].size(); j++)
        {
            printf("%f ", beamline[i][j]);
        }
        printf("\n");
    }
}

void CLimioptic::AddMatrix(int number, double *mat, double length)
{
    // Eine allgemeine 6x6-Matrix zu Beamline hinzufuegen
    // 'mat' enthaelt die 36 Parameter
    int i;
    vector<double> ele;
    ele.clear();
    ele.push_back(1.0);
    ele.push_back((double)(number));
    for (i = 0; i < 36; i++)
    {
        ele.push_back(mat[i]);
    }
    ele.push_back(length);
    beamline.push_back(ele);
}

void CLimioptic::AddDrift(int number, double gamma2, double length)
{
    // Siehe ApplyDrift fuer Dokumentation der Parameter

    vector<double> drift;
    drift.clear();
    drift.push_back(2.0);
    drift.push_back((double)(number));
    drift.push_back(gamma2);
    drift.push_back(length);
    beamline.push_back(drift);
    if (length < 0.) { spotsize = 9999999999999999.; }
}

void CLimioptic::AddBeamProfile(int index)
{
    vector<double> bp;
    bp.clear();
    bp.push_back(10.0);
    bp.push_back(1.);
    bp.push_back((double)index);
    beamline.push_back(bp);
}

void CLimioptic::AddWaist()
{
    vector<double> w;
    w.clear();
    w.push_back(21.0);
    w.push_back(1.);
    beamline.push_back(w);
}

void CLimioptic::AddSlit(double x, double dx, double y, double dy, int output)
{
    vector<double> slit;
    slit.clear();
    slit.push_back(9.0);
    slit.push_back(1.);
    slit.push_back(x);
    slit.push_back(dx);
    slit.push_back(y);
    slit.push_back(dy);
    slit.push_back(output);
    beamline.push_back(slit);
}

void CLimioptic::AddModifyEmittance(double factor1, double factor2)
{
    vector<double> mod;
    mod.clear();
    mod.push_back(11.0);
    mod.push_back(1.);
    mod.push_back(factor1);
    mod.push_back(factor2);
    beamline.push_back(mod);
}

void CLimioptic::ChangeBeamParameters(double dk, double dm, double strag_k, double strag_m, double strag_x, double strag_y, double strag_dx, double strag_dy, double percentage)
{
    vector<double> param;
    param.clear();
    param.push_back(12.0);
    param.push_back(1.);
    param.push_back(dk);
    param.push_back(dm);
    param.push_back(strag_k);
    param.push_back(strag_m);
    param.push_back(strag_x);
    param.push_back(strag_y);
    param.push_back(strag_dx);
    param.push_back(strag_dy);
    param.push_back(percentage);
    beamline.push_back(param);
}

void CLimioptic::ChangeBeamParameters2(double dk, double dm, double strag_k, double strag_m)
{
    vector<double> param;
    param.clear();
    param.push_back(112.0);
    param.push_back(1.);
    param.push_back(dk);
    param.push_back(dm);
    param.push_back(strag_k);
    param.push_back(strag_m);
    beamline.push_back(param);
}

void CLimioptic::AddThinLens(int number, double fx, double fy, double length)
{
    // Siehe ApplyThinLens fuer Dokumentation

    vector<double> ThinLens;
    ThinLens.clear();
    ThinLens.push_back(3.0);
    ThinLens.push_back((double)(number));
    ThinLens.push_back(fx);
    ThinLens.push_back(fy);
    ThinLens.push_back(length);
    beamline.push_back(ThinLens);
}

void CLimioptic::AddQuadrupolRadFoc(int number, double gamma2, double k, double l)
{
    vector<double> qrf;
    qrf.clear();
    qrf.push_back(4.0);
    qrf.push_back((double)(number));
    qrf.push_back(gamma2);
    qrf.push_back(k);
    qrf.push_back(l);
    beamline.push_back(qrf);
}

void CLimioptic::AddQuadrupolAxFoc(int number, double gamma2, double k, double l)
{
    vector<double> qaf;
    qaf.clear();
    qaf.push_back(5.0);
    qaf.push_back((double)(number));
    qaf.push_back(gamma2);
    qaf.push_back(k);
    qaf.push_back(l);
    beamline.push_back(qaf);
}

void CLimioptic::AddAMSQuadrupolRadFoc(int num, double gamma2, double kx, double ky, double l)
{
    vector<double> qrf;
    qrf.clear();
    qrf.push_back(14.0);
    qrf.push_back((double)(num));
    qrf.push_back(gamma2);
    qrf.push_back(kx);
    qrf.push_back(ky);
    qrf.push_back(l);
    beamline.push_back(qrf);
}

void CLimioptic::AddAMSQuadrupolAxFoc(int num, double gamma2, double kx, double ky, double l)
{
    vector<double> qaf;
    qaf.clear();
    qaf.push_back(15.0);
    qaf.push_back((double)(num));
    qaf.push_back(gamma2);
    qaf.push_back(kx);
    qaf.push_back(ky);
    qaf.push_back(l);
    beamline.push_back(qaf);
}

void CLimioptic::AddESD(int number, double gamma2, double alpha, double rho0, double r0,
                        double beta0, double korrektur)
{
    vector<double> esd;
    esd.clear();
    esd.push_back(6.0);
    esd.push_back((double)(number));
    esd.push_back(gamma2);
    esd.push_back(alpha);
    esd.push_back(rho0);
    esd.push_back(r0);
    esd.push_back(beta0);
    esd.push_back(korrektur);
    beamline.push_back(esd);
}

void CLimioptic::AddEdgeFocusing(int number, double r, double beta, double betaeff)
{
    vector<double> ef;
    ef.clear();
    ef.push_back(7.0);
    ef.push_back((double)(number));
    ef.push_back(r);
    ef.push_back(beta);
    ef.push_back(betaeff);
    beamline.push_back(ef);
}

void CLimioptic::AddEdgeFocusingY(int number, double r, double beta, double betaeff)
{
    vector<double> ef;
    ef.clear();
    ef.push_back(17.0);
    ef.push_back((double)(number));
    ef.push_back(r);
    ef.push_back(beta);
    ef.push_back(betaeff);
    beamline.push_back(ef);
}

void CLimioptic::AddHomDeflectingMagnet(int number, double gamma2, double r, double alpha, double korrektur)
{
    vector<double> hdm;
    hdm.clear();
    hdm.push_back(8.0);
    hdm.push_back((double)(number));
    hdm.push_back(gamma2);
    hdm.push_back(r);
    hdm.push_back(alpha);
    hdm.push_back(korrektur);
    beamline.push_back(hdm);
}

void CLimioptic::AddHomDeflectingMagnetY(int number, double gamma2, double r, double alpha, double korrektur)
{
    vector<double> hdm;
    hdm.clear();
    hdm.push_back(19.0);
    hdm.push_back((double)(number));
    hdm.push_back(gamma2);
    hdm.push_back(r);
    hdm.push_back(alpha);
    hdm.push_back(korrektur);
    beamline.push_back(hdm);
}

void CLimioptic::AddInhomDeflectingMagnet(int number, double rho, double phi, double n1)
{
    vector<double> idm;
    idm.clear();
    idm.push_back(18.0);
    idm.push_back((double)(number));
    idm.push_back(rho);
    idm.push_back(phi);
    idm.push_back(n1);
    beamline.push_back(idm);
}

void CLimioptic::ClearTrajectories()
{
    trajectories.clear();
}

void CLimioptic::PrintTrajectories()
{
    unsigned int i, j;
    int k;

    for (i = 0; i < particles.size(); i = i + particlesize)
    {
        for (j = 0; j < trajectories.size(); j = j + particles.size())
        {
            for (k = 0; k < particlesize; k++)
            {
                printf("%f ", trajectories[j + i + k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}


void CLimioptic::CalculateTrajectories()
{
    unsigned int i, j;
    int n, nelements, ibeamline, itraj, type;

    int tsize;

    ClearTrajectories();
    nelements = 0;
    for (i = 0; i < beamline.size(); i++)
    {
        nelements = nelements + beamline[i][1];
    }

    // Am Anfang stehen immer die unveraenderten Teilchen, also +1
    nelements = nelements + 1;
    tsize     = nelements * particles.size();
    trajectories.assign(tsize, 0.);

    // printf("tsize = %d\n",tsize);

    itraj = 0;     // Ab jetzt ist 'itraj' der Index fuer 'trajectories'
    n     = 0;     // 'n' ist die Anzahl der bereits berechneten Schritte
    // Wir sind fertig, wenn n == nelements

    // Am Anfang stehen die unveraenderten Teilchen
    for (j = 0; j < particles.size(); j++)
    {
        trajectories[itraj] = particles[j];
        itraj++ ;
    }
    n++;

    ibeamline = 0;           // Index fuer 'beamline'
    while (n < nelements)    // Solange noch ionenoptische Elemente uebrig sind ..
    {
        type = beamline[ibeamline][0];
        switch (type)
        {
        case 1:
            ApplyMatrix(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  &beamline[ibeamline].front() + 2);
            break;
        case 2:
            ApplyDrift(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3]);
            break;
        case 3:
            ApplyThinLens(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4]);
            break;
        case 4:
            ApplyQuadrupolRadFoc(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4]);
            break;
        case 5:
            ApplyQuadrupolAxFoc(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4]);
            break;
        case 6:
            ApplyESD(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4],
                  beamline[ibeamline][5], beamline[ibeamline][6], beamline[ibeamline][7]);
            break;
        case 7:
            ApplyEdgeFocusing(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4]);
            break;
        case 17:
            ApplyEdgeFocusingY(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4]);
            break;
        case 8:
            ApplyHomDeflectingMagnet(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5]);
            break;
        case 9:
            ApplySlit(&trajectories.front() + itraj,
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5], beamline[ibeamline][6]);
            break;
        case 10:
            ApplyBeamProfile(&trajectories.front() + itraj,
                  beamline[ibeamline][2]);
            break;
        case 21:
            ApplyWaist(&trajectories.front() + itraj);
            break;
        case 14:
            ApplyAMSQuadrupolRadFoc(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5]);
            break;
        case 11:
            ApplyModifyEmittance(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3]);
            break;
        case 12:
            ApplyChangeBeamParameters(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5],
                  beamline[ibeamline][6], beamline[ibeamline][7], beamline[ibeamline][8], beamline[ibeamline][9], beamline[ibeamline][10]);
            break;
        case 13:
            ApplyChangeBeamParameters2(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5]);
            break;
        case 15:
            ApplyAMSQuadrupolAxFoc(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5]);
            break;
        case 18:
            ApplyInhomDeflectingMagnet(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4]);
            break;
        case 19:
            ApplyHomDeflectingMagnetY(&trajectories.front() + itraj, (int)beamline[ibeamline][1],
                  beamline[ibeamline][2], beamline[ibeamline][3], beamline[ibeamline][4], beamline[ibeamline][5]);
            break;
        default:
            printf("Unknown ion optic element! Aborting ..\n");
            exit(0);
            break;
        }

        n     += beamline[ibeamline][1];
        itraj += beamline[ibeamline][1] * particles.size();
        ibeamline++;
    }
}

int CLimioptic::GetTrajectoriesSize()
{
    return (int)(trajectories.size());
}

double CLimioptic::GetSpotSize()
{
    //cout << "spotsize return = " << spotsize << endl;
    if (spotsize != spotsize)
    {
        return 9999999999999999.;
    }
    else
    {
        return spotsize;
    }
}

double CLimioptic::GetSigmaX() { return x_verteilung; }
double CLimioptic::GetSigmaY() { return y_verteilung; }

void CLimioptic::GetTrajectories(double *dst)
{
    int i, d;
    double *src;

    d   = (int)(trajectories.size());
    src = &trajectories.front();
    for (i = 0; i < d; i++)
    {
        dst[i] = src[i];
    }
}

int CLimioptic::GetTrajectorySize()
{
    return (int)(trajectories.size() / particles.size());
}

void CLimioptic::GetTrajectory(int iparticle, int iproperty, double *dst)
{
    // iproperty = 0,..,particlesize-1
    // iparticle = 0,..,NumberOfParticles-1
    unsigned int i, count;

    count = 0;
    for (i = 0; i < trajectories.size(); i = i + particles.size())
    {
        dst[count] = trajectories[i + particlesize * iparticle + iproperty];
        count++;
    }
}

void CLimioptic::ApplyMatrix(double *p, int nmat, double *mat)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum    = elesize / particlesize;  // Anzahl der Teilchen

    // printf("%f %f %f\n",mat[6],mat[7],mat[8]);

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = mat[0]  * p[i + 0 - elesize] + mat[1]  * p[i + 1 - elesize] + mat[2]  * p[i + 2 - elesize] +
                 mat[3]  * p[i + 3 - elesize] + mat[4]  * p[i + 4 - elesize] + mat[5]  * p[i + 5 - elesize];
            p1 = mat[6]  * p[i + 0 - elesize] + mat[7]  * p[i + 1 - elesize] + mat[8]  * p[i + 2 - elesize] +
                 mat[9]  * p[i + 3 - elesize] + mat[10] * p[i + 4 - elesize] + mat[11] * p[i + 5 - elesize];
            p2 = mat[12] * p[i + 0 - elesize] + mat[13] * p[i + 1 - elesize] + mat[14] * p[i + 2 - elesize] +
                 mat[15] * p[i + 3 - elesize] + mat[16] * p[i + 4 - elesize] + mat[17] * p[i + 5 - elesize];
            p3 = mat[18] * p[i + 0 - elesize] + mat[19] * p[i + 1 - elesize] + mat[20] * p[i + 2 - elesize] +
                 mat[21] * p[i + 3 - elesize] + mat[22] * p[i + 4 - elesize] + mat[23] * p[i + 5 - elesize];
            p4 = mat[24] * p[i + 0 - elesize] + mat[25] * p[i + 1 - elesize] + mat[26] * p[i + 2 - elesize] +
                 mat[27] * p[i + 3 - elesize] + mat[28] * p[i + 4 - elesize] + mat[29] * p[i + 5 - elesize];
            p5 = mat[30] * p[i + 0 - elesize] + mat[31] * p[i + 1 - elesize] + mat[32] * p[i + 2 - elesize] +
                 mat[33] * p[i + 3 - elesize] + mat[34] * p[i + 4 - elesize] + mat[35] * p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;

            p[i + 6] = p[i + 6 - elesize] + mat[36];

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                // Index des ionenoptischen Elements raufzaehlen
                p[i + 7] = p[i + 7] + 1.0;
            }

            i += particlesize;
        }
    }
}

void CLimioptic::ApplyDrift(double *p, int nmat, double gamma2, double length)
{
    // Die Transfermarix lautet:
    //  |  1 length 0    0   0       0        |
    //  |  0    1   0    0   0       0        |
    //  |  0    0   1 length 0       0        |
    //  |  0    0   0    1   0       0        |
    //  |  0    0   0    0   1 length/gamma^2 |
    //  |  0    0   0    0   0       1        |
    // Hierbei ist 'length' die Laenge der Driftstrecke in m
    // Fuer den Lorenzfaktor 'gamma2' gilt:
    // gamma2 = \gamma^2 = E/mc^2  mit  E = T + mc^2
    // Fuer langsame Teilchen gilt also die Naeherung gamma2=1
    // 'nmat' ist die Anzahl dieser Matrizen, die nacheinander angewendet werden.
    // '*p' zeigt auf das zu fuellende Array, das die Teilchenparameter aufnimmt.
    // Die Teilchenparameter des Schrittes vorher liegen im gleichen Array VOR dem
    // zu fuellenden Bereich (die muessen also mit negativem Index angesprochen werden).

    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum    = elesize / particlesize;  // Anzahl der Teilchen
    length  = length / nmat;

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] + length * p[i + 1 - elesize];
            p1 = p[i + 1 - elesize];
            p2 = p[i + 2 - elesize] + length * p[i + 3 - elesize];
            p3 = p[i + 3 - elesize];
            p4 = p[i + 4 - elesize]; // '+length/gamma2*p[i+5-elesize];' entfernt
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + length;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }

            if (imat == 0)
            {
                // Index des ionenoptischen Elements raufzaehlen
                p[i + 7] = p[i + 7] + 1.0;
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyBeamProfile(double *p, double index)
{
    int pnum, i, j, elesize, ip, durch;
    double xplus, yplus, aplus, bplus;


    elesize = particles.size();
    pnum    = elesize / particlesize;  // Anzahl der Teilchen

    xplus = 0;  // summe der x-quadrate
    yplus = 0;  // summe der y-quadrate
    aplus = 0;  // summe der dx-quadrate
    bplus = 0;  // summe der dy-quadrate

    //std::ofstream datei;
    //datei.open(string("particles") + SSTR(index) + string(".dat"), ios::out);
    FILE * datei;
    datei = fopen("particles.dat", "w");

    i = 0;
    for (ip = 0; ip < pnum; ip++)
    {
        xplus += (p[i + 0 - elesize] * p[i + 0 - elesize]);
        aplus += (p[i + 1 - elesize] * p[i + 1 - elesize]);
        yplus += (p[i + 2 - elesize] * p[i + 2 - elesize]);
        bplus += (p[i + 3 - elesize] * p[i + 3 - elesize]);


        if (!((p[i + 0 - elesize] == 0.) && (p[i + 1 - elesize] == 0.) && (p[i + 2 - elesize] == 0.) && (p[i + 3 - elesize] == 0.)))
        {
            fprintf(datei, "%f ", (p[i + 0 - elesize]));
            fprintf(datei, "%f ", (p[i + 1 - elesize]));
            fprintf(datei, "%f ", (p[i + 2 - elesize]));
            fprintf(datei, "%f ", (p[i + 3 - elesize]));
            fprintf(datei, "%f ", (p[i + 4 - elesize]));
            fprintf(datei, "%f\n", (p[i + 5 - elesize]));
        }


        for (j = 0; j < 8; j++)
        {
            p[i + j] = p[i + j - elesize];
        }

        p[i + 7] = p[i + 7] + 1.0; // Index des ionenoptischen Elements raufzaehlen
        i = i + particlesize;
    }

    //datei.close();
    fclose(datei);

    //pnum -= nichtdurch; // die durch einen schlitz abgefangenen partikel sollen nicht mitgezaehlt werden
    durch = (double)(pnum - nottransmitted);

    printf("transmission (beamprofile) =\t%f (%i/%i particles)\t@%f m\n", 1. - (double)nottransmitted / pnum, durch, pnum, p[i + 6 - particlesize]);
    printf("SigmaX=\t%f\tSigmaA=\t%f\tEmittanzX=\t%f\n", sqrt(xplus / (durch - 1)), sqrt(aplus / (durch - 1)), sqrt(xplus / (durch - 1)) * sqrt(aplus / (durch - 1)));
    printf("SigmaY=\t%f\tSigmaB=\t%f\tEmittanzY=\t%f\n", sqrt(yplus / (durch - 1)), sqrt(bplus / (durch - 1)), sqrt(yplus / (durch - 1)) * sqrt(bplus / (durch - 1)));

    x_verteilung = sqrt(xplus / (durch - 1));
    y_verteilung = sqrt(yplus / (durch - 1));

    /*
    datei.open("beamprofile.dat", ios::out | ios::app);

    datei << nottransmitted / pnum << " ";      //transmission
    datei << sqrt(xplus / (durch - 1)) << " ";  //strahlbreite x
    datei << sqrt(aplus / (durch - 1)) << " ";  //divergenz a
    datei << sqrt(xplus / (durch - 1)) * sqrt(aplus / (durch - 1)) << " "; //emittanz x
    datei << sqrt(yplus / (durch - 1)) << " ";  //strahlbreite y
    datei << sqrt(bplus / (durch - 1)) << " ";  //divergenz b
    datei << sqrt(yplus / (durch - 1)) * sqrt(bplus / (durch - 1)) << "\n"; //emittanz y
    datei.close();
    */
}

void CLimioptic::ApplyWaist(double *p)
{
    int pnum, i, j, elesize, ip;
    double xplus, yplus;

    xplus = 0.;
    yplus = 0.;

    elesize = particles.size();
    pnum    = elesize / particlesize;  // Anzahl der Teilchen

    i = 0;
    for (ip = 0; ip < pnum; ip++)
    {
        xplus += p[i + 0 - elesize] * p[i + 0 - elesize];
        yplus += p[i + 2 - elesize] * p[i + 2 - elesize];

        for (j = 0; j < 8; j++)
        {
            p[i + j] = p[i + j - elesize];
        }

        p[i + 7] = p[i + 7] + 1.0; // Index des ionenoptischen Elements raufzaehlen
        i = i + particlesize;
    }

    spotsize += (xplus + yplus);
}

void CLimioptic::ApplySlit(double *p, double x, double dx, double y, double dy, int output)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;
    int nmat = 1, durch = 0, nichtdurch = 0;

    double xplus = 0.;
    double yplus = 0.;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {

            if ((p[i + 0 - elesize] <= x + dx / 2.) &&
                    (p[i + 0 - elesize] >= x - dx / 2.) &&
                    (p[i + 2 - elesize] <= y + dy / 2.) &&
                    (p[i + 2 - elesize] >= y - dy / 2.))
            {
                p0 = p[i + 0 - elesize];
                p1 = p[i + 1 - elesize];
                p2 = p[i + 2 - elesize];
                p3 = p[i + 3 - elesize];
                p4 = p[i + 4 - elesize];
                p5 = p[i + 5 - elesize];

                durch++;

                // Teilchen, die schon vorher beblockt wurden, sollen nicht gezeahlt werden
                if ((p0 == 0) && (p1 == 0) && (p2 == 0) && (p3 == 0))
                {
                    durch--;
                }
            }
            else
            {
                xplus += p[i + 0 - elesize] * p[i + 0 - elesize];
                yplus += p[i + 2 - elesize] * p[i + 2 - elesize];
                p0 = 0.;
                p1 = 0.;
                p2 = 0.;
                p3 = 0.;
                p4 = 0.;
                p5 = 0.;
                nichtdurch++;
                nottransmitted++;
            }

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize];

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;   // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }

    spotsize += (xplus + yplus);

    if (output == 1)
    {
        printf("transmission (slit) =\t%f", (double) durch / (double) (durch + nichtdurch));
        printf(",  total: %f", (1. - (double) nottransmitted / pnum));
        printf("\t@ %f m\n", p[i + 6 - particlesize]);
    }

    //fstream datei( "slit.dat",ios::out|ios::app);
    //datei << (double)durch/(durch+nichtdurch)<<" ";
    //datei.close();
}

void CLimioptic::ApplyThinLens(double *p, int nmat, double fx, double fy, double length)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize];
            p1 = -1 / fx * p[i + 0 - elesize] + p[i + 1 - elesize];
            p2 = p[i + 2 - elesize];
            p3 = -1 / fy * p[i + 2 - elesize] + p[i + 3 - elesize];
            p4 = p[i + 4 - elesize];
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + length;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}


void CLimioptic::ApplyModifyEmittance(double *p, int nmat, double factor1, double factor2)
{
    /*
    Die Emittanz veraendern (z.B. durch den Stripping-Prozess)
    */
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * factor1;
            p1 = p[i + 1 - elesize] * factor2;
            p2 = p[i + 2 - elesize] * factor1;
            p3 = p[i + 3 - elesize] * factor2;
            p4 = p[i + 4 - elesize];
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize];

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}


void CLimioptic::ApplyChangeBeamParameters(double *p, int nmat, double dk, double dm,
        double strag_k, double strag_m, double strag_x, double strag_y, double strag_dx, double strag_dy, double percentage)
{
    /*
    dk, dm aendern. Zb bei Folie. strag = stragling
    */
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    std::default_random_engine generator;
    std::normal_distribution<double> k_distribution(dk, strag_k);
    std::normal_distribution<double> k0_distribution(0., strag_k);
    std::normal_distribution<double> m_distribution(dm, strag_m);
    std::normal_distribution<double> x_distribution(0., strag_x);
    std::normal_distribution<double> y_distribution(0., strag_y);
    std::normal_distribution<double> dx_distribution(0., strag_dx);
    std::normal_distribution<double> dy_distribution(0., strag_dy);
    std::uniform_real_distribution<double> uf_distribution(0., 1.);

    elesize = particles.size();
    pnum    = elesize / particlesize;  // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            if ((p[i + 0 - elesize] == 0.) && (p[i + 1 - elesize] == 0.) && (p[i + 2 - elesize] == 0.) && (p[i + 3 - elesize] == 0.))
            {
                p0 = 0.;
                p1 = 0.;
                p2 = 0.;
                p3 = 0.;
                p4 = 0.;
                p5 = 0.;
            }
            else
            {
                p0 = p[i + 0 - elesize] + x_distribution(generator);
                p1 = p[i + 1 - elesize] + dx_distribution(generator);
                p2 = p[i + 2 - elesize] + y_distribution(generator);
                p3 = p[i + 3 - elesize] + dy_distribution(generator);
                if (uf_distribution(generator) > (1. - percentage))
                {
                    p4 = p[i + 4 - elesize] + k_distribution(generator);
                }
                else
                {
                    p4 = p[i + 4 - elesize] + k0_distribution(generator);
                }
                p5 = p[i + 5 - elesize] + m_distribution(generator);
            }

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize];

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }

            // Index des ionenoptischen Elements raufzaehlen
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;
            }

            i = i + particlesize;
        }
    }
}


void CLimioptic::ApplyChangeBeamParameters2(double *p, int nmat, double dk, double dm, double strag_k, double strag_m)
{
    /*
    dk, dm aendern. Zb bei Folie. strag = stragling
    */
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    std::default_random_engine generator;
    std::normal_distribution<double> k_distribution(dk * 1000., strag_k * sqrt(dk * dk - 2.*dk + 2) * 1000.);
    //std::normal_distribution<double> m_distribution(dm, strag_m);

    elesize = particles.size();
    pnum    = elesize / particlesize;  // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            if ((p[i + 0 - elesize] == 0.) && (p[i + 1 - elesize] == 0.) && (p[i + 2 - elesize] == 0.) && (p[i + 3 - elesize] == 0.))
            {
                p0 = 0.;
                p1 = 0.;
                p2 = 0.;
                p3 = 0.;
                p4 = 0.;
                p5 = 0.;
            }
            else
            {
                p0 = p[i + 0 - elesize];
                p1 = p[i + 1 - elesize];
                p2 = p[i + 2 - elesize];
                p3 = p[i + 3 - elesize];
                p4 = k_distribution(generator);
                p5 = dm; //m_distribution(generator);
            }

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize];

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }

            // Index des ionenoptischen Elements raufzaehlen
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;
            }

            i = i + particlesize;
        }
    }
}


void CLimioptic::ApplyQuadrupolRadFoc(double *p, int nmat, double gamma2, double k, double l)
{
    //   Die Transportmatrix erster Ordnung eines magnetischen Quadrupols ist gleich
    //   der eines elektrostatischen Quadrupols. Der Parameter k haengt im einen Fall
    //   von B ab und im anderen Fall von E. Die Laenge l ist gleich der effektiven Laenge
    //   (siehe im Hinterberger Seite 220).

    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    l = l / nmat;

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cos(sqrt(k) * l) + p[i + 1 - elesize] * sin(sqrt(k) * l) / sqrt(k);
            p1 = p[i + 0 - elesize] * (-sqrt(k) * sin(sqrt(k) * l)) + p[i + 1 - elesize] * cos(sqrt(k) * l);
            p2 = p[i + 2 - elesize] * cosh(sqrt(k) * l) + p[i + 3 - elesize] * sinh(sqrt(k) * l) / sqrt(k);
            p3 = p[i + 2 - elesize] * sqrt(k) * sinh(sqrt(k) * l) + p[i + 3 - elesize] * cosh(sqrt(k) * l);
            p4 = p[i + 4 - elesize]; // entfernt: '+l/gamma2*p[i+5-elesize];'
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyAMSQuadrupolRadFoc(double *p, int nmat, double gamma2, double kx, double ky, double l)
{
    //   Die Transportmatrix erster Ordnung eines magnetischen Quadrupols ist gleich
    //   der eines elektrostatischen Quadrupols. Der Parameter k haengt im einen Fall
    //   von B ab und im anderen Fall von E. Die Laenge l ist gleich der effektiven Laenge
    //   (siehe im Hinterberger Seite 220).

    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    l = l / nmat;

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cos(sqrt(kx) * l) + p[i + 1 - elesize] * sin(sqrt(kx) * l) / sqrt(kx);
            p1 = p[i + 0 - elesize] * (-sqrt(kx) * sin(sqrt(kx) * l)) + p[i + 1 - elesize] * cos(sqrt(kx) * l);
            p2 = p[i + 2 - elesize] * cosh(sqrt(ky) * l) + p[i + 3 - elesize] * sinh(sqrt(ky) * l) / sqrt(ky);
            p3 = p[i + 2 - elesize] * sqrt(ky) * sinh(sqrt(ky) * l) + p[i + 3 - elesize] * cosh(sqrt(ky) * l);
            p4 = p[i + 4 - elesize]; // entfernt: '+l/gamma2*p[i+5-elesize];'
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyQuadrupolAxFoc(double *p, int nmat, double gamma2, double k, double l)
{
    //   Die Transportmatrix erster Ordnung eines magnetischen Quadrupols ist gleich
    //   der eines elektrostatischen Quadrupols. Der Parameter k haengt im einen Fall
    //   von B ab und im anderen Fall von E. Die Laenge l ist gleich der effektiven Laenge
    //   (siehe im Hinterberger Seite 220).

    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    l = l / nmat;

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cosh(sqrt(k) * l) + p[i + 1 - elesize] * sinh(sqrt(k) * l) / sqrt(k);
            p1 = p[i + 0 - elesize] * sqrt(k) * sinh(sqrt(k) * l) + p[i + 1 - elesize] * cosh(sqrt(k) * l);
            p2 = p[i + 2 - elesize] * cos(sqrt(k) * l) + p[i + 3 - elesize] * sin(sqrt(k) * l) / sqrt(k);
            p3 = p[i + 2 - elesize] * (-sqrt(k) * sin(sqrt(k) * l)) + p[i + 3 - elesize] * cos(sqrt(k) * l);
            p4 = p[i + 4 - elesize]; // entfernt: '+l/gamma2*p[i+5-elesize];'
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyAMSQuadrupolAxFoc(double *p, int nmat, double gamma2, double kx, double ky, double l)
{

    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    //std::cout << kx << "    " << ky << "\n";

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cosh(sqrt(kx) * l) + p[i + 1 - elesize] * sinh(sqrt(kx) * l) / sqrt(kx);
            p1 = p[i + 0 - elesize] * sqrt(kx) * sinh(sqrt(kx) * l) + p[i + 1 - elesize] * cosh(sqrt(kx) * l);
            p2 = p[i + 2 - elesize] * cos(sqrt(ky) * l) + p[i + 3 - elesize] * sin(sqrt(ky) * l) / sqrt(ky);
            p3 = p[i + 2 - elesize] * (-sqrt(ky) * sin(sqrt(ky) * l)) + p[i + 3 - elesize] * cos(sqrt(ky) * l);
            p4 = p[i + 4 - elesize];
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyESD(double *p, int nmat, double gamma2, double alpha, double rho0,
                          double r0, double beta0, double korrektur)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;
    double l, ne, kx, ky;

    alpha = alpha / nmat;

    korrektur = (korrektur - 1) * 1000.;

    elesize   = particles.size();
    pnum      = elesize / particlesize;  // Anzahl der Teilchen

    l  = alpha * rho0;  // Laenge der Sollbahn im Elektrostaten
    ne = 1. + rho0 / r0; // elektrischer Feldindex
    kx = (3. - ne - beta0 * beta0) / (rho0 * rho0);
    ky = (ne - 1.) / (rho0 * rho0);

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cos(sqrt(kx) * l) +
                 p[i + 1 - elesize] * sin(sqrt(kx) * l) / sqrt(kx) +
                 (p[i + 4 - elesize] + korrektur) * (1 - cos(sqrt(kx) * l)) / (rho0 * kx); // entfernt: '*(2.0-beta0*beta0)/(rho0*kx)*(1.0-cos(sqrt(kx)*l));'
            p1 = p[i + 0 - elesize] * (-sqrt(kx) * sin(sqrt(kx) * l)) +
                 p[i + 1 - elesize] * cos(sqrt(kx) * l) +
                 (p[i + 4 - elesize] + korrektur) * sin(sqrt(kx) * l) / sqrt(kx) / rho0; //(2.0-beta0*beta0)/(rho0*sqrt(kx))*sin(sqrt(kx)*l);
            p2 = p[i + 2 - elesize] * cos(sqrt(ky) * l) +
                 p[i + 3 - elesize] * sin(sqrt(ky) * l) / sqrt(ky);
            p3 = p[i + 2 - elesize] * (-sqrt(ky) * sin(sqrt(ky) * l)) +
                 p[i + 3 - elesize] * cos(sqrt(ky) * l);
            p4 = p[i + 4 - elesize]; //weg: 'p[i+0-elesize]*(-sin(sqrt(kx)*l)/(rho0*sqrt(kx)))+p[i+1-elesize]*(-(1.0-cos(sqrt(kx)*l))/(rho0*kx))+p[i+4-elesize]+p[i+5-elesize]*(l/gamma2-(2.0-beta0*beta0)/(rho0*rho0*kx)*(l-sin(sqrt(kx)*l)/sqrt(kx)));'
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyEdgeFocusing(double *p, int nmat, double r, double beta, double betaeff)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize];
            p1 = p[i + 0 - elesize] * tan(beta) / r + p[i + 1 - elesize];
            p2 = p[i + 2 - elesize];
            p3 = p[i + 2 - elesize] * (-tan(betaeff) / r) + p[i + 3 - elesize];
            p4 = p[i + 4 - elesize];
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + 0.0; // Die Kante hat die Laenge Null

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyEdgeFocusingY(double *p, int nmat, double r, double beta, double betaeff)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize];
            p3 = p[i + 2 - elesize] * tan(beta) / r
               + p[i + 3 - elesize];
            p2 = p[i + 2 - elesize];
            p1 = p[i + 0 - elesize] * (-tan(betaeff) / r)
               + p[i + 1 - elesize];
            p4 = p[i + 4 - elesize];
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + 0.0; // Die Kante hat die Laenge Null

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyHomDeflectingMagnet(double *p, int nmat, double gamma2, double r,
        double alpha, double korrektur)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;
    double l;

    alpha = alpha / nmat;

    l = r * alpha;
    korrektur = ((korrektur * korrektur) - 1.) * 1000.;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cos(alpha) +
                 p[i + 1 - elesize] * r * sin(alpha) +
                 p[i + 4 - elesize] * r * (1. - cos(alpha)) / 2 +
                 (p[i + 5 - elesize] + korrektur) * r * (1. - cos(alpha)) / 2; // angepasst
            p1 = p[i + 0 - elesize] * (-sin(alpha) / r) +
                 p[i + 1 - elesize] * cos(alpha) +
                 p[i + 4 - elesize] * sin(alpha) / 2 +
                 (p[i + 5 - elesize] + korrektur) * sin(alpha) / 2; // angepasst
            p2 = p[i + 2 - elesize] +
                 p[i + 3 - elesize] * r * alpha;
            p3 = p[i + 3 - elesize];
            p4 = p[i + 4 - elesize]; //weg:  p[i+0-elesize]*(-sin(alpha))+p[i+1-elesize]*(-r*(1.0-cos(alpha)))+p[i+4-elesize]+p[i+5-elesize]*(r*alpha/gamma2-r*(alpha-sin(alpha)));
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyHomDeflectingMagnetY(double *p, int nmat, double gamma2, double r,
        double alpha, double korrektur)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;
    double l;

    alpha = alpha / nmat;

    l = r * alpha;
    korrektur = ((korrektur * korrektur) - 1.) * 1000.;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p2 = p[i + 2 - elesize] * cos(alpha) +
                 p[i + 3 - elesize] * r * sin(alpha) +
                 p[i + 4 - elesize] * r * (1. - cos(alpha)) / 2 +
                 (p[i + 5 - elesize] + korrektur) * r * (1. - cos(alpha)) / 2; // angepasst
            p3 = p[i + 2 - elesize] * (-sin(alpha) / r) +
                 p[i + 3 - elesize] * cos(alpha) +
                 p[i + 4 - elesize] * sin(alpha) / 2 +
                 (p[i + 5 - elesize] + korrektur) * sin(alpha) / 2; // angepasst
            p0 = p[i + 0 - elesize] +
                 p[i + 1 - elesize] * r * alpha;
            p1 = p[i + 1 - elesize];
            p4 = p[i + 4 - elesize]; //weg:  p[i+0-elesize]*(-sin(alpha))+p[i+1-elesize]*(-r*(1.0-cos(alpha)))+p[i+4-elesize]+p[i+5-elesize]*(r*alpha/gamma2-r*(alpha-sin(alpha)));
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + l;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}

void CLimioptic::ApplyInhomDeflectingMagnet(double *p, int nmat, double rho, double phi, double n1)
{
    int pnum, imat, i, j, elesize, ip;
    double p0, p1, p2, p3, p4, p5;
    double w, kx, ky, cx, cy, sx, sy, dx, Nk, Nm;

    phi = phi / nmat;

    elesize = particles.size();
    pnum = elesize / particlesize; // Anzahl der Teilchen

    w    =   rho * phi;
    kx   =   sqrt(1 - n1) / rho / rho;
    ky   =   sqrt(n1) / rho;
    cx   =   cos(kx * w);
    cy   =   cos(ky * w);
    sx   =   sin(kx * w) / kx;
    sy   =   sin(ky * w) / ky;
    dx   =   (1 - cx / rho / kx / kx);
    Nk   =   0.5;
    Nm =  0.5;

    //std::cout << "n=\t" << n1 << "\n";


    i = 0;
    for (imat = 0; imat < nmat; imat++)
    {
        for (ip = 0; ip < pnum; ip++)
        {
            p0 = p[i + 0 - elesize] * cx + p[i + 1 - elesize] * sx + p[i + 4 - elesize] * dx * Nk + p[i + 5 - elesize] * dx * Nm;
            p1 = p[i + 0 - elesize] * (-sx) * kx * kx + p[i + 1 - elesize] * cx + p[i + 4 - elesize] * (sx / rho) * Nk + p[i + 5 - elesize] * (sx / rho) * Nm;
            p2 = p[i + 2 - elesize] * cy + p[i + 3 - elesize] * sy;
            p3 = p[i + 2 - elesize] * (-sy * ky * ky) + p[i + 3 - elesize] * cy;
            p4 = p[i + 4 - elesize]; // weg:  p[i+0-elesize]*(-sin(alpha))+p[i+1-elesize]*(-r*(1.0-cos(alpha)))+p[i+4-elesize]+p[i+5-elesize]*(r*alpha/gamma2-r*(alpha-sin(alpha)));
            p5 = p[i + 5 - elesize];

            p[i + 0] = p0;
            p[i + 1] = p1;
            p[i + 2] = p2;
            p[i + 3] = p3;
            p[i + 4] = p4;
            p[i + 5] = p5;
            p[i + 6] = p[i + 6 - elesize] + w;

            // Kopiere die restlichen Teilcheneigenschaften (falls vorhanden)
            for (j = 7; j < particlesize; j++)
            {
                p[i + j] = p[i + j - elesize];
            }
            if (imat == 0)
            {
                p[i + 7] = p[i + 7] + 1.0;    // Index des ionenoptischen Elements raufzaehlen
            }

            i = i + particlesize;
        }
    }
}


/////////////////////////////////////////////////////////////////////

