// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPBetheBlochModel.cc                                //
//  Description: Compute energy loss for FCP particle vs. charge.     //
//                                                                    //
//  Adapted from G4BetheBlochModel by Samir Banik.                    //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "FCPBetheBlochModel.hh"
#include "G4Electron.hh"
#include "G4EmCorrections.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4double FCPBetheBlochModel::xgi[]={ 0.0199, 0.1017, 0.2372, 0.4083, 0.5917,
                                      0.7628, 0.8983, 0.9801 };
                                      
G4double FCPBetheBlochModel::wgi[]={ 0.0506, 0.1112, 0.1569, 0.1813, 0.1813,
                                      0.1569, 0.1112, 0.0506 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

FCPBetheBlochModel::FCPBetheBlochModel(const G4ParticleDefinition* p,
                                         const G4String& nam)
  : G4VEmModel(nam),
  particle(0),
  limitKinEnergy(100.*keV),
  logLimitKinEnergy(G4Log(limitKinEnergy)),
  twoln10(2.0*G4Log(10.0)),
  bg2lim(0.0169),
  taulim(8.4146e-3),
  alphaprime(fine_structure_const/twopi)
{
  theElectron = G4Electron::Electron();
  corr = G4LossTableManager::Instance()->EmCorrections();
  fParticleChange = 0;

  // initial initialisation of memeber should be overwritten
  // by SetParticle
  mass = massSquare = ratio = 1.0;

  if(p) { SetParticle(p); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FCPBetheBlochModel::~FCPBetheBlochModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBetheBlochModel::MinEnergyCut(const G4ParticleDefinition*,
                                           const G4MaterialCutsCouple* couple)
{
  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double FCPBetheBlochModel::MaxSecondaryEnergy(const G4ParticleDefinition*,
                                                 G4double kinEnergy) 
{
  G4double tau  = kinEnergy/mass;
  G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCPBetheBlochModel::Initialise(const G4ParticleDefinition* p,
                                     const G4DataVector&)
{
  if(p) { SetParticle(p); }
  if(!fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBetheBlochModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxKinEnergy)
{
  G4double cross = 0.0;
  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  G4double maxEnergy = min(tmax,maxKinEnergy);
  //####################Get PDG Charge###################
  G4double charge  =  p->GetPDGCharge()/eplus;
  G4double charge2 = charge*charge; 
  if(cutEnergy < maxEnergy) {

    G4double totEnergy = kineticEnergy + mass;
    G4double energy2   = totEnergy*totEnergy;
    G4double beta2     = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;

    cross = 1.0/cutEnergy - 1.0/maxEnergy - beta2*G4Log(maxEnergy/cutEnergy)/tmax
          + 0.5*(maxEnergy - cutEnergy)/energy2;

    // radiative corrections of R. Kokoulin
    if (maxEnergy > limitKinEnergy) {

      G4double logtmax = G4Log(maxEnergy);
      G4double logtmin = G4Log(max(cutEnergy,limitKinEnergy));
      G4double logstep = logtmax - logtmin;
      G4double dcross  = 0.0;

      for (G4int ll=0; ll<8; ll++)
      {
        G4double ep = G4Exp(logtmin + xgi[ll]*logstep);
        G4double a1 = G4Log(1.0 + 2.0*ep/electron_mass_c2);
        G4double a3 = G4Log(4.0*totEnergy*(totEnergy - ep)/massSquare);
        dcross += wgi[ll]*(1.0/ep - beta2/tmax + 0.5*ep/energy2)*a1*(a3 - a1);
      }

      cross += dcross*logstep*alphaprime;
    }
    cross *= charge2*twopi_mc2_rcl2/beta2; // Charge correction applied
    //cross *= twopi_mc2_rcl2/beta2;
    
  }

  //  G4cout << "tmin= " << cutEnergy << " tmax= " << tmax
  //         << " cross= " << cross << G4endl;
  
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBetheBlochModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double cross = Z*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBetheBlochModel::CrossSectionPerVolume(
                                           const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBetheBlochModel::ComputeDEDXPerVolume(const G4Material* material,
                                                  const G4ParticleDefinition* p,
                                                  G4double kineticEnergy,
                                                  G4double cut)
{
  G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tau   = kineticEnergy/mass;
  G4double cutEnergy = min(cut,tmax);
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);
  //############## Get Charge of the particle #########
  G4double charge  =  p->GetPDGCharge()/eplus;
  G4double charge2 = charge*charge;
  //#################################################
  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc2 = eexc*eexc;

  G4double eDensity = material->GetElectronDensity();

  G4double dedx = G4Log(2.0*electron_mass_c2*bg2*cutEnergy/eexc2)
                 -(1.0 + cutEnergy/tmax)*beta2;

  G4double totEnergy = kineticEnergy + mass;
  G4double del = 0.5*cutEnergy/totEnergy;
  dedx += del*del;

  // density correction
  G4double x = G4Log(bg2)/twoln10;
  dedx -= material->GetIonisation()->DensityCorrection(x);

  // shell correction
  dedx -= 2.0*corr->ShellCorrection(p,material,kineticEnergy);

  // now compute the total ionization loss

  if (dedx < 0.0) dedx = 0.0 ;

  // radiative corrections of R. Kokoulin
  if (cutEnergy > limitKinEnergy) {

    G4double logtmax = G4Log(cutEnergy);
    G4double logstep = logtmax - logLimitKinEnergy;
    G4double dloss = 0.0;
    G4double ftot2= 0.5/(totEnergy*totEnergy);

    for (G4int ll=0; ll<8; ll++)
    {
      G4double ep = G4Exp(logLimitKinEnergy + xgi[ll]*logstep);
      G4double a1 = G4Log(1.0 + 2.0*ep/electron_mass_c2);
      G4double a3 = G4Log(4.0*totEnergy*(totEnergy - ep)/massSquare);
      dloss += wgi[ll]*(1.0 - beta2*ep/tmax + ep*ep*ftot2)*a1*(a3 - a1);
    }
    dedx += dloss*logstep*alphaprime;
  }

  //dedx *= twopi_mc2_rcl2*eDensity/beta2;

  dedx *= charge2*twopi_mc2_rcl2*eDensity/beta2;// Correction for charge

  //High order corrections
  dedx += corr->HighOrderCorrections(p,material,kineticEnergy,cutEnergy);

  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCPBetheBlochModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
                                            const G4MaterialCutsCouple*,
                                            const G4DynamicParticle* dp,
                                            G4double minKinEnergy,
                                            G4double maxEnergy)
{
  G4double tmax = MaxSecondaryKinEnergy(dp);
  G4double maxKinEnergy = min(maxEnergy,tmax);
  if(minKinEnergy >= maxKinEnergy) { return; }

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double totEnergy     = kineticEnergy + mass;
  G4double etot2         = totEnergy*totEnergy;
  G4double beta2         = kineticEnergy*(kineticEnergy + 2.0*mass)/etot2;
 
  G4double grej  = 1.;
  if(tmax > limitKinEnergy) {
    G4double a0    = G4Log(2.*totEnergy/mass);
    grej  += alphaprime*a0*a0;
  }

  G4double deltaKinEnergy, f;

  // sampling follows ...
  do {
    G4double q = G4UniformRand();
    deltaKinEnergy = minKinEnergy*maxKinEnergy
                    /(minKinEnergy*(1.0 - q) + maxKinEnergy*q);


    f = 1.0 - beta2*deltaKinEnergy/tmax 
            + 0.5*deltaKinEnergy*deltaKinEnergy/etot2;

    if(deltaKinEnergy > limitKinEnergy) {
      G4double a1 = G4Log(1.0 + 2.0*deltaKinEnergy/electron_mass_c2);
      G4double a3 = G4Log(4.0*totEnergy*(totEnergy - deltaKinEnergy)/massSquare);
      f *= (1. + alphaprime*a1*(a3 - a1));
    }

    if(f > grej) {
        G4cout << "FCPBetheBlochModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for edelta= " << deltaKinEnergy
               << " tmin= " << minKinEnergy << " max= " << maxKinEnergy
               << G4endl;
    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( grej*G4UniformRand() > f );

  G4double deltaMomentum =
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double totalMomentum = totEnergy*sqrt(beta2);
  G4double cost = deltaKinEnergy * (totEnergy + electron_mass_c2) /
                                   (deltaMomentum * totalMomentum);

  G4double sint = sqrt(1.0 - cost*cost);

  G4double phi = twopi * G4UniformRand() ;

  G4ThreeVector deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  G4ThreeVector direction = dp->GetMomentumDirection();
  deltaDirection.rotateUz(direction);

  // primary change
  kineticEnergy -= deltaKinEnergy;
  G4ThreeVector dir = totalMomentum*direction - deltaMomentum*deltaDirection;
  direction = dir.unit();
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(direction);

  // create G4DynamicParticle object for delta ray
  G4DynamicParticle* delta = new G4DynamicParticle(theElectron,
                                                 deltaDirection,deltaKinEnergy);
  vdp->push_back(delta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
