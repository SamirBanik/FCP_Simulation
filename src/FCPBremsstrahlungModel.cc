// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPBremsstrahlungModel.cc                            //
//  Description: Compute Brems energy loss for FCP                    //
//                                                                    //
//  Adapted from G4MuBremsstrahlungModel by Samir Banik.              //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "FCPBremsstrahlungModel.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4Exp.hh"
#include "G4Gamma.hh"
#include "G4Log.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

const G4double FCPBremsstrahlungModel::xgi[] = 
  {0.03377,0.16940,0.38069,0.61931,0.83060,0.96623};
const G4double FCPBremsstrahlungModel::wgi[] = 
  {0.08566,0.18038,0.23396,0.23396,0.18038,0.08566};
G4double FCPBremsstrahlungModel::fDN[] = {0.0};

FCPBremsstrahlungModel::FCPBremsstrahlungModel(const G4ParticleDefinition* p,
                                                 const G4String& nam)
  : G4VEmModel(nam),
    particle(0),
    sqrte(sqrt(G4Exp(1.))),
    bh(202.4),
    bh1(446.),
    btf(183.),
    btf1(1429.),
    fParticleChange(0),
    lowestKinEnergy(1.0*GeV),
    minThreshold(50*eV)
{
  theGamma = G4Gamma::Gamma();
  nist = G4NistManager::Instance();

  lowestKinEnergy = 100.*eV;//1.*GeV;Changed for test on Sep. 25, 2019  

  mass = rmass = cc = coeff = 1.0;

  if(0.0 == fDN[1]) {
    for(G4int i=1; i<93; ++i) {
      G4double dn = 1.54*nist->GetA27(i);
      fDN[i] = dn;
      if(1 < i) {
        fDN[i] /= std::pow(dn, 1./G4double(i));
      }
    }
  }

  if(p) { SetParticle(p); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FCPBremsstrahlungModel::~FCPBremsstrahlungModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition*,
                                               const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::MinPrimaryEnergy(const G4Material*,
                                                   const G4ParticleDefinition*,
                                                   G4double cut)
{
  return std::max(lowestKinEnergy,cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCPBremsstrahlungModel::Initialise(const G4ParticleDefinition* p,
                                         const G4DataVector& cuts)
{
  if(p) { SetParticle(p); }

  // define pointer to G4ParticleChange
  if(!fParticleChange) { fParticleChange = GetParticleChangeForLoss(); }

  if(IsMaster() && p == particle && lowestKinEnergy < HighEnergyLimit()) { 
    InitialiseElementSelectors(p, cuts); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPBremsstrahlungModel::InitialiseLocal(const G4ParticleDefinition* p,
                                              G4VEmModel* masterModel)
{
  if(p == particle && lowestKinEnergy < HighEnergyLimit()) {
    SetElementSelectors(masterModel->GetElementSelectors());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::ComputeDEDXPerVolume(
                                              const G4Material* material,
                                              const G4ParticleDefinition*,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy)
{
  G4double dedx = 0.0;
  if (kineticEnergy <= lowestKinEnergy) { return dedx; }

  G4double tmax = kineticEnergy;
  G4double cut  = std::min(cutEnergy,tmax);
  //G4cout<<"#### Cut => "<<cut<<"\t"<<"Minimum Threshold => "<<minThreshold<<G4endl;
  if(cut < minThreshold) { cut = minThreshold; }
  //cut = minThreshold;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector =
    material->GetAtomicNumDensityVector();

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double loss = 
      ComputFCPBremLoss((*theElementVector)[i]->GetZ(), kineticEnergy, cut);

    dedx += loss*theAtomicNumDensityVector[i];
  }
 
  if(dedx < 0.) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::ComputFCPBremLoss(G4double Z,
                                                   G4double tkin, G4double cut)
{
  G4double totalEnergy = mass + tkin;
  static const G4double ak1 = 0.05;
  static const G4int    k2=5;
  G4double loss = 0.;

  G4double vcut = cut/totalEnergy;
  G4double vmax = tkin/totalEnergy;

  G4double aaa = 0.;
  G4double bbb = vcut;
  if(vcut>vmax) { bbb = vmax; }
  G4int kkk = (G4int)((bbb-aaa)/ak1)+k2;
  if(kkk < 1) { kkk = 1; }

  G4double hhh=(bbb-aaa)/G4double(kkk);

  G4double aa = aaa;
  for(G4int l=0; l<kkk; l++)
  {
    for(G4int i=0; i<6; i++)
    {
      G4double ep = (aa + xgi[i]*hhh)*totalEnergy;
      loss += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    aa += hhh;
  }

  loss *=hhh*totalEnergy ;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::ComputeMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double cut)
{
  G4double totalEnergy = tkin + mass;
  static const G4double ak1 = 2.3;
  static const G4int    k2  = 4;
  G4double cross = 0.;

  if(cut >= tkin) return cross;

  G4double vcut = cut/totalEnergy;
  G4double vmax = tkin/totalEnergy;

  G4double aaa = G4Log(vcut);
  G4double bbb = G4Log(vmax);
  G4int    kkk = (G4int)((bbb-aaa)/ak1)+k2 ;
  if(kkk < 1) { kkk = 1; }

  G4double hhh = (bbb-aaa)/G4double(kkk);

  G4double aa = aaa;

  for(G4int l=0; l<kkk; l++)
  {
    for(G4int i=0; i<6; i++)
    {
      G4double ep = G4Exp(aa + xgi[i]*hhh)*totalEnergy;
      cross += ep*wgi[i]*ComputeDMicroscopicCrossSection(tkin, Z, ep);
    }
    aa += hhh;
  }

  cross *=hhh;

  //G4cout << "BR e= " << tkin<< "  cross= " << cross/barn << G4endl;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::ComputeDMicroscopicCrossSection(
                                           G4double tkin,
                                           G4double Z,
                                           G4double gammaEnergy)
//  differential cross section
{
  G4double dxsection = 0.;

  if(gammaEnergy > tkin) { return dxsection; }

  G4double E = tkin + mass ;
  G4double v = gammaEnergy/E ;
  G4double delta = 0.5*mass*mass*v/(E-gammaEnergy) ;
  G4double rab0  = delta*sqrte ;

  G4int iz = G4lrint(Z);
  if(iz < 1) { iz = 1; }
  else if(iz > 92) { iz = 92; }

  G4double z13 = 1.0/nist->GetZ13(iz);
  G4double dnstar = fDN[iz];

  G4double b,b1;

  if(1 == iz) {
    b  = bh;
    b1 = bh1;
  } else {
    b  = btf;
    b1 = btf1;
  }

  // nucleus contribution logarithm
  G4double rab1=b*z13;
  G4double fn=G4Log(rab1/(dnstar*(electron_mass_c2+rab0*rab1))*
              (mass+delta*(dnstar*sqrte-2.))) ;
  if(fn <0.) { fn = 0.; }
  // electron contribution logarithm
  G4double epmax1=E/(1.+0.5*mass*rmass/E) ;
  G4double fe=0.;
  if(gammaEnergy<epmax1)
  {
    G4double rab2=b1*z13*z13 ;
    fe=G4Log(rab2*mass/((1.+delta*rmass/(electron_mass_c2*sqrte))*
                              (electron_mass_c2+rab0*rab2))) ;
    if(fe<0.) { fe=0.; }
  }

  dxsection = coeff*(1.-v*(1. - 0.75*v))*Z*(fn*Z + fe)/gammaEnergy;

  return dxsection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FCPBremsstrahlungModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition*,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double cross = 0.0;
  if (kineticEnergy <= lowestKinEnergy) return cross;
  G4double tmax = std::min(maxEnergy, kineticEnergy);
  G4double cut  = std::min(cutEnergy, kineticEnergy);
  if(cut < minThreshold) cut = minThreshold;
  if (cut >= tmax) return cross;

  cross = ComputeMicroscopicCrossSection (kineticEnergy, Z, cut);
  if(tmax < kineticEnergy) {
    cross -= ComputeMicroscopicCrossSection(kineticEnergy, Z, tmax);
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCPBremsstrahlungModel::SampleSecondaries(
                              std::vector<G4DynamicParticle*>* vdp,
                              const G4MaterialCutsCouple* couple,
                              const G4DynamicParticle* dp,
                              G4double minEnergy,
                              G4double maxEnergy)
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  // check against insufficient energy
  G4double tmax = std::min(kineticEnergy, maxEnergy);
  G4double tmin = std::min(kineticEnergy, minEnergy);
  if(tmin < minThreshold) tmin = minThreshold;
  if(tmin >= tmax) return;

  // ===== sampling of energy transfer ======

  G4ParticleMomentum partDirection = dp->GetMomentumDirection();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple,particle,kineticEnergy);
  G4double Z = anElement->GetZ();

  G4double totalEnergy   = kineticEnergy + mass;
  G4double totalMomentum = sqrt(kineticEnergy*(kineticEnergy + 2.0*mass));

  G4double func1 = tmin*
    ComputeDMicroscopicCrossSection(kineticEnergy,Z,tmin);

  G4double lnepksi, epksi;
  G4double func2;

  G4double xmin = G4Log(tmin/MeV);
  G4double xmax = G4Log(kineticEnergy/tmin);

  do {
    lnepksi = xmin + G4UniformRand()*xmax;
    epksi   = MeV*G4Exp(lnepksi);
    func2   = epksi*ComputeDMicroscopicCrossSection(kineticEnergy,Z,epksi);

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while(func2 < func1*G4UniformRand());

  G4double gEnergy = epksi;

  // ===== sample angle =====

  G4double gam  = totalEnergy/mass;
  G4double rmax = gam*std::min(1.0, totalEnergy/gEnergy - 1.0);
  G4double rmax2= rmax*rmax;
  G4double x = G4UniformRand()*rmax2/(1.0 + rmax2);

  G4double theta = sqrt(x/(1.0 - x))/gam;
  G4double sint  = sin(theta);
  G4double phi   = twopi * G4UniformRand() ;
  G4double dirx  = sint*cos(phi), diry = sint*sin(phi), dirz = cos(theta) ;

  G4ThreeVector gDirection(dirx, diry, dirz);
  gDirection.rotateUz(partDirection);

  partDirection *= totalMomentum;
  partDirection -= gEnergy*gDirection;
  partDirection = partDirection.unit();

  // primary change
  kineticEnergy -= gEnergy;
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(partDirection);

  // save secondary
  G4DynamicParticle* aGamma = 
    new G4DynamicParticle(theGamma,gDirection,gEnergy);
  vdp->push_back(aGamma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
