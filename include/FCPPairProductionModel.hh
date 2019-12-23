#ifndef FCPPairProductionModel_h
#define FCPPairProductionModel_h 1

#include "G4VEmModel.hh"
#include "G4NistManager.hh"
#include "G4ElementData.hh"
#include "G4Physics2DVector.hh"
#include <vector>

class G4Element;
class G4ParticleChangeForLoss;
class G4ParticleChangeForGamma;

class FCPPairProductionModel : public G4VEmModel
{
public:

  FCPPairProductionModel(const G4ParticleDefinition* p = 0,
                          const G4String& nam = "FCPPairProd");

  virtual ~FCPPairProductionModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel* masterModel);
			
  virtual G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy);
				 
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
                                const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  virtual G4double MinPrimaryEnergy(const G4Material*,
                                    const G4ParticleDefinition*,
                                    G4double);

  inline void SetLowestKineticEnergy(G4double e);

  inline void SetParticle(const G4ParticleDefinition*);

protected:

  G4double ComputMuPairLoss(G4double Z, G4double tkin, G4double cut,
                            G4double tmax);

  G4double ComputeMicroscopicCrossSection(G4double tkin,
                                          G4double Z,
                                          G4double cut);

  virtual G4double ComputeDMicroscopicCrossSection(G4double tkin,
						   G4double Z,
						   G4double pairEnergy);

  inline G4double MaxSecondaryEnergyForElement(G4double kineticEnergy,
					       G4double Z);

private:

  void MakeSamplingTables();

  void DataCorrupted(G4int Z, G4double logTkin);

  inline G4double FindScaledEnergy(G4int Z, G4double rand, G4double logTkin,
				   G4double yymin, G4double yymax); 

  // hide assignment operator
  FCPPairProductionModel & operator=(const  FCPPairProductionModel &right);
  FCPPairProductionModel(const  FCPPairProductionModel&);

protected:

  const G4ParticleDefinition* particle;
  G4NistManager*              nist;

  G4double factorForCross;
  G4double sqrte;
  G4double particleMass;
  G4double z13;
  G4double z23;
  G4double lnZ;
  G4int    currentZ;

  static const G4double xgi[8],wgi[8];

private:

  G4ParticleDefinition*       theElectron;
  G4ParticleDefinition*       thePositron;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double minPairEnergy;
  G4double lowestKinEnergy;

  G4int nzdat;

  // gamma energy bins
  G4int    nYBinPerDecade;
  size_t   nbiny;
  size_t   nbine;
  G4double ymin;
  G4double dy;
  G4double emin;
  G4double emax;

  static const G4int zdat[5];
  static const G4double adat[5];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void FCPPairProductionModel::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void FCPPairProductionModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!particle) {
    particle = p;
    particleMass = particle->GetPDGMass();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
FCPPairProductionModel::MaxSecondaryEnergyForElement(G4double kineticEnergy,
						      G4double ZZ)
{
  G4int Z = G4lrint(ZZ);
  if(Z != currentZ) {
    currentZ = Z;
    z13 = nist->GetZ13(Z);
    z23 = z13*z13;
    lnZ = nist->GetLOGZ(Z);
  }
  return kineticEnergy + particleMass*(1.0 - 0.75*sqrte*z13);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
FCPPairProductionModel::FindScaledEnergy(G4int Z, G4double rand,
					  G4double logTkin,
					  G4double yymin, G4double yymax)
{
  G4double res = yymin;
  G4Physics2DVector* pv = fElementData->GetElement2DData(Z);
  if(!pv) { 
    DataCorrupted(Z, logTkin); 
  } else {
    G4double pmin = pv->Value(yymin, logTkin);
    G4double pmax = pv->Value(yymax, logTkin);
    G4double p0   = pv->Value(0.0, logTkin);
    if(p0 <= 0.0) { DataCorrupted(Z, logTkin); }
    else { res = pv->FindLinearX((pmin + rand*(pmax - pmin))/p0, logTkin); }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
