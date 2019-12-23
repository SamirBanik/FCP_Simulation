#ifndef FCPBremsstrahlungModel_h
#define FCPBremsstrahlungModel_h 1
// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPBremsstrahlungModel.cc                            //
//  Description: Compute energy loss for FCP particle vs. charge.     //
//                                                                    //
//  Adapted from G4MuBremsstrahlungModel by Samir Banik.              //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"

class G4Element;
class G4NistManager;
class G4ParticleChangeForLoss;


class FCPBremsstrahlungModel : public G4VEmModel
{

public:

  FCPBremsstrahlungModel(const G4ParticleDefinition* p = 0,
                          const G4String& nam = "FCPBrem");

  virtual ~FCPBremsstrahlungModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel* masterModel);

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
				const G4MaterialCutsCouple*);
			      
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

  inline void SetLowestKineticEnergy(G4double e);

  virtual G4double MinPrimaryEnergy(const G4Material*,
                                    const G4ParticleDefinition*, G4double);

protected:

  G4double ComputFCPBremLoss(G4double Z, G4double tkin, G4double cut);
  
  G4double ComputeMicroscopicCrossSection(G4double tkin,
					  G4double Z,
					  G4double cut);

  virtual G4double ComputeDMicroscopicCrossSection(G4double tkin,
						   G4double Z,
						   G4double gammaEnergy);

  inline void SetParticle(const G4ParticleDefinition*);

private:

  // hide assignment operator
  FCPBremsstrahlungModel & operator=(const  FCPBremsstrahlungModel &right);
  FCPBremsstrahlungModel(const  FCPBremsstrahlungModel&);

protected:

  const G4ParticleDefinition* particle;
  G4NistManager* nist;
  G4double mass;
  G4double charge;// Charge Declared : Samir
  G4double rmass;
  G4double cc;
  G4double coeff;
  G4double sqrte;
  G4double bh;
  G4double bh1;
  G4double btf;
  G4double btf1;

private:

  G4ParticleDefinition*       theGamma;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double lowestKinEnergy;
  G4double minThreshold;

  static const G4double xgi[6];
  static const G4double wgi[6];

  static G4double fDN[93];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void FCPBremsstrahlungModel::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void FCPBremsstrahlungModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!particle) {
    particle = p;
    mass = particle->GetPDGMass();
    charge = particle->GetPDGCharge();///eplus;
    rmass=mass/CLHEP::electron_mass_c2 ;
    cc=CLHEP::classic_electr_radius/rmass ;
    coeff= 16.*charge*charge*CLHEP::fine_structure_const*cc*cc/3. ;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
