#ifndef FCPBetheBlochModel_h
#define FCPBetheBlochModel_h 1
// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPBetheBlochModel.hh                                //
//  Description: Compute energy loss for FCP particle vs. charge.     //
//                                                                    //
//  Adapted from G4MuBetheBlochModel by Samir Banik.                  //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"

class G4ParticleChangeForLoss;
class G4EmCorrections;

class FCPBetheBlochModel : public G4VEmModel
{

public:

  FCPBetheBlochModel(const G4ParticleDefinition* p = 0,
                      const G4String& nam = "FCPBetheBloch");

  virtual ~FCPBetheBlochModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
				const G4MaterialCutsCouple*);
			
  virtual G4double ComputeCrossSectionPerElectron(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy);
				 
  virtual G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy);
				 				 
  virtual G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
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

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				      G4double kinEnergy);

private:

  inline void SetParticle(const G4ParticleDefinition* p);

  // hide assignment operator
  FCPBetheBlochModel & operator=(const  FCPBetheBlochModel &right);
  FCPBetheBlochModel(const  FCPBetheBlochModel&);

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theElectron;
  G4ParticleChangeForLoss*    fParticleChange;
  G4EmCorrections*            corr;

  G4double limitKinEnergy;
  G4double logLimitKinEnergy;
  G4double mass;
  G4double massSquare;
  G4double ratio;
  G4double twoln10;
  G4double bg2lim;
  G4double taulim;
  G4double alphaprime;
  static G4double xgi[8],wgi[8];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void FCPBetheBlochModel::SetParticle(const G4ParticleDefinition* p)
{
  if (!particle) {
    particle = p;
    mass = particle->GetPDGMass();
    massSquare = mass*mass;
    ratio = CLHEP::electron_mass_c2/mass;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
