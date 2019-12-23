// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPPairProduction.cc                                 //
//  Description: e+e- pair production by FCPs                         //
//                                                                    //
//  Adapted from G4MuPairProduction by Samir Banik.                       //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef FCPPairProduction_h
#define FCPPairProduction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4VEnergyLossProcess.hh"

class FCPPairProduction : public G4VEnergyLossProcess
{
public:

  FCPPairProduction(const G4String& processName = "FCPPairProd");

  virtual ~FCPPairProduction();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, G4double cut);

  virtual void PrintInfo();

  inline void SetLowestKineticEnergy(G4double e);

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

private:

  FCPPairProduction & operator=(const FCPPairProduction &right);
  FCPPairProduction(const FCPPairProduction&);

protected:

  const G4ParticleDefinition* theParticle;
  G4double                    lowestKinEnergy;
  G4bool                      isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void FCPPairProduction::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
