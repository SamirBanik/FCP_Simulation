#ifndef FCPBremsstrahlung_h
#define FCPBremsstrahlung_h 1
// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPBremsstrahlung.cc                                 //
//  Description: Radiation process for FCP particle.                  //
//                                                                    //
//  Adapted from G4MuBremsstrahlung by Samir Banik.                   //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VEnergyLossProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;

class FCPBremsstrahlung : public G4VEnergyLossProcess
{
public:

  FCPBremsstrahlung(const G4String& processName = "FCPBrems");

  virtual ~FCPBremsstrahlung();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, 
				    G4double cut);

  virtual void PrintInfo();

  inline void SetLowestKineticEnergy(G4double e);

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

private:

  FCPBremsstrahlung & operator=(const FCPBremsstrahlung &right);
  FCPBremsstrahlung(const FCPBremsstrahlung&);

protected:

  G4double  lowestKinEnergy;
  G4bool    isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void FCPBremsstrahlung::SetLowestKineticEnergy(G4double e) 
{
  lowestKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
