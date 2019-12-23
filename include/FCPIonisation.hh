#ifndef FCPIonisation_h
#define FCPIonisation_h 1
// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPIonisation.cc                                     //
//  Description: Ionisation and delta-ray production for FCP particle //
//                                                                    //
//  Adapted from G4MuIonisation by Samir Banik.                       //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4VEnergyLossProcess.hh"

class G4Material;
class G4ParticleDefinition;


class FCPIonisation : public G4VEnergyLossProcess
{

public:

  FCPIonisation(const G4String& name = "FCPIoni");

  virtual ~FCPIonisation();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
				    const G4Material*, G4double cut);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*);

private:

  // hide assignment operator
  FCPIonisation & operator=(const FCPIonisation &right);
  FCPIonisation(const FCPIonisation&);

  G4double    mass;
  G4double    ratio;

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;
  G4bool                      isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
