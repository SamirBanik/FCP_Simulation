// $Id: 
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPMultipleScattering.hh                             //
//                                                                    //
//  Adapted from G4MuMultipleScattering by Samir Banik.               //
//                                                                    //
//  Author:      Samir Banik                                          //
//  Date:        14-th June 2018                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////
#ifndef FCPMultipleScattering_h
#define FCPMultipleScattering_h 1

#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class FCPMultipleScattering : public G4VMultipleScattering
{
public:    // with description

  FCPMultipleScattering(const G4String& processName="FCPmsc");

  virtual ~FCPMultipleScattering();

  // returns true for charged particles, false otherwise
  G4bool IsApplicable (const G4ParticleDefinition& p);

  // Print few lines of informations about the process: validity range,
  void PrintInfo();

protected:

  // This function initialise models
  void InitialiseProcess(const G4ParticleDefinition*);

private:        // data members

  //  hide assignment operator
  FCPMultipleScattering & operator=(const  FCPMultipleScattering &right);
  FCPMultipleScattering(const  FCPMultipleScattering&);

  G4bool   isInitialized;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
