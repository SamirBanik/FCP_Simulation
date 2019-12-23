/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        G4FCPMinus.hh                                  /////
//      Description: charge, mass can be set as arguments           /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

#ifndef G4FCPMinus_h
#define G4FCPMinus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// ######################################################################
// ###                         FCPMINUS                               ###
// ######################################################################

class G4FCPMinus : public G4ParticleDefinition
{
private:
  static G4FCPMinus* theInstance;
  G4FCPMinus(G4double charge, G4double mass);
  ~G4FCPMinus(){;}
  
public:
  static G4FCPMinus* Definition();
  static G4FCPMinus* MakeFCPMinus(G4double charge=0.01, G4double mass=100*MeV);
  static G4FCPMinus* FCPMinusDefinition();
  static G4FCPMinus* FCPMinus();
};

#endif








