/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        G4FCPPlus.hh                                   /////
//      Description: charge, mass can be set as arguments           /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

#ifndef G4FCPPlus_h
#define G4FCPPlus_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

// ######################################################################
// ###                         FCPPLUS                               ###
// ######################################################################

class G4FCPPlus : public G4ParticleDefinition
{
private:
  static G4FCPPlus* theInstance;
  G4FCPPlus(G4double charge, G4double mass);
  ~G4FCPPlus(){;}
    
public:
  static G4FCPPlus* MakeFCPPlus(G4double charge=0.01, G4double mass=100*MeV);
  static G4FCPPlus* Definition();
  static G4FCPPlus* FCPPlusDefinition();
  static G4FCPPlus* FCPPlus();
};

#endif








