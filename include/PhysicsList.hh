/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        PhysicsList.hh                                 /////
//      Description: Breif definition of the PhysicsList class      /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "FCPPhysics.hh"

class PhysListMessenger;

class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList(G4int ver = 0);
  virtual ~PhysicsList();
  virtual void ConstructParticle();
  virtual void SetCuts();
  void AddFCPPhysics(G4double chg, G4double mass)
  {
    lPhys = new FCPPhysics(chg,mass);
  }
protected:
  PhysListMessenger* messenger;
  
private:
  FCPPhysics *lPhys;
  G4int verbose;

};
#endif
