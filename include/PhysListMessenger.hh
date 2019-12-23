#ifndef PhysListMessenger_hh
#define PhysListMessenger_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        PhysListMessenger.hh                                 //
//  Description: Brief definition of the PhysListMessenger class      //
//                                                                    //
//  Date:        1 July 2015                                          //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class PhysicsList;

class PhysListMessenger : public G4UImessenger{
public:
  PhysListMessenger(PhysicsList* thePhysicsList);
  ~PhysListMessenger();

  virtual void SetNewValue(G4UIcommand* cmd, G4String newValue);

private:
  PhysicsList* thePhysicsList;

  G4UIcommand* LIPCmd;	// FCP physics, charge and mass
};

#endif	
