///////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------
//      File:        ActionInitialization.hh                              /////
//      Description: Breif definition of the ActionInitialization class   /////
//      21 December 2019, Samir Banik                                     /////
// ****************************************************************************
// ---------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization();
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
