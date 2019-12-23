////////////////////////////////////////////////////////////////////////
//  $Id:                                                              //
//  File:        DetectorConstruction.hh                              //
//  Description: Brief definition of the DetectorConstruction class   //
//                                                                    //
//  Author:      Samir Banik                                          //
//  Date:        18-th Sep 2018                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////    

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Region.hh"
#include "TString.h"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  
  G4VPhysicalVolume* GetScoringFiducial() const { return fScoringFiducial; }
  G4VPhysicalVolume* GetWorldPhysical() const { return fWorldPhysical; }
  
protected:
  G4VPhysicalVolume*  fScoringFiducial;
  G4VPhysicalVolume* fWorldPhysical;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

