////////////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------------
//      File:        DetectorConstruction.cc                               /////
//      Description: An implementation of the DetectorConstruction class   /////
//      21 December 2019, Samir Banik                                      /////
// *****************************************************************************
// ----------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringFiducial(0),
  fWorldPhysical(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   //
  G4bool checkOverlaps = true;
  G4double DetThickness = 2.5*cm;
  G4double DetDiameter  = 7.6*cm;
  G4double world_sizeXY = 1.5*DetDiameter;
  G4double world_sizeZ  = 1.5*DetDiameter;

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  //%%%%%%%%%%%% World %%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Ge");

  //############################################################################
  //#########              The detector Geometry                    ############
  //############################################################################
  G4VSolid* solidDetector
    = new G4Tubs("Detector",            // its name
		 0.*cm, 0.5*DetDiameter, 0.5*DetThickness, 0.*deg, 360*deg); // Its radius and thickness
  
  G4LogicalVolume* logicDetector
    = new G4LogicalVolume(
			  solidDetector,        // its solid
			  shape2_mat, // its material
			  "Detector");          // its name
  
  G4VPhysicalVolume* physicalDetector
    = new G4PVPlacement(
			0,                // no rotation
			G4ThreeVector(0*cm,0*cm,0*cm), // its position
			logicDetector,       // its logical volume
			"Detector",           // its name
			logicWorld,          // its mother  volume
			false,            // no boolean operation
			0,                // copy number
			checkOverlaps);  // checking overlaps
  G4VisAttributes * DetectorVisAtt = new G4VisAttributes(cyan);
  DetectorVisAtt->SetForceWireframe(false);
  logicDetector->SetVisAttributes(DetectorVisAtt);


  //############################################################################
  //############ Returning the scoring volume
  fScoringFiducial = physicalDetector;
  fWorldPhysical = physWorld;
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
