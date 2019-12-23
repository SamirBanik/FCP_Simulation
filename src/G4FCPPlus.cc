///////////////////////////////////////////////////////////////////////// 
// ----------------------------------------------------------------------
//      File:        G4FCPPlus.cc                                   /////
//      Description: charge, mass can be set as arguments           /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

#include "G4FCPPlus.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                          FCPPLUS                              ###
// ######################################################################
G4FCPPlus* G4FCPPlus::theInstance = 0;

G4FCPPlus* G4FCPPlus::MakeFCPPlus(G4double charge, G4double mass)
{
  if (theInstance) return theInstance;
  const G4String name = "FCP+";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
 
  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  anInstance = new G4ParticleDefinition(
					name, mass, 0,  +fabs(charge)*eplus,
					1,               0,                0,          
					0,               0,                0,             
					"lepton",        -1,                0,        -90,
					true,             0,             NULL,
					false,           "FCP" 
					);

  // Bohr Magnetron
  G4double muB =  0.5*eplus*hbar_Planck/(anInstance->GetPDGMass()/c_squared) ;
  anInstance->SetPDGMagneticMoment( muB * 1.0011659209);
  }
  theInstance = reinterpret_cast<G4FCPPlus*>(anInstance);
  return theInstance;
}

G4FCPPlus*  G4FCPPlus::Definition()
{
  return theInstance;
}

G4FCPPlus*  G4FCPPlus::FCPPlusDefinition()
{
  return theInstance;
}

G4FCPPlus*  G4FCPPlus::FCPPlus()
{
  return theInstance;
}

