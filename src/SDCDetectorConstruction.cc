//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: SDCDetectorConstruction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file SDCDetectorConstruction.cc
/// \brief Implementation of the SDCDetectorConstruction class

#include "SDCDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDCDetectorConstruction::SDCDetectorConstruction()
: G4VUserDetectorConstruction(),fCheckOverlaps(true),
worldMaterial(0),physicalWorld(0), logicalWorld(0)

{
    // Default dimensions and positions
    // World
    worldSizeX = 50 *cm;
    worldSizeY = 50 *cm;
    worldSizeZ = 50 *cm;
    
    
    //Defaults for Boro
    BoroPositionX = 0.0 *cm;//-4.0 *cm;
    BoroPositionY = 0.0 *cm;
    BoroPositionZ = 0.0 *cm;//-6.0 *cm;
    
    BoroDimensionX = 15 *cm;
    BoroExternalRadius = 15 *mm;
    BoroInternalRadius = 0.001 *mm;
    
    
    // Defaults for the First Detection Plane
    FirstDetectionPlaneSizeX = 100 *um;
    FirstDetectionPlaneSizeY = worldSizeY;
    FirstDetectionPlaneSizeZ = worldSizeZ;
    
    FirstDetectionPlanePositionX = -6.0 *cm;
    FirstDetectionPlanePositionY = 0.0 *m;
    FirstDetectionPlanePositionZ = 0.0 *m;
    
    
    // ***** Ge crystal *****
    GermaniumInternalRadius = 0.0 *mm;
    GermaniumExternalRadius = 79.0 *mm;
    GermaniumDimensionX = 71.0 *mm;
    
    GermaniumPositionX = -10.0 *cm;
    GermaniumPositionY = 0.0 *cm;
    GermaniumPositionZ = 0.0 *cm;
    
    //Colour //
    green = new G4VisAttributes( G4Colour(0, 1, 0));
    green -> SetVisibility(true);
    
    blue = new G4VisAttributes(G4Colour(0, 0, 1));
    blue -> SetVisibility(true);
    
    red = new G4VisAttributes(G4Colour(1, 0, 0));
    red -> SetVisibility(true);
    
    magenta = new G4VisAttributes(G4Colour(1, 0, 1));
    magenta -> SetVisibility(true);
    
    gray = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    gray -> SetVisibility(true);
    
    yellow = new G4VisAttributes(G4Colour(1, 1, 0));
    yellow -> SetVisibility(true);
    
    black = new G4VisAttributes(G4Colour(0, 0, 0));
    black -> SetVisibility(true);
    
    cyan = new G4VisAttributes(G4Colour(0, 1, 1));
    cyan -> SetVisibility(true);

  // **Material definition**
}



SDCDetectorConstruction::~SDCDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//////////////////////////////////////////////////////////////////////////////////

void SDCDetectorConstruction::ComputeWordGeometricalParameters()
{
}

//////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* SDCDetectorConstruction::Construct()
{
    G4cout << "SONO DENTRO Construct()" << G4endl;
    
    // Call the method where the materials are defined
    DefineMaterials();
    
    // Define volumes
    return ConstructWorld();
    
}
void SDCDetectorConstruction::DefineMaterials()
{
    //-------- NIST Materials ----------------------------------------------------
    //  Material Information imported from NIST database.
    
    G4NistManager* pNISTManager = G4NistManager::Instance();
    air  = pNISTManager -> FindOrBuildMaterial("G4_AIR");
    water  = pNISTManager -> FindOrBuildMaterial("G4_WATER");
    vacuum = pNISTManager -> FindOrBuildMaterial("G4_Galactic");
    tantalum = pNISTManager -> FindOrBuildMaterial("G4_Ta");
    aluminium = pNISTManager -> FindOrBuildMaterial("G4_Al");
    kapton = pNISTManager -> FindOrBuildMaterial("G4_KAPTON");
    PMMA = pNISTManager -> FindOrBuildMaterial("G4_PLEXIGLASS");
    glass = pNISTManager -> FindOrBuildMaterial("G4_GLASS_PLATE");
    germanium   = pNISTManager->FindOrBuildMaterial("G4_Ge");
    
    
    G4UnitDefinition::BuildUnitsTable();
    Isotopo_B10 = new G4Isotope(name="Isotopo_B10", iz=5, n=5, a=10.01*CLHEP::g/CLHEP::mole);
    Isotopo_B11 = new G4Isotope(name="Isotopo_B11", iz=5, n=6, a=11.00*CLHEP::g/CLHEP::mole);
    
    Elemento_B10  = new G4Element(name= "Elemento_B10", symbol="B10", ncomponents=1);
    Elemento_B11  = new G4Element(name= "Elemento_B11", symbol="B10", ncomponents=1);
    Elemento_B10->AddIsotope(Isotopo_B10, abundance=100*perCent);
    Elemento_B11->AddIsotope(Isotopo_B11, abundance=100*perCent);
    
    Materiale_B10 = new G4Material("Materiale_B10", density=1.15*g/cm3, ncomp=1);
    Materiale_B10 -> AddElement(Elemento_B10, 1);
    Materiale_B11 = new G4Material("Materiale_Boro11", density=1.15*g/cm3, ncomp=1);
    Materiale_B11-> AddElement(Elemento_B11, 1);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* SDCDetectorConstruction::ConstructWorld()
{  
 
    solidWorld = new G4Box("solidWorld",
                           worldSizeX/2,
                           worldSizeY/2,
                           worldSizeZ/2);
    
    logicalWorld = new G4LogicalVolume(solidWorld,
                                       vacuum,
                                       "logicalWorld",
                                       0,
                                       0,
                                       0);
    
    
    physicalWorld = new G4PVPlacement(0,
                                    G4ThreeVector(),
                                       "physicalWorld",
                                       logicalWorld,
                                       0,
                                       false,
                                       0,
                                       fCheckOverlaps);
    
    // Visible attributes of the World
    gray -> SetForceWireframe(true);
    logicalWorld -> SetVisAttributes(gray);
    

    BoroSolidVolume = new G4Tubs("BoroSolidVolume",
                                 BoroInternalRadius/2,
                                 BoroExternalRadius/2,
                                 BoroDimensionX/2,
                                 0.*deg,
                                 360.*deg);
    
    
    BoroLogicalVolume = new G4LogicalVolume(BoroSolidVolume,
                                            tantalum,
                                            "BoroLogicalVolume",
                                            0,
                                            0,
                                            0);
    
    
    BoroPhysicalVolume = new G4PVPlacement (0,
                                            G4ThreeVector(BoroPositionX, BoroPositionY,  BoroPositionZ),
                                            BoroLogicalVolume,
                                            "BoroPhysicalVolume",
                                            logicalWorld,
                                            false,
                                            0,
                                            fCheckOverlaps);
    
    yellow-> SetForceSolid(true);
    BoroLogicalVolume -> SetVisAttributes(yellow);
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    G4double phi = 90. *deg;
    // Matrix definition for a 90 deg rotation with respect to Y axis
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    //Crystal solid
    GermaniumSolidVolume = new G4Tubs("GermaniumSolidVolume",
                                        GermaniumInternalRadius,
                                        GermaniumExternalRadius/2.,
                                        GermaniumDimensionX/2.,
                                        0,
                                        CLHEP::twopi);
    
    GermaniumLogicalVolume = new G4LogicalVolume(GermaniumSolidVolume,
                                                 germanium,
                                                 "GermaniumLogicalVolume",
                                                 0,
                                                 0,
                                                 0);
    
    
    GermaniumPhysicalVolume = new G4PVPlacement (G4Transform3D(rm,
                                                               G4ThreeVector(GermaniumPositionX, GermaniumPositionY,  GermaniumPositionZ)),
                                            GermaniumLogicalVolume,
                                            "GermaniumPhysicalVolume",
                                            logicalWorld,
                                            false,
                                            0,
                                            fCheckOverlaps);
    
     return physicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void SDCDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
 
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("GermaniumPhysicalVolume");
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("GermaniumLogicalVolume",cryst);
  
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
