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
// $Id: SDCDetectorConstruction.hh 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file SDCDetectorConstruction.hh
/// @brief Definition of the SDCDetectorConstruction class (Mandatory)

#ifndef SDCDetectorConstruction_h
#define SDCDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4RotationMatrix.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4SubtractionSolid.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"


class G4LogicalVolume;
class G4Material;
class G4UImanager;

/// Detector construction class to define materials (with their physical properties) and detector geometry.
///

class SDCDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    SDCDetectorConstruction();
    ~SDCDetectorConstruction();

  public:

    virtual G4VPhysicalVolume* Construct();
    void SetWorldMaterial(G4String);
    virtual void ConstructSDandField();
public:
    
    const
    G4VPhysicalVolume* GetWorld()      {return physicalWorld;};
    G4double          GetSize()       {return worldSizeX;};
    G4Material*       GetMaterial()   {return worldMaterial;};
    
    
      void               PrintWorldParameters();
    
    
  private:
    
    // Materials
    
    G4String name, symbol;             // a=mass of a mole;
    G4double a, z, density;            // z=mean number of protons;
    G4int iz, n, ncomp;                       //iz=nb of protons  in an isotope;
    
    G4int ncomponents, natoms;
    G4double abundance, fractionmass;
    G4double temperature, pressure;
    

    G4Material*         water;
    G4Material*         air;
    G4Material*         vacuum;
    G4Material*         tantalum;
    G4Material*         aluminium;
    G4Material*         kapton;
    G4Material*         PMMA;
    G4Material*         glass;
    G4Material*         germanium;
    G4Material*         lead;
    
    G4Isotope* Isotopo_B11;
    G4Isotope* Isotopo_B10;
    G4Element* Elemento_B10;
    G4Element* Elemento_B11;
    
    G4Material* Materiale_B10;
    G4Material* Materiale_B11;
    G4Material* DetectionPlaneMaterial;

    G4bool  fCheckOverlaps;
    
    
    G4Material*         worldMaterial;
    
    G4double    worldSizeX, worldSizeY, worldSizeZ;
    
    // Definition of the position and dimension

    G4double DetectionPlaneSizeX, DetectionPlaneSizeY, DetectionPlaneSizeZ;
    G4double DetectionPlanePositionX, DetectionPlanePositionY, DetectionPlanePositionZ;
    G4double FrontShieldSizeX, FrontShieldSizeY, FrontShieldSizeZ;
    G4double FrontShieldPositionX, FrontShieldPositionY, FrontShieldPositionZ;
    
    G4double FrontShieldSlitSizeX, FrontShieldSlitSizeY, FrontShieldSlitSizeZ;
    G4double FrontShieldSlitPositionX, FrontShieldSlitPositionY, FrontShieldSlitPositionZ;
    
    G4double BoroPositionX, BoroPositionY, BoroPositionZ;
    G4double BoroDimensionX, BoroInternalRadius, BoroExternalRadius;
    
    G4double GermaniumInternalRadius;
    G4double GermaniumExternalRadius;
    G4double GermaniumDimensionX;
    
    G4double GermaniumPositionX;
    G4double GermaniumPositionY;
    G4double GermaniumPositionZ;
    
    
    // Definition of the 'solid volumes' variables
    G4Box*  solidWorld;
    G4Box*  DetectionPlaneSolidVolume;
    G4Box*  FrontShieldSolidVolume;
G4Box*  FrontShieldSlitSolidVolume;
    
    G4Sphere* sampleSolid;
    G4Tubs* BoroSolidVolume;
    
    G4Tubs* GermaniumSolidVolume;
    
    
    // Definition of the 'physical volumes' variables
    G4VPhysicalVolume*  physicalWorld;
    G4VPhysicalVolume*  DetectionPlanePhysicalVolume;
    G4VPhysicalVolume*  FrontShieldPhysicalVolume;
    G4VPhysicalVolume*  FrontShieldSlitPhysicalVolume;

    //
    G4VPhysicalVolume* BoroPhysicalVolume;
    G4LogicalVolume * GermaniumLogicalVolume;
    G4VPhysicalVolume *GermaniumPhysicalVolume;
    
    // Definition of the 'logical volumes' variables
    G4LogicalVolume*    logicalWorld;
    G4LogicalVolume*    FirstLogicalDetectionPlane;
    G4LogicalVolume*    BoroLogicalVolume;
    G4LogicalVolume*    DetectionPlaneLogicalVolume;
    G4LogicalVolume*    FrontShieldLogicalVolume;
    G4LogicalVolume*    FrontShieldSlitLogicalVolume;
    
    //G4MultiFunctionalDetector* Germanium_Detector;
    G4VPrimitiveScorer* PrimitiveScorer;
    
private:
    G4VisAttributes *green;
    G4VisAttributes *blue;
    G4VisAttributes *red;
    G4VisAttributes *magenta;
    G4VisAttributes *yellow;
    G4VisAttributes *cyan;
    G4VisAttributes *gray;
    G4VisAttributes *black;
    
    void DefineMaterials();
   void ComputeWordGeometricalParameters();
    
    G4VPhysicalVolume* ConstructWorld();
    
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

