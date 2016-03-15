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
    G4Isotope* Isotopo_B11;
    G4Isotope* Isotopo_B10;
    G4Element* Elemento_B10;
    G4Element* Elemento_B11;
    G4Material* Materiale_B10;
    G4Material* Materiale_B11;

    G4bool  fCheckOverlaps;
    
    
    G4Material*         worldMaterial;
    
    G4double    worldSizeX, worldSizeY, worldSizeZ;
    
    G4double FirstDetectionPlaneSizeX, FirstDetectionPlaneSizeY, FirstDetectionPlaneSizeZ;
    G4double FirstDetectionPlanePositionX, FirstDetectionPlanePositionY, FirstDetectionPlanePositionZ;
    
    // Definition of the position and dimension
    G4double BoroPositionX, BoroPositionY, BoroPositionZ;
    G4double BoroDimensionX, BoroDimensionY, BoroDimensionZ;
    
    
    // Definition of the parameters
    G4double crystalDiameter, crystalHeight;
    G4double holeDiameter, holeHeight;
    G4double innerRadius;
    G4double sampleRadius;
    
    // Definition of the 'solid volumes' variables
    G4Box*  solidWorld;
    G4Box*  FirstSolidDetectionPlane;

    G4Sphere* sampleSolid;
    G4Tubs* BoroSolidVolume;
    G4Tubs* Germanium_SolidVolume1;
    G4Tubs* Germanium_SolidVolume2;
    G4SubtractionSolid* Germanium_SolidVolume;

    
    // Definition of the 'physical volumes' variables
    G4VPhysicalVolume*  physicalWorld;
    G4VPhysicalVolume*  FirstPhysicalDetectionPlane;
    
    //
    G4VPhysicalVolume* BoroPhysicalVolume;
    G4VPhysicalVolume* Germanium_PhysicalVolume;
    //
    //

    
    
    // Definition of the 'logical volumes' variables
    G4LogicalVolume*    logicalWorld;
    G4LogicalVolume*    FirstLogicalDetectionPlane;
    G4LogicalVolume*    BoroLogicalVolume;
    G4LogicalVolume*    Germanium_LogicalVolume;
    
    
    G4MultiFunctionalDetector* Germanium_Detector;
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
    
    G4RotationMatrix rotationMatrix;
    G4Transform3D BoroTransform;
    G4ThreeVector BoroPosition;
    
    void DefineMaterials();
   void ComputeWordGeometricalParameters();
    
    G4VPhysicalVolume* ConstructWorld();
    
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

