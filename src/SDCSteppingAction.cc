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
// $Id: SDCSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file SDCSteppingAction.cc
/// \brief Implementation of the SDCSteppingAction class

#include "SDCSteppingAction.hh"

#include "G4Track.hh"
#include "G4Ions.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDCSteppingAction::SDCSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SDCSteppingAction::~SDCSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SDCSteppingAction::UserSteppingAction(const G4Step* step)
{
    
    
    G4StepPoint* PreStep = step->GetPreStepPoint();
    G4StepPoint* PostStep = step->GetPostStepPoint();
    
    G4double PreStepX = PreStep->GetPosition().x();
    G4double PreStepY = PreStep->GetPosition().y();
    G4double PreStepZ = PreStep->GetPosition().z();
    G4double parentID = step->GetTrack()->GetParentID();
    G4double trackID = step->GetTrack()->GetTrackID();
    
    G4double PostStepX = PostStep->GetPosition().x();
    G4double PostStepY = PostStep->GetPosition().y();
    G4double PostStepZ = PostStep->GetPosition().z();
    
    G4int eventNum = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
    G4double eKin = step -> GetPreStepPoint() -> GetKineticEnergy();
    G4double PosX = step->GetTrack()->GetPosition().x();
    G4double PosY = step->GetTrack()->GetPosition().y();
    G4double PosZ = step->GetTrack()->GetPosition().z();
    G4String material = step -> GetTrack() -> GetMaterial() -> GetName();
    const G4LogicalVolume* VolumeVertex = step -> GetTrack() -> GetLogicalVolumeAtVertex();
    G4String VolumeAtVertex = VolumeVertex -> GetName();
    
    G4String volume =  step->GetTrack()->GetVolume()->GetName();

    G4Track* theTrack = step->GetTrack();
    G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();
    
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4String processName = process->GetProcessName();
    
     if (volume == "DetectionPlanePhysicalVolume" && step->GetTrack()->GetDefinition()->GetParticleName() == "gamma")
     {
     std::ofstream WriteDataIn("GammaDistributionDaSteppingAction.txt", std::ios::app);
     WriteDataIn
     <<   eventNum         << '\t' << "   "
     <<   particleName     << '\t' << "   "
     <<   parentID         << '\t' << "   "
     <<   trackID         << '\t' << "   "
     <<   PosX/CLHEP::mm   << '\t' << "   "
     <<   PosY/CLHEP::mm   << '\t' << "   "
     <<   PosZ/CLHEP::mm   << '\t' << "   "
     <<   eKin/MeV           << '\t' << "   "
     <<   VolumeAtVertex         << '\t' << "   "
     <<   G4endl;
     }
    
    
    /*
    if((step->GetTrack()->GetDefinition()->GetParticleName() == "gamma"))
    {
        std::ofstream WriteDataIn("Gamma", std::ios::app);
        WriteDataIn
        <<   eventNum         << '\t' << "   "
        <<   particleName     << '\t' << "   "
        <<   parentID         << '\t' << "   "
        <<   PosX/CLHEP::mm   << '\t' << "   "
        <<   PosY/CLHEP::mm   << '\t' << "   "
        <<   PosZ/CLHEP::mm   << '\t' << "   "
        <<   eKin/CLHEP::megaelectronvolt             << '\t' << "   "
        <<   material         << '\t' << "   "
        <<   processName      << '\t' << "   "
        <<   VolumeAtVertex   << '\t' << "   "
        <<   G4endl;
        //theTrack -> SetTrackStatus(fStopAndKill);
        theTrack -> SetTrackStatus(fKillTrackAndSecondaries);

    }
     */
    
    /*
    
  //This Stepping action kills long-lived nuclei (they do not decay)
  G4String particleType = step->GetTrack()->GetDefinition()
    ->GetParticleType();
    
  if (particleType == "nucleus" && step->GetTrack()->GetParentID()>0)
    {
      G4double energy = step->GetTrack()->GetKineticEnergy();
      if (energy < 0.1*keV)
        {
          G4Ions* ion = (G4Ions*) step->GetTrack()->GetDefinition();
          G4double lifetime = ion->GetPDGLifeTime();
          G4double excitationEnergy = ion->GetExcitationEnergy();
          //stable and excited nuclei --> track them as usual
          if (lifetime < 0 || excitationEnergy > 0) return;
          if (lifetime > 1.0*microsecond) //kill long-lived nuclei
            {              
              step->GetTrack()->SetTrackStatus(fStopAndKill);
            }
          //stable nuclei are unaffected
        }
    }
     */
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
