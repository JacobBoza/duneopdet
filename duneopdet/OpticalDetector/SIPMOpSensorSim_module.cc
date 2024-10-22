//=========================================================
// SIPMOpSensorSim_module.cc
// This module produces detected photons (creating OpDetDivRec)
// from photon detectors taking SimPhotons as input.
// Applies quantum efficiency and cross-talk
//
// Gleb Sinev, Duke, 2015
// Anne Christensen, CSU, 2019
// Alex Himmel, FNAL, 2021
//=========================================================

#ifndef SIPMOpSensorSim_h
#define SIPMOpSensorSim_h

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h" //vitor
#include "art_root_io/TFileDirectory.h"



// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes

#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larana/OpticalDetector/OpDetResponseInterface.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "duneopdet/OpticalDetector/AlgoSSPLeadingEdge.h"
#include "dunecore/DuneObj/OpDetDivRec.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"


// CLHEP includes

#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"

// C++ includes

#include <vector>
#include <map>
#include <cmath>
#include <memory>

// ROOT includes

#include "TTree.h"

namespace opdet {

  class SIPMOpSensorSim : public art::EDProducer{

  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>  InputTag            { Name("InputTag"),            Comment("Input tag for SimPhotons") };
      fhicl::Atom<double>         QuantumEfficiency   { Name("QuantumEfficiency"),   Comment("Probabilityof recording a photon") };
      fhicl::Atom<double>         DarkNoiseRate       { Name("DarkNoiseRate"),       Comment("Rate in Hz") };
      fhicl::Atom<double>         CrossTalk           { Name("CrossTalk"),           Comment("Cross talk (1->2 PE) probability") };
      fhicl::Atom<double>         Correction          { Name("Correction"),          Comment("Adjust the amount of total light. Kept seprate from QE for clarity."), 1};
      fhicl::OptionalAtom<double> LateLightCorrection { Name("LateLightCorrection"), Comment("Adjust the amount of late light")};
      fhicl::OptionalAtom<double> LateLightBoundary   { Name("LateLightBoundary"),   Comment("Boundary that defines late light (ns)") };

    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit SIPMOpSensorSim(Parameters const & config);
    void produce(art::Event&) override;

  private:

    // The parameters read from the FHiCL file
    art::InputTag fInputTag;            // Input tag for OpDet collection
    std::vector < std::string > fInputModule;
    art::ProductToken< std::vector<sim::SimPhotons> > fInputToken;
    art::Handle< std::vector< sim::SimPhotons > > PhotonHandle;




    double        fQE;
    double        fDarkNoiseRate;	      // In Hz
    double        fCrossTalk;           // Probability of SiPM producing 2 PE signal

    double        fTimeBegin;           // Earliest and latest possible times for reading out data
    double        fTimeEnd;             // Used for defining the dark noise range

    bool          fCorrectLateLight;    // Do we apply a late light correction?
    double        fLateLightCorrection; // How much to correct the late light
    double        fLateLightBoundary;   // What is the boundary which defines late light?

    // Random number engines
    CLHEP::HepRandomEngine& fSIPMEngine;
    CLHEP::RandExponential  fRandExponential;
    CLHEP::RandFlat         fRandFlat;
    CLHEP::RandPoissonQ     fRandPoissPhot;
  
    // Produce waveform on one of the optical detectors
    // These cannot be const due to the random number generators
    void PhotonsToPE(sim::SimPhotons const& simPh, sim::OpDetDivRec& dr_plusnoise, art::Event& evt);
    void AddDarkNoise(sim::OpDetDivRec &);
    unsigned short CrossTalk();

    void SetBeginEndTimes();

  };
}

#endif

namespace opdet {

  DEFINE_ART_MODULE(SIPMOpSensorSim)

}

namespace opdet {

  //---------------------------------------------------------------------------
  // Constructor
  SIPMOpSensorSim::SIPMOpSensorSim(Parameters const & config)
    : art::EDProducer{config}
    , fInputTag{     config().InputTag()}
    , 
fInputToken{   consumes< std::vector<sim::SimPhotons> >(fInputTag) } //compiles with this one
    , fDarkNoiseRate{config().DarkNoiseRate()}
    , fCrossTalk{    config().CrossTalk()}
    , fSIPMEngine(
        art::ServiceHandle<rndm::NuRandomService>()->registerAndSeedEngine(
          createEngine(0, "HepJamesRandom", "sipm"),
          "HepJamesRandom",
          "sipm", 
          config.get_PSet(), 
          "SeedSiPM"))
    , fRandExponential(fSIPMEngine)
    , fRandFlat(fSIPMEngine)
    , fRandPoissPhot(fSIPMEngine)
  {

    // This module produces (infrastructure piece)
    produces< std::vector< sim::OpDetDivRec > >();


    if (fDarkNoiseRate < 0.0) {
      throw art::Exception(art::errors::Configuration) 
        << "fDarkNoiseRate: " << fDarkNoiseRate << '\n'
        << "Dark noise rate should be non-negative!\n";
    }


    double tempQE = config().QuantumEfficiency() * config().Correction();
    if (tempQE < 0 || tempQE > 1) {
      throw art::Exception(art::errors::Configuration) 
        << "QuantumEfficiency in SIPMSensorSim: " << config().QuantumEfficiency() << "\n"
        << "                   with correction: " << config().Correction() << "\n"
        << "       leading to total efficiency: " << tempQE << "\n"
        << "It should be between 0 and 1!\n";
    }

    // Correct out the prescaling applied during simulation
    auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
    fQE = tempQE / LarProp->ScintPreScale();

    if (fQE > 1.0001 ) {
      throw art::Exception(art::errors::Configuration)
        << "QuantumEfficiency in SIPMSensorSim: " << config().QuantumEfficiency() << "\n"
        << "                   with correction: " << config().Correction() << "\n"
        << "       leading to total efficiency: " << tempQE << "\n"
        << "is too large.\n"
        << "It is larger than the prescaling applied during simulation, " << LarProp->ScintPreScale() << ".\n"
        << "Final QE must be equal to or smaller than the QE applied at simulation time.\n";
    }

    //This will be handled on the WaveformDigitizerSim_module.cc 
    // Check for non-trivial channel mapping which is not supported

    // Set time ranges if needed for dark noise
    if (fDarkNoiseRate > 0) SetBeginEndTimes();

    // Check for optional parameters for adjusting late light.
    // Apply adjustments if both are given.
    // Throw an error if only one is specified.
    bool haveCorrection = config().LateLightCorrection(fLateLightCorrection);
    bool haveBoundary   = config().LateLightBoundary(fLateLightBoundary);
    fCorrectLateLight = haveBoundary && haveCorrection;

    if (haveBoundary != haveCorrection) {
      throw art::Exception(art::errors::Configuration)
        << "Must specify both LateLightCorrection and LateLightBoundary in order to apply correction.\n"
        << "Leave both out to not adjust late light.\n";
    }

    if (fCorrectLateLight && fLateLightCorrection*tempQE > LarProp->ScintPreScale())
    {
      throw art::Exception(art::errors::Configuration)
        << "QuantumEfficiency in SIPMSensorSim: " << config().QuantumEfficiency() << "\n"
        << "                   with correction: " << config().Correction() << "\n"
        << "         and late light correction: " << fLateLightCorrection << "\n"
        << "       leading to total efficiency: " << tempQE * fLateLightCorrection << "\n"
        << "is too large.\n"
        << "It is larger than the prescaling applied during simulation, " << LarProp->ScintPreScale() << ".\n"
        << "Final QE must be equal to or smaller than the QE applied at simulation time.\n";
    }
  }


  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::produce(art::Event& event)
  {
    // A pointer that will store produced OpDetDivRec
    auto OpDetDivRecPtr = std::make_unique< std::vector< sim::OpDetDivRec > >();

    // Get OpDetBacktrackerRecord from the event
    auto const & simph_handle = event.getValidHandle(fInputToken);

    // For every optical detector:
    for (auto const& simph : *simph_handle) {
	int opDet = simph.OpChannel();
      auto DivRecPlusNoise = sim::OpDetDivRec(opDet);

      PhotonsToPE(simph, DivRecPlusNoise, event);
    // sim::OnePhoton photon;
    // double time = photon.Time;
      // Generate dark noise
      if (fDarkNoiseRate > 0.0) AddDarkNoise(DivRecPlusNoise);

      // AddAfterPulsing();

      OpDetDivRecPtr->emplace_back(DivRecPlusNoise);
    }

    event.put(std::move(OpDetDivRecPtr));
  }

  //---------------------------------------------------------------------------


 void SIPMOpSensorSim::PhotonsToPE(sim::SimPhotons const& simph, sim::OpDetDivRec& dr_plusnoise, art::Event& evt)
  {

    // Don't do anything without any records
    auto photon_handles = evt.getMany<std::vector<sim::SimPhotons>>();   
    art::ServiceHandle<opdet::OpDetResponseInterface const> odresponse;
    if (photon_handles.size() == 0)
        throw art::Exception(art::errors::ProductNotFound)
          << "sim SimPhotons retrieved and you requested them.";

    for (auto const& ph_handle : photon_handles) {
          // Do some checking before we proceed
          if (!ph_handle.isValid()) continue;

          if ((*ph_handle).size() > 0) {
              for (auto const& itOpDet : (*ph_handle)) {
              //Reset Counters
              int fCountOpDetDetected = 0;
              //Reset t0 for visible light
            //  fT0_vis = 999.;
      //Need to take all detected photons and add the detected photons to OpDetDivRec         
              //Get data from HitCollection entry
              int fOpChannel = itOpDet.OpChannel();
              const sim::SimPhotons& TheHit = itOpDet;
                  
          for (const sim::OnePhoton& Phot : TheHit) {
		if (odresponse->detected(fOpChannel, Phot)) ++fCountOpDetDetected;
        }
      
         int nphot = fRandPoissPhot.fire(fQE * (double)fCountOpDetDetected);
        for(int truePh=0; truePh<nphot; ++truePh) {
	  const sim::OnePhoton& Phot = TheHit[nphot];
          // Determine actual PE with Cross Talk
          unsigned int PE = 1+CrossTalk();
          for(unsigned int i = 0; i < PE; i++) {
            // Add to collection
            dr_plusnoise.AddPhoton(simph.OpChannel(), // Channel
                                   Phot.MotherTrackID,    // TrackID
                                   Phot.Time);          // Time
          }
       }
      }
     }
    }
  }

  //---------------------------------------------------------------------------
  unsigned short SIPMOpSensorSim::CrossTalk()
  {
    // Use recursion to allow cross talk to create cross talk
    if      (fCrossTalk <= 0.0)                 return 0;
    else if (fRandFlat.fire(1.0) > fCrossTalk)  return 0;
    else                                        return 1+CrossTalk();
  }

  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::AddDarkNoise(sim::OpDetDivRec & dr_plusnoise)
  {

    double tau = 1.0/fDarkNoiseRate*1.e6; // Typical time between Dark counts in us
                                          // Multiply by 10^6 since fDarkNoiseRate is in Hz

    double darkNoiseTime = fRandExponential.fire(tau) + fTimeBegin;
    while (darkNoiseTime < fTimeEnd) {
	    //Currently using all zeros as an indicator of Dark Noise for the trackID.
      int PE = 1+CrossTalk();
      for(int j = 0; j < PE; j++) {
	     
   dr_plusnoise.AddPhoton(dr_plusnoise.OpDetNum(), 0, darkNoiseTime);
      }

      // Find next time to simulate a single PE pulse
      darkNoiseTime += fRandExponential.fire(tau);
    }
  }

  //---------------------------------------------------------------------------
  void SIPMOpSensorSim::SetBeginEndTimes()
  {

    // Get the window start/end time from the detector properties services
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    auto const geometry  = art::ServiceHandle< geo::Geometry >();

    double maxDrift = 0.0;
    for (geo::TPCGeo const& tpc : geometry->Iterate<geo::TPCGeo>())
      if (maxDrift < tpc.DriftDistance()) maxDrift = tpc.DriftDistance();

    // Start at -1 drift window
    fTimeBegin = -1*maxDrift/detProp.DriftVelocity();

    // Take the TPC readout window size and convert
    // to us with the electronics clock frequency
    fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();
  }


} // end opdet namespace 
