#include "Analysis/CalAnalyzer.h"
#include "SimCore/Event/SimParticle.h"
#include "DetDescr/HcalID.h"
#include "Framework/NtupleManager.h"
#include "Event/HcalHit.h"
#include "Event/EcalHit.h"
#include "SimCore/Event/SimCalorimeterHit.h"

namespace ldmx {

    void CalAnalyzer::configure(framework::config::Parameters& ps) {
     

        //These values are hard-coded for the v12 detector
      
         HcalStartZ_ = ps.getParameter<double>("HcalStartZ_"); //Z position of the start of the back Hcal for the v12 detector

      	//EcalStartZ_ = ps.getParameter<double>("EcalStartZ_"); 
      	//PEthresh_ = ps.getParameter<double>("PEthresh_"); //PE threshold for vetoing
  
	return;

    }


  //count events
    int i=0;
    

    void CalAnalyzer::analyze(const framework::Event& event) {
      
        //Grab the SimParticle Map
        auto particle_map{event.getMap<int, ldmx::SimParticle>("SimParticles")};
	
		    i=i+1;


        //Loop over all SimParticles
        for (auto const& it : particle_map) {
          SimParticle p = it.second;
          int pdgid = p.getPdgID();
          histograms_.fill("SimParticlePdgID", pdgid);



          /*Only select the neutron for now
          //Note other neutrons can be created in the event, currently I don't distinguish between them for now
          if(pdgid == 2112){
            std::vector<double> vert = p.getVertex();
            std::vector<double> endvert = p.getEndPoint();
            std::vector<double> mom = p.getMomentum();
            //std::vector<double> endmom = p.getEndMomentum();
            //Fill SimParticle Histograms
            histograms_.fill("SimParticleEnergy", p.getEnergy());
            histograms_.fill("SimParticleX", vert[0]);
            histograms_.fill("SimParticleY", vert[1]);
            histograms_.fill("SimParticleZ", vert[2]);
            histograms_.fill("SimParticleEndX", endvert[0]);
            histograms_.fill("SimParticleEndY", endvert[1]);
            histograms_.fill("SimParticleEndZ", endvert[2]);
            histograms_.fill("SimParticlePx", mom[0]);
            histograms_.fill("SimParticlePy", mom[1]);
            histograms_.fill("SimParticlePz", mom[2]);
            //histograms_.fill("SimParticleEndPx", endmom[0]);
            //histograms_.fill("SimParticleEndPy", endmom[1]);
            //histograms_.fill("SimParticleEndPz", endmom[2]);
	    }*/



        }
	


        //Grab the Hcal Sim Hits
	   auto hcalSimHits{event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits")};
        double hcal_simesum = 0.;
        for (auto hit : hcalSimHits ) {
            //Grab the hit information
	    int pdgid = hit.getID();
            float energy = hit.getEdep();
            float x = hit.getPosition()[0];
            float y = hit.getPosition()[1];
            float z = hit.getPosition()[2];

            //Fill the Ecal Sim Hits Histograms
	    histograms_.fill("hcalsimenergy", energy);
	    // histograms_.fill("hcalsimx", x);
	    // histograms_.fill("hcalsimy", y);
            histograms_.fill("hcalsimz", z);
	    histograms_.fill("hcalsimx_hcalsimy", x, y);
            histograms_.fill("hcalsimenergy_hcalsimz", z, energy);

            hcal_simesum = hcal_simesum + energy;
          }
	  histograms_.fill("hcalsimenergysum", hcal_simesum);




	  
        //Grab the Ecal Sim Hits
         auto ecalSimHits{event.getCollection<ldmx::SimCalorimeterHit>("EcalSimHits")};
        double ecal_simesum = 0.;
	double z1 =0.;
        for (auto hit : ecalSimHits ) {
            //Grab the hit information
	    int pdgid = hit.getID();
            float energy = hit.getEdep();
            float x = hit.getPosition()[0];
            float y = hit.getPosition()[1];
            float z = hit.getPosition()[2];

            //Fill the Ecal Sim Hits Histograms
	    histograms_.fill("ecalsimenergy", energy);
	    // histograms_.fill("ecalsimx", x);
	    // histograms_.fill("ecalsimy", y);
            histograms_.fill("ecalsimz", z);
	    histograms_.fill("ecalsimx_ecalsimy", x, y);
	    histograms_.fill("ecalsimenergy_ecalsimz", z, energy);

            ecal_simesum = ecal_simesum + energy;
	    z1=z;

	    if (energy > 1){
	    std::cout<<energy<<"   ";
	    }

          }
	  histograms_.fill("ecalsimenergysum", ecal_simesum);
	  //	  histograms_.fill("ecalsimenergysum_ecalsimz", ecal_simesum, z1); 

	  	  

        //Grab the Hcal Reconstructed Hits
        std::vector<HcalHit> hcalHits = event.getCollection<HcalHit>("HcalRecHits");
	float hminZ = 9999; //Minimum z is the minimum amount of material required to veto event
        double hcal_esum = 0.;
        for (const HcalHit &hit : hcalHits ) {
            //Grab the hit information
            HcalID detID(hit.getID());
            int section = detID.getSection(); //Front, back, side, etc.
            int strip = detID.getStrip();
            int layer = detID.getLayerID();
            int PE = hit.getPE();
            int minPE = hit.getMinPE();
            float energy = hit.getEnergy();
            float time = hit.getTime();
            float amplitude = hit.getAmplitude();
            float x = hit.getXPos();
            float y = hit.getYPos();
            float z = hit.getZPos();

            //Fill the Hcal Reconstructed Hits Histograms
            //histograms_.fill("section", section);
            //histograms_.fill("layer", layer);
            //histograms_.fill("strip", strip);
            //histograms_.fill("hcalamplitude", amplitude);
            histograms_.fill("hcalenergy", energy);
            //histograms_.fill("hcaltime", time);
            //histograms_.fill("hcalx", x);
            //histograms_.fill("hcaly", y);
            histograms_.fill("hcalz", z);
            //histograms_.fill("PE", PE);
	    histograms_.fill("hcalenergy_hcalz", z, energy);
	    histograms_.fill("hcalx_hcaly", x, y);
            hcal_esum = hcal_esum + energy;
	    //   if(z < hminZ && PE >= PEthresh_ && z >= HcalStartZ_){
	    //   hminZ = z; }
          }
          histograms_.fill("hcalenergysum", hcal_esum);
	  // histograms_.fill("hcalminZ", hminZ);

	  //	  histograms_.fill("hcalenergysum_hcalz", hcal_esum, z); 




	



                //Grab the Ecal Reconstructed Hits
	  std::vector<EcalHit> ecalHits = event.getCollection<EcalHit>("EcalRecHits");




	  float eminZ = 9999; //Minimum z is the minimum amount of material required to veto event
          double ecal_esum = 0.;
          for (const EcalHit &hit : ecalHits ) {
              //Grab the hit information

              float energy = hit.getEnergy();
              float time = hit.getTime();
              float amplitude = hit.getAmplitude();
              float x = hit.getXPos();
              float y = hit.getYPos();
              float z = hit.getZPos();

              //Fill the Ecal Reconstructed Hits Histograms
              //histograms_.fill("ecalamplitude", amplitude);
              histograms_.fill("ecalenergy", energy);
              //histograms_.fill("ecaltime", time);
              //histograms_.fill("ecalx", x);
              //histograms_.fill("ecaly", y);
              histograms_.fill("ecalz", z);
	      histograms_.fill("ecalenergy_ecalz", z, energy);
	      histograms_.fill("ecalx_ecaly", x, y);

	      if (energy > 1){
		std::cout<<energy<<" | ";
		}
	      



              ecal_esum = ecal_esum + energy;
	      //  if(z < eminZ && PE >= PEthresh_ && z >= HcalStartZ_){
	      // eminZ = z; }

            }
	    histograms_.fill("ecalenergysum", ecal_esum); 
	    //  histograms_.fill("ecalminZ", eminZ);

	    //	    std::cout<<ecal_esum<<"   ";



	    
	  
  
	    //   std::vector<double> z_arr;
       
        double min_z = 0;
        double max_z = 1000;
        double n_bin = 100;






	



        for (int i = 0; i < n_bin; i++){
            int z = min_z + i*(max_z - min_z) / (n_bin - 1);
	  	double ecal_esum_z = 0.;


            for (const EcalHit &hit : ecalHits ) {
                float energy = hit.getEnergy();
                float z_pos = hit.getZPos();

                if(z_pos <= z){
                    ecal_esum_z = ecal_esum_z + energy;
              }

            }
	     histograms_.fill("ecalenergysum_ecalz_upstream", z, ecal_esum_z);

	}

 

        for (int i = 0; i < n_bin; i++){
            int z = min_z + i*(max_z - min_z) / (n_bin - 1);
	  	double ecal_esum_z = 0.;

            for (const EcalHit &hit : ecalHits ) {
                float energy = hit.getEnergy();
                float z_pos = hit.getZPos();

                if(z <= z_pos){
                    ecal_esum_z = ecal_esum_z + energy;
              }

            }
	     histograms_.fill("ecalenergysum_ecalz_downstream", z, ecal_esum_z);

	   
	}

	



	
	




        double max_z2 = 5000;



	 for (int i = 0; i < n_bin; i++){
            int z = min_z + i*(max_z2 - min_z) / (n_bin - 1);
	  	double hcal_esum_z = 0.;

            for (const HcalHit &hit : hcalHits ) {
                float energy = hit.getEnergy();
                float z_pos = hit.getZPos();

                if(z_pos <= z){
                    hcal_esum_z = hcal_esum_z + energy;
              }

            }
	     histograms_.fill("hcalenergysum_hcalz_upstream", z, hcal_esum_z);

	   
	}


	  for (int i = 0; i < n_bin; i++){
            int z = min_z + i*(max_z2 - min_z) / (n_bin - 1);
	  	double hcal_esum_z = 0.;

            for (const HcalHit &hit : hcalHits ) {
                float energy = hit.getEnergy();
                float z_pos = hit.getZPos();

                if(z <= z_pos){
                    hcal_esum_z = hcal_esum_z + energy;
              }

            }
	    histograms_.fill("hcalenergysum_hcalz_downstream", z, hcal_esum_z);

	   
	}

