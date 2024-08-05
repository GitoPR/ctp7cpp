#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <hls_math.h>


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "cicada.h"
#include "algo_unpacked.h"
#include "UCTSummaryCard.hpp"
#include "PU_LUT.h"
#include "calo_out_coordinates.h"
#include "superregion.h"
#include "bitonicSort64.h"

const uint16_t NRegionsPerLink = 7; // Bits 8-21, 22-39, 40-55,..., 104-119, keeping ranges (7, 0) and (127, 120) unused
const uint16_t MaxRegions = N_CH_IN * NRegionsPerLink;

 /*
  * algo_unpacked interface exposes fully unpacked input and output link data.
  * This version assumes use of 10G 8b10b links, and thus providing 192bits/BX/link.
  *
  * !!! N.B.: Do NOT use the first and the last bytes of link_in/link_out (i.e. link_in/link_out[].range(7,0)
  * and link_in/link_out[].range(127, 120) as these fields are reserved for transmission of 8b10b input/output 
  * link alignment markers and CRC checksum word
  *
  * The remaining 112 bits are available for algorithm use.
  *
  * !!! N.B. 2: make sure to assign every bit of link_out[] data. First byte should be assigned zero.
  */
typedef struct {

  ap_uint<5> iphi;
  ap_uint<5> ieta;
  
 } my_region_t;

const my_region_t region_lut[252] =

  {
    {0,0}, 
    {0,1}, 
    {0,2},
    {0,3},
    {0,4},
    {0,5},
    {0,6},
    {0,7},
    {0,8},
    {0,9},
    {0,10},
    {0,11},
    {0,12},
    {0,13},
    {0,0},
    {1,1},
    {1,2},
    {1,3},
    {1,4},
    {1,5},
    {1,6},
    {1,7},
    {1,8},
    {1,9},
    {1,10},
    {1,11},
    {1,12},
    {1,13},
    {2,0},
    {2,1},
    {2,2},
    {2,3},
    {2,4},
    {2,5},
    {2,6},
    {2,7},
    {2,8},
    {2,9},
    {2,10},
    {2,11},
    {2,12},
    {2,13},
    {3,0},
    {3,1},
    {3,2},
    {3,3},
    {3,4},
    {3,5},
    {3,6},
    {3,7},
    {3,8},
    {3,9},
    {3,10},
    {3,11},
    {3,12},
    {3,13},
    {4,0},
    {4,1},
    {4,2},
    {4,3},
    {4,4},
    {4,5},
    {4,6},
    {4,7},
    {4,8},
    {4,9},
    {4,10},
    {4,11},
    {4,12},
    {4,13},
    {5,0},
    {5,1},
    {5,2},
    {5,3},
    {5,4},
    {5,5},
    {5,6},
    {5,7},
    {5,8},
    {5,9},
    {5,10},
    {5,11},
    {5,12},
    {5,13},
    {6,0},
    {6,1},
    {6,2},
    {6,3},
    {6,4},
    {6,5},
    {6,6},
    {6,7},
    {6,8},
    {6,9},
    {6,10},
    {6,11},
    {6,12},
    {6,13},
    {7,0},
    {7,1},
    {7,2},
    {7,3},
    {7,4},
    {7,5},
    {7,6},
    {7,7},
    {7,8},
    {7,9},
    {7,10},
    {7,11},
    {7,12},
    {7,13},
    {8,0},
    {8,1},
    {8,2},
    {8,3},
    {8,4},
    {8,5},
    {8,6},
    {8,7},
    {8,8},
    {8,9},
    {8,10},
    {8,11},
    {8,12},
    {8,13},
    {9,0},
    {9,1},
    {9,2},
    {9,3},
    {9,4},
    {9,5},
    {9,6},
    {9,7},
    {9,8},
    {9,9},
    {9,10},
    {9,11},
    {9,12},
    {9,13},
    {10,0},
    {10,1},
    {10,2},
    {10,3},
    {10,4},
    {10,5},
    {10,6},
    {10,7},
    {10,8},
    {10,9},
    {10,10},
    {10,11},
    {10,12},
    {10,13},
    {11,0},
    {11,1},
    {11,2},
    {11,3},
    {11,4},
    {11,5},
    {11,6},
    {11,7},
    {11,8},
    {11,9},
    {11,10},
    {11,11},
    {11,12},
    {11,13},
    {12,0},
    {12,1},
    {12,2},
    {12,3},
    {12,4},
    {12,5},
    {12,6},
    {12,7},
    {12,8},
    {12,9},
    {12,10},
    {12,11},
    {12,12},
    {12,13},
    {13,0},
    {13,1},
    {13,2},
    {13,3},
    {13,4},
    {13,5},
    {13,6},
    {13,7},
    {13,8},
    {13,9},
    {13,10},
    {13,11},
    {13,12},
    {13,13},
    {14,0},
    {14,1},
    {14,2},
    {14,3},
    {14,4},
    {14,5},
    {14,6},
    {14,7},
    {14,8},
    {14,9},
    {14,10},
    {14,11},
    {14,12},
    {14,13},
    {15,0},
    {15,1},
    {15,2},
    {15,3},
    {15,4},
    {15,5},
    {15,6},
    {15,7},
    {15,8},
    {15,9},
    {15,10},
    {15,11},
    {15,12},
    {15,13},
    {16,0},
    {16,1},
    {16,2},
    {16,3},
    {16,4},
    {16,5},
    {16,6},
    {16,7},
    {16,8},
    {16,9},
    {16,10},
    {16,11},
    {16,12},
    {16,13},
    {17,0},
    {17,1},
    {17,2},
    {17,3},
    {17,4},
    {17,5},
    {17,6},
    {17,7},
    {17,8},
    {17,9},
    {17,10},
    {17,11},
    {17,12},
    {17,13},
  };

// for converting from ieta iphi to iregion index
class iregion_lut {
public:
  constexpr static int rows = 14;
  constexpr static int cols = 18;
  ap_uint<32> iRegion_lut[rows][cols];

  iregion_lut() {
    int count = 1;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        iRegion_lut[i][j] = count++;
      }
    }
  }
};



typedef ap_uint<8> loop; //loop type(i guess)
//Jet class

class jetInfo{
public:
  ap_uint<12> seedEnergy;
  ap_uint<12> energy;
  ap_uint<5> phiMax;
  ap_uint<5> etaMax;

  jetInfo(){
    seedEnergy = 0;
    energy = 0;
    phiMax = 0;
    etaMax = 0;
  }

  jetInfo& operator=(const jetInfo& rhs){
    seedEnergy = rhs.seedEnergy;
    energy = rhs.energy;
    phiMax = rhs.phiMax;
    etaMax = rhs.etaMax;
    return *this;
  }
};


//setting up max et region finder

static constexpr int nSTPhi = 18;
static constexpr int nSTEta = 14;

namespace gctobj {

  class towerMax {
  public:
    ap_uint<12> energy;
    ap_uint<7> phi;
    ap_uint<7> eta;
    ap_uint<5> ieta;
    ap_uint<5> iphi;
    ap_uint<7> towerEta;
    ap_uint<7> towerPhi;

    towerMax() {
      energy = 0;
      phi = 0;
      eta = 0;
      ieta = 0; 
      iphi= 0; 
      towerEta = 0; 
      towerPhi = 0; 
    }
  };

  typedef struct {
    int et;
    int eta;
    int phi;
    int towerEta;
    int towerPhi;
    int ieta;
    int iphi; 
  } GCTsupertower_t;


  typedef struct {
    GCTsupertower_t cr[nSTPhi];
  } etaStrip_t;

  typedef struct {
    GCTsupertower_t pk[nSTEta];
  } etaStripPeak_t;

  inline GCTsupertower_t bestOf2(const GCTsupertower_t& calotp0, const GCTsupertower_t& calotp1) {
    GCTsupertower_t x;
    x = (calotp0.et > calotp1.et) ? calotp0 : calotp1;
    return x;
  }


  inline GCTsupertower_t getPeakBin18N(const etaStrip_t& etaStrip) {
    GCTsupertower_t best01 = bestOf2(etaStrip.cr[0], etaStrip.cr[1]);
    GCTsupertower_t best23 = bestOf2(etaStrip.cr[2], etaStrip.cr[3]);
    GCTsupertower_t best45 = bestOf2(etaStrip.cr[4], etaStrip.cr[5]);
    GCTsupertower_t best67 = bestOf2(etaStrip.cr[6], etaStrip.cr[7]);
    GCTsupertower_t best89 = bestOf2(etaStrip.cr[8], etaStrip.cr[9]);
    GCTsupertower_t best1011 = bestOf2(etaStrip.cr[10], etaStrip.cr[11]);
    GCTsupertower_t best1213 = bestOf2(etaStrip.cr[12], etaStrip.cr[13]);
    GCTsupertower_t best1415 = bestOf2(etaStrip.cr[14], etaStrip.cr[15]);
    GCTsupertower_t best1617 = bestOf2(etaStrip.cr[16], etaStrip.cr[17]);



    GCTsupertower_t best0123 = bestOf2(best01, best23);
    GCTsupertower_t best4567 = bestOf2(best45, best67);
    GCTsupertower_t best891011 = bestOf2(best89, best1011);
    GCTsupertower_t best12131415 = bestOf2(best1213, best1415);

    GCTsupertower_t best0to7 = bestOf2(best0123, best4567);
    GCTsupertower_t best8to15 = bestOf2(best12131415, best891011);

    GCTsupertower_t best8to17 = bestOf2(best8to15, best1617);

    GCTsupertower_t bestOf18 = bestOf2(best0to7, best8to17);

    return bestOf18;
  }

  inline towerMax getPeakBin14N(const etaStripPeak_t& etaStrip) {
    towerMax x;

    GCTsupertower_t best01 = bestOf2(etaStrip.pk[0], etaStrip.pk[1]);
    GCTsupertower_t best23 = bestOf2(etaStrip.pk[2], etaStrip.pk[3]);
    GCTsupertower_t best45 = bestOf2(etaStrip.pk[4], etaStrip.pk[5]);
    GCTsupertower_t best67 = bestOf2(etaStrip.pk[6], etaStrip.pk[7]);
    GCTsupertower_t best89 = bestOf2(etaStrip.pk[8], etaStrip.pk[9]);
    GCTsupertower_t best1011 = bestOf2(etaStrip.pk[10], etaStrip.pk[11]);
    GCTsupertower_t best1213 = bestOf2(etaStrip.pk[12], etaStrip.pk[13]);

    GCTsupertower_t best0123 = bestOf2(best01, best23);
    GCTsupertower_t best4567 = bestOf2(best45, best67);
    GCTsupertower_t best891011 = bestOf2(best89, best1011);

    GCTsupertower_t best8to13 = bestOf2(best891011, best1213);
    GCTsupertower_t best0to7 = bestOf2(best0123, best4567);

    GCTsupertower_t bestOf14 = bestOf2(best0to7, best8to13);

    x.energy = bestOf14.et;
    x.towerPhi = bestOf14.towerPhi;
    x.towerEta = bestOf14.towerEta;
    x.ieta = bestOf14.ieta;
    x.iphi = bestOf14.iphi; 
    
    return x;
  }
  inline towerMax getTowerMax(GCTsupertower_t temp[nSTEta][nSTPhi]) {
    etaStripPeak_t etaStripPeak;

    for (int i = 0; i < nSTEta; i++) {
      etaStrip_t test;
      for (int j = 0; j < nSTPhi; j++) {
        test.cr[j] = temp[i][j];
      }
      etaStripPeak.pk[i] = getPeakBin18N(test);
    }

    towerMax peakIn14;
    peakIn14 = getPeakBin14N(etaStripPeak);
    return peakIn14;
  }


} // namespace gctobj


//apparently calo out coordinates header file has initialized this array already.



// Functions for 3x3 region


jetInfo getJetValues(gctobj::GCTsupertower_t tempX[nSTEta][nSTPhi], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){
#pragma HLS ARRAY_PARTITION variable=tempX complete dim=0
#pragma HLS latency min=6

  ap_uint<12> temp[nSTEta+2][nSTPhi+2] ;
#pragma HLS ARRAY_PARTITION variable=temp complete dim=0

  ap_uint<12> eta_slice[3] ;
#pragma HLS ARRAY_PARTITION variable=eta_slice complete dim=0

  jetInfo jet_tmp;


  for(loop i=0; i<nSTEta+2; i++){
#pragma HLS UNROLL
    for(loop k=0; k<nSTPhi+2; k++){
#pragma HLS UNROLL
      temp[i][k] = 0 ;
    }
  }

  for(loop i=0; i<nSTEta; i++){
#pragma HLS UNROLL
    
    //    std::cout<< "tempX[i][17].et : " << tempX[i][17].et << "\n"<< "tempX[i][0].et : " << tempX[i][0].et << "\n"<<std::endl; 
    
    temp[i+1][0] = tempX[i][17].et;
    temp[i+1][19]= tempX[i][0].et;

    for(loop k=0; k<nSTPhi; k++){
#pragma HLS UNROLL
      temp[i+1][k+1] = tempX[i][k].et ;
      
    }
  }


  ap_uint<5> seed_eta1,  seed_phi1 ;

  seed_eta1 = seed_eta ; //to start from corner                                                                                                                                                     
  seed_phi1 = seed_phi ; //to start from corner                                                                                                                                                     
  ap_uint<12> tmp1, tmp2, tmp3 ;

  for(loop j=0; j<nSTEta; j++){
    for(loop k=0; k<nSTPhi; k++){
#pragma HLS UNROLL
      if(j== seed_eta1 && k == seed_phi1){
	//std::cout << "seed_eta1 : " << j << "seed_phi1 : " << k << std::endl; 
	for(loop m=0; m<3 ; m++){
#pragma HLS UNROLL
	  tmp1 = temp[j+m][k] ;
	  tmp2 = temp[j+m][k+1] ;
	  tmp3 = temp[j+m][k+2] ;
	  eta_slice[m] = tmp1 + tmp2 + tmp3 ; // Sum the energies of 3 3x1 adjacent slices to make the 3x3. 
	}
      }
    }
  }

  jet_tmp.energy=eta_slice[0] + eta_slice[1] + eta_slice[2];
  /*  std::cout << "eta_slice[0] : " << eta_slice[0] << std::endl; 
  std::cout << "eta_slice[1] : " << eta_slice[1] << std::endl;
  std::cout << "eta_slice[2] : " << eta_slice[2] << std::endl;
  std::cout << " jet_tmp.energy : " <<  jet_tmp.energy << std::endl; */
  

  for(loop i=0; i<nSTEta; i++){
    if(i+1>=seed_eta && i<=seed_eta+1){
      for(loop k=0; k<nSTPhi; k++){
#pragma HLS UNROLL
	if(k+1>=seed_phi && k<=seed_phi+1)  tempX[i][k].et = 0 ; // set the 3x3 energies to 0. 
      }
    }
  }
  
  

  return jet_tmp ;
} //end of the getJetValues function



void algo_unpacked(ap_uint<128> link_in[N_CH_IN], ap_uint<192> link_out[N_CH_OUT])
{

// !!! Retain these 4 #pragma directives below in your algo_unpacked implementation !!!
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=4
#pragma HLS INTERFACE ap_ctrl_hs port=return
   

    // null algo specific pragma: avoid fully combinatorial algo by specifying min latency
    // otherwise algorithm clock input (ap_clk) gets optimized away
#pragma HLS latency min=4

//  for (int lnk = 0; lnk < N_CH_OUT ; lnk++) {
//#pragma HLS UNROLL
////  pass-through "algo"
//       link_out[lnk].range(7,0) = 0;            // reserved for 8b10b control word
//       link_out[lnk].range(127,120) = 0;        // reserved for CRC word
//       link_out[lnk].range(119, 8) = link_in[lnk].range(119, 8) ;   // input to output pass-through
//    }
  
        ap_uint<192> tmp_link_out[N_CH_OUT];
#pragma HLS ARRAY_PARTITION variable=tmp_link_out    complete dim=0
        for (int idx = 0; idx < N_CH_OUT; idx++){
#pragma HLS UNROLL
                tmp_link_out[idx] = 0;
        }

        input_t et_calo_ad[N_INPUT_1_1];
#pragma HLS ARRAY_RESHAPE variable=et_calo_ad complete dim=0
        result_t cicada_out[N_LAYER_10];
#pragma HLS ARRAY_PARTITION variable=cicada_out complete dim=0
        region_t centr_region[NR_CNTR_REG];
#pragma HLS ARRAY_PARTITION variable=centr_region complete dim=1
	
	/*	std::cout << "nSTEta : " << nSTEta << std::endl; 
		std::cout<<"nSTPhi : " << nSTPhi << std::endl;  */ 

	//gctobj::GCTsupertower_t temp[nSTEta][nSTPhi];
        regionLoop: for(int iRegion = 0; iRegion < NR_CNTR_REG; iRegion++) {
#pragma HLS UNROLL
                if(iRegion > MaxRegions) {
                        fprintf(stderr, "Too many regions - aborting");
                        exit(1);
                }
                int link_idx = iRegion / NRegionsPerLink;
                int bitLo = ((iRegion - link_idx * NRegionsPerLink) % NRegionsPerLink) * 16 + 8;
                int bitHi = bitLo + 15;
                uint16_t region_raw = link_in[link_idx].range(bitHi, bitLo);
                et_calo_ad[iRegion] = (region_raw & 0x3FF >> 0);   // 10 bits
                centr_region[iRegion].et = (region_raw & 0x3FF >> 0);   // 10 bits
                centr_region[iRegion].eg_veto = (region_raw & 0x7FF) >> 10;   // 1 bit
                centr_region[iRegion].tau_veto = (region_raw & 0xFFF) >> 11;   // 1 bit
                centr_region[iRegion].rloc_phi = (region_raw & 0x3FFF) >> 12;   // 2 bit
                centr_region[iRegion].rloc_eta = (region_raw & 0xFFFF) >> 14;   // 2 bit
		
		
		//int ieta = iRegion % 14  ; 
		//int iphi = iRegion / 14 ; 
		
		//temp[ieta][iphi].eta =  ieta ; // 4*ieta + centr_region[iRegion].rloc_eta; 
		//temp[ieta][iphi].phi = iphi ;  //4*iphi + centr_region[iRegion].rloc_phi;
		//temp[ieta][iphi].et =  centr_region[iRegion].et;
		
		//std::cout << "centr_region[iRegion].rloc_eta : " << centr_region[iRegion].rloc_eta << std::endl; 
		//std::cout << "IRegion : " << iRegion << std::endl;
		//std::cout << "centr_region[iRegion].rloc_phi : " << centr_region[iRegion].rloc_phi << "\n"<< std::endl; */
		
		
		/*std::cout << "(region_raw & 0x3FF >> 0) : " << (region_raw & 0x3FF >> 0) << std::endl; 
		std::cout << "centr_region[iRegion].et : " << centr_region[iRegion].et << std::endl; 
		
		std::cout << "temp[ieta][iphi].eta : " << temp[ieta][iphi].eta << std::endl; 
		//std::cout << "temp[ieta][iphi].phi : " << temp[ieta][iphi].phi << std::endl;*/ 
		//std::cout << "temp[" << ieta << "][" << iphi<< "].et : " << temp[ieta][iphi].et << "\n"<< std::endl;
		//temp[centr_Region.rloc_eta]
        }
	//for(loop i = 0 ; i < nSTEta ; i ++) {std::cout << "temp["<<i<< "][17].et : " << temp[i][17].et <<std::endl;  }
	
	//gctobj::towerMax maxTower  = gctobj::getTowerMax(temp);
	       
	  
	  /*std::cout << "algorithm max et = " << maxTower.energy <<std::endl;
	  std::cout << "algorithm eta = " << maxTower.eta << std::endl;
	  std::cout <<"algorithm phi = " << maxTower.phi << "\n" << std::endl; */
	  
	/* 	jetInfo test_jet; 
	  
	  test_jet.etaMax = maxTower.eta; 
	  test_jet.phiMax = maxTower.phi;
	  
	  jetInfo tmp_jet;

	  
	  tmp_jet = getJetValues(temp,maxTower.eta, maxTower.phi);
	  test_jet.energy = tmp_jet.energy; 
	  
	  std::cout << "test_jet.energy : " << test_jet.energy << std::endl; 
	  std::cout << "test_jet.phi : " << test_jet.phiMax << std::endl; 
	  std::cout << "test_jet.eta : "<< test_jet.etaMax <<  "\n"<<std:: endl; */

	  
        // Anomlay detection algorithm
        cicada(et_calo_ad, cicada_out);

////////////////////////////////////////////////////////////
        // Objets from input
        ap_uint<10> et_calo[NR_CNTR_REG];
        ap_uint<10> pu_sub_et_calo[NR_CNTR_REG];
        ap_uint<10> et_3by3_calo[NR_CNTR_REG];
        ap_uint<10> et_3by3_cntr[NR_CNTR_REG];

        ap_uint<10> et_jet_boosted[NR_SCNTR_REG];
        ap_uint<9> rIdx_boostedjet[NR_SCNTR_REG];

        ap_uint<32> so_in_jet_boosted[64];
        ap_uint<32> so_out_jet_boosted[64];

        ap_uint<NR_CNTR_REG> tmp = 0;
        ap_uint<PUM_LEVEL_BITSIZE> pum_level;
        ap_uint<5> pum_bin;

        region_t centr_region_pu_sub[NR_CNTR_REG];

///////////////////////////////////////////////////////////

        algo_config_t algo_config;

        algo_config.egamma_IsoFact =  0.3;
        algo_config.egamma_seed = 5;
        algo_config.jet_seed = 10;
        algo_config.pum_thr = 0;
        algo_config.tau_IsoFact =  0.3;
        algo_config.tau_seed = 10;

///////////////////////////////////////////////////////////

#pragma HLS INTERFACE ap_none port=algo_config

#pragma HLS ARRAY_RESHAPE variable=centr_region_pu_sub complete dim=1

#pragma HLS ARRAY_RESHAPE variable=so_in_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_jet_boosted complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_calo complete dim=0
#pragma HLS ARRAY_RESHAPE variable=pu_sub_et_calo complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=rIdx_boostedjet complete dim=0

////////////////////////////////////////////////////////////////////////
        //  "pum bin" calculation
        for (int i = 0; i < NR_CNTR_REG; i++)
        {
#pragma HLS UNROLL
                if (centr_region[i].et > algo_config.pum_thr)
                {
                        tmp.set_bit((i), true);
                }
                else
                {
                        tmp.set_bit((i), false);
                }
        }

        // Count number of ones in tmp variable
        pum_level = popcount(tmp);
        pum_bin = pum_level / 14;

////////////////////////////////////////////////////////////
         // Unpack calo ET values in et_calo array
        for (int idx = 0; idx < NR_CNTR_REG; idx++)
        {
#pragma HLS UNROLL
                et_calo[idx] = centr_region[idx].et;
        }

////////////////////////////////////////////////////////////
        // Calculate pile-up subtracted ET values: pu_sub_et_calo
        pu_lut_cntr(pum_bin, et_calo, pu_sub_et_calo);

////////////////////////////////////////////////////////////
        et_3by3(pu_sub_et_calo, et_3by3_calo);
	
	gctobj::GCTsupertower_t temp[nSTEta][nSTPhi]; //create temp region matrix for finding max jet et. 
        for (int idx = 0; idx < NR_CNTR_REG; idx++)
        {
#pragma HLS UNROLL
                centr_region_pu_sub[idx].et = pu_sub_et_calo[idx];
                centr_region_pu_sub[idx].eg_veto = centr_region[idx].eg_veto;
                centr_region_pu_sub[idx].tau_veto = centr_region[idx].tau_veto;
                centr_region_pu_sub[idx].rloc_eta = centr_region[idx].rloc_eta;
                centr_region_pu_sub[idx].rloc_phi = centr_region[idx].rloc_phi;
		
		int ieta = region_lut[idx].ieta  ; // make a look up table converting idx to ieta and iphi
                int iphi = region_lut[idx].iphi ;

		
		calo_coor_t calo_coor_event = calo_coor[idx];
		// std::cout << "calo_coor_event.side : " << calo_coor_event.side << std::endl;

		int towerEta;
		int towerPhi = calo_coor_event.iphi + centr_region_pu_sub[idx].rloc_phi;
		if (calo_coor_event.side > 0){
		  towerEta = -(calo_coor_event.ieta + centr_region_pu_sub[idx].rloc_eta) ;
		  
		}
		else{
		  towerEta = calo_coor_event.ieta + centr_region_pu_sub[idx].rloc_eta;
		  
		}

		
		temp[ieta][iphi].ieta =  ieta ; // 4*ieta + centr_region[iRegion].rloc_eta;
                temp[ieta][iphi].iphi = iphi ;  //4*iphi + centr_region[iRegion].rloc_phi;
                temp[ieta][iphi].et =  centr_region_pu_sub[idx].et;
		temp[ieta][iphi].towerEta = towerEta;
		temp[ieta][iphi].towerPhi = towerPhi;
		//temp[ieta][iphi].eta = getUCTTowerEta(towerEta);
		//temp[ieta][iphi].phi = getUCTTowerPhi(towerPhi);
		
		  }


	
	for(int i = 0 ; i< 2 ; i++ ) {
	gctobj::towerMax maxTower  = gctobj::getTowerMax(temp);

        //std::cout << "maxTower.energy : " << maxTower.energy << "\n" << std::endl;

	// testing out iregiont_t

	/*iregion_t iRegion;

	for (int k = 0; k < 18; k++) {
	  for(int j = 0; j< 14; j ++){

	    std::cout << "iRegion[" << j << "][" << k << "] :" << iRegion[j][k] << std::endl;  
	  }
	  }*/
	jetInfo tmp_jet;
	
	
	tmp_jet = getJetValues(temp,maxTower.ieta, maxTower.iphi);
	
	
        tmp_jet.etaMax = maxTower.towerEta;
	tmp_jet.phiMax = maxTower.towerPhi;
	
	//std::cout << "test_jet.energy : " << tmp_jet.energy << std::endl;
	//std::cout << "test_jet.phi : " << tmp_jet.phiMax << std::endl;
	//std::cout << "test_jet.eta : "<< tmp_jet.etaMax <<  "\n"<<std:: endl;
	}	



	  
////////////////////////////////////////////////////////////
        // Prepare algorithm results
        boostedjet(algo_config.jet_seed, centr_region_pu_sub, et_3by3_calo, et_jet_boosted, rIdx_boostedjet);

        for (int idx = 0; idx < NR_SCNTR_REG; idx++)
        {
#pragma HLS UNROLL
	  ap_uint<9> idx_jet_in = rIdx_boostedjet[idx];
                so_in_jet_boosted[idx].range(9, 0) = et_jet_boosted[idx];
                so_in_jet_boosted[idx].range(18, 10) = idx_jet_in;
                so_in_jet_boosted[idx].range(25, 19) = centr_region[idx_jet_in].rloc_phi;
                so_in_jet_boosted[idx].range(31, 26) = centr_region[idx_jet_in].rloc_eta;
        }
        so_in_jet_boosted[63] = 0;

        // Sorting objects
        bitonicSort64(so_in_jet_boosted, so_out_jet_boosted);

        // Assign the algorithm outputs
        tmp_link_out[0].range(31, 28) = cicada_out[0].range(15, 12);
        tmp_link_out[0].range(63, 60) = cicada_out[0].range(11, 8);
        tmp_link_out[0].range(95, 92) = cicada_out[0].range(7, 4);
        tmp_link_out[0].range(127, 124) = cicada_out[0].range(3, 0);


        // printing the leading boosted jet information from pattern based algorithm
        int word = 32;
        for (int idx = 0; idx < 6; idx++) {
#pragma HLS UNROLL
            /* Boosted jets
               output scheme: 6x32-bits
               For the format the interface document states that the current jet
               collection is 8 bits phi then 8 bits in eta (7 bits position then 1 bit
               for +/- eta) then 11 bits et and 1 bit for flag.
            */
            ap_uint<9> idx_srt;
            idx_srt = so_out_jet_boosted[idx].range(18, 10);

            int bLoET = idx*word;
            int bHiET = bLoET + 10;
            tmp_link_out[0].range(bHiET, bLoET) = so_out_jet_boosted[idx].range(9, 0);

            int bLoEta = bHiET + 1;
            int bHiEta = bLoEta + 7;
            int ieta = 0x003F & calo_coor[idx_srt].ieta + so_out_jet_boosted[idx].range(31, 26);
            ap_uint<1> isNegativeSide = calo_coor[idx_srt].side;
            tmp_link_out[0].range(bHiEta, bLoEta) = ieta_lut[isNegativeSide][ieta];

            int bLoPhi = bHiEta + 1;
            int bHiPhi = bLoPhi + 7;
            int test1 = 0x007F & calo_coor[idx_srt].iphi + so_out_jet_boosted[idx].range(25, 19);
            int iphi = !signbit(test1 - 72) ? (0x007F & test1 - 0x0048) : (0x007F & test1);
            tmp_link_out[0].range(bHiPhi, bLoPhi) = iphi_lut[iphi];
	    //  std::cout<<"jet number: "<<idx<<"\t"<<"jet et: "<<so_out_jet_boosted[idx].range(9, 0)<<"\t"<<"eta: "<<ieta<<"\t"<<"phi: "<<iphi<<std::endl;
	    
       }
	//std::cout << "End of run" << "\n" << std::endl; 
	
        for(int i = 0; i < N_CH_OUT; i++){
#pragma HLS unroll
                link_out[i] = tmp_link_out[i];
        }
}

////////////////////////////////////////////////////////////
// count number of ones in bitString
ap_uint<8> popcount(ap_uint<NR_CNTR_REG> bitString)
{
#pragma HLS PIPELINE II=4
        ap_uint<9> popcnt = 0;

        for (ap_uint<9> b = 0; b < NR_CNTR_REG; b++)
        {
#pragma HLS unroll
                popcnt += ((bitString >> b) & 1);
        }
        return popcnt;
}

