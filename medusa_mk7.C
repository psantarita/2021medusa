#include "TFile.h"      //Write root files
#include "TTree.h"      // Make trees
#include "TObject.h"    // Write objects
#include "TCanvas.h"    // TCanvas
#include "TH1.h"        //  1D histograms
#include "TH2.h"        //  2D histograms
#include "TStopwatch.h"  // stopwatch for time
#include "TChain.h"     // Tchain to process several root files at the time
#include "TMath.h"      //TMath::Sort
#include "TGraph.h"
#include "TGraph2D.h"

#include "TF1.h"
#include <sstream>      // for ostringstream
#include <stdlib.h>     // srand, rand
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <math.h>       // sqrt
#include <string>       //string
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>

#include <vector>

#include "TGeoPolygon.h"

#include "TObjArray.h"
#include "TMath.h"
#include "TGeoShape.h"
#include "TGeoManager.h"
#include "TVirtualGeoPainter.h"
using namespace std;

//Colours for format and visuals
#define RESET   "\033[0m"
#define BLUE    "\033[34m"	/* Blue */
#define WHITE   "\033[37m"	/* White */
#define GREEN   "\033[32m"	/* Green */
#define YELLOW  "\033[33m"	/* Yellow */
#define RED     "\033[31m"	/* Red */


Double_t position(Int_t); // calculates the position using the adc number
Double_t positionDE(Int_t); // calculates the position using the adc number for DE detectos
Double_t positionc(Int_t); // calculates the position using the adc number for chunky
Double_t p_alpha(Double_t); // alpha momentum from energy meassured in DSSD
Double_t angle_theta(Double_t); // angle between Z axis(reaction target) and X position (DSSD)
Double_t angle_theta_ch(Double_t); // angle between Z axis(reaction target) and X position (DSSD)
float scattering(float energy, float theta,float phi, float mass,float Eb);
Double_t thetaLaBr(Int_t);
Double_t phiLaBr(Int_t);
Double_t energyloss(Double_t,Double_t,Double_t);
Double_t energyloss_table(Double_t aEnergy_i,
  Double_t atheta,
  Double_t mat_thick,
  Double_t aenergyEL[],
  Double_t arangemg[],
  Double_t apar_indexvsE[],
  Double_t apar_indexvsR[]);
Double_t cat_x(float energy, float theta,float phi);
  // to execute this progam type in root: .x medusa.C++
  class single_event {
    public:
      Double_t labrE;
      Double_t dssdE;
      Double_t tdcT;
      Double_t theta;
      Double_t phi;
      Double_t catx;
      Double_t caty;
      Int_t run;
      single_event(Double_t labrE_i,
                   Double_t dssdE_i,
                   Double_t tdcT_i,
                   Double_t theta_i,
                   Double_t phi_i,
                   Double_t catx_i,
                   Double_t caty_i,
                   Int_t run_i
                  ){

        labrE =labrE_i;
        dssdE =dssdE_i;
        tdcT= tdcT_i;
        phi = phi_i;
        theta =theta_i;
        catx =  catx_i,
        caty = caty_i;
        run = run_i;
      }
  };
  Int_t pnpoly(int nvert, Double_t *vertx, Double_t *verty, Double_t testx, Double_t testy)
  {
    Int_t i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
      if ( ((verty[i]>=testy) != (verty[j]>=testy)) &&
       (testx <= (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
         c = !c;
    }
    return c;
  }

  void medusa_mk7(){
    TFile f("medusa43.root","RECREATE"); // opens the root file
    f.Print();
    //calibration file read

    ofstream myfiletxt;
    myfiletxt.open ("gammaerrs.log",ios::out);
    myfiletxt<< "calibcounter" << "\t" << "eadc[k]" << "\t" << "prerealampl" << "\t"  << "realampl" <<endl;
    Double_t calibp0_0[1000]={0}; //allocate space for the values from calibration
    Double_t calibp1_0[1000]={0};
    Double_t calibp2_0[1000]={0};
    Double_t calibp3_0[1000]={0};

/////////////////////////calibration routine :////////////////////////
    Double_t clb_array[30][300][5]={0};
    //23 calibration files
    //223 channels
    //4 coeficients
cout <<  "start finding the calib files" << endl;
    for (Int_t i = 0; i < 26; i++) {
      Double_t d2=1,d3=0,d4=0,d5=0; // clean variables
      Int_t d1,fread_dump_0; // fread_dump is a variable to dump the fread return value, this is only to run the progam without waring
      FILE *out_0; // file declaration
      auto filenameklb = Form("calib_backup/clbfilestotal/calibtotal_%i.txt",i);
      out_0=fopen(filenameklb,"r"); // opens the calib file with the name typed
      if (out_0!=NULL) { // only works for actual files
        cout << GREEN<< "running... "<< i << " batch calibration file found: " << filenameklb << RESET <<  endl;
        while(!feof(out_0))  { // reads until is over
          fread_dump_0=fscanf(out_0,"%i %le %le %le %le",&d1,&d2,&d3,&d4,&d5);
          clb_array[i][d1][0]=d2; // power 1   // power 2
          clb_array[i][d1][1]=d3; /// powe 0   // power 1
          clb_array[i][d1][2]=d4;/// power 2   // power 0
          clb_array[i][d1][3]=d5; /// power 3  // power 3
          if (i==21) {
            cout << d1 << "\t" << clb_array[i][d1][0] << "\t"
                               << clb_array[i][d1][1] << "\t"
                               << clb_array[i][d1][2] << "\t"
                               << clb_array[i][d1][3] << "\t"
                               << endl;

          }
        }
      } else {
        cout <<RED<< "calibration file num: "<< i << " not found " <<  filenameklb << RESET<< endl;
      }
    }
cout <<  "done finding the calib files" << endl;
/// TDC calibration  need to improve this, but this will do the trick for now
    Double_t tdc_coef_p1[10]={},tdc_coef_p0[10]={};
    tdc_coef_p1[0] = 0.33017435 ;
    tdc_coef_p1[1] = 0.32570295 ;
    tdc_coef_p1[2] = 0.32370842 ;
    tdc_coef_p1[3] = 0.32752085 ;
    tdc_coef_p1[4] = 0.32502424 ;
    tdc_coef_p1[5] = 0.32200503 ;
    tdc_coef_p1[6] = 0.3297731  ;
    tdc_coef_p1[7] = 0.33459657 ;
    tdc_coef_p1[8] = 0.32591252 ;
    tdc_coef_p1[9] = 0.32245771 ;

    tdc_coef_p0[0] = -221.36628031;
    tdc_coef_p0[1] = -220.47835101;
    tdc_coef_p0[2] = -220.05450636;
    tdc_coef_p0[3] = -219.65443781;
    tdc_coef_p0[4] = -221.23500122;
    tdc_coef_p0[5] = -220.31178139;
    tdc_coef_p0[6] = -220.65946449;
    tdc_coef_p0[7] = -219.6821548 ;
    tdc_coef_p0[8] = -221.0706936 ;
    tdc_coef_p0[9] = -220.79030484;
///////////////////////// end of calibration /////////////////////////////////
///////////Energy Calibration for MEDUSA//////////////////
// Double_t d2_0=1,d3_0=0,d4_0=0,d5_0=0; // clean variables
// Int_t d1_0,fread_dump_0; // fread_dump is a variable to dump the fread return value, this is only to run the progam without waring
// FILE *out_0; // file declaration
// string mystringcalib_0;
// // ENERGY  AND TIME CALIBRATION DSSDs, initial labr3, TDCs
// mystringcalib_0 = "rawdata/3alpha23.txt";  /// calibration from beam time
// out_0=fopen(mystringcalib_0.c_str(),"r"); // opens the calib file with the name typed
// if (out_0!=NULL) { // only works for actual files
//   cout <<  "found DSSDs energy calibration file" << endl;
//   while(!feof(out_0))  { // reads until is over
//     fread_dump_0=fscanf(out_0,"%i %le %le %le %le",&d1_0,&d2_0,&d3_0,&d4_0,&d5_0);
//     calibp0_0[d1_0]=d2_0; // order 1
//     calibp1_0[d1_0]=d3_0; // order 0
//     calibp2_0[d1_0]=d4_0; // order 2
//     calibp3_0[d1_0]=d5_0; // order 3
//     //cout << d1_0 << "\t" <<calibp0_0[d1_0] << "\t" << calibp1_0[d1_0] << "\t" << calibp2_0[d1_0]<< "\t" << calibp3_0[d1_0] << endl;
//   }
// } else {
//   cout << "missing DSSDs energy calibration file" << endl;
// }
//////////////////////////////////////////////////////////
///POLYGONS///

Double_t x_t[29],y_t[29];
x_t[0]=26.2865;      y_t[0]=8.48897;
x_t[1]=10.6387;      y_t[1]=7.03363;
x_t[2]=10.1369;      y_t[2]=6.75248;
x_t[3]=11.1861;      y_t[3]=6.52095;
x_t[4]=13.0566;      y_t[4]=6.47133;
x_t[5]=15.292     ; y_t[5]=6.7194;
x_t[6]=22.7737;      y_t[6]=7.34785;
x_t[7]=26.7883;      y_t[7]=7.66207;
x_t[8]=31.031;      y_t[8]=7.92668;
x_t[9]=35.5018;      y_t[9]=8.20783;
x_t[10]=39.8814;      y_t[10]=8.50551;
x_t[11]=44.9909;      y_t[11]=8.63782;
x_t[12]=49.5073;      y_t[12]=8.68743;
x_t[13]=51.3321;      y_t[13]=8.65436;
x_t[14]=53.0657;      y_t[14]=8.52205;
x_t[15]=55.0274;      y_t[15]=8.38975;
x_t[16]=54.9818;      y_t[16]=8.52205;
x_t[17]=53.385;      y_t[17]=9.03473;
x_t[18]=51.8339;      y_t[18]=9.24972;
x_t[19]=46.0858;      y_t[19]=9.7624;
x_t[20]=44.1697;      y_t[20]=9.91125;
x_t[21]=40.885;      y_t[21]=10.1428;
x_t[22]=40.0639;      y_t[22]=10.0766;
x_t[23]=37.5091;      y_t[23]=9.72933;
x_t[24]=33.9051;      y_t[24]=9.26626;
x_t[25]=30.8942;      y_t[25]=8.91896;
x_t[26]=27.427;      y_t[26]=8.60474;
x_t[27]=26.7427;      y_t[27]=8.53859;
x_t[28]=26.2865;      y_t[28]=8.48897;


  /////cUT FOR 6.1
  Double_t x_61[23],y_61[23];
  x_61[0]=20.5839; y_61[0]=7.92668;
  x_61[1]=15.885; y_61[1]=7.57938;
  x_61[2]=11.6423; y_61[2]=7.33131;
  x_61[3]=10; y_61[3]=7.03363;
  x_61[4]=9.90876; y_61[4]=6.76902;
  x_61[5]=10.5931; y_61[5]=6.57056;
  x_61[6]=12.4635; y_61[6]=6.66979;
  x_61[7]=28.3394; y_61[7]=7.56284;
  x_61[8]=36.6423; y_61[8]=7.9763;
  x_61[9]=45.4927; y_61[9]=8.25744;
  x_61[10]=50.6934; y_61[10]=8.27398;
  x_61[11]=55.6661; y_61[11]=8.14168;
  x_61[12]=57.1259; y_61[12]=8.04245;
  x_61[13]=54.3431; y_61[13]=9.10088;
  x_61[14]=53.6131; y_61[14]=9.24972;
  x_61[15]=45.0821; y_61[15]=10.0932;
  x_61[16]=44.2609; y_61[16]=9.94432;
  x_61[17]=42.2993; y_61[17]=9.66318;
  x_61[18]=37.8285; y_61[18]=9.20011;
  x_61[19]=31.6697; y_61[19]=8.68743;
  x_61[20]=25.0547; y_61[20]=8.2409;
  x_61[21]=20.4927; y_61[21]=7.94322;
  x_61[22]=20.5839; y_61[22]=7.92668;


  /////cUT FOR 6.9
  Double_t x_69[27],y_69[27];
  x_69[0]=20.3558; y_69[0]=8.68743;
  x_69[1]=14.1515; y_69[1]=8.2409;
  x_69[2]=12.6916; y_69[2]=8.12514;
  x_69[3]=11.5967; y_69[3]=7.81092;
  x_69[4]=11.4142; y_69[4]=7.61246;
  x_69[5]=12.0529; y_69[5]=7.48015;
  x_69[6]=14.2883; y_69[6]=7.61246;
  x_69[7]=28.2938; y_69[7]=8.38975;
  x_69[8]=38.5584; y_69[8]=8.81974;
  x_69[9]=42.938; y_69[9]=8.90243;
  x_69[10]=45.4927; y_69[10]=8.96858;
  x_69[11]=50.1004; y_69[11]=8.86935;
  x_69[12]=53.6131; y_69[12]=8.70397;
  x_69[13]=55.5748; y_69[13]=8.57166;
  x_69[14]=55.5292; y_69[14]=8.70397;
  x_69[15]=54.2062; y_69[15]=9.16703;
  x_69[16]=51.6515; y_69[16]=9.48126;
  x_69[17]=47.0894; y_69[17]=9.91125;
  x_69[18]=42.3449; y_69[18]=10.2585;
  x_69[19]=40.5657; y_69[19]=10.3082;
  x_69[20]=39.2427; y_69[20]=10.242;
  x_69[21]=35.5474; y_69[21]=9.79548;
  x_69[22]=31.6241; y_69[22]=9.48126;
  x_69[23]=26.7427; y_69[23]=9.11742;
  x_69[24]=22.4088; y_69[24]=8.85281;
  x_69[25]=20.2646; y_69[25]=8.63782;
  x_69[26]=20.3558; y_69[26]=8.68743;

  ///// 12C ground state cut:
  Double_t x_GS[18],y_GS[18];
x_GS[0]=54.4343; y_GS[0]=6.14057;
x_GS[1]=42.8011; y_GS[1]=4.32139;
x_GS[2]=30.5748; y_GS[2]=3.16373;
x_GS[3]=17.4818; y_GS[3]=2.07222;
x_GS[4]=9.63504; y_GS[4]=1.50992;
x_GS[5]=7.08029; y_GS[5]=1.04686;
x_GS[6]=7.71898; y_GS[6]=0.550717;
x_GS[7]=12.281; y_GS[7]=0.51764;
x_GS[8]=24.7354; y_GS[8]=1.27839;
x_GS[9]=41.615; y_GS[9]=2.83297;
x_GS[10]=47.2263; y_GS[10]=3.42834;
x_GS[11]=67.8011; y_GS[11]=4.371;
x_GS[12]=66.7974; y_GS[12]=5.0656;
x_GS[13]=62.8285; y_GS[13]=5.99173;
x_GS[14]=59.3157; y_GS[14]=6.4548;
x_GS[15]=56.9891; y_GS[15]=6.4548;
x_GS[16]=56.9891; y_GS[16]=6.4548;
x_GS[17]=54.4343; y_GS[17]=6.14057;



  Double_t x_44[26],y_44[26];
  x_44[0]=18.4854; y_44[0]=6.57056;
  x_44[1]=12.281; y_44[1]=6.00827;
  x_44[2]=8.85949; y_44[2]=5.69405;
  x_44[3]=7.90146; y_44[3]=5.2806;
  x_44[4]=8.12956; y_44[4]=5.01599;
  x_44[5]=9.58942; y_44[5]=4.91676;
  x_44[6]=10.5474; y_44[6]=4.99945;
  x_44[7]=28.385; y_44[7]=6.42172;
  x_44[8]=37.5547; y_44[8]=7.08324;
  x_44[9]=47.135; y_44[9]=7.64553;
  x_44[10]=49.6442; y_44[10]=7.69515;
  x_44[11]=52.7464; y_44[11]=7.74476;
  x_44[12]=57.3996; y_44[12]=7.67861;
  x_44[13]=58.2664; y_44[13]=7.61246;
  x_44[14]=58.5401; y_44[14]=7.66207;
  x_44[15]=54.2974; y_44[15]=9.10088;
  x_44[16]=52.792; y_44[16]=9.39857;
  x_44[17]=48.0018; y_44[17]=9.84509;
  x_44[18]=46.542; y_44[18]=9.54741;
  x_44[19]=44.7172; y_44[19]=9.11742;
  x_44[20]=40.9307; y_44[20]=8.67089;
  x_44[21]=36.4142; y_44[21]=8.15821;
  x_44[22]=29.2062; y_44[22]=7.49669;
  x_44[23]=21.0858; y_44[23]=6.75248;
  x_44[24]=19.5347; y_44[24]=6.65325;
  x_44[25]=18.4854; y_44[25]=6.57056;


int arrSize61 = sizeof(x_61)/sizeof(x_61[0]);
int arrSize44 = sizeof(x_44)/sizeof(x_44[0]);
int arrSize69 = sizeof(x_69)/sizeof(x_69[0]);
int arrSizeGS = sizeof(x_GS)/sizeof(x_GS[0]);
int arrSizet = sizeof(x_t)/sizeof(x_t[0]);
cout << pnpoly(arrSizet,x_t,y_t,20,7)<< endl;
cout << pnpoly(arrSize69,x_69,y_69,20,7)<< endl;
cout << pnpoly(arrSize61,x_61,y_61,20,7)<< endl;

///////
    TStopwatch StopWatch; //stopwatch to keep on track of efficiency
    StopWatch.Start(); // start of the stopwatch
    TChain chain("data");
    long int n_entries=0;
    /// eloss C
    ///// energyloss C
    FILE *out_eloss;
    Double_t energyEL[300003]={0}; // order 1
    Double_t rangemg[300003] ={0}; // order 3
    Double_t indexarr[300003]={0};
    string mystringeloss;
    mystringeloss= "alphainCatima.txt";
    Int_t index=0;
    Int_t fread_dumpe=0;
    Double_t e1=0,e2=0;
    out_eloss=fopen(mystringeloss.c_str(),"r"); // opens the calib file with the name typed
    if (out_eloss!=NULL) { // only works for actual files
      cout <<  "Energy_C loss file found" << endl;
      while(!feof(out_eloss))  { // reads until is over
        fread_dumpe=fscanf(out_eloss,"%le %le",&e1,&e2);
        indexarr[index]=index;
        energyEL[index]=e1; // order 1
        rangemg[index] =e2; // order 3
        index++;
        //cout << index<< " , " << e1 << " , "<< e2  << endl;
        //printf("%le %le \n",e1, e2 );
      }
    } else {
      cout << "Energy_C loss file found file not found" << endl;
    }
    ///// energyloss Au
    FILE *out_eloss_Au ;
    Double_t energyEL_Au[15003]={0}; // order 1
    Double_t dedxmm_Au[15003]  ={0}; // order 0
    Double_t dedxmg_Au[15003]  ={0}; // order 2
    Double_t rangemm_Au[15003] ={0}; // order 3
    Double_t rangemg_Au[15003] ={0}; // order 3
    Double_t indexarr_Au[15003] ={0};
    string mystringeloss_Au;
    mystringeloss_Au= "alphainAu.txt";
    Int_t index_Au=0;
    Int_t fread_dumpe_Au=0;
    Double_t e1_Au=0,e2_Au=0,e3_Au=0,e4_Au=0,e5_Au=0;
    out_eloss_Au=fopen(mystringeloss_Au.c_str(),"r"); // opens the calib file with the name typed
    if (out_eloss_Au!=NULL) { // only works for actual files
      cout <<  "Energy_Au loss file found" << endl;
      while(!feof(out_eloss_Au))  { // reads until is over
        fread_dumpe_Au=fscanf(out_eloss_Au,"%le %le %le %le %le",&e1_Au,&e2_Au,&e3_Au,&e4_Au,&e5_Au);
        indexarr_Au[index_Au]=index_Au;
        energyEL_Au[index_Au]=e1_Au;
        dedxmm_Au[index_Au]  =e2_Au;
        dedxmg_Au[index_Au]  =e3_Au;
        rangemm_Au[index_Au] =e4_Au;
        rangemg_Au[index_Au] =e5_Au;
        index_Au++;
        //cout << index_Au << " , " << e1_Au << " , "<< e2_Au << " , "<< e3_Au << " , "<< e4_Au << " , "<< e5_Au << endl;
      }
    } else {
      cout << "Energy_Au loss file found file not found" << endl;
    }
    //cout << index << endl;
    //line for C
    Double_t par_indexvsE[6]={0};
    Double_t par_indexvsR[9]={0};
    //cout << index << endl;
    TGraph *indexvsE = new TGraph(index,energyEL,indexarr);
    TGraph *indexvsR = new TGraph(index,rangemg,indexarr);
    TF1 *line_indexvsE  = new TF1("line_indexvsE","pol1",0,index);
    TF1 *line_indexvsR  = new TF1("line_indexvsR","pol9",0,rangemg[index]);
    indexvsE->Fit("line_indexvsE","Q");
    indexvsR->Fit("line_indexvsR","Q");
    line_indexvsE->GetParameters(&par_indexvsE[0]);
    line_indexvsR->GetParameters(&par_indexvsR[0]);
    //line for Au
    Double_t par_indexvsE_Au[6]={0};
    Double_t par_indexvsR_Au[9]={0};
    TGraph *indexvsE_Au = new TGraph(index_Au,energyEL_Au,indexarr_Au);
    TGraph *indexvsR_Au = new TGraph(index_Au,rangemg_Au,indexarr_Au);
    TF1 *line_indexvsE_Au  = new TF1("line_indexvsE_Au","pol1",0,index_Au);
    TF1 *line_indexvsR_Au  = new TF1("line_indexvsR_Au","pol9",0,rangemg_Au[index_Au]);
    indexvsE_Au->Fit("line_indexvsE_Au","Q");
    indexvsR_Au->Fit("line_indexvsR_Au","Q");
    line_indexvsE_Au->GetParameters(&par_indexvsE_Au[0]);
    line_indexvsR_Au->GetParameters(&par_indexvsR_Au[0]);

    // for (Int_t  p = 0; p < 6; p++) {
    //   cout << par_indexvsE[p] << endl;
    //   cout << par_indexvsR[p] << endl;
    // }

    indexvsE->Write();
    indexvsR->Write();
    indexvsE_Au->Write();
    indexvsR_Au->Write();

    // chain.Add("rawdata/R20_0.root");
    //chain.Add("rawdata/R310_0.root");
    // chain.Add("rawdata/R295_0.root");

    //// end of testing single files nomsaying?

    long long eventsuntil[26]={0};
    TString Filelist[26] = {
                 "rawdata/09aug.root",// 0
                 "rawdata/10aug.root",// 0
                 "rawdata/11aug.root",// 0
                 "rawdata/12aug.root",//1
                 "rawdata/16aug.root",//1
                 "rawdata/17aug.root",//1
                 "rawdata/18aug.root",//1
                 "rawdata/23aug.root", // 2
                 "rawdata/24aug.root", // 2
                 "rawdata/25aug0.root", // 2
                 "rawdata/25aug1.root", // 2
                 "rawdata/26aug0.root", // 2
                 "rawdata/26aug1.root", // 2
                 "rawdata/27aug.root",//3
                 "rawdata/31aug.root",//3
                 "rawdata/01sep.root",//3
                 "rawdata/07sep.root",//4
                 "rawdata/08sep.root",//4
                 "rawdata/10sep.root",//4
                 "rawdata/13sep.root",//4
                 "rawdata/14sep.root",//5
                 "rawdata/15sep.root",//5
                 "rawdata/20sep.root",//5
                 "rawdata/21sep.root",//5
                 "rawdata/22sep.root", //6
                 "rawdata/24sep.root" //6
               };



    long long numnum = sizeof(eventsuntil)/sizeof(*eventsuntil);
    long int eventcounterc = 0;
    long int eventcountercarr[50] = {0};
    for (Int_t i = 0; i < numnum; i++) {
        chain.Add(Filelist[i]);
        eventsuntil[i]=chain.GetEntries();
        cout <<  eventsuntil[i]-eventcounterc << endl;
        eventcountercarr[i] = eventcounterc=+eventsuntil[i];
        // cout << "AAAAAAAA: " <<eventcounterc<< endl;
    }

    cout <<  endl;
    auto otherarr = chain.GetListOfFiles();
    otherarr->Print();

    Int_t numofchainelements = otherarr->GetEntries();
    cout<< "number of link in the chain: " <<  numofchainelements << endl;

    Int_t eadc[600]; // allocate the space for the adc
    Int_t eampl[600]; // allocate the space for the ampl

    Int_t emult=0;// start Multiplicity variables
    n_entries = chain.GetEntries(); // number of entries
   cout << "Total Events in all chains: " << n_entries << endl;
    // printf("i = %lld\n", n_entries);
    Int_t t=0;

    chain.SetBranchAddress("emult", &emult); // this is to direct the read variables into a virutal branch
    chain.SetBranchAddress("eadc", eadc);// this is without the & because eadc is an array
    chain.SetBranchAddress("eampl",eampl);// this is without the & because eampl is an array

    //Tree set up to store the most IMportant DATa, angle of scatter, labrt data, etc.

    TTree *imdat = new TTree("imdat","Tree data");

    Double_t scat_angphi=0,
             scat_ene=0,scat_ang=0,labr3t_ene=0,Ex_scat_12Cpluslabr=0,
             Ex_scat_16Opluslabr=0,TDC_tim=0,
             Ex_scat_12C=0,Ex_scat_16O=0;

    //catania singlehit
    // Double_t cathit[2]={0};
    // imdat->Branch("cathit",&cathit,"cathit[2]/D");


    // imdat->Branch("scat_ene",&scat_ene,"scat_ene/D");
    // imdat->Branch("scat_ang",&scat_ang,"scat_ang/D");
    // imdat->Branch("scat_angphi",&scat_angphi,"scat_angphi/D");
    // imdat->Branch("labr3t_ene",&labr3t_ene,"labr3t_ene/D");
    // imdat->Branch("TDC_tim",&TDC_tim,"TDC_tim/D");

    auto event = new single_event(0.,0.,0.,0.,0.,0.,0.,0);
    imdat->Branch("EventBranch",&event);
    // Double_t labr3_me[5];
    // imdat->Branch("labr3_me", labr3_me, "labr3_me[5]/D");

    // imdat->Branch("polyID", polyID, "polyID[5]/I");

    Int_t labr3_hitcounter=0;
    // imdat->Branch("labr3_hitcounter",&labr3_hitcounter, "labr3_hitcounter/I");
    Int_t tdc_hitcounter=0;
    // imdat->Branch("tdc_hitcounter",&tdc_hitcounter, "tdc_hitcounter/I");
    // Int_t reactionID=0;
    // imdat->Branch("reactionID",&reactionID, "reactionID/I");
    // const int maxnumofhits = 20;
    //
    // Double_t labr3_me[maxnumofhits];
    // Double_t tdc_mt[maxnumofhits];
    //
    // imdat->Branch("labr3_me", labr3_me, "labr3_me[labr3_hitcounter]/D"); // Branch with the adc number
    // imdat->Branch("tdc_mt",tdc_mt,"tdc_mt[tdc_hitcounter]/D"); // amplitude number

    // branhes for double hits in detectors



    // branches for TDC and LaBr3, these need to be of variable lenght
    // Int_t mult_tdc=0,mult_labr3=0; //variables of multiplicty
    // Double_t tdc_branch[10]={0},labr3_branch[10]={0};
    //
    // imdat->Branch("mult_labr3", &mult_labr3, "mult_labr3/I");  // Multiplicity of the event
    // imdat->Branch("mult_tdc", &mult_tdc, "mult_tdc/I");  // Multiplicity of the event
    // imdat->Branch("tdc_branch", tdc_branch, "tdc_branch[mult_tdc]/F"); // Branch with the adc number
    // imdat->Branch("labr3_branch",labr3_branch,"labr3_branch[mult_labr3]/F"); // amplitude number

    //tac histogram
    // charged particles left
    TH2* mapleftDEcm          = new TH2D("mapleftDEcm", "mapleftDEmm;mm;mm",    400,-60,60,400,-60,60);
    TH2* maprightDEcm         = new TH2D("maprightDEcm", "maprightDEmm;mm;mm",    400,-60,60,400,-60,60);
    TH2* mapleft              = new TH2D("mapleft", "Map left;Angle(rad);Angle(rad)",    200,-1.5,1.5,200,-1.5,1.5);
    TH2* mapleftDE            = new TH2D("mapleftDE", "Map left DE;Angle(rad);Angle(rad)",   200,-1.5,1.5,200,-1.5,1.5);
    TH2* maprightDE            = new TH2D("maprightDE", "Map left DE;Angle(rad);Angle(rad)",   200,-1.5,1.5,200,-1.5,1.5);


    TH2* mapthetaR = new TH2D("mapthetaR","mapthetaR",600,-1.5,1.5,600,-1.5 ,1.5);
    TH2* mapphiR= new TH2D("mapphiR","mapphiR",600,-3,3,600,-3,3);

    TH2* mapthetaL = new TH2D("mapthetaL","mapthetaL",600,0,1.5,600,0,1.5);
    TH2* mapphiL= new TH2D("mapphiL","mapphiL",600,-3,3,600,-3,3);

    TH2* xEDER = new TH2D("xEDER","xEDER",600,-30,30,600,-30,30);
    TH2* xEDEL = new TH2D("xEDEL","xEDEL",600,-30,30,600,-30,30);
    TH2* yEDER = new TH2D("yEDER","yEDER",600,-30,30,600,-30,30);
    TH2* yEDEL = new TH2D("yEDEL","yEDEL",600,-30,30,600,-30,30);

    TH2* maprightEmm            = new TH2D("maprightEmm", "maprightEmm;Angle(rad);Angle(rad)",200,-50,50,200,-50,50);


    TH1I* adchits     = new TH1I("adchits", "adchits;counts;ADC number",160,-10,150);

    TH1* Ex_C12_scat_single   = new TH1D("Ex_C12_scat_single", "Ex_C12_scat_single;Energy(MeV);Counts", 1000, -10, 10);
    TH1* Ex_C12_scat          = new TH1D("Ex_C12_scat", "Ex_C12_scat;Energy(MeV);Counts", 1000, -10, 10);
    TH2* hitmap               = new TH2D("hitmap", "hitmap;Angle(rad);Angle(rad)", 200,-1.5,1.5, 200, -0.6, 0.6);

    TH2* EDEL                  = new TH2D("EDEL", "E-DE;E;DE",   400,0,15,400,0,15);
    TH2* EDER                  = new TH2D("EDER", "E-DE;E;DE",   400,0,15,400,0,15);

    TH2* mapchunky            = new TH2D("mapchunky", "mapchunky;Angle(rad);Angle(rad)",    200,0.186,1.34,200, -0.6, 0.6 );
    TH2* scatanglechunky      = new TH2D("scatanglechunky","scatanglechunky;Angle(rad);Angle(rad)",400,-1.7,1.7,400,1.4,15);
    TH2* hitmapsingle         = new TH2D("hitmapsingle", "hitmapsingle;Angle(rad);Angle(rad)", 200,-1.7,1.7, 200, -0.6, 0.6);
    TH2* scatangle            = new TH2D("scatangle","scatangle;Energy(MeV);Angle(rad)",400,-1.7,1.7,400,1.4,15);
    TH2* scatangleall         = new TH2D("scatangleall","scatangleal;Energy(MeV);Angle(rad)",400,-1.7,1.7,400,1.4,15);
    TH2* mapright             = new TH2D("mapright", "mapright;Angle(rad);Angle(rad)",   200,-1.5,1.5,200,-1.5,1.5  );
    TH2D* ADCmap               = new TH2D("ADCmap","ADCmap;ADC number;Energy(MeV)",226,0,226,800,0,16);// maps of adc

    TH2* mapleftEcm           = new TH2D("mapleftEcm", "mapleftDEmm;mm;mm",    400,-60,60,400,-60,60);
    TH2* maprightEcm          = new TH2D("maprightEcm", "mapleftDEmm;mm;mm",    400,-60,60,400,-60,60);

    TH2* mapchunkymm         = new TH2D("mapchunkymm", "mapchunkymm;mm;mm",    400,-60,60,400,-60,60);

    TH2D* RvsDL               = new TH2D("RvsDL","Right vs D_Left ;Energy(MeV);Energy(MeV)",200,0,12,200,0,12);// maps of adc

    TH1D* labr_GS_t     = new TH1D("labr_GS_t", "labr_GS_total;Energy(MeV);Counts",1000,0.,8);
    TH1D* labr_44_t     = new TH1D("labr_44_t", "labr_44_total;Energy(MeV);Counts",1000,0.,8);
    TH1D* labr_61_t     = new TH1D("labr_61_t", "labr_61_total;Energy(MeV);Counts",1000,0.,8);
    TH1D* labr_69_t     = new TH1D("labr_69_t", "labr_69_total;Energy(MeV);Counts",1000,0.,8);
    TH1D* labr_t_t      = new TH1D("labr_t_t", "labr_t_t;Energy(MeV);Counts",1000,0.,8);

    TH2* catania_single              = new TH2D("catania_single", "catania_single;Energy(MeV);Energy(MeV)",500,0,70,500,0,12);
    // CATANIA histograms
    TH2* catania44               = new TH2D("catania44", "catania44;Energy(MeV);Energy(MeV)",500,0,70,500,0,12);
    TH2* catania61               = new TH2D("catania61", "catania61;Energy(MeV);Energy(MeV)",500,0,70,500,0,12);
    TH2* catania69               = new TH2D("catania69", "catania69;Energy(MeV);Energy(MeV)",500,0,70,500,0,12);
    TH2* cataniat                = new TH2D("cataniat", "cataniat;Energy(MeV);Energy(MeV)",500,0,70,500,0,12);
    TH2* cataniaGS               = new TH2D("cataniaGS", "cataniaGS;Energy(MeV);Energy(MeV)",500,0,70,500,0,12);
    //hitogramas for gamma spectra
    TH1D *ceadc[300]; // hitogram array declaration for DSSD
    TH1D *adcraw[300]; // hitogram array declaration for DSSD

    TH1D *labr44[numofchainelements][10]; // hitogram array declaration for labr
    TH1D *labr61[numofchainelements][10]; // hitogram array declaration for labr
    TH1D *labr69[numofchainelements][10]; // hitogram array declaration for labr
    TH1D *labrGS[numofchainelements][10]; // hitogram array declaration for labr

    TH1D *labr44c[numofchainelements][10]; // hitogram array declaration for labr
    TH1D *labr61c[numofchainelements][10]; // hitogram array declaration for labr
    TH1D *labr69c[numofchainelements][10]; // hitogram array declaration for labr
    TH1D *labrGSc[numofchainelements][10]; // hitogram array declaration for labr
    //
    // TH1D *labr44d[numofchainelements][10]; // hitogram array declaration for labr
    // TH1D *labr61d[numofchainelements][10]; // hitogram array declaration for labr
    // TH1D *labr69d[numofchainelements][10]; // hitogram array declaration for labr
    // TH1D *labrGSd[numofchainelements][10]; // hitogram array declaration for labr

    TH1D* labrEh          = new TH1D("labrEh", "labrE;Energy(MeV);Counts", 500, 0, 10);

    for (Int_t i = 0; i < numofchainelements; i++) {
      for(Int_t w=0;w<10;w++) {
        std::ostringstream name,name44,name61,name69,nameGS; // string to call them DSSD_(number)+"root file name
        //Raw data evolution in channels
        name44<<"labr44_"<<i<<"_"<<w; // add the number of the detector to the string NEED TO KEEP THIS LINE AND THE NEXT ONE DUE TO CALIBRATION
        labr44[i][w] = new TH1D(name44.str().c_str(),name44.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration
        // Initial calibrtation
        name61<<"labr61_"<<i<<"_"<<w; // add the number of the detector to the string
        labr61[i][w] = new TH1D(name61.str().c_str(),name61.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration
        //Calibrated with time as well
        name69<<"labr69_"<<i<<"_"<<w; // add the number of the detector to the string
        labr69[i][w] = new TH1D(name69.str().c_str(),name69.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration

        nameGS<<"labrGS_"<<i<<"_"<<w; // add the number of the detector to the string
        labrGS[i][w] = new TH1D(nameGS.str().c_str(),nameGS.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration

        std::ostringstream namec,name44c,name61c,name69c,nameGSc; // string to call them DSSD_(number)+"root file name
        //Raw data evolution in channels
        name44c<<"labr44c_"<<i<<"_"<<w; // add the number of the detector to the string NEED TO KEEP THIS LINE AND THE NEXT ONE DUE TO CALIBRATION
        labr44c[i][w] = new TH1D(name44c.str().c_str(),name44c.str().c_str(),4000,0.5,10);// name.str().c_str() converts the string to be use in the histogram declaration
        // Initial calibrtation
        name61c<<"labr61c_"<<i<<"_"<<w; // add the number of the detector to the string
        labr61c[i][w] = new TH1D(name61c.str().c_str(),name61c.str().c_str(),4000,0.5,10);// name.str().c_str() converts the string to be use in the histogram declaration
        //Calibrated with time as well
        name69c<<"labr69c_"<<i<<"_"<<w; // add the number of the detector to the string
        labr69c[i][w] = new TH1D(name69c.str().c_str(),name69c.str().c_str(),4000,0.5,10);// name.str().c_str() converts the string to be use in the histogram declaration
        nameGSc<<"labrGSc_"<<i<<"_"<<w; // add the number of the detector to the string
        labrGSc[i][w] = new TH1D(nameGSc.str().c_str(),nameGSc.str().c_str(),4000,0.5,10);// name.str().c_str() converts the string to be use in the histogram declaration

        // std::ostringstream name44d,name61d,name69d,nameGSd; // string to call them DSSD_(number)+"root file name
        // name44d<<"labr44d_"<<i<<"_"<<w; // add the number of the detector to the string NEED TO KEEP THIS LINE AND THE NEXT ONE DUE TO CALIBRATION
        // labr44d[i][w] = new TH1D(name44d.str().c_str(),name44d.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration
        // // Initial calibrtation
        // name61d<<"labr61d_"<<i<<"_"<<w; // add the number of the detector to the string
        // labr61d[i][w] = new TH1D(name61d.str().c_str(),name61d.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration
        // //Calibrated with time as well
        // name69d<<"labr69d_"<<i<<"_"<<w; // add the number of the detector to the string
        // labr69d[i][w] = new TH1D(name69d.str().c_str(),name69d.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration
        //
        // nameGSd<<"labrGSd_"<<i<<"_"<<w; // add the number of the detector to the string
        // labrGSd[i][w] = new TH1D(nameGSd.str().c_str(),nameGSd.str().c_str(),4096,0,4096);// name.str().c_str() converts the string to be use in the histogram declaration

      }
    }
    cout << "histo check " << endl;
// TH1D *nametestlabrh[10];
// for (Int_t i = 0; i < 10; i++) {
//   std::ostringstream nametestlabr;
//   nametestlabr<<"nametestlabr_"<<i; // add the number of the detector to the string NEED TO KEEP THIS LINE AND THE NEXT ONE DUE TO CALIBRATION
//   nametestlabrh[i] = new TH1D(nametestlabr.str().c_str(),nametestlabr.str().c_str(),1000,0.5,8);// name.str().c_str() converts the string to be use in the histogram declaration
//   // Initial calibrtation
// }

    // TH1D *tdc[10];
    // for (Int_t i = 0; i < 10; i++) {
    //   std::ostringstream nametdc;
    //   nametdc<<"TDC"<<i;
    //   tdc[i] = new TH1D(nametdc.str().c_str(),nametdc.str().c_str(),1200,0,1200);// name.str().c_str() converts the string to be use in the histogram declaration
    //   tdc[i]->GetXaxis()->SetTitle("ns");
    //   tdc[i]->GetYaxis()->SetTitle("Counts");
    // }

    Double_t Ethreshold=0.1;
    Double_t Eb= 12.22;



    ////
    // Energy at the start of the target 12.984; for 0.5 mg E= 3.2375*4 , for 0.4, 3.20365*4 ,, for 0.39 3.20468 for 0.43 3.20054, for 0.45 3.19847, for 0.52 3.19122, for 0.51 3.19225, for 0.2 3.22427
    // FOR 0.3  3.21397 for 1.39 3.09981
    Double_t Pb=sqrt(2*(4.00260)*(Eb));
    Double_t MO16= 15.9949;
    Double_t PO=sqrt(2*(MO16)*(3.195));
    Double_t Mc=12;
    Double_t Ma=4.002603;
    Double_t Qv = 0; // scattering carbon alpha
    Double_t Mal = 26.981539; //Al mass
    Double_t Mbe8 = 8.005305;
    //global counters
    Int_t numadcs=0;
    // Useful things
    Double_t pi=4.*atan(1.);
    Double_t dtr = pi/180;

    Int_t g=0; //variable for the multihistograms
    Int_t calibcounter = 0;
    Long64_t chaincounter = 0;
    // n_entries = eventcountercarr[0]/2;
    cout << "Total Events: " << n_entries << endl;
    auto date = TDatime();
    cout << "start: ";
    date.Print();
    long int a_counter=0;

    for (long int i =  0; i <= n_entries; i++) { // over all entries of the chain
      a_counter++;



      // if (chaincounter<10) {
      //   cout << chaincounter << endl;
      // }
      //Long64_t chaincounter = chain.LoadTree(i);
      chain.GetEntry(i); // this one directs the attention to the "i" entry
      // if (i<5) {
      //   // chain.Show(i);
      //   // cout << " calib counter: " <<  calibcounter
      //   //      << " event " << i
      //   //      << " chaincounter " << chaincounter<< endl;
      // }
      /// keep track of hits on detectors
      Int_t hitondl=0,
            hitonl=0,
            hitondr=0,
            hitonr=0,
            hitanode=0,
            hitdynod=0,
            hitonchunky=0,
            hitontdc=0,
            hitonlabr=0,
            labrrawec=0;

      Int_t hitonlabrd=0,
            labrrawecd=0;
      ///counters for positions and angles
      //left DE

      Int_t labr[15];
      Int_t labrrawe[15];

      Int_t labrd[15];
      Int_t labrrawed[15];

      Double_t xdlmm[16]={0};
      Double_t xdl[16]={0};
      Double_t xdle[16]={0};
      Int_t s_xdle[16]={0};
      Int_t xdlc=0,xdlec=0;

      Double_t ydlmm[16]={0};
      Double_t ydl[16]={0};
      Double_t ydle[16]={0};
      Int_t s_ydle[16]={0};
      Int_t ydlc=0,ydlec=0;
      //left E
      Double_t xlmm[16]={0};
      Double_t xl[16]={0};
      Double_t xle[16]={0};
      Int_t s_xle[16]={0};
      Int_t xlc=0,xlec=0;

      Double_t ylmm[16]={0};
      Double_t yl[16]={0};
      Double_t yle[16]={0};
      Int_t s_yle[16]={0};
      Int_t ylc=0,ylec=0;

      //right DE
      Double_t xdrmm[16]={0};
      Double_t xdr[16]={0};
      Double_t xdre[16]={0};
      Int_t s_xdre[16]={0};
      Int_t xdrc=0,xdrec=0;

      Double_t ydrmm[16]={0};
      Double_t ydr[16]={0};
      Double_t ydre[16]={0};
      Int_t s_ydre[16]={0};
      Int_t ydrc=0,ydrec=0;
      //right E
      Double_t xrmm[16]={0};
      Double_t xr[16]={0};
      Double_t xre[16]={0};
      Int_t s_xre[16]={0};
      Int_t xrc=0,xrec=0;

      Double_t yrmm[16]={0};
      Double_t yr[16]={0};
      Double_t yre[16]={0};
      Int_t s_yre[16]={0};
      Int_t yrc=0,yrec=0;
      ///chunky
      Double_t xchmm[16]={0};
      Double_t xch[16]={0};
      Double_t xche[16]={0};
      Int_t s_xche[16]={0};
      Int_t xchc=0,xchec=0;

      Double_t ychmm[16]={0};
      Double_t ych[16]={0};
      Double_t yche[16]={0};
      Int_t s_yche[16]={0};
      Int_t ychc=0,ychec=0;
      ///Labr3
      // Double_t labr[10]={0};
      Double_t labre[10]={0};
      Int_t labrc=0,labrec=0;
      // Double_t labrd[10]={0};
      Double_t labred[10]={0};
      Int_t labrcd=0,labrecd=0;
      /// Energies, angles:
      Double_t E_dl[256]={0};
      Double_t tt_dl[256]={0};
      Double_t phit_dl[256]={0};

      Double_t E_dr[256]={0};
      Double_t tt_dr[256]={0};
      Double_t phit_dr[256]={0};

      Double_t  E_r[256]={0};
      Double_t  E_l[256]={0};
      Double_t tt_l[256]={0};
      Double_t phit_l[256]={0};

      Double_t tt_r[256]={0};
      Double_t phit_r[256]={0};


      Double_t  E_ch[256]={0};
      Double_t  tdc_ch[10]={0};
      Double_t  tdc_t[10]={0};


      for (Int_t k = 0; k < emult; k++) {  //start fisrt loop over multiplicty

        if (1==1) {// turned off channels.

        Double_t realampl=0; // clean the real amplitude
        Double_t prerealampl=0;
        Int_t realchan = eadc[k];

        if (realchan<128) { // Medusa
          prerealampl = eampl[k];
          realampl = (prerealampl*clb_array[calibcounter][eadc[k]][0] ) +
                     (clb_array[calibcounter][eadc[k]][1]);

        // realampl = (pow(prerealampl,3)*clb_array[calibcounter-1][eadc[k]][3] ) +
        //           (pow(prerealampl,2)*clb_array[calibcounter-1][eadc[k]][2] ) +
        //           (prerealampl*clb_array[calibcounter-1][eadc[k]][0] ) +
        //           (clb_array[calibcounter-1][eadc[k]][1]);


        }
        if (realchan>=128  && realchan<128+10) { // Labr3
          prerealampl = eampl[k];
          // if (calibcounter==0) {
          //   cout << calibcounter << " " << eadc[k] << " "
          //                        << " " << clb_array[calibcounter][eadc[k]][0]
          //                        << " " << clb_array[calibcounter][eadc[k]][1]
          //                        << " " << clb_array[calibcounter][eadc[k]][2] << endl;
          // }

          realampl =  (pow(prerealampl,2)*clb_array[calibcounter][eadc[k]][2] ) +
                      (prerealampl*clb_array[calibcounter][eadc[k]][0] ) +
                      (clb_array[calibcounter][eadc[k]][1] );
          // nametestlabrh[realchan-128]->Fill(realampl);

          // if (realampl<0.) {
          //   cout << "red flag " << endl;
          // }
          hitonlabr++;
          labrrawe[labrrawec]=eampl[k];
          labr[labrc]=eadc[k];
          labre[labrec]=realampl;
          labrc++;labrec++;labrrawec++;
          ADCmap->Fill( eadc[k],realampl);
          if (realampl<0.) {
            // cout<<"RED FLAG!!"<<endl;
            // cout << "calib counter: " << calibcounter << endl;
            // cout << "chanel: " << eadc[k] << endl;
            // cout << "End of RED FLAG" << endl;
            myfiletxt<< calibcounter << "\t" << eadc[k] << "\t" << prerealampl << "\t"  << realampl <<endl;
          }

        //  cout << prerealampl <<  "eampl[k] " <<  eadc[k] << " " << clb_array[calibcounter-1][eadc[k]][0] <<" " <<
        //                            clb_array[calibcounter-1][eadc[k]][1] <<" " <<
        //                            clb_array[calibcounter-1][eadc[k]][2] <<" " <<
        //                            clb_array[calibcounter-1][eadc[k]][3] << endl;
        // cout <<    realampl <<  endl;
        // if (realchan== 130) {
        //   labrEh->Fill(realampl);
        // }

        }
        // if (realchan>=144  && realchan<144+10) { // Labr3
        //   prerealampl = eampl[k];
        //   hitonlabrd++;
        //   labrrawed[labrrawecd]=eampl[k];
        //   labrd[labrcd]=eadc[k]-144;
        //   realampl = 0;
        //   labred[labrecd]=realampl;
        //   // cout <<  labrd[labrcd] << " , "<< labrrawed[labrrawecd] <<endl;
        //   labrcd++;labrecd++;labrrawecd++;
        //
        //
        // }
        if (realchan>=160  && realchan<191) { // chunks
          realampl = 0;
        }
        if (realchan>200) { // chunks
          realampl = 0;
        }
        ADCmap->Fill(realchan,realampl);
        //left side
        if (eadc[k]>=0 && eadc[k]<32) {hitondl++;}
        if (eadc[k]>=32 && eadc[k]<64) {hitonl++;}
        //right side
        if (eadc[k]>=64 && eadc[k]<95) {hitondr++;}
        if (eadc[k]>=95 && eadc[k]<95+32) {hitonr++;}
        //gamma detectors, anode(signal) and dynode(time)
        if (eadc[k]>=128 && eadc[k]<=137) {hitanode++;}
        if (eadc[k]>=144 && eadc[k]<=154) {hitdynod++;}
        /// chunky

        ///tdc
        if (eadc[k] >= 512 &&  eadc[k] <= 512+9) {
          realampl = (eampl[k]*tdc_coef_p1[eadc[k]-512]) + tdc_coef_p1[eadc[k]-512];
          tdc_ch[hitontdc]=eadc[k];
          tdc_t[hitontdc]=realampl;
          hitontdc++;

        }
        /// just to now the max number of multiplicty
        if (eadc[k]>numadcs) {numadcs=eadc[k];}
        ///calibrate the gamma detectors
        // if ( eadc[k] >= 128 &&   eadc[k] <= 128+10 ) {
        //
        //
        // }

        // Positions, angles of dssds

        ///left DE
        if (eadc[k]>=0 && eadc[k]<16) {
          xdlmm[xdlc]=-position(realchan);
          xdl[xdlc]=angle_theta(-position(realchan))-(45*dtr);
          xdle[xdlc]=realampl;
          xdlc++;
        }
        if (eadc[k]>=16 && eadc[k]<32) {
          ydlmm[ydlc]=-position(realchan);
          ydl[ydlc]=angle_theta(-position(realchan))+(0*dtr);
          ydle[ydlc]=realampl;
          ydlc++;
        }
        // left E
        if (eadc[k]>=32 && eadc[k]<32+16) {
          xlmm[xlc]=-position(realchan);
          xl[xlc]=angle_theta(-position(realchan))-(45*dtr);
          xle[xlc]=realampl;
          xlc++;
        }
        if (eadc[k]>=32+16 && eadc[k]<32+16+16) {
          ylmm[ylc]=position(realchan);
          yl[ylc]=angle_theta(-position(realchan))+(0*dtr);
          yle[ylc]=realampl;
          ylc++;
        }
        ///right DE
        if (eadc[k]>=64 && eadc[k]<64+16) {
          xdrmm[xdrc]=-position(realchan);
          xdr[xdrc]=angle_theta(-position(realchan))+(45*dtr);
          xdre[xdrc]=realampl;
          xdrc++;
        }
        if (eadc[k]>=64+16 && eadc[k]<64+16+16) {
          ydrmm[ydrc]=position(realchan);
          ydr[ydrc]=angle_theta(position(realchan))+(0*dtr);
          ydre[ydrc]=realampl;
          ydrc++;
        }
        ///right E
        if (eadc[k]>=96 && eadc[k]<96+16) {
          xrmm[xrc]=-position(realchan);
          // cout << angle_theta(position(realchan))+(45*dtr) << endl;
          xr[xrc]=angle_theta(position(realchan))+(45*dtr);
          // cout<< xr[xrc]<<endl;
          xre[xrc]=realampl;
          xrc++;
        }
        if (eadc[k]>=96+16 && eadc[k]<96+16+16) {
          yrmm[yrc]=position(realchan);
          yr[yrc]=angle_theta(-position(realchan))+(0*dtr);
          yre[yrc]=realampl;
          yrc++;
        }
        //cube position and angle assingment
        ///chunky

      }

    } //  end fisrt loop over multiplicty

      //////// left DE
      if (xdlc>1) {
        TMath::Sort(xdlc,xdle,s_xdle,kTRUE);
      }
      if (ydlc>1) {
        TMath::Sort(ydlc,ydle,s_ydle,kTRUE);
      }
      //// left E
      if (xlc>1) {
        TMath::Sort(xlc,xle,s_xle,kTRUE);
      }
      if (ylc>1) {
        TMath::Sort(ylc,yle,s_yle,kTRUE);
      }
      //////// right DE
      if (xdrc>1) {
        TMath::Sort(xdrc,xdre,s_xdre,kTRUE);
      }
      if (ydrc>1) {
        TMath::Sort(ydrc,ydre,s_ydre,kTRUE);
      }
      //// right E
      if (xrc>1) {
        TMath::Sort(xrc,xre,s_xre,kTRUE);
      }
      if (yrc>1) {
        TMath::Sort(yrc,yre,s_yre,kTRUE);
      }

      if (xlc > 0 && xlc==ylc && ((xlc+ylc)%2)==0 ) {

        for (int j = 0; j < xlc;j++) {
          if (abs(xle[s_xle[j]]-yle[s_yle[j]])<0.2 &&  xle[s_xle[j]]> 0.3 && yle[s_yle[j]]> 0.3) {
            mapleftEcm->Fill(xlmm[s_xle[j]],ylmm[s_yle[j]]);
            E_l[j]=(xle[s_xle[j]]+yle[s_yle[j]])/2;
            Double_t ttb     = acos(cos(xl[s_xle[j]])*cos(yl[s_yle[j]]));
            Double_t phi_tb  = acos(sin(yl[s_yle[j]])/sin(ttb));////////new phi angle
            tt_l[j]=ttb;
            phit_l[j]=phi_tb;
        }
      }
    }
      if (xrc > 0 && xrc==yrc && ((xrc+yrc)%2)==0 ) {
        for (int j = 0; j < xrc;j++) {
          if (abs(xre[s_xre[j]]-yre[s_yre[j]])<0.2 &&  xre[s_xre[j]]> 0.3 && yre[s_yre[j]]> 0.3) {
            maprightEcm->Fill(xrmm[s_xre[j]],yrmm[s_yre[j]]);
            E_r[j]=(xre[s_xre[j]]+yre[s_yre[j]])/2;
            Double_t tta     = acos(cos(xr[s_xre[j]])*cos(yr[s_yre[j]]));
            Double_t phi_ta  = acos(sin(yr[s_yre[j]])/sin(tta));////////new phi angle
            tt_r[j]=tta;
            phit_r[j]=phi_ta;
          }
        }
      }

      if (xdlc > 0 && xdlc==ydlc && ((xdlc+ydlc)%2)==0 ) {
        for (int j = 0; j < xdlc;j++) {
          if (abs(xdle[s_xdle[j]]-ydle[s_ydle[j]])<0.2 &&  xdle[s_xdle[j]]> 1.1 && ydle[s_ydle[j]]> 1.1  ) {
            Double_t messE  =((xdle[s_xdle[j]]+ydle[s_ydle[j]])/2);
            Double_t tt     = acos(cos(xdl[s_xdle[j]])*cos(ydl[s_ydle[j]]));
            Double_t phi_t  = acos(sin(ydl[s_ydle[j]])/sin(tt));////////new phi angle
            // E_dl[j]=energyloss_table(messE,tt,0.251,energyEL,rangemg,par_indexvsE,par_indexvsR);
            E_dl[j]=messE;
            tt_dl[j]=tt;
            phit_dl[j]=phi_t-(pi/2.);
            mapleftDE->Fill(xdl[s_xdle[j]],ydl[s_ydle[j]]);
            hitmap->Fill(xdl[s_xdle[j]],ydl[s_ydle[j]]);
            mapleftDEcm->Fill(xdlmm[s_xdle[j]],ydlmm[s_ydle[j]]);
          }
        }

      }
      if (xdrc > 0 && xdrc==ydrc && ((xdrc+ydrc)%2)==0 ) {
        for (int j = 0; j < xdrc;j++) {
          if (abs(xdre[s_xdre[j]]-ydre[s_ydre[j]])<0.2 &&  xdre[s_xdre[j]]> 1.1 && ydre[s_ydre[j]]> 1.1 ) {
            Double_t messE  =((xdre[s_xdre[j]]+ydre[s_ydre[j]])/2);
            Double_t tt     = -acos(cos(xdr[s_xdre[j]])*cos(ydr[s_ydre[j]]));
            Double_t phi_t  = acos(sin(ydr[s_ydre[j]])/sin(tt));////////new phi angle
            //E_dr[j]=energyloss_table(messE,tt,0.251,energyEL,rangemg,par_indexvsE,par_indexvsR);
            E_dr[j]=messE;
            tt_dr[j]=tt;
            phit_dr[j]=phi_t-(pi/2.);
            maprightDE->Fill(xdr[s_xdre[j]],ydr[s_ydre[j]]);
            hitmap->Fill(xdr[s_xdre[j]],ydr[s_ydre[j]]);
            // scatangleall->Fill(-tt_dr[j],E_dr[j]);
            maprightDEcm->Fill(xdrmm[s_xdre[j]],ydrmm[s_ydre[j]]);
          }
        }

      }

      if (ydlc==1 && xdlc==1 && E_dl[0]>1.2  && xdrc==0 && ydrc==0  && xrc==0 && yrc==0) {
        scat_ene=energyloss_table(E_dl[0]+E_l[0],tt_dl[0],0.251,energyEL,rangemg,par_indexvsE,par_indexvsR);;
        scat_ang= -acos(cos(xdl[s_xdle[0]])*cos(ydl[s_ydle[0]]));
        scat_angphi=phit_dl[0];
        catania_single->Fill(cat_x(scat_ene,scat_ang,scat_angphi),Eb- scat_ene);
        Ex_C12_scat_single->Fill(scattering(scat_ene,scat_ang,scat_angphi,12,Eb));
        //cathit[0]=cat_x(scat_ene,scat_ang,scat_angphi);
        //cathit[1]=Eb-E_dl[0];
        //imdat->Fill();
      }
      if (ydrc==1 && xdrc==1 && E_dr[0]>1.2 && xdlc==0 && xdlec==0 && xlc==0 && xlec==0 ) {
        scat_ene= energyloss_table(E_dr[0]+E_r[0],tt_dr[0],0.251,energyEL,rangemg,par_indexvsE,par_indexvsR);
        scat_ang= acos(cos(xdr[s_xdre[0]])*cos(ydr[s_ydre[0]]));
        scat_angphi=phit_dr[0];
        catania_single->Fill(cat_x(scat_ene,scat_ang,scat_angphi),Eb- scat_ene);
        Ex_C12_scat_single->Fill(scattering(scat_ene,scat_ang,scat_angphi,12,Eb));
        //cathit[0]=cat_x(scat_ene,scat_ang,scat_angphi);
        //cathit[1]=Eb-E_dr[0];
        //imdat->Fill();
      }

if (hitontdc==1 && hitonlabr==1 && labre[0]>0.2 ) {
// if (hitontdc==1 && hitonlabr==1 && labre[0]>1.) {

  if (ydlc==1 && xdlc==1 && E_dl[0]>1.2 && xdrc==0 && ydrc==0  && xrc==0 && yrc==0) {
    scat_ene=energyloss_table(E_dl[0]+E_l[0],tt_dl[0],0.251,energyEL,rangemg,par_indexvsE,par_indexvsR);;
    scat_ang= tt_dl[0];
    scat_angphi = phit_dl[0];
    scatangle->Fill(scat_ang,scat_ene);
    hitmapsingle->Fill(xdl[s_xdle[0]],ydl[s_ydle[0]]); /// single hit map
    labr3t_ene=labre[0];
    TDC_tim=tdc_t[0];
    Double_t catx = cat_x(scat_ene,scat_ang,scat_angphi);
    Double_t caty= Eb-scat_ene;
    single_event *event1 = new single_event(labr3t_ene,scat_ene,TDC_tim,scat_ang,scat_angphi,catx,caty,calibcounter);
    event = event1;
    if (pnpoly(arrSizeGS,x_GS,y_GS,catx,caty)==1) {
      labrGS[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labrGSc[calibcounter][labr[0]-128]->Fill(labre[0]);
      // labrGSd[calibcounter][labrd[0]]->Fill(labrrawed[0]);
      labr_GS_t->Fill(labre[0]);
      cataniaGS->Fill(catx,caty);
    }
    if (pnpoly(arrSize61,x_61,y_61,catx,caty)==1 ) {
      labr61[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labr61c[calibcounter][labr[0]-128]->Fill(labre[0]);
      labr_61_t->Fill(labre[0]);
      catania61->Fill(catx,caty);
    }
    if (pnpoly(arrSize44,x_44,y_44,catx,caty)==1) {
      labr44[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labr44c[calibcounter][labr[0]-128]->Fill(labre[0]);
      labr_44_t->Fill(labre[0]);
      catania44->Fill(catx,caty);
    }
    if (pnpoly(arrSize69,x_69,y_69,catx,caty)==1  ) {
      labr69[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labr69c[calibcounter][labr[0]-128]->Fill(labre[0]);
      labr_69_t->Fill(labre[0]);
      catania69->Fill(catx,caty);
    }
    if (pnpoly(arrSizet,x_t,y_t,catx,caty)==1) {
      labr_t_t->Fill(labre[0]);
      cataniat->Fill(catx,caty);
    }
    imdat->Fill();
  }

  if (ydrc==1 && xdrc==1 && E_dr[0]>1.2 && xdlc==0 && ydlc==0 && xlc==0 && ylc==0) {
    scat_ene= energyloss_table(E_dr[0]+E_r[0],tt_dr[0],0.251,energyEL,rangemg,par_indexvsE,par_indexvsR);
    scat_ang= tt_dr[0];
    scat_angphi = phit_dr[0];
    scatangle->Fill(scat_ang,scat_ene);
    hitmapsingle->Fill(xdr[s_xdre[0]],ydr[s_ydre[0]]); /// single hit map
    labr3t_ene=labre[0];
    TDC_tim=tdc_t[0];
    Double_t catx = cat_x(scat_ene,scat_ang,scat_angphi);
    Double_t caty= Eb-scat_ene;
    single_event *event1 = new single_event(labr3t_ene,scat_ene,TDC_tim,scat_ang,scat_angphi,catx,caty,calibcounter);
    event = event1;
    if (pnpoly(arrSizeGS,x_GS,y_GS,catx,caty)==1) {
      labrGS[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labrGSc[calibcounter][labr[0]-128]->Fill(labre[0]);
      // labrGSd[calibcounter][labrd[0]]->Fill(labrrawed[0]);
      labr_GS_t->Fill(labre[0]);
      cataniaGS->Fill(catx,caty);
    }
    if (pnpoly(arrSize61,x_61,y_61,catx,caty)==1 ) {
      labr61[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labr61c[calibcounter][labr[0]-128]->Fill(labre[0]);
      labr_61_t->Fill(labre[0]);
      catania61->Fill(catx,caty);
    }
    if (pnpoly(arrSize44,x_44,y_44,catx,caty)==1) {
      labr44[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labr44c[calibcounter][labr[0]-128]->Fill(labre[0]);
      labr_44_t->Fill(labre[0]);
      catania44->Fill(catx,caty);
    }
    if (pnpoly(arrSize69,x_69,y_69,catx,caty)==1  ) {
      labr69[calibcounter][labr[0]-128]->Fill(labrrawe[0]);
      labr69c[calibcounter][labr[0]-128]->Fill(labre[0]);
      labr_69_t->Fill(labre[0]);
      catania69->Fill(catx,caty);

    }
    if (pnpoly(arrSizet,x_t,y_t,catx,caty)==1) {
      labr_t_t->Fill(labre[0]);
      cataniat->Fill(catx,caty);
    }

    imdat->Fill(); //filling the tree

  }
}


        if ( (hitondl==2 && hitonl==2 ) && abs(xdle[s_xdle[0]]-ydle[s_ydle[0]])<0.2 &&
         xdle[s_xdle[0]]> 1.2 && ydle[s_ydle[0]]> 1.2 &&
         abs(xle[s_xle[0]]-yle[s_yle[0]])<0.3 &&  xle[s_xle[0]]> 0.2 && yle[s_yle[0]]> 0.2) {
           mapthetaL->Fill(tt_l[0],tt_dl[0]);
           mapphiL->Fill(phit_l[0],phit_dl[0]);
           EDEL->Fill((xle[s_xle[0]]+yle[s_yle[0]])/2,((xdle[s_xdle[0]]+ydle[s_ydle[0]])/2));
           yEDEL->Fill(ylmm[s_yle[0]],ydlmm[s_ydle[0]]);
           xEDEL->Fill(xlmm[s_xle[0]],xdlmm[s_xdle[0]]);
           // if(abs(xlmm[s_xle[0]]-xdlmm[s_xdle[0]])>5.){ cout << "alert " << abs(xlmm[s_xle[0]]-xdlmm[s_xdle[0]]) << endl;}

        }

      //lpace holder for R
      if ( (hitondr==2 && hitonr==2 ) && abs(xdre[s_xdre[0]]-ydre[s_ydre[0]])<0.2 &&
       xdre[s_xdre[0]]> 1.2 && ydre[s_ydre[0]]> 1.2 &&
       abs(xre[s_xre[0]]-yre[s_yre[0]])<0.3 &&  xre[s_xre[0]]> 0.2 && yre[s_yre[0]]> 0.2) {
           mapthetaR->Fill(tt_r[0],tt_dr[0]);
           mapphiR->Fill(phit_r[0],phit_dr[0]);
           EDER->Fill((xre[s_xre[0]]+yre[s_yre[0]])/2,((xdre[s_xdre[0]]+ydre[s_ydre[0]])/2));
           yEDER->Fill(yrmm[s_yre[0]],ydrmm[s_ydre[0]]);
           xEDER->Fill(xrmm[s_xre[0]],xdrmm[s_xdre[0]]);
      }

      if (i==eventsuntil[calibcounter]) {
        // cout << "calibcoutner " << calibcounter <<endl;
        calibcounter++;
      }



            if (i % int(n_entries/1000) == 0) {
            cout<< YELLOW << "\rgoing good: " << ((i*1.0)/(n_entries*1.))*100.  << " % of events processed" << RESET <<flush;
            }



          }//end of loop over all events


          f.Write("",TObject::kOverwrite); // DO NOT forget to write the root file with dem new histograms
          f.Print();
          StopWatch.Stop(); // Measure time for efficiency
          cout << endl <<  "Total Real Time: "<< StopWatch.RealTime() << "s" << endl;
          cout << "total events processed: " << a_counter << endl;
        }// end program

        Double_t position(Int_t realchan){
          Int_t imap_frontl[16]={8,9,10,11,12,13,14,15,7,6,5,4,3,2,1,0};
          Int_t imap_backl[16]={7,6,5,4,3,2,1,0,8,9,10,11,12,13,14,15};

          Int_t imap_frontr[16]={8,9,10,11,12,13,14,15,7,6,5,4,3,2,1,0};
          Int_t imap_backr[16]={7,6,5,4,3,2,1,0,8,9,10,11,12,13,14,15};

          Int_t pre_pos=0;
          Double_t pos=0.;
          Double_t random=gRandom->Uniform(0,1);

          if (realchan>=0  && realchan<16){pre_pos=imap_frontl[realchan];     pos = ((50./16.)*((1.*pre_pos)+random))-25.+6.34;}//front leftDE
          if (realchan>=16 && realchan<32){pre_pos=imap_backl[realchan-16];  pos = ((50./16.)*((1.*pre_pos)+random))-25.;}//back left DE

          if (realchan>=32 && realchan<48){pre_pos=imap_frontl[realchan-32];  pos = ((50./16.)*((1.*pre_pos)+random))-25.+6.34;} //fornt left E
          if (realchan>=48 && realchan<64){pre_pos=imap_backl[realchan-48];  pos = ((50./16.)*((1.*pre_pos)+random))-25.;} //back left E

          if (realchan>=64  && realchan<64+16){pre_pos=imap_frontr[realchan-64];     pos = ((50./16.)*((1.*pre_pos)+random))-25-6.34;}// front rigth DE
          if (realchan>=64+16 && realchan<64+32){pre_pos=imap_frontr[realchan-64-16];  pos = ((50./16.)*((1.*pre_pos)+random))-25;} // back right DE

          if (realchan>=96 && realchan<96+16){pre_pos=imap_frontr[realchan-96];  pos = ((50./16.)*((1.*pre_pos)+random))-25.-6.34;} // front right E
          if (realchan>=96+16 && realchan<96+32){pre_pos=imap_frontr[realchan-96-16];  pos = ((50./16.)*((1.*pre_pos)+random))-25.;} // back right E
          //cout <<realchan <<" " << pos << endl;
          return pos ;

        }
        Double_t positionc(Int_t realchan){
          Double_t pre_pos=0;
          Double_t pos=0;
          Double_t random=gRandom->Uniform(-0.5,0.5);

          Int_t imap_front[16]={8,9,10,11,12,13,14,15,7,6,5,4,3,2,1,0};
          Int_t imap_back[16]={7,6,5,4,3,2,1,0,8,9,10,11,12,13,14,15};

          if (realchan>=160  && realchan<160+16){pre_pos=imap_front[realchan-160];     pos = ((50./16.)*(pre_pos+random))-25.;}
          if (realchan>=160+16 && realchan<160+16+16){pre_pos=imap_back[realchan-16-160];  pos = ((50./16.)*(pre_pos+random))-25.;}
          //cout << pos << endl;
          return pos ;
        }



        Double_t p_alpha(Double_t Energy){
          Double_t pA = sqrt(2.*(4.002603)*Energy);
          return pA;
        }


        Double_t angle_theta(Double_t y){
          Double_t thetay=atan(y/(45.));
          return thetay;
        }

        Double_t angle_theta_ch(Double_t y){
          Double_t thetay=atan(y/(123.5));
          return thetay;
        }




        Double_t energyloss(Double_t Energy_i, Double_t theta,Double_t mat_thick){

          Double_t thickness = mat_thick/cos(theta);

          Double_t Range_i=(0.120167)+
          (0.726836*pow(Energy_i,1))+
          (-0.155274*pow(Energy_i,2))+
          (0.110193*pow(Energy_i,3))+
          (-0.0285561*pow(Energy_i,4))+
          (0.00452181*pow(Energy_i,5))+
          (-0.0004435*pow(Energy_i,6))+
          (2.62466e-05*pow(Energy_i,7))+
          (-8.58008e-07*pow(Energy_i,8))+
          (1.18918e-08*pow(Energy_i,9));

          Double_t Range_f = Range_i + thickness;

          Double_t Energy_f= (-0.274762)+
          (1.85389*pow(Range_f,1))+
          (-0.279884*pow(Range_f,2))+
          (0.0398548*pow(Range_f,3))+
          (-0.00386446*pow(Range_f,4))+
          (0.000245444*pow(Range_f,5))+
          (-9.99131e-06*pow(Range_f,6))+
          (2.49485e-07*pow(Range_f,7))+
          (-3.45484e-09*pow(Range_f,8))+
          (2.01367e-11*pow(Range_f,9));

          return Energy_f;
        }
        Double_t cat_x(float energy, float theta,float phi){
            Double_t Eb= 12.2;
            Double_t Pb=sqrt(2*(4.00260)*(Eb));
            Double_t Ma= 4.001506;
            Double_t Mc = 12.0;

            float pA = sqrt(2*(4.002603)*energy);

            Float_t Ax=pA*sin(theta)*cos(phi);
            Float_t Ay=pA*sin(theta)*sin(phi);
            Float_t Az=pA*cos(theta);

            Float_t Bx=0;
            Float_t By=0;
            Float_t Bz=Pb;

            Float_t Cx=Bx-Ax;
            Float_t Cy=By-Ay;
            Float_t Cz=Bz-Az;

            Float_t E12C = ((Cx*Cx)+(Cy*Cy)+(Cz*Cz))/(2);
            //cout <<  E12C << endl;
            //Float_t Ea = ((Ax*Ax)+(Ay*Ay)+(Az*Az))/(2*Ma);
            return E12C;
          }

          float scattering(float energy, float theta,float phi, float mass, float Eb)
          {


            Double_t Pb=sqrt(2*(4.00260)*(Eb));
            Double_t MO16= mass;

            float pA = sqrt(2*(4.002603)*energy);

            Float_t Ax=pA*sin(theta)*cos(phi);
            Float_t Ay=pA*sin(theta)*sin(phi);
            Float_t Az=pA*cos(theta);

            Float_t Bx=0;
            Float_t By=0;
            Float_t Bz=Pb;

            Float_t Ox=Bx-Ax;
            Float_t Oy=By-Ay;
            Float_t Oz=Bz-Az;

            Float_t EO16 = ((Ox*Ox)+(Oy*Oy)+(Oz*Oz))/(2*MO16);

            Float_t Ex= Eb-EO16-energy;
            return Ex;
          }
        Double_t energyloss_table(Double_t aEnergy_i,
          Double_t atheta,
          Double_t mat_thick,
          Double_t aenergyEL[],
          Double_t arangemg[],
          Double_t apar_indexvsE[],
          Double_t apar_indexvsR[]){
            Double_t actualenergy = 0;
            //find the right index, we have a linear relation between E and indexe
            //this relation is given by Index = E (apar_indexvsE[1])+(apar_indexvsE[0]
            if (aEnergy_i<13.) {
              Double_t thickness =  mat_thick/cos(atheta);
              Int_t index = int((aEnergy_i*apar_indexvsE[1])+(apar_indexvsE[0]));
              ///intermediate value for the range
              //initial range
              Double_t av_range = (arangemg[index]+arangemg[index+1])/2;
              //  cout << "range i " <<  rangemg[index] << " av range " << av_range << " range i+1 " <<  rangemg[index+1] << endl;
              Double_t ac_thickness = av_range + thickness;
              Int_t j=index;
              do {
                j++;
              } while(arangemg[j]<ac_thickness);
              Int_t newindex= j;

              actualenergy = (aenergyEL[newindex]+aenergyEL[newindex+1])/2;
            }


            return actualenergy;

            //cout <<"actual thickness: " <<ac_thickness << endl;
            // we need to find the closest two ranges in the table of ranges.
            // to do that we grab ac_thickness and convert it into index base using the pol
            // Double_t ac_thickness_preindex=0;
            // for (Int_t g = 0; g < 9; g++) {ac_thickness_preindex+=apar_indexvsR[g]*pow(ac_thickness,g);}
            // Int_t ac_thickness_index = int(ac_thickness_preindex);
            //cout<<"index of actual thickness: " << ac_thickness_index<< endl;
            //cout<<"range "<< rangemg[ac_thickness_index] << " actual " <<  ac_thickness<<" range+1 "<< rangemg[ac_thickness_index+1] << endl;
            // now he have a headstart to look for the index in the range so we need to compare the actual range with the table values
            // but we have a center value
            // Int_t acac_thickness_index=0;
            // // cout <<  ac_thickness_index -10 << endl;
            // // cout <<  ac_thickness_index +10 << endl;
            // Int_t lowball=0,highball=0;
            // Int_t lowball_temp=0,highball_temp=0;
            // Double_t limit=0,nulllimit=0;

            //cout << "start" << endl;

            // if (ac_thickness_index+100<15000) {
            //   for (Int_t j = ac_thickness_index -30 ; j <= ac_thickness_index+30; j++) {arangemg_temp[lowball]=arangemg[j];lowball++;}
            // }
            //
            // for (Int_t d = 0; d < 70; d++) {
            //   if (arangemg_temp[70]<ac_thickness) {
            //     lowball_temp++;
            //   }
            // }
            // acac_thickness_index=lowball_temp+ac_thickness_index;
            // //cout << acac_thickness_index << endl;
            // //cout << "lowball: " <<lowball << endl;
            // //cout << "end" << endl;
            //
            // //cout<<"acrange "<< arangemg[lowball] << " actual " <<  ac_thickness<<" acrange+1 "<< arangemg[lowball+1] << endl;
            // // now we have the index of the sandwich energy

            //Double_t actualenergy = (aenergyEL[acac_thickness_index]+aenergyEL[acac_thickness_index+1])/2;
            // cout<< "Ener mess: " << aEnergy_i <<endl;
            // cout <<"thickness: " << ac_thickness <<endl;
            // cout << "Ori Ener: " << actualenergy<< endl;


            //cout << "actual energy: "<< actualenergy << endl;
            //delete aenergyEL,arangemg,apar_indexvsE,apar_indexvsR;
            //this give us a start point to look in the array and interpolate
            //between two adjacent Energies
            // for (Int_t r = 0; r < 15002; r++) {
            //   cout << energyEL[r] << endl;
            // }
            // cout <<"index "<<  index << endl;
            // cout <<"input "<<  Energy_i << endl;
            // cout <<"table energy "<< energyEL[index] << endl;
            // cout <<"table energy "<< energyEL[index+1] << endl;


          }


          /*
          how to draw

          TFile *_file0 = TFile::Open("R175_0.root")
          TCanvas *c1 = new TCanvas("c1","multipads",900,700);
          c1->Divide(2,0);
          c1->cd(1);
          mapleft->SetStats(0);
          mapleft->Draw("COLZ");
          c1->cd(2);
          mapright->Draw("COLZ");
          mapright->SetStats(0);



          TFile *_file0 = TFile::Open("tchain.root")
          alphaenergyl->Draw()
          TH1 *back = alphaenergyl->ShowBackground()
          TH1 *alphaenergylco = alphaenergyl->Clone();
          alphaenergylco->Add(back,-1)
          alphaenergylco->GetYaxis()->SetRangeUser(0,50)
          alphaenergylco->Draw()
          */
