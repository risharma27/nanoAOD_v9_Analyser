
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=0){
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events"); //"Events"
  //Declare an instance of our code class
  nano9Ana m_selec;
  
   if(sample==0){
    //Add one file to chain. This is the input file.
    chain->Add("inputs/DYJetsToLL_M-50.root");
    //Set Names of outputfiles
    hstfilename = "hst_output/hst_DY.root";
    sumfilename = "sum_output/sum_DY.txt";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }
 
  
  if(sample==1){
    chain->Add("inputs/TTTo2L2Nu.root");
    hstfilename = "hst_output/hst_tt.root";
    sumfilename = "sum_output/sum_tt.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }

   if(sample==2){
    chain->Add("inputs/HiggsToZZTo4L_M-125.root");
    hstfilename = "hst_output/hst_Higgs.root";
    sumfilename = "sum_output/sum_Higgs.txt";
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }

 
  
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);

}
