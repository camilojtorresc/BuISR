// https://root.cern/doc/master/classROOT_1_1RDataFrame.html

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "TTree.h"

using namespace std;
using namespace ROOT; // RDataFrame's namespace

void RDataFrameExample()
{
    // -------------- Creating an RDataFrame --------------

    // single file -- all constructors are equivalent

    const char *fdata1 = "../../BuHI/Histos/finaltree_Bujk_AOD_HI2016_sample1.root";
    const char *fdata2 = "../../BuHI/Histos/finaltree_Bujk_AOD_HI2016_sample2.root";

    TFile *f1 = TFile::Open(fdata1);
    TTree *t1 = f1->Get<TTree>("Butree");
 
    RDataFrame d1("Butree", fdata1);
    RDataFrame d2("Butree", f1); // same as TTreeReader    
    RDataFrame d3(*t1);

    auto c1 = d1.Count();
    auto c2 = d2.Count();
    auto c3 = d3.Count();

    cout << "read single file 1: " << *c1 << "\n";
    cout << "read single file 2: " << *c2 << "\n";
    cout << "read single file 3: " << *c3 << "\n\n";

    // multiple files -- all constructors are equivalent
    TChain chain("Butree");
    chain.Add(fdata1);
    chain.Add(fdata2);
    RDataFrame d4(chain);
    
    RDataFrame d5("Butree", {fdata1, fdata2});

    vector<string> files = {fdata1, fdata2};
    RDataFrame d6("Butree", files);

    RDataFrame d7("Butree", "../BuHI/Histos/finaltree_Bujk_AOD_HI2016_sample*.root"); // the glob is passed as-is to TChain's constructor
    

    auto c4 = d4.Count();
    auto c5 = d5.Count();
    auto c6 = d6.Count();
    auto c7 = d7.Count();

    cout << "read multiple file 4: " << *c4 << "\n";
    cout << "read multiple file 5: " << *c5 << "\n";
    cout << "read multiple file 6: " << *c6 << "\n";
    cout << "read multiple file 7: " << *c7 << "\n\n";

    // Create an empty RDataFrame

    RDataFrame dem(10); // a RDF with 10 entries (and no columns/branches, for now)
    dem.Foreach([] { static int i = 0; cout << i++ << " "; }); // silly example usage: count to ten

    cout << "\n\n";

    // input csv files as RDataFrame

    //auto df = ROOT::RDF::MakeCsvDataFrame("input.csv");
    // use df as usual

    // -------------- Filling a histogram --------------

    // Fill a TH1D with the "Bupt" branch
    auto h1 = d1.Histo1D("Bupt");

    TCanvas* canv = new TCanvas("Bupt","Bupt",50,50,800,600);
    canv->cd();
    canv->Draw();
    h1->Draw();

    canv->SaveAs("Bupt.png");

    cout << "\n";

    // -------------- Applying a filter --------------
    auto c8 = d1.Filter("Bupt > 20.0").Count(); // computations booked, not run

    cout << "Number of events with pT higher than 20 GeV: " << *c8 << "\n"; // computations run here, upon first access to the result

    // Defining a cut
    auto ptCut = [](double x) { return x >= 10. && x <= 40.; }; // a C++11 lambda function checking "x > 4"
    auto c9 = d1.Filter(ptCut, {"Bupt"}).Count();

    cout << "Number of events with pT[10,40]: " << *c9 << "\n\n";

    // -------------- Defining custom columns --------------

    auto absDif = [](int x, double y) { return abs(x-y); }; //Define a cut
    auto zMean = d1.Define("NtrkCorrDiff", absDif, {"NtrkQ","NtrkQCorr"}).Mean("NtrkCorrDiff");
    std::cout << "Mean difference between Ntrk and NtrkCorr: " << *zMean << std::endl;

    // Create a new root file with two colums from an empty RDF
    RDataFrame d( *(d1.Count()) ); // an RDF that will generate entries (currently empty)

    int x= -1;
    auto d_with_columns = d.Define("x", [&x] { return x; })
                        .Define("xx", [&x] { return x*x; });

    d_with_columns.Snapshot("myNewTree", "newfile.root");

    // Create a new root file filtered with a new column

    auto new_d = d1.Filter("Bupt >= 5 && Bupt <= 40")
                .Filter("abs(Bueta)<1.5")
                .Define("Bupt_weight", "Bupt*WeightD");
                
    
    new_d.Snapshot("Butree","finaltree_sample1.root");

    
    // -------------- Running on a range of entries --------------

    // Here we store a dataframe that loops over only the first 30 entries in a variable
    auto d30 = d.Range(30);
    auto c10 = d30.Count();

    // This is how you pick all entries from 15 onwards
    auto d15on = d.Range(15, 0);
    auto c11 = d15on.Count();

    // We can specify a stride too, in this case we pick an event every 3
    auto d15each3 = d.Range(0, 15, 3); // from 0 to 15 every 3 events
    auto c12 = d15each3.Count();

    cout << "\nNumber of entries for first 30: " << *c10 << "; from 15 onwards: " << *c11 << "; event every 3: " << *c12 << "\n\n";

    // -------------- Executing multiple actions in the same event loop --------------
    auto h2 = d1.Filter("Bupt > 10").Histo1D("Bupt");
    auto h3 = d1.Histo1D("Bupt");  

    TCanvas* canv2 = new TCanvas("Bupt","Bupt",50,50,800,600);
    canv2->cd();
    canv2->Draw();

    h3->SetLineColor(2);

    h2->Draw();       // event loop is run once here
    h3->Draw("SAME"); // no need to run the event loop again  

    canv2->SaveAs("Bupt_multiple.png");

    // It is therefore good practice to declare all your transformations and actions before accessing their results, 
    // allowing RDataFrame to run the loop once and produce all results in one go.

    cout << "\n";

    // -------------- Going parallel --------------

    // Multi core tasks, dividing events fairly between cores
    ROOT::EnableImplicitMT();

    // -------------- Working with collections and object selections --------------

    // RDataFrame reads collections as the special type ROOT::RVec.
    // RVec is a container similar to std::vector (and can be used just like a std::vector) 
    // but it also offers a rich interface to operate on the array elements in a vectorised fashion, 
    // similarly to Python's NumPy arrays.

    // h is filled with all the elements of `good_pts`, for each event
    //auto h4 = d1.Define("good_pts", [](const RVec &Bupt) { return Bupt[Bupt > 7 && Bupt < 40]; })
    //       .Histo1D("good_pts");

    
    

    return;
}