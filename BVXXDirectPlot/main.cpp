// 2025.4.5 kobayashi

#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>

#include <VxxReader.h>
#include <KeywordArgs.h>
#include <Template.h>
#include <netscan_data_types_ui.h>
#include <magic_enum.hpp>

#include <TROOT.h>
#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TMath.h>
#include <TLatex.h>
#include <TCut.h>
#include <TLegend.h>
#include <TColor.h>

#include "logon.hpp"
#include "MyPalette.hpp"
#include "ShowProgress.hpp"

void position(TCanvas *c1, TTree *tree, Int_t *AreaParam);
void position_projection(TCanvas *c1, TTree *tree, Int_t *AreaParam, size_t entries);
void angle(TCanvas *c1, TTree *tree, Double_t angle_range, Double_t angle_resolution);
void angle_projection(TCanvas *c1, TTree *tree, Double_t angle_range);
void d_angle(TCanvas *c1, TTree *tree);
void d_angle_SN(TCanvas *c1, TTree *tree, TString da_cutX, TString da_cutY);
void d_angle_rl1(TCanvas *c1, TTree *tree);
void d_angle_rl2(TCanvas *c1, TTree *tree);
void ph_vph(TCanvas *c1, TTree *tree, Int_t i, Double_t interval);
void ranking1(TCanvas *c1, TTree *tree, Int_t vph_standard);
void ranking2(TCanvas *c1, TTree *tree, Int_t vph_standard);
void ranking3(TCanvas *c1, TTree *tree, Int_t vph_standard);
void ranking4(TCanvas *c1, TTree *tree, Int_t vph_standard);
void ranking5(TCanvas *c1, TTree *tree, Int_t vph_standard);
void ranking6(TCanvas *c1, TTree *tree, Int_t vph_standard);
void distribution_map(TCanvas *c1, TTree *tree);
void distribution_map_re(TCanvas *c1, TTree *tree);
void distribution_xy(TCanvas *c1, TTree *tree);
void ph_vph_allangle(TCanvas *c1, TTree *tree);
void ph_vph_smallangle(TCanvas *c1, TTree *tree);

int main(int argc, char* argv[]) {
    // Check if the correct number of arguments is provided
    if (argc < 3) {
        std::cerr << "\n Usage: " << argv[0] << " input_bvxx_path pl\n" << std::endl;
        return 1;
    }

    // Get the BVXX file path from command line arguments
    std::string bvxxpath = argv[1];
    int pl = argv[2] ? std::stoi(argv[2]) : 0; // Default plate number is 0

    // 一通り完成したら引数で受け取れるように変更する
    std::string Palette = argv[3] ? argv[3] : "kBird"; // Default palette is kBird  // argparseに変える
    Double_t angle_range = 6.0;
    Double_t angle_resolution = 0.1;
    Double_t cut_slope = 0.015;
    Double_t cut_intercept = 0.01;
    Double_t dlat_range = 0.05;
    Double_t drad_range = 1.0;

    TString da_cutX = Form("(%f * (ax < 0 ? -ax : ax) + %f)", cut_slope, cut_intercept);
    TString da_cutY = Form("(%f * (ay < 0 ? -ay : ay) + %f)", cut_slope, cut_intercept);

    // Measure time taken for the process
    TStopwatch t;

    // Read the BVXX file
    std::cout << "\nReading BVXX file... " << std::endl;
    vxx::BvxxReader br;
    std::vector<vxx::base_track_t> btvec = br.ReadAll(bvxxpath, pl, 0);
    if (btvec.empty()) {
        std::cerr << "\n Error: Failed to read any data from the input file." << std::endl;
        return 1;
    }
    Double_t elapsedtime_read = t.CpuTime();
    std::cout << "Reading BVXX file completed. - Elapsed " << elapsedtime_read << " [s] (CPU)" << std::endl;

    // Display information about the data
    size_t entries = btvec.size();
    std::cout << "\n Input file : " << bvxxpath << std::endl;
    std::cout << " # of BT    : " << entries << " tracks" << std::endl;

	// Progress bar
	int page = 0;
	const int total = 35; // total pages

    // Start plotting
    TDatime starttime;
    Int_t year = starttime.GetYear();
    Int_t month = starttime.GetMonth();
    Int_t day = starttime.GetDay();
    Int_t hour = starttime.GetHour();
    Int_t minute = starttime.GetMinute();
    Int_t second = starttime.GetSecond();
	TString StartTime = Form("%d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, minute, second);
	std::cout << " Plot start : " << StartTime << std::endl;
	ShowProgress(page, static_cast<Double_t>(page) / total);

    // Don't show ROOT information messages
    gErrorIgnoreLevel = kError;

    // Make TTree from BVXX file
    // TString tempfile = bvxxpath + "_" + Form("%d", static_cast<int>(time(nullptr))) + ".root";
    // TFile* file = new TFile(tempfile, "recreate");
    t.Start();
    TTree* tree = new TTree("tree", "");
    tree->SetDirectory(0); // Do not attach to any directory

    // Create branches for the tree
    const UInt_t NumberOfImager = 72;
    UInt_t ShotID1, ViewID1, ImagerID1, ShotID2, ViewID2, ImagerID2;
    Double_t x, y, ax, ay, ph1, ph2, ax1, ay1, ax2, ay2, z1, z2;
    const std::vector<std::pair<std::string, void*>> intBranches = {
        {"ShotID1", &ShotID1}, {"ViewID1", &ViewID1}, {"ImagerID1", &ImagerID1}, 
        {"ShotID2", &ShotID2}, {"ViewID2", &ViewID2}, {"ImagerID2", &ImagerID2}
    };
    for (const auto& branch : intBranches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/I").c_str());
    }
    const std::vector<std::pair<std::string, void*>> doubleBranches = {
        {"x", &x}, {"y", &y}, {"ax", &ax}, {"ay", &ay}, 
        {"ph1", &ph1}, {"ph2", &ph2}, {"ax1", &ax1}, {"ay1", &ay1}, 
        {"ax2", &ax2}, {"ay2", &ay2}, {"z1", &z1}, {"z2", &z2}
    };
    for (const auto& branch : doubleBranches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/D").c_str());
    }

    // Fill the tree with data
    for (const auto& bt : btvec) {
        ShotID1 = vxx::hts_shot_id(bt.m[0].col, bt.m[0].row);
        ShotID2 = vxx::hts_shot_id(bt.m[1].col, bt.m[1].row);
        ViewID1 = ShotID1 / NumberOfImager;
        ViewID2 = ShotID2 / NumberOfImager;
        ImagerID1 = ShotID1 % NumberOfImager;
        ImagerID2 = ShotID2 % NumberOfImager;
        x = bt.x;
        y = bt.y;
        ax = bt.ax;
        ay = bt.ay;
        ph1 = bt.m[0].ph;
        ph2 = bt.m[1].ph;
        ax1 = bt.m[0].ax;
        ay1 = bt.m[0].ay;
        ax2 = bt.m[1].ax;
        ay2 = bt.m[1].ay;
        z1 = bt.m[0].z;
        z2 = bt.m[1].z;

        tree->Fill();
    }

    // Set up style
    logon();

    // Set up color palette
    if (auto PaletteOpt = magic_enum::enum_cast<MyPalette::Palette>(Palette)) {
        MyPalette::SetPalette(static_cast<int>(*PaletteOpt));
    } else if (auto PaletteOpt = magic_enum::enum_cast<EColorPalette>(Palette)) {
        MyPalette::SetPalette(static_cast<int>(*PaletteOpt));
    } else {
        std::cerr << "\n Error: " << Palette << " is an unknown palette." << std::endl;
        return 1;
    }
    Float_t r1, g1, b1, r2, g2, b2, r3, g3, b3; // Define some colors
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.15))->GetRGB(r1, g1, b1);
    gROOT->GetColor(90)->SetRGB(r1, g1, b1);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.5))->GetRGB(r2, g2, b2);
    gROOT->GetColor(92)->SetRGB(r2, g2, b2);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.85))->GetRGB(r3, g3, b3);
    gROOT->GetColor(91)->SetRGB(r3, g3, b3);

    // Make canvas and PDF file
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((bvxxpath + ".pdf[").c_str());

    Int_t MinX = tree->GetMinimum("x");
    Int_t MaxX = tree->GetMaximum("x");
    Int_t MinY = tree->GetMinimum("y");
    Int_t MaxY = tree->GetMaximum("y");
    Int_t RangeX = MaxX - MinX;
    Int_t RangeY = MaxY - MinY;
    Int_t LowX, UpX, LowY, UpY, bin;
    if (RangeX >= RangeY)
    {
        LowX = MinX - 10000;
        UpX = MaxX + 10000;
        LowY = MinY - (RangeX - RangeY + 20000) * 0.5;
        UpY = MaxY + (RangeX - RangeY + 20000) * 0.5;
        bin = (RangeX + 20000) * 0.001;
    } else {
        LowX = MinX - (RangeY - RangeX + 20000) * 0.5;
        UpX = MaxX + (RangeY - RangeX + 20000) * 0.5;
        LowY = MinY - 10000;
        UpY = MaxY + 10000;
        bin = (RangeY + 20000) * 0.001;
    }
    Int_t AreaParam[11] = {bin, LowX, UpX, LowY, UpY, MinX, MaxX, MinY, MaxY, RangeX, RangeY};

    position(c1, tree, AreaParam);
    c1->Print((bvxxpath + ".pdf").c_str()); c1->Clear();
	ShowProgress(page, static_cast<Double_t>(page) / total);

    position_projection(c1, tree, AreaParam, entries);
    c1->Print((bvxxpath + ".pdf").c_str()); c1->Clear();
    gDirectory->Delete("Position*");
	ShowProgress(page, static_cast<Double_t>(page) / total);

    angle(c1, tree, angle_range, angle_resolution);
    c1->Print((bvxxpath + ".pdf").c_str()); c1->Clear();
    ShowProgress(page, static_cast<Double_t>(page) / total);

    angle_projection(c1, tree, angle_range);
    c1->Print((bvxxpath + ".pdf").c_str()); c1->Clear();
    gDirectory->Delete("Angle*");
    ShowProgress(page, static_cast<Double_t>(page) / total);

    d_angle(c1, tree);
    c1->Print((bvxxpath + ".pdf").c_str()); c1->Clear();
    gDirectory->Delete("*a*");
	ShowProgress(page, static_cast<Double_t>(page) / total);

    d_angle_SN(c1, tree, da_cutX, da_cutY);
    c1->Print((bvxxpath + ".pdf").c_str()); c1->Clear();
    gDirectory->Delete("*a*");
	ShowProgress(page, static_cast<Double_t>(page) / total);

    // Close PDF file
    c1->Print((bvxxpath + ".pdf]").c_str());
	if (page < total) page = total;
    ShowProgress(page, 1.0);

    // End plotting
    TDatime endtime;
    year = endtime.GetYear();
    month = endtime.GetMonth();
    day = endtime.GetDay();
    hour = endtime.GetHour();
    minute = endtime.GetMinute();
    second = endtime.GetSecond();
	TString EndTime = Form("%d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, minute, second);
    Double_t elapsedtime = t.CpuTime();
	std::cout << "\n Plot end   : " << EndTime << " - Elapsed " << elapsedtime << " [s] (CPU)" << std::endl;
    std::cout << " Output     : " << bvxxpath << ".pdf" << std::endl;


    delete c1;
    delete tree;

    return 0;
}

void position(TCanvas *c1, TTree *tree, Int_t *AreaParam) {
    Int_t bin  = AreaParam[0];
    Int_t LowX = AreaParam[1];
    Int_t UpX  = AreaParam[2];
    Int_t LowY = AreaParam[3];
    Int_t UpY  = AreaParam[4];

    gStyle->SetOptStat("e");
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(1.2, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    TH2D* Position = new TH2D(
        "Position", "Position;x [mm];y [mm];/mm^{2}", bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001
    );
    tree->Draw("y*0.001:x*0.001 >> Position", "", "colz");
}

void position_projection(TCanvas *c1, TTree *tree, Int_t *AreaParam, size_t entries)
{
    Int_t bin    = AreaParam[0];
    Int_t LowX   = AreaParam[1];
    Int_t UpX    = AreaParam[2];
    Int_t LowY   = AreaParam[3];
    Int_t UpY    = AreaParam[4];
    Int_t MinX   = AreaParam[5];
    Int_t MaxX   = AreaParam[6];
    Int_t MinY   = AreaParam[7];
    Int_t MaxY   = AreaParam[8];
    Int_t RangeX = AreaParam[9];
    Int_t RangeY = AreaParam[10];

    c1->Divide(2, 2);
    c1->GetPad(1)->SetRightMargin(0.235);
	c1->GetPad(1)->SetLeftMargin(0.23);
    c1->GetPad(2)->SetRightMargin(0.3);
	c1->GetPad(2)->SetLeftMargin(0.165);
    c1->GetPad(3)->SetRightMargin(0.235);
	c1->GetPad(3)->SetLeftMargin(0.23);
    c1->GetPad(4)->SetRightMargin(0.3);
	c1->GetPad(4)->SetLeftMargin(0.165);

	TH2D* Position = (TH2D*)gDirectory->Get("Position");

    c1->cd(1);
    Position->Draw("colz");

    c1->cd(2);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    Position->SetFillColor(91);
    Position->ProjectionY()->Draw("hbar");
    Position->ProjectionY()->SetTitle("");
    Position->ProjectionY()->SetStats(0);

    c1->cd(3);
    gStyle->SetStatX(0.765);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    Position->SetFillColor(90);
    Position->ProjectionX()->Draw("bar");
    Position->ProjectionX()->SetTitle("");
    Position->ProjectionX()->SetStats(0);

    c1->cd(4);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    gStyle->SetTitleOffset(1.5, "y");
    gStyle->SetOptStat("");
    TH1D* track_density = new TH1D(
        "track_density", ";Track Density [/mm^{2}];Frequency", 10000, 0, 100000
    );
    Int_t Xbins = ((TH2D*)Position)->GetNbinsX();
    Int_t Ybins = ((TH2D*)Position)->GetNbinsY();
    Double_t max_density = 0.0;
    for (int xBin = 0; xBin < Xbins; ++xBin) {
        for (int yBin = 0; yBin < Ybins; ++yBin) {
            Double_t density = ((TH2D*)Position)->GetBinContent(xBin + 1, yBin + 1);
            if (density > 0.0) track_density->Fill(density);
            if (max_density < density) max_density = density;
        }
    }
    track_density->GetXaxis()->SetRangeUser(0.0, max_density);
    track_density->SetFillColorAlpha(92, 0.7);
    track_density->Draw();

	Int_t density_entries = track_density->GetEntries();
    Double_t density_mean = track_density->GetMean();
    Double_t density_stddev = track_density->GetStdDev();
    TLegend* density_lg = new TLegend(0.67, 0.7, 0.9, 0.9);
    density_lg->SetName("density_lg");
    density_lg->SetFillStyle(0);
    density_lg->SetBorderSize(0);
    density_lg->SetTextSize(0.04);
    density_lg->AddEntry(track_density, "Track Density [/mm^{2}]", "");
    density_lg->AddEntry(track_density, Form("%d areas", density_entries), "");
    density_lg->AddEntry(track_density, Form("Mean   %.2f", density_mean), "");
    density_lg->AddEntry(track_density, Form("Std Dev   %.2f", density_stddev), "");
    density_lg->Draw();
}

void angle(TCanvas *c1, TTree *tree, Double_t angle_range, Double_t angle_resolution)
{
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetOptStat("e");
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(0.9, "y");
    gStyle->SetTitleOffset(1.6, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    UInt_t angle_bin = 2 / angle_resolution * angle_range;

    TString angtitle = Form("Angle;tan#it{#theta}_{x};tan#it{#theta}_{y};/(%g rad)^{2}", angle_resolution);
    TH2D* Angle = new TH2D(
        "Angle", angtitle, angle_bin, -angle_range, angle_range, angle_bin, -angle_range, angle_range
    );
    tree->Draw("ay:ax >> Angle", "", "colz");
}

void angle_projection(TCanvas *c1, TTree *tree, Double_t angle_range)
{
    c1->Divide(2, 2);
    c1->GetPad(1)->SetRightMargin(0.235);
	c1->GetPad(1)->SetLeftMargin(0.23);
    c1->GetPad(2)->SetRightMargin(0.3);
	c1->GetPad(2)->SetLeftMargin(0.165);
    c1->GetPad(3)->SetRightMargin(0.235);
	c1->GetPad(3)->SetLeftMargin(0.23);
    c1->GetPad(4)->SetRightMargin(0.3);
	c1->GetPad(4)->SetLeftMargin(0.165);

	TH2D* Angle = (TH2D*)gDirectory->Get("Angle");

    c1->cd(1);
    Angle->Draw("colz");

    c1->cd(2);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    Angle->SetFillColor(91);
    Angle->ProjectionY()->Draw("hbar");
    Angle->ProjectionY()->SetTitle("");
    Angle->ProjectionY()->SetStats(0);

    c1->cd(3);
    gStyle->SetStatX(0.765);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    Angle->SetFillColor(90);
    Angle->ProjectionX()->Draw("bar");
    Angle->ProjectionX()->SetTitle("");
    Angle->ProjectionX()->SetStats(0);

    c1->cd(4);
    gStyle->SetOptStat("e");
    gStyle->SetStatFormat("8.6f");
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.5, "y");
    TString ang1Dtitle = ";#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}};Frequency";
	TH1D* Angle_1D = new TH1D(
        "Angle_1D", ang1Dtitle, angle_range * 10, 0.0, angle_range
    );
    Angle_1D->SetFillColorAlpha(92, 0.7);
	tree->Draw("sqrt(ay*ay+ax*ax)>>Angle_1D");
}

void d_angle(TCanvas *c1, TTree *tree)
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.3, "y");

    c1->Divide(2, 2);
    c1->GetPad(1)->SetRightMargin(0.13);
	c1->GetPad(1)->SetLeftMargin(0.13);
    c1->GetPad(2)->SetRightMargin(0.13);
	c1->GetPad(2)->SetLeftMargin(0.13);
    c1->GetPad(3)->SetRightMargin(0.13);
	c1->GetPad(3)->SetLeftMargin(0.13);
    c1->GetPad(4)->SetRightMargin(0.13);
	c1->GetPad(4)->SetLeftMargin(0.13);

    // Helper lambda to create 2D histograms
    auto createHistogram = [&](const char* name, const char* title, const char* drawExpr) {
        TH2D* sig_hist = new TH2D(name, title, 100, -2.0, 2.0, 100, -0.1, 0.1);
        tree->Draw((std::string(drawExpr) + " >> " + name).c_str(), "", "goff");
        return sig_hist;
    };

    // Plot 1: ax - ax1
    c1->cd(1);
    TH2D* axdax1 = createHistogram(
        "axdax1",
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x1} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x1}",
        "ax-ax1:ax"
    );
    axdax1->Draw("colz");

    // Plot 2: ay - ay1
    c1->cd(2);
    TH2D* ayday1 = createHistogram(
        "ayday1",
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y1} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y1}",
        "ay-ay1:ay"
    );
    ayday1->Draw("colz");

    // Plot 3: ax - ax2
    c1->cd(3);
    TH2D* axdax2 = createHistogram(
        "axdax2",
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x2} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x2}",
        "ax-ax2:ax"
    );
    axdax2->Draw("colz");

    // Plot 4: ay - ay2
    c1->cd(4);
    TH2D* ayday2 = createHistogram(
        "ayday2",
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y2} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y2}",
        "ay-ay2:ay"
    );
    ayday2->Draw("colz");
}

void d_angle_SN(TCanvas *c1, TTree *tree, TString da_cutX, TString da_cutY)
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.3, "y");
    gStyle->SetTitleOffset(1.15, "z");

    c1->Divide(2, 2);
    c1->GetPad(1)->SetRightMargin(0.13);
	c1->GetPad(1)->SetLeftMargin(0.13);
    c1->GetPad(2)->SetRightMargin(0.13);
	c1->GetPad(2)->SetLeftMargin(0.13);
    c1->GetPad(3)->SetRightMargin(0.13);
	c1->GetPad(3)->SetLeftMargin(0.13);
    c1->GetPad(4)->SetRightMargin(0.13);
	c1->GetPad(4)->SetLeftMargin(0.13);

    // Helper lambda to create and normalize 2D histograms
    auto createAndNormalizeHistogram = [&](const char* name, const char* title, const char* drawExpr, const TCut& cut) {
        TH2D* sig_hist = new TH2D(name, title, 100, -2.0, 2.0, 100, -0.1, 0.1);
        TH2D* all_hist = new TH2D("all_hist", "", 100, -2.0, 2.0, 100, -0.1, 0.1);

        tree->Draw((std::string(drawExpr) + " >> " + name).c_str(), cut, "goff");
        tree->Draw((std::string(drawExpr) + " >> all_hist").c_str(), "", "goff");

        for (Int_t xBin = 1; xBin <= 100; ++xBin) {
            for (Int_t yBin = 1; yBin <= 100; ++yBin) {
                Double_t entries = all_hist->GetBinContent(xBin, yBin);
                Double_t sig_like = sig_hist->GetBinContent(xBin, yBin);
                if (entries != sig_like) {
                    sig_hist->SetBinContent(xBin, yBin, sig_like / (entries - sig_like));
                } else {
                    sig_hist->SetBinContent(xBin, yBin, 0.0); // Avoid division by zero
                }
            }
        }

        delete all_hist; // Clean up temporary histogram
        return sig_hist;
    };

    // Plot 1: ax - ax1
    c1->cd(1);
    TCut cut_temp = Form(
        "(ax-ax2)*(ax-ax2)<%s*%s&&(ay-ay2)*(ay-ay2)<%s*%s&&(ay-ay1)*(ay-ay1)<%s*%s",
        da_cutX.Data(), da_cutX.Data(),
        da_cutY.Data(), da_cutY.Data(),
        da_cutY.Data(), da_cutY.Data()
    );
    TH2D* axdax1 = createAndNormalizeHistogram(
        "axdax1",
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x1} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x1};S/N",
        "ax-ax1:ax",
        cut_temp
    );
    axdax1->Draw("colz");

    // Plot 2: ay - ay1
    c1->cd(2);
    cut_temp = Form(
        "(ax-ax1)*(ax-ax1)<%s*%s&&(ay-ay2)*(ay-ay2)<%s*%s&&(ax-ax2)*(ax-ax2)<%s*%s",
        da_cutX.Data(), da_cutX.Data(),
        da_cutY.Data(), da_cutY.Data(),
        da_cutX.Data(), da_cutX.Data()
    );
    TH2D* ayday1 = createAndNormalizeHistogram(
        "ayday1",
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y1} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y1};S/N",
        "ay-ay1:ay",
        cut_temp
    );
    ayday1->Draw("colz");

    // Plot 3: ax - ax2
    c1->cd(3);
    cut_temp = Form(
        "(ax-ax1)*(ax-ax1)<%s*%s&&(ay-ay2)*(ay-ay2)<%s*%s&&(ay-ay1)*(ay-ay1)<%s*%s",
        da_cutX.Data(), da_cutX.Data(),
        da_cutY.Data(), da_cutY.Data(),
        da_cutY.Data(), da_cutY.Data()
    );
    TH2D* axdax2 = createAndNormalizeHistogram(
        "axdax2",
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x2} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x2};S/N",
        "ax-ax2:ax",
        cut_temp
    );
    axdax2->Draw("colz");

    // Plot 4: ay - ay2
    c1->cd(4);
    cut_temp = Form(
        "(ax-ax1)*(ax-ax1)<%s*%s&&(ay-ay1)*(ay-ay1)<%s*%s&&(ax-ax2)*(ax-ax2)<%s*%s",
        da_cutX.Data(), da_cutX.Data(),
        da_cutY.Data(), da_cutY.Data(),
        da_cutX.Data(), da_cutX.Data()
    );
    TH2D* ayday2 = createAndNormalizeHistogram(
        "ayday2",
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y2} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y2};S/N",
        "ay-ay2:ay",
        cut_temp
    );
    ayday2->Draw("colz");
}
