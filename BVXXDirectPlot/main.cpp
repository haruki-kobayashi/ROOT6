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
#include <TLatex.h>
#include <TCut.h>
#include <TLegend.h>
#include <TColor.h>

#include "logon.hpp"
#include "MyPalette.hpp"
#include "ShowProgress.hpp"

void position(TCanvas *c1, TTree *tree, int *AreaParam, double *TDRange);
void position_projection(TCanvas *c1, TTree *tree, int *AreaParam, size_t entries, double *TDRange);
void angle(TCanvas *c1, TTree *tree, double angle_max, double angle_resolution);
void angle_projection(TCanvas *c1, TTree *tree, double angle_max);
void d_angle(TCanvas *c1, TTree *tree);
void d_angle_Ncut(TCanvas *c1, TTree *tree, TString da_cutX, TString da_cutY, uint8_t da_cutPH);
void d_angle_rl(TCanvas *c1, TTree *tree, double angle_max, double dlat_range, double drad_range, int face);
void ph_vph(TCanvas *c1, TTree *tree, int vph_range, int i, double interval);
void ranking1(TCanvas *c1, TTree *tree, int vph_standard);
void ranking2(TCanvas *c1, TTree *tree, int vph_standard);
void ranking3(TCanvas *c1, TTree *tree, int vph_standard);
void ranking4(TCanvas *c1, TTree *tree, int vph_standard);
void ranking5(TCanvas *c1, TTree *tree, int vph_standard);
void ranking6(TCanvas *c1, TTree *tree, int vph_standard);
void distribution_map(TCanvas *c1, TTree *tree);
void distribution_map_re(TCanvas *c1, TTree *tree);
void distribution_xy(TCanvas *c1, TTree *tree);
void ph_vph_allangle(TCanvas *c1, TTree *tree);
void ph_vph_smallangle(TCanvas *c1, TTree *tree);

int main(int argc, char* argv[]) {
    // Check if the correct number of arguments is provided
    if (argc < 3) {
        std::cerr << "\n Usage: " << argv[0] << " bvxx_file pl\n" << std::endl;
        return 1;
    }

    // Parse command line arguments
    std::string bvxxfile = argv[1];
    int pl = argv[2] ? std::stoi(argv[2]) : 0; // Default plate number is 0

    // 一通り完成したらargparseで受け取れるように変更する
    std::string Palette = argv[3] ? argv[3] : "kBird"; // Default palette is kBird
    std::string outputfile = argv[4] ? argv[4] : bvxxfile;
    std::string output = (outputfile.size() > 4 && outputfile.substr(outputfile.size() - 4) == ".pdf") 
                         ? outputfile 
                         : outputfile + ".pdf";
    double TDRange[2] = {0.0, 0.0};
    double angle_max = 6.0;
    double angle_resolution = 0.1;
    double da_cut_slope = 0.08;
    double da_cut_intercept = 0.02;
    uint8_t da_cutPH = 9;
    double dlat_range = 0.05;
    double drad_range = 1.0;
    int vph_range = 100;

    TString da_cutX = Form("(%f * (ax < 0 ? -ax : ax) + %f)", da_cut_slope, da_cut_intercept);
    TString da_cutY = Form("(%f * (ay < 0 ? -ay : ay) + %f)", da_cut_slope, da_cut_intercept);

    // Measure time taken for the process
    TStopwatch t;

    // Don't show ROOT information messages
    gErrorIgnoreLevel = kError;

    // Create TTree
    std::cout << "\nReading BVXX file and filling the TTree... " << std::endl;
    TTree* tree = new TTree("tree", "");
    TTree* subtree = new TTree("subtree", "");

    // Create branches for the TTree
    const uint8_t NumberOfImager = 72;
    uint8_t ph1, ph2;
    uint32_t ShotID1, ViewID1, ImagerID1, ShotID2, ViewID2, ImagerID2, vph1, vph2;
    double x, y, ax, ay, ax1, ay1, ax2, ay2, dax1, day1, dax2, day2, dx, dy, dz, tan;
    // tree
    const std::vector<std::pair<std::string, void*>> uint8Branches = {
        {"ph1", &ph1}, {"ph2", &ph2}
    };
    for (const auto& branch : uint8Branches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/b").c_str());
    }
    const std::vector<std::pair<std::string, void*>> uint32Branches = {
        {"ShotID1", &ShotID1}, {"ViewID1", &ViewID1}, {"ImagerID1", &ImagerID1}, 
        {"ShotID2", &ShotID2}, {"ViewID2", &ViewID2}, {"ImagerID2", &ImagerID2}, 
        {"vph1", &vph1}, {"vph2", &vph2}
    };
    for (const auto& branch : uint32Branches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/I").c_str());
    }
    const std::vector<std::pair<std::string, void*>> doubleBranches = {
        {"x", &x}, {"y", &y}, {"ax", &ax}, {"ay", &ay}, 
        {"ax1", &ax1}, {"ay1", &ay1}, {"ax2", &ax2}, {"ay2", &ay2}, 
        {"dax1", &dax1}, {"day1", &day1}, {"dax2", &dax2}, {"day2", &day2}, {"tan", &tan}
    };
    for (const auto& branch : doubleBranches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/D").c_str());
    }
    // subtree
    const std::vector<std::pair<std::string, void*>> doubleBranchesSub = {
        {"x", &x}, {"y", &y}, {"dx", &dx}, {"dy", &dy}, {"dz", &dz}
    };
    for (const auto& branch : doubleBranchesSub) {
        subtree->Branch(branch.first.c_str(), branch.second, (branch.first + "/D").c_str());
    }

    // Read the BVXX file and fill the TTree
    vxx::BvxxReader br;
    if (br.Begin(bvxxfile, pl, 0))
    {
        vxx::HashEntry h;
        vxx::base_track_t b;

        while (br.NextHashEntry(h))
        {
            while (br.NextBaseTrack(b))
            {
                ShotID1 = vxx::hts_shot_id(b.m[0].col, b.m[0].row);
                ShotID2 = vxx::hts_shot_id(b.m[1].col, b.m[1].row);
                ViewID1 = ShotID1 / NumberOfImager;
                ViewID2 = ShotID2 / NumberOfImager;
                ImagerID1 = ShotID1 % NumberOfImager;
                ImagerID2 = ShotID2 % NumberOfImager;
                x = b.x;
                y = b.y;
                ax = b.ax;
                ay = b.ay;
                ph1 = static_cast<uint8_t>(b.m[0].ph * 0.0001);
                ph2 = static_cast<uint8_t>(b.m[1].ph * 0.0001);
                vph1 = static_cast<uint32_t>(b.m[0].ph % 10000);
                vph2 = static_cast<uint32_t>(b.m[1].ph % 10000);
                ax1 = b.m[0].ax;
                ay1 = b.m[0].ay;
                ax2 = b.m[1].ax;
                ay2 = b.m[1].ay;

                dax1 = ax - ax1;
                day1 = ay - ay1;
                dax2 = ax - ax2;
                day2 = ay - ay2;
                dz = b.m[1].z - b.m[0].z;
                dx = ax * dz - (ax2 * dz * 0.5) - (ax1 * dz * 0.5); // At the center of base layer
                dy = ay * dz - (ay2 * dz * 0.5) - (ay1 * dz * 0.5); // At the center of base layer
                tan = sqrt(ax * ax + ay * ay);

                tree->Fill();
                subtree->Fill();
            }
        }
        br.End();
    }

    double elapsedtime_read = t.CpuTime();
    std::cout << "TTree created. - Elapsed " << elapsedtime_read << " [s] (CPU)" << std::endl;

    // Display information
    size_t entries = tree->GetEntriesFast();
    std::cout << "\n Input file : " << bvxxfile << std::endl;
    std::cout << " # of BT    : " << entries << " tracks" << std::endl;

    // Start plotting
    TDatime starttime;
    int year = starttime.GetYear();
    int month = starttime.GetMonth();
    int day = starttime.GetDay();
    int hour = starttime.GetHour();
    int minute = starttime.GetMinute();
    int second = starttime.GetSecond();
	TString StartTime = Form("%d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, minute, second);
    t.Start();
	std::cout << " Plot start : " << StartTime << std::endl;

	// Progress bar
	int page = 0;
	int total = 50; // total pages
	ShowProgress(page, static_cast<double>(page) / total);

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

    // Create canvas and PDF file
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());

    int MinX = tree->GetMinimum("x");
    int MaxX = tree->GetMaximum("x");
    int MinY = tree->GetMinimum("y");
    int MaxY = tree->GetMaximum("y");
    int RangeX = MaxX - MinX;
    int RangeY = MaxY - MinY;
    int LowX, UpX, LowY, UpY, bin;
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
    int AreaParam[11] = {bin, LowX, UpX, LowY, UpY, MinX, MaxX, MinY, MaxY, RangeX, RangeY};

    position(c1, tree, AreaParam, TDRange);
    c1->Print(output.c_str()); c1->Clear();
	ShowProgress(page, static_cast<double>(page) / total);

    position_projection(c1, tree, AreaParam, entries, TDRange);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("Position*");
	ShowProgress(page, static_cast<double>(page) / total);

    angle(c1, tree, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    ShowProgress(page, static_cast<double>(page) / total);

    angle_projection(c1, tree, angle_max);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("Angle*");
    ShowProgress(page, static_cast<double>(page) / total);

    d_angle(c1, tree);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
	ShowProgress(page, static_cast<double>(page) / total);

    d_angle_Ncut(c1, tree, da_cutX, da_cutY, da_cutPH);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
	ShowProgress(page, static_cast<double>(page) / total);

    d_angle_rl(c1, tree, angle_max, dlat_range, drad_range, 1);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
    ShowProgress(page, static_cast<double>(page) / total);

    d_angle_rl(c1, tree, angle_max, dlat_range, drad_range, 2);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
    ShowProgress(page, static_cast<double>(page) / total);

    for (int i = 0; i < 10; ++i) // 0.0-0.1 ~ 0.9-1.0
    {
        ph_vph(c1, tree, vph_range, i, 0.1);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*ph*");
        ShowProgress(page, static_cast<double>(page) / total);
    }

    int phvph_loop = static_cast<int>((angle_max - 0.1) * 2) + 2;
    for (int i = 2; i < phvph_loop; ++i) // 1.0-1.1 ~
    {
        ph_vph(c1, tree, vph_range, i, 0.5);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*ph*");
        ShowProgress(page, static_cast<double>(page) / total);
    }

    // Close PDF file
    c1->Print((output + "]").c_str());
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
    double elapsedtime = t.CpuTime();
	std::cout << "\n Plot end   : " << EndTime << " - Elapsed " << elapsedtime << " [s] (CPU)" << std::endl;
    std::cout << " Output     : " << output << std::endl;

    // delete c1;
    // delete tree;

    return 0;
}

void position(TCanvas *c1, TTree *tree, int *AreaParam, double *TDRange) {
    int bin  = AreaParam[0];
    int LowX = AreaParam[1];
    int UpX  = AreaParam[2];
    int LowY = AreaParam[3];
    int UpY  = AreaParam[4];

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
    if (TDRange[1] > 0.0) Position->GetZaxis()->SetRangeUser(TDRange[0], TDRange[1]);
    tree->Draw("y*0.001:x*0.001 >> Position", "", "colz");
}

void position_projection(TCanvas *c1, TTree *tree, int *AreaParam, size_t entries, double *TDRange)
{
    int bin    = AreaParam[0];
    int LowX   = AreaParam[1];
    int UpX    = AreaParam[2];
    int LowY   = AreaParam[3];
    int UpY    = AreaParam[4];
    int MinX   = AreaParam[5];
    int MaxX   = AreaParam[6];
    int MinY   = AreaParam[7];
    int MaxY   = AreaParam[8];
    int RangeX = AreaParam[9];
    int RangeY = AreaParam[10];

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
    int Xbins = ((TH2D*)Position)->GetNbinsX();
    int Ybins = ((TH2D*)Position)->GetNbinsY();
    double max_density = 0.0;
    for (int xBin = 0; xBin < Xbins; ++xBin) {
        for (int yBin = 0; yBin < Ybins; ++yBin) {
            double density = ((TH2D*)Position)->GetBinContent(xBin + 1, yBin + 1);
            if (density > 0.0) track_density->Fill(density);
            if (max_density < density) max_density = density;
        }
    }
    if (TDRange[1] == 0.0) {
        TDRange[0] = 0.0;
        TDRange[1] = max_density;
    }
    track_density->GetXaxis()->SetRangeUser(TDRange[0], TDRange[1]);
    track_density->SetFillColorAlpha(92, 0.7);
    track_density->Draw();

	int density_entries = track_density->GetEntries();
    double density_mean = track_density->GetMean();
    double density_stddev = track_density->GetStdDev();
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

void angle(TCanvas *c1, TTree *tree, double angle_max, double angle_resolution)
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

    uint32_t angle_bin = 2 / angle_resolution * angle_max;

    TString angtitle = Form("Angle;tan#it{#theta}_{x};tan#it{#theta}_{y};/(%g rad)^{2}", angle_resolution);
    TH2D* Angle = new TH2D(
        "Angle", angtitle, angle_bin, -angle_max, angle_max, angle_bin, -angle_max, angle_max
    );
    tree->Draw("ay:ax >> Angle", "", "colz");
}

void angle_projection(TCanvas *c1, TTree *tree, double angle_max)
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
        "Angle_1D", ang1Dtitle, angle_max * 10, 0.0, angle_max
    );
    Angle_1D->SetFillColorAlpha(92, 0.7);
	tree->Draw("tan>>Angle_1D");
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
        TH2D* hist = new TH2D(name, title, 100, -2.0, 2.0, 100, -0.1, 0.1);
        tree->Draw((std::string(drawExpr) + " >> " + name).c_str(), "", "goff");
        return hist;
    };

    // Plot 1: ax - ax1
    c1->cd(1);
    TH2D* axdax1 = createHistogram(
        "axdax1", 
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x1} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x1}", 
        "dax1:ax"
    );
    axdax1->Draw("colz");

    // Plot 2: ay - ay1
    c1->cd(2);
    TH2D* ayday1 = createHistogram(
        "ayday1", 
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y1} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y1}", 
        "day1:ay"
    );
    ayday1->Draw("colz");

    // Plot 3: ax - ax2
    c1->cd(3);
    TH2D* axdax2 = createHistogram(
        "axdax2", 
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x2} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x2}", 
        "dax2:ax"
    );
    axdax2->Draw("colz");

    // Plot 4: ay - ay2
    c1->cd(4);
    TH2D* ayday2 = createHistogram(
        "ayday2", 
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y2} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y2}", 
        "day2:ay"
    );
    ayday2->Draw("colz");
}

void d_angle_Ncut(TCanvas *c1, TTree *tree, TString da_cutX, TString da_cutY, uint8_t da_cutPH)
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

    // Helper lambda to create 2D histograms and apply noise cuts
    auto createNoiseCutHistogram = [&](const char* name, const char* title, const char* drawExpr, const TCut& cut) {
        TH2D* hist = new TH2D(name, title, 100, -2.0, 2.0, 100, -0.1, 0.1);
        tree->Draw((std::string(drawExpr) + " >> " + name).c_str(), cut, "goff");
        return hist;
    };

    // Plot 1: ax - ax1
    c1->cd(1);
    TCut cut_temp = Form(
        "dax2*dax2<%s*%s&&day2*day2<%s*%s&&ph1>%d&&ph2>%d", 
        da_cutX.Data(), da_cutX.Data(), 
        da_cutY.Data(), da_cutY.Data(), 
        da_cutPH, da_cutPH
    );
    TH2D* axdax1 = createNoiseCutHistogram(
        "axdax1", 
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x1} : tan#it{#theta}_{x} (Noise cut);tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x1}", 
        "dax1:ax", 
        cut_temp
    );
    axdax1->Draw("colz");

    // Plot 2: ay - ay1
    c1->cd(2);
    cut_temp = Form(
        "dax2*dax2<%s*%s&&day2*day2<%s*%s&&ph1>%d&&ph2>%d", 
        da_cutX.Data(), da_cutX.Data(), 
        da_cutY.Data(), da_cutY.Data(), 
        da_cutPH, da_cutPH
    );
    TH2D* ayday1 = createNoiseCutHistogram(
        "ayday1", 
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y1} : tan#it{#theta}_{y} (Noise cut);tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y1}", 
        "day1:ay", 
        cut_temp
    );
    ayday1->Draw("colz");

    // Plot 3: ax - ax2
    c1->cd(3);
    cut_temp = Form(
        "dax1*dax1<%s*%s&&day1*day1<%s*%s&&ph1>%d&&ph2>%d", 
        da_cutX.Data(), da_cutX.Data(), 
        da_cutY.Data(), da_cutY.Data(), 
        da_cutPH, da_cutPH
    );
    TH2D* axdax2 = createNoiseCutHistogram(
        "axdax2", 
        "tan#it{#theta}_{x} #minus tan#it{#theta}_{x2} : tan#it{#theta}_{x} (Noise cut);tan#it{#theta}_{x};tan#it{#theta}_{x} #minus tan#it{#theta}_{x2}", 
        "dax2:ax", 
        cut_temp
    );
    axdax2->Draw("colz");

    // Plot 4: ay - ay2
    c1->cd(4);
    cut_temp = Form(
        "dax1*dax1<%s*%s&&day1*day1<%s*%s&&ph1>%d&&ph2>%d", 
        da_cutX.Data(), da_cutX.Data(), 
        da_cutY.Data(), da_cutY.Data(), 
        da_cutPH, da_cutPH
    );
    TH2D* ayday2 = createNoiseCutHistogram(
        "ayday2", 
        "tan#it{#theta}_{y} #minus tan#it{#theta}_{y2} : tan#it{#theta}_{y} (Noise cut);tan#it{#theta}_{y};tan#it{#theta}_{y} #minus tan#it{#theta}_{y2}", 
        "day2:ay", 
        cut_temp
    );
    ayday2->Draw("colz");
}

void d_angle_rl(TCanvas *c1, TTree *tree, double angle_max, double dlat_range, double drad_range, int face)
{
    gStyle->SetOptStat("");
    gStyle->SetStatY(0.97);
    gStyle->SetStatX(0.95);
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.25);
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.0, "y");

    c1->Divide(1, 2);

    TString suffix = (face == 1) ? "1" : "2";

    c1->cd(1);
    TString lattitle = Form("#Deltalateral %s;tan#it{#theta};#frac{tan#it{#theta}_{x%s} #times tan#it{#theta}_{y} #plus tan#it{#theta}_{y%s} #times tan#it{#theta}_{x}}{#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}}", suffix.Data(), suffix.Data(), suffix.Data());
    TH2D* lat = new TH2D("lat", lattitle, 200, 0.0, angle_max, 200, -dlat_range, dlat_range);
    tree->Draw(Form("(ax%s*ay-ay%s*ax)/tan:tan>>lat", suffix.Data(), suffix.Data()), "", "colz");

    c1->cd(2);
    TString radtitle = Form("#Deltaradial %s;tan#it{#theta};#frac{tan#it{#theta}_{x%s} #times tan#it{#theta}_{x} #plus tan#it{#theta}_{y%s} #times tan#it{#theta}_{y}}{#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}} #minus #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}", suffix.Data(), suffix.Data(), suffix.Data());
    TH2D* rad = new TH2D("rad", radtitle, 200, 0.0, angle_max, 200, -drad_range, drad_range);
    tree->Draw(Form("(ax%s*ax+ay%s*ay)/tan-tan:tan>>rad", suffix.Data(), suffix.Data()), "", "colz");
}

void ph_vph(TCanvas *c1, TTree *tree, int vph_range, int i, double interval)
{
    gStyle->SetOptStat("em");
    gStyle->SetStatFormat("6.1f");
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.4);
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.25);
    gStyle->SetTitleOffset(1.1, "x");

    // Set angular range
    double tan_low = i * interval;
    double tan_up = i * interval + 0.1;
    TString range = Form("tan>=%.1f&&tan<%.1f", tan_low, tan_up);

    c1->Divide(3, 2);

    // Helper lambda to create and draw histograms
    auto createAndDrawHistogram = [&](int pad, const char* name, const char* title, const char* drawExpr, int bins, double xMin, double xMax) {
        c1->cd(pad);
        TString histTitle = Form("%s (%.1f #leq tan#it{#theta} < %.1f);%s;", title, tan_low, tan_up, title);
        TH1D* hist = new TH1D(name, histTitle, bins, xMin, xMax);
        hist->SetFillColorAlpha(92, 0.7);
        tree->Draw(Form("%s >> %s", drawExpr, name), range);
    };

    // Create and draw histograms
    createAndDrawHistogram(1, "phsum", "PHsum", "ph1+ph2", 32, 0.5, 32.5);
    createAndDrawHistogram(2, "ph1", "PH1", "ph1", 16, 0.5, 16.5);
    createAndDrawHistogram(3, "ph2", "PH2", "ph2", 16, 0.5, 16.5);

    gStyle->SetOptStat("em");
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.25);

    createAndDrawHistogram(4, "vphsum", "VPHsum", "vph1+vph2", 2 * vph_range, 0, 2 * vph_range);
    createAndDrawHistogram(5, "vph1", "VPH1", "vph1", vph_range, 0, vph_range);
    createAndDrawHistogram(6, "vph2", "VPH2", "vph2", vph_range, 0, vph_range);
}
