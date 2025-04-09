// 2025.4.5 kobayashi

#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <variant>

#include <VxxReader.h>
#include <KeywordArgs.h>
#include <Template.h>
#include <netscan_data_types_ui.h>

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
#include <TCut.h>
#include <TLegend.h>
#include <TColor.h>

#include "logon.hpp"
#include "MyPalette.hpp"
#include "ShowProgress.hpp"

struct RankingParams {
    double tan_low;
    double tan_up;
    uint32_t range_max;
    double lin_max;
};

void position(TCanvas *c1, TTree *tree, const double *AreaParam, const double *TDRange) noexcept;
void position_projection(TCanvas *c1, TTree *tree, const size_t entries, const double *TDRange) noexcept;
void angle(TCanvas *c1, TTree *tree, const double angle_max, const double angle_resolution) noexcept;
void angle_projection(TCanvas *c1, TTree *tree, const double angle_max, const double angle_resolution) noexcept;
void d_angle(TCanvas *c1, TTree *tree) noexcept;
void d_angle_Ncut(TCanvas *c1, TTree *tree, const TString da_cutX, const TString da_cutY, const uint8_t da_cutPH) noexcept;
void d_angle_rl(TCanvas *c1, TTree *tree, const double angle_max, const double dlat_range, const double drad_range, const int face) noexcept;
void ph_vph(TCanvas *c1, TTree *tree, const uint32_t vph_range, const uint8_t i, const double interval) noexcept;
void ranking(TCanvas *c1, TTree *tree, const double lat_range, const RankingParams (&params)[3]) noexcept;
void distribution_map(TCanvas *c1, TTree *subtree, const double *AreaParam, const double bin_dmap, const float pitch, const double markersize) noexcept;
void distribution_map_re(TCanvas *c1);
void distribution_xy(TCanvas *c1, TF1 *gaus);
void phvph_2D(TCanvas *c1, TTree *tree, const uint32_t vph_range, const double angle_max, const double angle_resolution) noexcept;

int main(int argc, char* argv[])
{
    // Check if the correct number of arguments is provided
    if (argc < 3) {
        std::cerr << "\n Usage: " << argv[0] << " bvxx_file pl\n" << std::endl;
        return 1;
    }

    // Parse command line arguments
    const std::string bvxxfile = argv[1];
    const int pl = argv[2] ? std::stoi(argv[2]) : 0; // Default plate number is 0

    // 一通り完成したらargparseで受け取れるように変更する
    const std::variant<int, std::string, MyPalette::Palette, EColorPalette> Palette = argv[3] ? argv[3] : "kBird"; // Default palette is kBird
    const std::string outputfile = argv[4] ? argv[4] : bvxxfile;
    const std::string output = (outputfile.size() > 4 && outputfile.substr(outputfile.size() - 4) == ".pdf") 
                                ? outputfile 
                                : outputfile + ".pdf";
    const double TDRange[2] = {0.0, 0.0};
    const double angle_max = 6.0;
    const double angle_resolution = 0.1;
    const double da_cut_slope = 0.08;
    const double da_cut_intercept = 0.02;
    const uint8_t da_cutPH = 9;
    const double dlat_range = 0.05;
    const double drad_range = 1.0;
    const uint32_t vph_range = 100;
    const uint32_t ranking_range_min = 30;
    const double ranking_lat_range = 0.05;

    const TString da_cutX = Form("(%f * (ax < 0 ? -ax : ax) + %f)", da_cut_slope, da_cut_intercept);
    const TString da_cutY = Form("(%f * (ay < 0 ? -ay : ay) + %f)", da_cut_slope, da_cut_intercept);

    // Measure time taken for the process
    TStopwatch t;

    // Don't show ROOT information messages
    gErrorIgnoreLevel = kError;

    // Create TTree
    std::cout << "\nReading BVXX file and filling the TTree... " << std::endl;
    TTree* tree = new TTree("tree", "");
    TTree* subtree = new TTree("subtree", "");

    // Create branches for the TTree
    constexpr uint8_t NumberOfImager = 72;
    uint8_t ph1, ph2;
    uint32_t ShotID1, ViewID1, ImagerID1, ShotID2, ViewID2, ImagerID2, vph1, vph2;
    double x, y, ax, ay, ax1, ay1, ax2, ay2, dax1, day1, dax2, day2, dx, dy, dz, tan, lin, linl;
    // tree
    const std::vector<std::pair<std::string, void*>> uint8Branches = {
        {"ph1", &ph1}, {"ph2", &ph2}
    };
    for (const auto& branch : uint8Branches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/b").c_str());
    }
    const std::vector<std::pair<std::string, void*>> uint32Branches = {
        {"ViewID1", &ViewID1}, {"ImagerID1", &ImagerID1}, 
        {"ViewID2", &ViewID2}, {"ImagerID2", &ImagerID2}, 
        {"vph1", &vph1}, {"vph2", &vph2}
    };
    for (const auto& branch : uint32Branches) {
        tree->Branch(branch.first.c_str(), branch.second, (branch.first + "/I").c_str());
    }
    const std::vector<std::pair<std::string, void*>> doubleBranches = {
        {"x", &x}, {"y", &y}, {"ax", &ax}, {"ay", &ay}, 
        {"ax1", &ax1}, {"ay1", &ay1}, {"ax2", &ax2}, {"ay2", &ay2}, 
        {"dax1", &dax1}, {"day1", &day1}, {"dax2", &dax2}, {"day2", &day2}, 
        {"tan", &tan}, {"lin", &lin}, {"linl", &linl}
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
                lin = sqrt(dax1*dax1 + day1*day1 + dax2*dax2 + day2*day2);
                linl = sqrt(((ax*ay1-ay*ax1)/tan)*((ax*ay1-ay*ax1)/tan)+((ax*ay2-ay*ax2)/tan)*((ax*ay2-ay*ax2)/tan));

                tree->Fill();

                if ((ph1+ph2)>24 && tan>1.0 && tan<1.1) {
                    subtree->Fill();
                }
            }
        }
        br.End();
    }

    const double elapsedtime_read = t.CpuTime();
    std::cout << "TTree created. - Elapsed " << elapsedtime_read << " [s] (CPU)" << std::endl;

    // Display information
    const size_t entries = tree->GetEntriesFast();
    std::cout << "\n Input file : " << bvxxfile << std::endl;
    std::cout << " # of BT    : " << entries << " tracks" << std::endl;

    // Start plotting
    const TDatime starttime;
    uint32_t year = starttime.GetYear();
    uint8_t month = starttime.GetMonth();
    uint8_t day = starttime.GetDay();
    uint8_t hour = starttime.GetHour();
    uint8_t minute = starttime.GetMinute();
    uint8_t second = starttime.GetSecond();
	const TString StartTime = Form("%d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, minute, second);
    t.Start();
	std::cout << " Plot start : " << StartTime << std::endl;

	// Progress bar
	int page = 0;
	const int total = 50; // total pages
	ShowProgress(page, static_cast<double>(page) / total);

    // Set up style
    logon();

    // Set up color palette
    MyPalette::SetPalette(Palette);
    Float_t r1, g1, b1, r2, g2, b2, r3, g3, b3; // Define some colors
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.15))->GetRGB(r1, g1, b1);
    gROOT->GetColor(90)->SetRGB(r1, g1, b1);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.85))->GetRGB(r3, g3, b3);
    gROOT->GetColor(91)->SetRGB(r3, g3, b3);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.5))->GetRGB(r2, g2, b2);
    gROOT->GetColor(92)->SetRGB(r2, g2, b2);
    gROOT->GetColor(93)->SetRGB(r2 * 0.6, g2 * 0.6, b2 * 0.6);
    gStyle->SetHistFillColor(92);
    gStyle->SetHistLineColor(93);

    // Create canvas and PDF file
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());

    const int MinX = tree->GetMinimum("x");
    const int MaxX = tree->GetMaximum("x");
    const int MinY = tree->GetMinimum("y");
    const int MaxY = tree->GetMaximum("y");
    const int RangeX = MaxX - MinX;
    const int RangeY = MaxY - MinY;
    double LowX, UpX, LowY, UpY, bin, bin_dmap;
    float pitch;
	double markersize;
    if (RangeX >= RangeY) {
        pitch = 5.0;  // 5.0mm pitch
        LowX = MinX - 10000;
        UpX = MaxX + 10000;
        LowY = MinY - (RangeX - RangeY + 20000) * 0.5;
        UpY = MaxY + (RangeX - RangeY + 20000) * 0.5;
        bin = (RangeX + 20000) * 0.001;
        bin_dmap = (RangeX + 20000) * 0.0002;
		markersize = 17400 * std::pow(static_cast<double>(RangeX), -0.85) + 0.02;
        if (RangeX < 100000) {
            pitch = 1.0;  // 1.0mm pitch
            bin_dmap *= 5;
            markersize *= 0.2;
        } else if (RangeX < 150000) {
            pitch = 2.5;  // 2.5mm pitch
            bin_dmap *= 2;
            markersize *= 0.5;
        }
    } else {
        pitch = 5.0;  // 5.0mm pitch
        LowX = MinX - (RangeY - RangeX + 20000) * 0.5;
        UpX = MaxX + (RangeY - RangeX + 20000) * 0.5;
        LowY = MinY - 10000;
        UpY = MaxY + 10000;
        bin = (RangeY + 20000) * 0.001;
        bin_dmap = (RangeY + 20000) * 0.0002;  // 5.0mm pitch
		markersize = 17400 * std::pow(static_cast<double>(RangeY), -0.85) + 0.02;
        if (RangeY < 100000) {
            pitch = 1.0;  // 1.0mm pitch
            bin_dmap *= 5;
            markersize *= 0.2;
        } else if (RangeY < 150000) {
            pitch = 2.5;  // 2.5mm pitch
            bin_dmap *= 2;
            markersize *= 0.5;
        }
    }
    const double AreaParam[5] = {bin, LowX, UpX, LowY, UpY};

    position(c1, tree, AreaParam, TDRange);
    c1->Print(output.c_str()); c1->Clear();
	ShowProgress(page, static_cast<double>(page) / total);

    position_projection(c1, tree, entries, TDRange);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("position_2D*");
	ShowProgress(page, static_cast<double>(page) / total);

    angle(c1, tree, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    ShowProgress(page, static_cast<double>(page) / total);

    angle_projection(c1, tree, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("angle_2D*");
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

    phvph_2D(c1, tree, vph_range, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*ph*");
    ShowProgress(page, static_cast<double>(page) / total);

    for (uint8_t i = 0; i < 10; ++i) // 0.0-0.1 ~ 0.9-1.0
    {
        ph_vph(c1, tree, vph_range, i, 0.1);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*ph*");
        ShowProgress(page, static_cast<double>(page) / total);
    }

    const uint8_t phvph_loop = static_cast<uint8_t>((angle_max - 0.1) * 2) + 2;
    for (uint8_t i = 2; i < phvph_loop; ++i) // 1.0-1.1 ~
    {
        ph_vph(c1, tree, vph_range, i, 0.5);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*ph*");
        ShowProgress(page, static_cast<double>(page) / total);
    }

    // Set VPH range of track ranking plot
    uint32_t range = 200;
    uint32_t range_min = 65;
    uint32_t range_max = 260;
    uint32_t vph_entries = 10000;
    uint32_t vph_mean = 0;
    uint32_t vph_sigma = 40;
    uint8_t cut_type = 0;
    TH1D* vph_temp = new TH1D("vph_temp", "vph_temp", 51, 15, 270);
    TCut cutvph_temp;

    do {
        range_min -= 5;
        if (range_min < 60) {
            vph_entries = vph_temp->GetEntries();
            vph_sigma = vph_temp->GetRMS();
            range = vph_sigma * 3;
            if (vph_entries < 1000) cut_type = 1;
        }

        range_max = range_min + range;

        if (cut_type == 0) {
            cutvph_temp = Form("(vph1+vph2)>%d && (vph1+vph2)<%d && tan>1.5 && tan<1.51 && lin>0.03 && lin<0.05", range_min, range_max);
        } else {
            cutvph_temp = Form("(vph1+vph2)>%d && (vph1+vph2)<%d && tan>1.5 && tan<1.51 && lin>0.08 && lin<0.10", range_min, range_max);
        }
        tree->Draw("(vph1+vph2)>>vph_temp", cutvph_temp, "goff");
        vph_mean = vph_temp->GetMean();
        if (range_min < 15) break;
    } while (vph_mean < range_min + 0.5 * vph_sigma);

    TF1* gaus = new TF1("gaus", "gaus", DBL_MIN, DBL_MAX);
    vph_temp->Fit("gaus", "q", "", range_min, vph_mean + vph_sigma); c1->Clear();
    uint32_t vph_standard = gaus->GetParameter(1) + 5;
    if (vph_standard < ranking_range_min) vph_standard = ranking_range_min;
    gDirectory->Delete("*temp");

    ranking(c1, tree, ranking_lat_range, {
        {0.0, 0.1, vph_standard + 120, 0.1}, 
        {0.1, 0.2, vph_standard + 70, 0.1}, 
        {0.2, 0.3, vph_standard + 45, 0.1}
    });
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("rank*");
	ShowProgress(page, static_cast<double>(page) / total);

    ranking(c1, tree, ranking_lat_range, {
        {0.3, 0.4, vph_standard + 25, 0.14}, 
        {0.4, 0.5, vph_standard + 20, 0.15}, 
        {0.5, 0.6, vph_standard + 15, 0.15}
    });
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("rank*");
	ShowProgress(page, static_cast<double>(page) / total);

    ranking(c1, tree, ranking_lat_range, {
        {0.6, 0.7, vph_standard + 10, 0.2}, 
        {0.7, 0.8, vph_standard, 0.2}, 
        {0.8, 0.9, vph_standard, 0.2}
    });
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("rank*");
	ShowProgress(page, static_cast<double>(page) / total);

    ranking(c1, tree, ranking_lat_range, {
        {0.9, 1.0, vph_standard, 0.25}, 
        {1.0, 1.1, vph_standard, 0.3}, 
        {1.1, 1.3, vph_standard, 0.35}
    });
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("rank*");
	ShowProgress(page, static_cast<double>(page) / total);

    ranking(c1, tree, ranking_lat_range, {
        {1.3, 1.5, vph_standard, 0.35}, 
        {1.5, 2.0, vph_standard, 0.5}, 
        {2.0, 2.5, vph_standard, 0.6}
    });
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("rank*");
	ShowProgress(page, static_cast<double>(page) / total);

    ranking(c1, tree, ranking_lat_range, {
        {2.5, 3.0, vph_standard, 0.7}, 
        {3.0, 4.0, vph_standard, 0.9}, 
        {4.0, 5.0, vph_standard, 1.0}
    });
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("rank*");
	ShowProgress(page, static_cast<double>(page) / total);

    distribution_map(c1, subtree, AreaParam, bin_dmap, pitch, markersize);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*temp");
    ShowProgress(page, static_cast<double>(page) / total);

    MyPalette::SetPalette(MyPalette::Palette::kRedWhiteBlue);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.5))->GetRGB(r2, g2, b2);
    gROOT->GetColor(92)->SetRGB(r2, g2, b2);
    gROOT->GetColor(93)->SetRGB(r2 * 0.6, g2 * 0.6, b2 * 0.6);

    distribution_map_re(c1);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("dz*");
	ShowProgress(page, static_cast<double>(page) / total);

    distribution_xy(c1, gaus);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("d*");
    delete gaus;
	ShowProgress(page, static_cast<double>(page) / total);

    MyPalette::SetPalette(Palette);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.5))->GetRGB(r2, g2, b2);
    gROOT->GetColor(92)->SetRGB(r2, g2, b2);
    gROOT->GetColor(93)->SetRGB(r2 * 0.6, g2 * 0.6, b2 * 0.6);

    // Close PDF file
    c1->Print((output + "]").c_str());
	if (page < total) page = total;
    ShowProgress(page, 1.0);

    // End plotting
    const TDatime endtime;
    year = endtime.GetYear();
    month = endtime.GetMonth();
    day = endtime.GetDay();
    hour = endtime.GetHour();
    minute = endtime.GetMinute();
    second = endtime.GetSecond();
	const TString EndTime = Form("%d-%02d-%02d %02d:%02d:%02d", year, month, day, hour, minute, second);
    double elapsedtime = t.CpuTime();
	std::cout << "\n Plot end   : " << EndTime << " - Elapsed " << elapsedtime << " [s] (CPU)" << std::endl;
    std::cout << " Output     : " << output << std::endl;

    // delete c1;
    // delete tree;

    return 0;
}

void position(TCanvas *c1, TTree *tree, const double *AreaParam, const double *TDRange) noexcept
{
    uint32_t bin = static_cast<uint32_t>(AreaParam[0]);
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

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

    TH2D* position_2D = new TH2D(
        "position_2D", "Position;x [mm];y [mm];/mm^{2}", bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001
    );
    if (TDRange[1] > 0.0) position_2D->GetZaxis()->SetRangeUser(TDRange[0], TDRange[1]);
    tree->Draw("y*0.001:x*0.001 >> position_2D", "", "colz");
}

void position_projection(TCanvas *c1, TTree *tree, const size_t entries, const double *TDRange) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

	TH2D* position_2D = (TH2D*)gDirectory->Get("position_2D");

    c1->cd(1);
    position_2D->Draw("colz");

    c1->cd(2);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    position_2D->SetFillColor(91);
    position_2D->ProjectionY()->Draw("hbar");
    position_2D->ProjectionY()->SetTitle("");
    position_2D->ProjectionY()->SetStats(0);

    c1->cd(3);
    gStyle->SetStatX(0.765);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    position_2D->SetFillColor(90);
    position_2D->ProjectionX()->Draw("bar");
    position_2D->ProjectionX()->SetTitle("");
    position_2D->ProjectionX()->SetStats(0);

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
    int Xbins = ((TH2D*)position_2D)->GetNbinsX();
    int Ybins = ((TH2D*)position_2D)->GetNbinsY();
    double min_density = 0.0;
    double max_density = 0.0;
    for (int xBin = 0; xBin < Xbins; ++xBin) {
        for (int yBin = 0; yBin < Ybins; ++yBin) {
            double density = ((TH2D*)position_2D)->GetBinContent(xBin + 1, yBin + 1);
            if (density > 0.0) track_density->Fill(density);
            if (max_density < density) max_density = density;
        }
    }
    if (TDRange[1] != 0.0) {
        min_density = TDRange[0];
        max_density = TDRange[1];
    }
    track_density->GetXaxis()->SetRangeUser(min_density, max_density);
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

void angle(TCanvas *c1, TTree *tree, const double angle_max, const double angle_resolution) noexcept
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
    TH2D* angle_2D = new TH2D(
        "angle_2D", angtitle, angle_bin, -angle_max, angle_max, angle_bin, -angle_max, angle_max
    );
    tree->Draw("ay:ax >> angle_2D", "", "colz");
}

void angle_projection(TCanvas *c1, TTree *tree, const double angle_max, const double angle_resolution) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

	TH2D* angle_2D = (TH2D*)gDirectory->Get("angle_2D");

    c1->cd(1);
    angle_2D->Draw("colz");

    c1->cd(2);
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    angle_2D->SetFillColor(91);
    angle_2D->ProjectionY()->Draw("hbar");
    angle_2D->ProjectionY()->SetTitle("");
    angle_2D->ProjectionY()->SetStats(0);

    c1->cd(3);
    gStyle->SetStatX(0.765);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
    angle_2D->SetFillColor(90);
    angle_2D->ProjectionX()->Draw("bar");
    angle_2D->ProjectionX()->SetTitle("");
    angle_2D->ProjectionX()->SetStats(0);

    c1->cd(4);
    gStyle->SetOptStat("e");
    gStyle->SetStatFormat("8.6f");
    gStyle->SetStatX(0.7);
    gStyle->SetStatY(0.97);
    gStyle->SetStatW(0.25);
    gStyle->SetStatH(0.17);
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.5, "y");

    uint32_t angle_bin = angle_max / angle_resolution;

    TString ang1Dtitle = ";#sqrt{tan^{2}#it{#theta}_{x}#plus tan^{2}#it{#theta}_{y}};Frequency";
	TH1D* angle_1D = new TH1D("angle_1D", ang1Dtitle, angle_bin, 0.0, angle_max);
    angle_1D->SetFillColorAlpha(92, 0.7);
	tree->Draw("tan>>angle_1D");
}

void d_angle(TCanvas *c1, TTree *tree) noexcept
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.3, "y");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.13);
        c1->GetPad(pad)->SetLeftMargin(0.13);
    }

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
        "tan#it{#theta}_{x}#minus tan#it{#theta}_{x1} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x}#minus tan#it{#theta}_{x1}", 
        "dax1:ax"
    );
    axdax1->Draw("colz");

    // Plot 2: ay - ay1
    c1->cd(2);
    TH2D* ayday1 = createHistogram(
        "ayday1", 
        "tan#it{#theta}_{y}#minus tan#it{#theta}_{y1} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y}#minus tan#it{#theta}_{y1}", 
        "day1:ay"
    );
    ayday1->Draw("colz");

    // Plot 3: ax - ax2
    c1->cd(3);
    TH2D* axdax2 = createHistogram(
        "axdax2", 
        "tan#it{#theta}_{x}#minus tan#it{#theta}_{x2} : tan#it{#theta}_{x};tan#it{#theta}_{x};tan#it{#theta}_{x}#minus tan#it{#theta}_{x2}", 
        "dax2:ax"
    );
    axdax2->Draw("colz");

    // Plot 4: ay - ay2
    c1->cd(4);
    TH2D* ayday2 = createHistogram(
        "ayday2", 
        "tan#it{#theta}_{y}#minus tan#it{#theta}_{y2} : tan#it{#theta}_{y};tan#it{#theta}_{y};tan#it{#theta}_{y}#minus tan#it{#theta}_{y2}", 
        "day2:ay"
    );
    ayday2->Draw("colz");
}

void d_angle_Ncut(TCanvas *c1, TTree *tree, const TString da_cutX, const TString da_cutY, const uint8_t da_cutPH) noexcept
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.3, "y");
    gStyle->SetTitleOffset(1.15, "z");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.13);
        c1->GetPad(pad)->SetLeftMargin(0.13);
    }

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
        "tan#it{#theta}_{x}#minus tan#it{#theta}_{x1} : tan#it{#theta}_{x} (Noise cut);tan#it{#theta}_{x};tan#it{#theta}_{x}#minus tan#it{#theta}_{x1}", 
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
        "tan#it{#theta}_{y}#minus tan#it{#theta}_{y1} : tan#it{#theta}_{y} (Noise cut);tan#it{#theta}_{y};tan#it{#theta}_{y}#minus tan#it{#theta}_{y1}", 
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
        "tan#it{#theta}_{x}#minus tan#it{#theta}_{x2} : tan#it{#theta}_{x} (Noise cut);tan#it{#theta}_{x};tan#it{#theta}_{x}#minus tan#it{#theta}_{x2}", 
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
        "tan#it{#theta}_{y}#minus tan#it{#theta}_{y2} : tan#it{#theta}_{y} (Noise cut);tan#it{#theta}_{y};tan#it{#theta}_{y}#minus tan#it{#theta}_{y2}", 
        "day2:ay", 
        cut_temp
    );
    ayday2->Draw("colz");
}

void d_angle_rl(TCanvas *c1, TTree *tree, const double angle_max, const double dlat_range, const double drad_range, const int face) noexcept
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
    TString lattitle = Form("#Deltalateral %s;tan#it{#theta};#frac{tan#it{#theta}_{x%s}#times tan#it{#theta}_{y}#plus tan#it{#theta}_{y%s}#times tan#it{#theta}_{x}}{#sqrt{tan^{2}#it{#theta}_{x}#plus tan^{2}#it{#theta}_{y}}}", suffix.Data(), suffix.Data(), suffix.Data());
    TH2D* lat = new TH2D("lat", lattitle, 200, 0.0, angle_max, 200, -dlat_range, dlat_range);
    tree->Draw(Form("(ax%s*ay-ay%s*ax)/tan:tan>>lat", suffix.Data(), suffix.Data()), "", "colz");

    c1->cd(2);
    TString radtitle = Form("#Deltaradial %s;tan#it{#theta};#frac{tan#it{#theta}_{x%s}#times tan#it{#theta}_{x}#plus tan#it{#theta}_{y%s}#times tan#it{#theta}_{y}}{#sqrt{tan^{2}#it{#theta}_{x}#plus tan^{2}#it{#theta}_{y}}}#minus #sqrt{tan^{2}#it{#theta}_{x}#plus tan^{2}#it{#theta}_{y}}", suffix.Data(), suffix.Data(), suffix.Data());
    TH2D* rad = new TH2D("rad", radtitle, 200, 0.0, angle_max, 200, -drad_range, drad_range);
    tree->Draw(Form("(ax%s*ax+ay%s*ay)/tan-tan:tan>>rad", suffix.Data(), suffix.Data()), "", "colz");
}

void phvph_2D(TCanvas *c1, TTree *tree, const uint32_t vph_range, const double angle_max, const double angle_resolution) noexcept
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.5, "y");

    c1->Divide(3, 2);
    for (int pad = 1; pad <= 6; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
    }

    uint32_t phvph_bin = static_cast<uint32_t>(angle_max / angle_resolution);

    auto createAndDraw2DHistogram = [&](int pad, const char* name, const char* title, const char* drawExpr, int yBins, double yMin, double yMax, int NDiv) {
        c1->cd(pad);
        TH2D* hist = new TH2D(name, title, phvph_bin, 0.0, angle_max, yBins, yMin, yMax);
        tree->Draw(Form("%s:tan >> %s", drawExpr, name), "", "colz");
        hist->GetYaxis()->SetNdivisions(NDiv);
    };

    createAndDraw2DHistogram(1, "phs0", "PHsum;tan#it{#theta};PHsum", "(ph1+ph2)", 27, 5.5, 32.5, 14);
    createAndDraw2DHistogram(2, "phs1", "PH1;tan#it{#theta};PH1", "ph1", 11, 5.5, 16.5, 11);
    createAndDraw2DHistogram(3, "phs2", "PH2;tan#it{#theta};PH2", "ph2", 11, 5.5, 16.5, 11);
    createAndDraw2DHistogram(4, "vphs0", "VPHsum;tan#it{#theta};VPHsum", "(vph1+vph2)", 2 * vph_range, 0.5, 2 * vph_range + 0.5, 10);
    createAndDraw2DHistogram(5, "vphs1", "VPH1;tan#it{#theta};VPH1", "vph1", vph_range, 0.5, vph_range + 0.5, 10);
    createAndDraw2DHistogram(6, "vphs2", "VPH2;tan#it{#theta};VPH2", "vph2", vph_range, 0.5, vph_range + 0.5, 10);
}

void ph_vph(TCanvas *c1, TTree *tree, const uint32_t vph_range, const uint8_t i, const double interval) noexcept
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
        TString histTitle = Form("%s (%.1f#leq tan#it{#theta} < %.1f);%s;", title, tan_low, tan_up, title);
        TH1D* hist = new TH1D(name, histTitle, bins, xMin, xMax);
        hist->SetFillColorAlpha(92, 0.7);
        tree->Draw(Form("%s >> %s", drawExpr, name), range);
    };

    // Create and draw histograms
    createAndDrawHistogram(1, "phsum", "PHsum", "ph1+ph2", 32, 0.5, 32.5);
    createAndDrawHistogram(2, "ph1", "PH1", "ph1", 16, 0.5, 16.5);
    createAndDrawHistogram(3, "ph2", "PH2", "ph2", 16, 0.5, 16.5);

    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);

    createAndDrawHistogram(4, "vphsum", "VPHsum", "vph1+vph2", 2 * vph_range, 0, 2 * vph_range);
    createAndDrawHistogram(5, "vph1", "VPH1", "vph1", vph_range, 0, vph_range);
    createAndDrawHistogram(6, "vph2", "VPH2", "vph2", vph_range, 0, vph_range);
}

void ranking(TCanvas *c1, TTree *tree, const double lat_range, const RankingParams (&params)[3]) noexcept
{
    gStyle->SetOptStat("e");
    gStyle->SetStatFormat("6.2f");
    gStyle->SetStatY(0.16);
    gStyle->SetStatX(0.85);
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.2);
    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.5, "y");

    c1->Divide(3, 2);
    for (int pad = 1; pad <= 6; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
    }

    TCut range;

    auto createAndDrawRank = [&](int pad, const RankingParams& param) {
        c1->cd(pad);
        range = Form("tan>=%.1f&&tan<%.1f", param.tan_low, param.tan_up);

        TString rank_title = Form("2D TRank (%.1f#leq tan#it{#theta} < %.1f);#sqrt{(dtan#it{#theta}_{x1})^{2}#plus (dtan#it{#theta}_{y1})^{2}#plus (dtan#it{#theta}_{x2})^{2}#plus (dtan#it{#theta}_{y2})^{2}};VPHsum", param.tan_low, param.tan_up);
        TH2D* rank = new TH2D("rank", rank_title, 50, 0.0, param.lin_max, param.range_max, 0, param.range_max);
        tree->Draw("(vph1+vph2):lin >> rank", range, "colz");

        c1->cd(pad + 3);
        TString rankl_title = Form("2D TRank Lateral (%.1f#leq tan#it{#theta} < %.1f);Angular difference (Lateral);VPHsum", param.tan_low, param.tan_up);
        TH2D* rankl = new TH2D("rankl", rankl_title, 50, 0.0, lat_range, param.range_max, 0, param.range_max);
        rankl->GetXaxis()->SetNdivisions(505);
        tree->Draw("(vph1+vph2):linl >> rankl", range, "colz");
    };

    for (int i = 0; i < 3; ++i) {
        createAndDrawRank(i + 1, params[i]);
    }
}

void distribution_map(TCanvas *c1, TTree *subtree, const double *AreaParam, const double bin_dmap, const float pitch, const double markersize) noexcept
{
    uint32_t bin = static_cast<uint32_t>(bin_dmap);
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    gStyle->SetOptStat("");
    gStyle->SetStatFormat(".4g");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    TString dz_title = Form("#Deltaz (1.0 < tan#it{#theta} < 1.1);x [mm];y [mm];Average of#Deltaz [#mum] at each %.1f#times %.1f mm^{2}", pitch, pitch);
    TString dz_1D_title = Form("#Deltaz at each %.1f#times %.1f mm^{2}     ;Average of#Deltaz [#mum] at each %.1f#times %.1f mm^{2};Frequency", pitch, pitch);
    TString dx_title = Form("#Deltax (1.0 < tan#it{#theta} < 1.1, interpolated);x [mm];y [mm];Average of#Deltax [#mum] in each %.1f#times %.1f mm^{2}", pitch, pitch);
    TString dx_1D_title = Form("#Deltax (1.0 < tan#it{#theta} < 1.1, interpolated);Average of#Deltax [#mum] in each %.1f#times %.1f mm^{2};Frequency", pitch, pitch);
    TString dy_title = Form("#Deltay (1.0 < tan#it{#theta} < 1.1, interpolated)     ;x [mm];y [mm];Average of#Deltay [#mum] in each %.1f#times %.1f mm^{2}", pitch, pitch);
    TString dy_1D_title = Form("#Deltay (1.0 < tan#it{#theta} < 1.1, interpolated)     ;Average of#Deltay [#mum] in each %.1f#times %.1f mm^{2};Frequency", pitch, pitch);

    TH2D* dz_2D = new TH2D("dz_2D", dz_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
    TH1D* dz_1D = new TH1D("dz_1D", dz_1D_title, 400, 150, 550);
    TH1D* dz_temp = new TH1D("dz_temp", "dz_temp", 400, 150, 550);
    TH2D* dx_2D = new TH2D("dx_2D", dx_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
    TH1D* dx_1D = new TH1D("dx_1D", dx_1D_title, 500, -100, 100);
    TH1D* dx_temp = new TH1D("dx_temp", "dx_temp", 100, -1000, 1000);
    TH2D* dy_2D = new TH2D("dy_2D", dy_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
    TH1D* dy_1D = new TH1D("dy_1D", dy_1D_title, 500, -100, 100);
    TH1D* dy_temp = new TH1D("dy_temp", "dy_temp", 100, -1000, 1000);

    TGraphErrors* white_frame = new TGraphErrors;
    white_frame->SetMarkerColor(0);
    white_frame->SetMarkerStyle(21);
    white_frame->SetMarkerSize(markersize);

    float pitch_half = pitch * 0.5;
    int nf = 0;

    for (int ix = 0; ix <= bin; ++ix) {
        for (int iy = 0; iy <= bin; ++iy) {
            double xcenter = dz_2D->GetXaxis()->GetBinCenter(ix);
            double ycenter = dz_2D->GetYaxis()->GetBinCenter(iy);
            TCut area = Form("(%f-x*0.001)^2<%f^2 && (%f-y*0.001)^2<%f^2", xcenter, pitch_half, ycenter, pitch_half);

            subtree->Draw("dz >> dz_temp", area, "goff");
            subtree->Draw("dx >> dx_temp", area, "goff");
            subtree->Draw("dy >> dy_temp", area, "goff");

            if (dz_temp->GetEntries() != 0) {
                double thickness = dz_temp->GetMean();
                dz_1D->Fill(thickness);
                dz_2D->SetBinContent(ix, iy, thickness);

                double dx_pitch = dx_temp->GetMean();
                dx_1D->Fill(dx_pitch);
                dx_2D->SetBinContent(ix, iy, dx_pitch);

                double dy_pitch = dy_temp->GetMean();
                dy_1D->Fill(dy_pitch);
                dy_2D->SetBinContent(ix, iy, dy_pitch);
            } else {
                white_frame->SetPoint(nf++, xcenter, ycenter);
            }
        }
    }

    auto setRange = [](TH1D* hist, TH2D* hist2D, double mean, double sigma) {
        double range = 5 * sigma;
        hist->GetXaxis()->SetRangeUser(mean - range, mean + range);
        hist2D->GetZaxis()->SetRangeUser(mean - range, mean + range);
    };

    setRange(dz_1D, dz_2D, dz_1D->GetMean(), dz_1D->GetRMS());
    setRange(dx_1D, dx_2D, 0.0, dx_1D->GetRMS());
    setRange(dy_1D, dy_2D, 0.0, dy_1D->GetRMS());

    c1->cd(1);
    dz_2D->Draw("colz");

    c1->cd(2);
    dz_1D->SetFillColorAlpha(92, 0.7);
    dz_1D->Draw();

    TLegend* dz_lg = new TLegend(0.68, 0.7, 0.9, 0.9);
    dz_lg->SetName("dz_lg");
    gDirectory->Add(dz_lg); // Add to gDirectory to retrieve later
    dz_lg->SetFillStyle(0);
    dz_lg->SetBorderSize(0);
    dz_lg->SetTextSize(0.04);
    dz_lg->AddEntry(dz_1D, Form("Areas      %d", static_cast<uint32_t>(dz_1D->GetEntries())), "");
    dz_lg->AddEntry(dz_1D, Form("Mean      %.2f [#mum]", dz_1D->GetMean()), "");
    dz_lg->AddEntry(dz_1D, Form("Std Dev   %.2f [#mum]", dz_1D->GetRMS()), "");
    dz_lg->Draw();

    c1->cd(3);
    dx_2D->Draw("colz");
    white_frame->Draw("same P");
    gPad->RedrawAxis();

    c1->cd(4);
    dy_2D->Draw("colz");
    white_frame->Draw("same P");
    gPad->RedrawAxis();
}

void distribution_map_re(TCanvas *c1)
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    // Retrieve objects from gDirectory
    TH2D* dz_2D = (TH2D*)gDirectory->Get("dz_2D");
    TH1D* dz_1D = (TH1D*)gDirectory->Get("dz_1D");
    TLegend* dz_lg = (TLegend*)gDirectory->Get("dz_lg");
    TH2D* dx_2D = (TH2D*)gDirectory->Get("dx_2D");
    TH2D* dy_2D = (TH2D*)gDirectory->Get("dy_2D");

    c1->cd(1);
    dz_2D->Draw("colz");

    c1->cd(2);
    dz_1D->SetFillColorAlpha(0, 0.7);
    dz_1D->Draw();
    dz_lg->Draw();

    c1->cd(3);
    dx_2D->Draw("colz");

    c1->cd(4);
    dy_2D->Draw("colz");
}

void distribution_xy(TCanvas *c1, TF1 *gaus)
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    // Retrieve objects from gDirectory
	TH2D* dx_2D = (TH2D*)gDirectory->Get("dx_2D");
	TH1D* dx_1D = (TH1D*)gDirectory->Get("dx_1D");
	TH2D* dy_2D = (TH2D*)gDirectory->Get("dy_2D");
	TH1D* dy_1D = (TH1D*)gDirectory->Get("dy_1D");

    c1->cd(1);
    dx_2D->Draw("colz");

    c1->cd(2);
    dy_2D->Draw("colz");

    c1->cd(3);
    dx_1D->SetFillColorAlpha(0, 0.7);
    dx_1D->Draw();
    double dxmean = dx_1D->GetMean();
    double dx_5sigma = 5 * dx_1D->GetRMS();
	double dx_range = TMath::Max(dx_5sigma - dxmean, dxmean + dx_5sigma);
    dx_1D->GetXaxis()->SetRangeUser(-dx_range, dx_range);
    dxmean = dx_1D->GetMean();
    dx_5sigma = 5 * dx_1D->GetRMS();
    dx_range = TMath::Max(dx_5sigma - dxmean, dxmean + dx_5sigma);
    dx_1D->GetXaxis()->SetRangeUser(-dx_range, dx_range);

    dx_1D->Fit("gaus", "q");
    double dx_fit_mean = gaus->GetParameter(1);
    double dx_fit_mean_err = gaus->GetParError(1);
    double dx_fit_sigma = gaus->GetParameter(2);
    double dx_fit_sigma_err = gaus->GetParError(2);

	int dx_entries = dx_1D->GetEntries();
	TLegend* dx_lg = new TLegend(0.74, 0.57, 0.96, 0.9);
	dx_lg->SetFillStyle(0);
	dx_lg->SetBorderSize(0);
	dx_lg->SetTextSize(0.04);
	dx_lg->AddEntry(dx_1D, Form("Areas   %d", dx_entries), "");
	dx_lg->AddEntry(dx_1D, Form("Mean   %.2g", dx_fit_mean), "");
	dx_lg->AddEntry(dx_1D, Form("    #pm%.2g", dx_fit_mean_err), "");
	dx_lg->AddEntry(dx_1D, Form("Sigma  %.2g", dx_fit_sigma), "");
	dx_lg->AddEntry(dx_1D, Form("    #pm%.2g", dx_fit_sigma_err), "");
	dx_lg->Draw();

	c1->cd(4);
    dy_1D->SetFillColorAlpha(0, 0.7);
    dy_1D->Draw();
    double dymean = dy_1D->GetMean();
    double dy_5sigma = 5 * dy_1D->GetRMS();
	double dy_range = TMath::Max(dy_5sigma - dymean, dymean + dy_5sigma);
	dy_1D->GetXaxis()->SetRangeUser(-dy_range, dy_range);
    dymean = dy_1D->GetMean();
    dy_5sigma = 5 * dy_1D->GetRMS();
    dy_range = TMath::Max(dy_5sigma - dymean, dymean + dy_5sigma);
	dy_1D->GetXaxis()->SetRangeUser(-dy_range, dy_range);

    dy_1D->Fit("gaus", "q");
    double dy_fit_mean = gaus->GetParameter(1);
    double dy_fit_mean_err = gaus->GetParError(1);
    double dy_fit_sigma = gaus->GetParameter(2);
    double dy_fit_sigma_err = gaus->GetParError(2);

	int dy_entries = dy_1D->GetEntries();
	TLegend* dy_lg = new TLegend(0.68, 0.57, 0.9, 0.9);
	dy_lg->SetFillStyle(0);
	dy_lg->SetBorderSize(0);
	dy_lg->SetTextSize(0.04);
	dy_lg->AddEntry(dy_1D, Form("Areas   %d", dy_entries), "");
	dy_lg->AddEntry(dy_1D, Form("Mean   %.2g", dy_fit_mean), "");
	dy_lg->AddEntry(dx_1D, Form("    #pm%.2g", dy_fit_mean_err), "");
	dy_lg->AddEntry(dy_1D, Form("Sigma  %.2g", dy_fit_sigma), "");
	dy_lg->AddEntry(dx_1D, Form("    #pm%.2g", dy_fit_sigma_err), "");
	dy_lg->Draw();
}
