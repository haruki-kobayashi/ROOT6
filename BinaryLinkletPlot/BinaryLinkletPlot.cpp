#define YAML_CPP_STATIC_DEFINE

#include <csignal>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TCut.h>
#include <TLegend.h>
#include <TColor.h>
#include <TError.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include <VxxReader/netscan_data_types_ui.h>
#include <argparse/argparse.hpp>
#include <yaml-cpp/yaml.h>

#include <MyHeader/Vec3.hpp>
#include <MyHeader/CalcTrackPair.hpp>
#include <MyHeader/NinjaFormat.hpp>
#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[]);
int ReadYaml(const std::string &param_file, std::unordered_map<std::string, bool> &on_plot,
    std::unordered_map<std::string, bool> &off_plot, argparse::ArgumentParser &parser,
    int &font_number, bool &hideGrid, std::string &palette_arg, int &NContours, bool &invertpalette,
    bool &negatepalette, std::vector<double> &cutX0, std::vector<double> &cutX1,
    std::vector<double> &cutY0, std::vector<double> &cutY1,std::vector<double> &cutAx0,
    std::vector<double> &cutAx1, std::vector<double> &cutAy0, std::vector<double> &cutAy1,
    std::vector<double> &cutTan0, std::vector<double> &cutTan1, std::vector<double> &cutPH0ang,
    std::vector<int> &cutPH0th, std::vector<double> &cutPH1ang, std::vector<int> &cutPH1th,
    std::vector<double> &cutVPH0ang, std::vector<int> &cutVPH0th, std::vector<double> &cutVPH1ang,
    std::vector<int> &cutVPH1th, double &angle_correction0, double &angle_correction1,
    std::vector<double> &TD_range, double &angle_max, double &angle_resolution,
    int &dxyz_cutPH, bool &failedFlag);
void ReadLinklet(std::string linkletfile, TTree* tree, TTree* subtree,
    uint32_t &pl0, uint32_t &pl1, size_t &entries, double &gap,
    const std::vector<double> &cutX0, const std::vector<double> &cutX1,
    const std::vector<double> &cutY0, const std::vector<double> &cutY1,
    const std::vector<double> &cutAx0, const std::vector<double> &cutAx1,
    const std::vector<double> &cutAy0, const std::vector<double> &cutAy1,
    const std::vector<double> &cutTan0, const std::vector<double> &cutTan1,
    const std::vector<double> &cutDar, const double cutDal,
    const std::vector<double> &cutDr, const double cutDl,
    const std::vector<double> &cutMd, const double cutOa,
    const std::vector<double> &cutPH0ang, const std::vector<int> &cutPH0th,
    const std::vector<double> &cutPH1ang, const std::vector<int> &cutPH1th,
    const std::vector<double> &cutVPH0ang, const std::vector<int> &cutVPH0th,
    const std::vector<double> &cutVPH1ang, const std::vector<int> &cutVPH1th,
    const double angle_correction0, const double angle_correction1, const int dxyz_cutPH);
void position(TCanvas *c1, TTree *tree, const size_t entries, uint32_t pl, uint8_t idx,
              const std::vector<double> &TD_range, const double *AreaParam) noexcept;
void angle(TCanvas *c1, TTree *tree, uint32_t pl, uint8_t idx,
    const double angle_max, const double angle_resolution) noexcept;
void difference_xy(TCanvas *c1, TTree *tree, double gap, uint32_t pl0, uint32_t pl1,
    const double angle_max, const double angle_resolution) noexcept;
void difference_rl(TCanvas *c1, TTree *tree, double gap, uint32_t pl0, uint32_t pl1,
    const double angle_max, const double angle_resolution) noexcept;
void difference_xyrl(TCanvas *c1, TTree *tree, double gap, uint32_t pl0, uint32_t pl1) noexcept;
void difference_1D(TCanvas *c1, TTree *tree, int type, const std::vector<double> &angle_list,
    int start_num, std::vector<double> &mean, std::vector<double> &sigma, std::vector<double> &sigma_err,
    bool &noEntries, double gap, uint32_t pl0, uint32_t pl1) noexcept;
void deviation(TCanvas *c1, TTree *tree, const double angle_max, const std::vector<double> &angle_list,
    const std::vector<std::vector<double>> &sigma, const std::vector<std::vector<double>> &sigma_err,
    double gap, uint8_t pl0, bool xy, bool dividedSqrt2) noexcept;
void md_oa(TCanvas *c1, TTree *tree, uint8_t pl0, const size_t entries,
    const double angle_max, const double angle_resolution) noexcept;
void dxdydz(TCanvas *c1, TTree *tree, const double *AreaParam) noexcept;
void dxdy(TCanvas *c1, TTree *tree, const double *AreaParam) noexcept;

namespace {
    // Ctrl+Cで終了したときの処理用にグローバル変数を定義
    TCanvas* global_c1 = nullptr;
    std::string global_output = "";

    // その他のグローバル変数
    bool global_darkmode = false;
}

// Ctrl+Cで終了したときの処理
void handleSIGINT(int) {
    std::cout << "\n*** Catched Ctrl+C: ";

    // PDFファイルの生成を完了させる
    if (global_c1 && !global_output.empty()) {
        global_c1->Print((global_output + "]").c_str());
        std::cout << "Output PDF file closed and saved. ***" << std::endl;
        std::cout << " Output     : " << global_output << std::endl;
    } else { // プロット開始前はそのまま終了
        std::cout << "Terminated. ***" << std::endl;
    }

    std::exit(0);
}

int main(int argc, char* argv[])
{
    // Ctrl+Cで終了したときの処理を設定
    std::signal(SIGINT, handleSIGINT);

    // 処理時間を計測
    TStopwatch sw;

    // argparseを使用して引数を解析
    std::cout << "\nInitializing..." << std::endl;
    argparse::ArgumentParser parser("BinaryLinkletPlot.exe", "1.0.0");
    parse_arguments(parser, argc, argv);

    // 引数を取得
    const auto linkletfile = parser.get<std::string>("input_binary_linklet");
    const auto output_arg = parser.get<std::string>("--output");
    const auto param_file = parser.get<std::string>("--param_file");
    std::unordered_map<std::string, bool> on_plot, off_plot;
    for (const auto& plot : parser.get<std::vector<std::string>>("--on_plot")) {
        on_plot[plot] = true;
    }
    for (const auto& plot : parser.get<std::vector<std::string>>("--off_plot")) {
        off_plot[plot] = true;
    }
    auto font_number = parser.get<int>("--font_number");
    auto hideGrid = parser.get<bool>("--hide_grid");
    global_darkmode = parser.get<bool>("--dark_mode");
    auto palette_arg = parser.get<std::string>("--palette");
    auto NContours = parser.get<int>("--contours");
    auto invertpalette = parser.get<bool>("--invert_palette");
    auto negatepalette = parser.get<bool>("--negate_palette");
    auto cutX0 = parser.get<std::vector<double>>("--cut_x0");
    auto cutX1 = parser.get<std::vector<double>>("--cut_x1");
    auto cutY0 = parser.get<std::vector<double>>("--cut_y0");
    auto cutY1 = parser.get<std::vector<double>>("--cut_y1");
    auto cutAx0 = parser.get<std::vector<double>>("--cut_ax0");
    auto cutAx1 = parser.get<std::vector<double>>("--cut_ax1");
    auto cutAy0 = parser.get<std::vector<double>>("--cut_ay0");
    auto cutAy1 = parser.get<std::vector<double>>("--cut_ay1");
    auto cutTan0 = parser.get<std::vector<double>>("--cut_tan0");
    auto cutTan1 = parser.get<std::vector<double>>("--cut_tan1");
    auto cutDar = parser.get<std::vector<double>>("--cut_dar");
    auto cutDal = parser.get<double>("--cut_dal");
    auto cutDr = parser.get<std::vector<double>>("--cut_dr");
    auto cutDl = parser.get<double>("--cut_dl");
    auto cutMd = parser.get<std::vector<double>>("--cut_md");
    auto cutOa = parser.get<double>("--cut_oa");
    auto cutPH0ang = parser.get<std::vector<double>>("--cut_ph0_angle");
    auto cutPH0th = parser.get<std::vector<int>>("--cut_ph0_threshold");
    auto cutPH1ang = parser.get<std::vector<double>>("--cut_ph1_angle");
    auto cutPH1th = parser.get<std::vector<int>>("--cut_ph1_threshold");
    auto cutVPH0ang = parser.get<std::vector<double>>("--cut_vph0_angle");
    auto cutVPH0th = parser.get<std::vector<int>>("--cut_vph0_threshold");
    auto cutVPH1ang = parser.get<std::vector<double>>("--cut_vph1_angle");
    auto cutVPH1th = parser.get<std::vector<int>>("--cut_vph1_threshold");
    auto angle_correction0 = parser.get<double>("--angle_correction0");
    auto angle_correction1 = parser.get<double>("--angle_correction1");
    auto TD_range = parser.get<std::vector<double>>("--track_density_range");
    auto angle_max = parser.get<double>("--angle_max");
    auto angle_resolution = parser.get<double>("--angle_resolution");
    auto dxyz_cutPH = parser.get<int>("--dxyz_cut_ph");

    // YAMLファイルから引数を取得
    if (parser.is_used("--param_file")) {
        bool failedFlag;
        int retVal = ReadYaml(param_file, on_plot, off_plot, parser, font_number,
            hideGrid, palette_arg, NContours, invertpalette, negatepalette,
            cutX0, cutX1, cutY0, cutY1, cutAx0, cutAx1, cutAy0, cutAy1,
            cutTan0, cutTan1, cutPH0ang, cutPH0th, cutPH1ang, cutPH1th,
            cutVPH0ang, cutVPH0th, cutVPH1ang, cutVPH1th, angle_correction0, angle_correction1,
            TD_range, angle_max, angle_resolution, dxyz_cutPH, failedFlag);
        if (failedFlag)
            return retVal;
    }

    // 出力ファイル名の設定
    const std::string output = (output_arg.empty())
        ? (linkletfile + ".pdf")
        : ((output_arg.size() > 4 && output_arg.substr(output_arg.size() - 4) == ".pdf")
            ? output_arg
            : output_arg + ".pdf");

    // 引数を利用する変数を設定
    const int font_code = 10 * font_number + 2;
    std::set<std::string> all_list = {
        "pos", "ang", "diff_xy", "diff_rl", "diff_xyrl", "md_oa",
        "diff_x", "diff_y", "diff_r", "diff_l", "diff_ax", "diff_ay", "diff_ar", "diff_al",
        "dev_xy", "dev_rl", "dev_xy2", "dev_rl2", "dxyz", "dxy"
    };
    if (on_plot.empty()) {
        for (const auto& plot : all_list) {
            on_plot[plot] = true;
        }
    }
    for (const auto& [plot, _] : on_plot) {
        if (off_plot.count(plot)) {
            on_plot[plot] = false; // off_plotに含まれるものはon_plotから除外
        }
    }
    std::set<std::string> plot_list;
    for (const auto& [plot, isEnabled] : on_plot) {
        if (isEnabled) {
            plot_list.insert(plot); // この時点でon_plotフラグがtrueになっているものをplot_listに追加
        }
    }
    if (plot_list.empty()) {
        std::cerr << "\nError: No plots are specified to be generated." << std::endl;
        return 1;
    }

    // PH, VPHカットが正しく設定されているか確認して設定する
    auto check_and_prepare_cut_args = [&](int idx, const char* ph_or_vph) -> bool {
        std::string ang_opt = (idx == 0)
            ? (std::string("-") + ph_or_vph + "0-ang")
            : (std::string("-") + ph_or_vph + "1-ang");
        std::string th_opt = (idx == 0)
            ? (std::string("-") + ph_or_vph + "0-th")
            : (std::string("-") + ph_or_vph + "1-th");
        std::vector<double>& ang = (std::string(ph_or_vph) == "ph")
            ? (idx == 0 ? cutPH0ang : cutPH1ang)
            : (idx == 0 ? cutVPH0ang : cutVPH1ang);
        std::vector<int>& th = (std::string(ph_or_vph) == "ph")
            ? (idx == 0 ? cutPH0th : cutPH1th)
            : (idx == 0 ? cutVPH0th : cutVPH1th);
        const char* label = (std::string(ph_or_vph) == "ph")
            ? (idx == 0 ? "PH0" : "PH1")
            : (idx == 0 ? "VPH0" : "VPH1");
        if (parser.is_used(ang_opt.c_str()) || parser.is_used(th_opt.c_str())) {
            if (!parser.is_used(ang_opt.c_str()) || !parser.is_used(th_opt.c_str())) {
                std::cerr << "\nError: " << label << " cut angles and thresholds must be used together." << std::endl;
                return false;
            }
            if (ang.empty() || th.empty()) {
                std::cerr << "\nError: " << label << " cut angles and thresholds must be specified." << std::endl;
                return false;
            }
            if (ang.size() != th.size()) {
                std::cerr << "\nError: The number of " << label
                    << " cut angles and thresholds must be the same." << std::endl;
                return false;
            }
        }
        ang.insert(ang.begin(), 0.0); // tanθ=0.0を先頭に追加
        th.push_back(0); // thresholdの末尾に0を追加
        return true;
    };
    if (!check_and_prepare_cut_args(0, "ph")) return 1;
    if (!check_and_prepare_cut_args(1, "ph")) return 1;
    if (!check_and_prepare_cut_args(0, "vph")) return 1;
    if (!check_and_prepare_cut_args(1, "vph")) return 1;

    // カラーパレットの設定
    std::variant<int, std::string> Palette;
    try {
        Palette = std::stoi(palette_arg); // 整数として解釈
    } catch (const std::invalid_argument&) {
        Palette = palette_arg; // 整数でなかったら文字列として扱う
    }
    std::visit([NContours](auto&& arg) { MyPalette::SetPalette(arg, NContours); }, Palette);
    if (invertpalette) MyPalette::InvertPalette();
    if (negatepalette) MyPalette::NegatePalette();
    // カラー定義
    float r1, g1, b1, r2, g2, b2, r3, g3, b3;
    int nColors = gStyle->GetNumberOfColors();
    gROOT->GetColor(gStyle->GetColorPalette(nColors * 0.15))->GetRGB(r1, g1, b1);
    gROOT->GetColor(90)->SetRGB(r1, g1, b1);
    gROOT->GetColor(gStyle->GetColorPalette(nColors * 0.85))->GetRGB(r3, g3, b3);
    gROOT->GetColor(91)->SetRGB(r3, g3, b3);
    gROOT->GetColor(gStyle->GetColorPalette(nColors * 0.5))->GetRGB(r2, g2, b2);
    gROOT->GetColor(92)->SetRGB(r2, g2, b2);
    gROOT->GetColor(93)->SetRGB(r2 * 0.6, g2 * 0.6, b2 * 0.6);
    gROOT->GetColor(94)->SetRGB( 200./255., 200./255., 200./255.);
    gROOT->GetColor(95)->SetRGB(  60./255.,  60./255.,  60./255.);
    // カラーの設定
    gStyle->SetHistFillColor(92); // ヒストグラム
    gStyle->SetHistLineColor(93); // ヒストグラムの枠
    gStyle->SetFuncColor(93);     // グラフ
    gStyle->SetGridColor(global_darkmode ? 95 : 94);       // グリッド
    gStyle->SetCanvasColor(global_darkmode ? 1 : 0);       // キャンバス(全体の背景)
    gStyle->SetPadColor(global_darkmode ? 1 : 0);          // Pad(グラフの背景はこれ)
    gStyle->SetStatColor(global_darkmode ? 1 : 0);         // 統計box
    gStyle->SetAxisColor(global_darkmode ? 0 : 1, "xyz");  // 軸
    gStyle->SetLabelColor(global_darkmode ? 0 : 1, "xyz"); // 軸ラベル(数値)
    gStyle->SetTitleColor(global_darkmode ? 0 : 1, "xyz"); // 軸タイトル
    gStyle->SetTitleTextColor(global_darkmode ? 0 : 1);    // メインのタイトル
    gStyle->SetFrameLineColor(global_darkmode ? 0 : 1);    // 描画エリアの枠
    gStyle->SetStatTextColor(global_darkmode ? 0 : 1);     // 統計box内のtext
    gStyle->SetLineColor(global_darkmode ? 0 : 1);         // 統計boxの枠など

    // スタイルの設定
    gStyle->SetOptStat("");               // 統計boxの表示項目
    gStyle->SetPadRightMargin(0.1);       // Pad右側のマージン
    gStyle->SetPadLeftMargin(0.1);        // Pad左側のマージン
    gStyle->SetPadTopMargin(0.1);         // Pad上側のマージン
    gStyle->SetPadBottomMargin(0.11);     // Pad下側のマージン
    gStyle->SetLabelOffset(0.008, "xyz"); // 軸ラベル(数値)と軸の距離
    gStyle->SetTitleOffset(1.1, "xyz");   // 軸titleと軸の距離
    gStyle->SetTitleY(0.985);             // タイトルのy方向の位置
    gStyle->SetHistFillStyle(1011); // ヒストグラム内部のパターン 1011=塗りつぶし
    gStyle->SetHistLineStyle(1);    // ヒストグラムの線種 1=直線
    gStyle->SetHistLineWidth(1);    // ヒストグラムの線幅 pixel
    gStyle->SetPadGridX(!hideGrid); // グリッドの表示
    gStyle->SetPadGridY(!hideGrid); // グリッドの表示
    gStyle->SetPadTickX(1);         // 上側x軸の目盛り表示
    gStyle->SetPadTickY(1);         // 右側y軸の目盛り表示
    gStyle->SetStatFont(font_code);         // 統計box内のフォント
    gStyle->SetLabelFont(font_code, "xyz"); // 軸ラベルのフォント
    gStyle->SetTitleFont(font_code, "xyz"); // 軸titleのフォント
    gStyle->SetTitleFont(font_code, "");    // titleのフォント
    gStyle->SetTextFont(font_code);         // textのフォント
    gStyle->SetLegendFont(font_code);       // 凡例のフォント
    gStyle->SetLabelSize(0.04, "xyz");      // 軸ラベルのサイズ
	gStyle->SetTitleSize(0.042, "xyz");     // 軸タイトルのサイズ

    // エラーメッセージ未満のROOTのメッセージを非表示に設定
    gErrorIgnoreLevel = kError;

    // TTreeの作成
    std::cout << "Starting to read binary linklet file..." << std::endl;
    TTree* tree = new TTree("tree", "");
    TTree* subtree = new TTree("subtree", "");

    uint32_t pl0 = 0, pl1 = 0;
    size_t entries = 0;
    double gap = 0.0;

    ReadLinklet(linkletfile, tree, subtree, pl0, pl1, entries, gap,
        cutX0, cutX1, cutY0, cutY1, cutAx0, cutAx1, cutAy0, cutAy1, cutTan0, cutTan1,
        cutDar, cutDal, cutDr, cutDl, cutMd, cutOa,
        cutPH0ang, cutPH0th, cutPH1ang, cutPH1th, cutVPH0ang, cutVPH0th, cutVPH1ang, cutVPH1th,
        angle_correction0, angle_correction1, dxyz_cutPH);

    double elapsed_time = sw.RealTime();
    double cpu_time = sw.CpuTime();
    std::cout << Form(
        "Linklet data successfully loaded. - Elapsed %.2f [s] (CPU: %.2f [s])",
        elapsed_time, cpu_time
    ) << std::endl;

    // 情報表示
    std::cout << "\n Input file   : " << linkletfile << std::endl;
    std::cout << Form(" Plate #      : %.3d & %.3d", pl0, pl1) << std::endl;
    std::cout << Form(" # of linklets: %d", entries) << std::endl;
    std::cout << Form(" Gap mean     : %.2f [um] (PL%.3d - PL%.3d)", gap, pl1, pl0) << std::endl;

    // プロット開始
    TDatime time_now;
    const char* Time_Now = Form(
        "%d-%02d-%02d %02d:%02d:%02d",
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(),
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    sw.Start();
	std::cout << " Plot start   : " << Time_Now << std::endl;

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());
    global_c1 = c1;
    global_output = output;

    // 角度範囲の設定 0.0-1.0: step 0.1, 1.0-: step 0.2
    std::vector<double> angle_list;
    angle_list.push_back(0.0);
    for (double a = 0.1; a <= std::min(1.0, angle_max); a = std::round((a + 0.1) * 1e6) / 1e6) {
        angle_list.push_back(a);
    }
    double a = (angle_list.back() < 1.0) ? 1.0 : angle_list.back();
    for (double b = a + 0.2; b <= angle_max + 1e-9; b = std::round((b + 0.2) * 1e6) / 1e6) {
        angle_list.push_back(b);
    }

    const size_t n_ranges = (angle_list.size() >= 2) ? (angle_list.size() - 1) : 0;
    const double diff_page = std::ceil(static_cast<double>(n_ranges) / 6.);

	// プログレスバーの初期化
    int page = 0;
    std::unordered_map<std::string, int> page_map = {
        {"pos", 2}, {"ang", 2},
        {"diff_xy", 1}, {"diff_rl", 1}, {"diff_xyrl", 1}, {"md_oa", 1},
        {"diff_x", diff_page}, {"diff_y", diff_page}, {"diff_r", diff_page}, {"diff_l", diff_page},
        {"diff_ax", diff_page}, {"diff_ay", diff_page}, {"diff_ar", diff_page}, {"diff_al", diff_page},
        {"dev_xy", 1}, {"dev_rl", 1}, {"dev_xy2", 1}, {"dev_rl2", 1}, {"dxyz", 1}, {"dxy", 1}
    };
    const int total = std::accumulate( // 合計ページ数
        plot_list.begin(),
        plot_list.end(),
        0, // 初期値
        [&page_map](int acc, const std::string& key) {
            auto it = page_map.find(key);
            return acc + (it != page_map.end() ? it->second : 0);
        }
    );
	MyUtil::ShowProgress(page, total);

    // データの座標の範囲を取得し、表示範囲とビンの数を決定する
    // フィルムの長辺の端から1cm外側までを最大表示範囲とし、縦横比を正しく保って表示する
    // 位置分布等のビン幅は1mm。dx, dy, dzのプロットのビン幅は計算時間短縮のためフィルムサイズによって変える
    const int MinX = std::min(tree->GetMinimum("x0"), tree->GetMinimum("x1"));
    const int MaxX = std::max(tree->GetMaximum("x0"), tree->GetMaximum("x1"));
    const int MinY = std::min(tree->GetMinimum("y0"), tree->GetMinimum("y1"));
    const int MaxY = std::max(tree->GetMaximum("y0"), tree->GetMaximum("y1"));
    const int RangeX = MaxX - MinX; // データ領域
    const int RangeY = MaxY - MinY; // データ領域
    double LowX, UpX, LowY, UpY, bin, bin_dxdydz, pitch; // 表示範囲とビンの数とdx, dy, dzのプロットのビン幅
    if (RangeX >= RangeY) {
        pitch = 5.0; // 5.0 mm pitch for dx, dy, dz plot
        LowX = MinX - 10000;
        UpX = MaxX + 10000;
        LowY = MinY - (RangeX - RangeY + 20000) * 0.5;
        UpY = MaxY + (RangeX - RangeY + 20000) * 0.5;
        bin = (RangeX + 20000) * 0.001;
        bin_dxdydz = (RangeX + 20000) * 0.0002;
        if (RangeX < 100000) {
            pitch = 1.0; // 1.0 mm pitch for dx, dy, dz plot
            bin_dxdydz *= 5;
        } else if (RangeX < 150000) {
            pitch = 2.5; // 2.5 mm pitch for dx, dy, dz plot
            bin_dxdydz *= 2;
        }
    } else {
        pitch = 5.0; // 5.0 mm pitch for dx, dy, dz plot
        LowX = MinX - (RangeY - RangeX + 20000) * 0.5;
        UpX = MaxX + (RangeY - RangeX + 20000) * 0.5;
        LowY = MinY - 10000;
        UpY = MaxY + 10000;
        bin = (RangeY + 20000) * 0.001;
        bin_dxdydz = (RangeY + 20000) * 0.0002; // 5.0mm pitch
        if (RangeY < 100000) {
            pitch = 1.0; // 1.0 mm pitch for dx, dy, dz plot
            bin_dxdydz *= 5;
        } else if (RangeY < 150000) {
            pitch = 2.5; // 2.5 mm pitch for dx, dy, dz plot
            bin_dxdydz *= 2;
        }
    }
    const double AreaParam[7] = {bin, LowX, UpX, LowY, UpY, bin_dxdydz, pitch};

    if (std::find(plot_list.begin(), plot_list.end(), "pos") != plot_list.end()) {
        position(c1, tree, entries, pl0, 0, TD_range, AreaParam);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("pos*");
        gDirectory->Delete("track_density");
        MyUtil::ShowProgress(page, total);

        position(c1, tree, entries, pl1, 1, TD_range, AreaParam);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("pos*");
        gDirectory->Delete("track_density");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "ang") != plot_list.end()) {
        angle(c1, tree, pl0, 0, angle_max, angle_resolution);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("ang*");
        MyUtil::ShowProgress(page, total);

        angle(c1, tree, pl1, 1, angle_max, angle_resolution);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("ang*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "diff_xy") != plot_list.end()) {
        difference_xy(c1, tree, gap, pl0, pl1, angle_max, angle_resolution);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*diff*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "diff_rl") != plot_list.end()) {
        difference_rl(c1, tree, gap, pl0, pl1, angle_max, angle_resolution);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*diff*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "diff_xyrl") != plot_list.end()) {
        difference_xyrl(c1, tree, gap, pl0, pl1);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*diff*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "md_oa") != plot_list.end()) {
        md_oa(c1, tree, pl0, entries, angle_max, angle_resolution);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("hist*");
        MyUtil::ShowProgress(page, total);
    }

    std::vector<std::vector<double>> mean(8, std::vector<double>(n_ranges, 0.0));
    std::vector<std::vector<double>> sigma(8, std::vector<double>(n_ranges, 0.0));
    std::vector<std::vector<double>> sigma_err(8, std::vector<double>(n_ranges, 0.0));

    bool noEntries;

	const std::vector<std::string> type_list = {
        "ax", "x", "ay", "y", "ar", "r", "al", "l"
    };
    for (int type = 0; type < 8; ++type)
	{
		if (std::find(plot_list.begin(), plot_list.end(), Form("diff_%s", type_list[type])) == plot_list.end())
            continue;
        for (int i = 0; i < static_cast<int>(diff_page); ++i)
        {
            std::vector<double> mean_vec = mean[type];
            std::vector<double> sigma_vec = sigma[type];
            std::vector<double> sigma_err_vec = sigma_err[type];
            difference_1D(c1, tree, type, angle_list, 6 * i,
                mean_vec, sigma_vec, sigma_err_vec, noEntries, gap, pl0, pl1);
            if (!mean_vec.empty()) mean[type] = mean_vec;
            if (!sigma_vec.empty()) sigma[type] = sigma_vec;
            if (!sigma_err_vec.empty()) sigma_err[type] = sigma_err_vec;
            c1->Print(output.c_str()); c1->Clear();
			gDirectory->Delete("hist*");
			gDirectory->Delete("lg*");
			MyUtil::ShowProgress(page, total);
			if (noEntries) {
                for (int j = i + 1; j < static_cast<int>(diff_page); ++j) {
                    MyUtil::ShowProgress(page, total);
                }
                break;
            }
		}
	}

    if (std::find(plot_list.begin(), plot_list.end(), "dev_xy") != plot_list.end()) {
        deviation(c1, tree, angle_max, angle_list, sigma, sigma_err, gap, pl0, true, false);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*dev*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "dev_rl") != plot_list.end()) {
        deviation(c1, tree, angle_max, angle_list, sigma, sigma_err, gap, pl0, false, false);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*dev*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "dev_xy2") != plot_list.end()) {
        deviation(c1, tree, angle_max, angle_list, sigma, sigma_err, gap, pl0, true, true);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*dev*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "dev_rl2") != plot_list.end()) {
        deviation(c1, tree, angle_max, angle_list, sigma, sigma_err, gap, pl0, false, true);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*dev*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "dxyz") != plot_list.end()) {
        dxdydz(c1, subtree, AreaParam);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*_temp");
        gDirectory->Delete("dz*");
        MyUtil::ShowProgress(page, total);
    }

    if (std::find(plot_list.begin(), plot_list.end(), "dxy") != plot_list.end()) {
        dxdy(c1, subtree, AreaParam);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("d*");
        MyUtil::ShowProgress(page, total);
    }

    gDirectory->Delete("subtree");

    // PDFファイルを閉じる
    c1->Print((output + "]").c_str());
	if (page < total) page = total;
    MyUtil::ShowProgress(page, total);

    // プロット終了
    time_now.Set();
    Time_Now = Form(
        "%d-%02d-%02d %02d:%02d:%02d",
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(),
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    elapsed_time = sw.RealTime();
    cpu_time = sw.CpuTime();
    std::cout << Form(
        "\n Plot end     : %s - Elapsed %.2f [s] (CPU: %.2f [s])",
        Time_Now, elapsed_time, cpu_time
    ) << std::endl;
    std::cout << " Output       : " << output << std::endl;

    delete c1;
    gDirectory->Delete("tree");

    return 0;
}

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[]) {
    parser.set_usage_max_line_width(80);
    parser.add_description("Tips: You can combine single-character arguments.\n"
        "      For example, \"-d -p kBird -i\" is equivalent to \"-dpi kBird\".");

    // 必須引数: binary_linkletファイル
    parser.add_argument("input_binary_linklet")
        .help("Path to the binary linklet file to be processed.")
        .required();
    parser.add_group("Optional arguments");
    // オプション引数: 出力ファイル名
    parser.add_argument("-o", "--output")
        .help("Output PDF file name. [default: input_binary_linklet.pdf]")
        .default_value(std::string());
    // オプション引数: パラメータファイル
    parser.add_argument("-param", "--param_file")
        .help("Path to the parameter YAML file.\n"
            "If not specified, the default parameters and command line will be used.\n"
            "The command line options will override the parameters in the YAML file.")
        .default_value(std::string());

    // オプション引数: その他の設定
    parser.add_argument("-on", "--on_plot")
        .help("Plot the selected plots. If not specified, all plots are selected.\n"
            "[pos, ang, diff_xy, diff_rl, diff_xyrl, md_oa,\n"
            " diff_x, diff_y, diff_r, diff_l,\n"
            " diff_ax, diff_ay, diff_ar, diff_al,\n"
            " dev_xy, dev_rl, dev_xy2, dev_rl2, dxyz, dxy]")
        .choices("pos", "ang", "diff_xy", "diff_rl", "diff_xyrl", "md_oa",
            "diff_x", "diff_y", "diff_r", "diff_l",
            "diff_ax", "diff_ay", "diff_ar", "diff_al",
            "dev_xy", "dev_rl", "dev_xy2", "dev_rl2", "dxyz", "dxy")
        .nargs(0, 20);
    parser.add_argument("-off", "--off_plot")
        .help("Do not plot the selected plots.\n"
            "[pos, ang, diff_xy, diff_rl, diff_xyrl, md_oa,\n"
            " diff_x, diff_y, diff_r, diff_l,\n"
            " diff_ax, diff_ay, diff_ar, diff_al,\n"
            " dev_xy, dev_rl, dev_xy2, dev_rl2, dxyz, dxy]")
        .choices("pos", "ang", "diff_xy", "diff_rl", "diff_xyrl", "md_oa",
            "diff_x", "diff_y", "diff_r", "diff_l",
            "diff_ax", "diff_ay", "diff_ar", "diff_al",
            "dev_xy", "dev_rl", "dev_xy2", "dev_rl2", "dxyz", "dxy")
        .nargs(0, 20);
    parser.add_argument("-f", "--font_number")
        .help("Font number (default: Helvetica).\n"
            "4, 13, 6, or 2 are recommended.\n"
            "Refer to https://root.cern.ch/doc/master/classTAttText.html\n"
            "for details.")
        .default_value(4)
        .scan<'i', int>();
    parser.add_argument("-g", "--hide_grid")
        .help("Hide grid. [default: false]")
        .flag();
    parser.add_argument("-d", "--dark_mode")
        .help("Dark mode. [default: false]")
        .flag();
    parser.add_argument("-p", "--palette")
        .help("Color palette to use.\n"
            "All ROOT palettes (https://root.cern.ch/doc/master/classTColor.html)\n"
            "and the following custom palettes are available:\n"
            " - kBirdDark, kBlueWhiteRed, kBlueBlackRed,\n"
            " - kGreenWhiteMagenta, kGreenBlackMagenta,\n"
            " - kLegacy (ROOT5's default).")
        .default_value(std::string("kBird"));
    parser.add_argument("-c", "--contours")
        .help("Number of contours in the color palette.")
        .default_value(256)
        .scan<'i', int>();
    parser.add_argument("-i", "--invert_palette")
        .help("Invert the color palette vertically. [default: false]")
        .flag();
    parser.add_argument("-n", "--negate_palette")
        .help("Negate the color palette. [default: false]")
        .flag();

    parser.add_group("Cut options");
    parser.add_argument("-x0", "--cut_x0")
        .help("Cut X0 range (μm).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-x1", "--cut_x1")
        .help("Cut X1 range (μm).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-y0", "--cut_y0")
        .help("Cut Y0 range (μm).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-y1", "--cut_y1")
        .help("Cut Y1 range (μm).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-ax0", "--cut_ax0")
        .help("Cut ax0 range (tanθ).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-ax1", "--cut_ax1")
        .help("Cut ax1 range (tanθ).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-ay0", "--cut_ay0")
        .help("Cut ay0 range (tanθ).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-ay1", "--cut_ay1")
        .help("Cut ay1 range (tanθ).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-tan0", "--cut_tan0")
        .help("Cut tan0 (tan0) range.\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-tan1", "--cut_tan1")
        .help("Cut tan1 (sqrt(ax1*ax1+ay1*ay1)) range.\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-dr", "--cut_dr")
        .help("Cut of radial distance (μm).\n"
            "First = slope, Second = intercept.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-dl", "--cut_dl")
        .help("Cut of lateral distance (μm).")
        .default_value(0.0)
        .scan<'g', double>();
    parser.add_argument("-dar", "--cut_dar")
        .help("Cut of radial angle difference (tanθ).\n"
            "First = slope, Second = intercept.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-dal", "--cut_dal")
        .help("Cut of lateral angle difference (tanθ).")
        .default_value(0.0)
        .scan<'g', double>();
    parser.add_argument("-md", "--cut_md")
        .help("Cut of minimum distance (μm).\n"
            "First = slope, Second = intercept.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-oa", "--cut_oa")
        .help("Cut of opening angle (rad).")
        .default_value(0.0)
        .scan<'g', double>();
    parser.add_argument("-ph0-ang", "--cut_ph0_angle")
        .help("Angle (tanθ) list of PH sum cut of PL0. Must be used with -ph0-th together.\n"
            "For example, \"-ph0-ang 0.1 0.2 0.5 -ph0-th 20 18 16\" defines cuts\n"
            "as [0.0,0.1): ≧20, [0.1,0.2): ≧18, [0.2,0.5): ≧16.")
        .default_value(std::vector<double>())
        .nargs(0, 100)
        .scan<'g', double>();
    parser.add_argument("-ph0-th", "--cut_ph0_threshold")
        .help("Threshold list of PH sum cut of PL0. Must be used with -ph0-ang together.\n"
            "For example, \"-ph0-ang 0.1 0.2 0.5 -ph0-th 20 18 16\" defines cuts\n"
            "as [0.0,0.1): ≧20, [0.1,0.2): ≧18, [0.2,0.5): ≧16.")
        .default_value(std::vector<int>())
        .nargs(0, 100)
        .scan<'d', int>();
    parser.add_argument("-ph1-ang", "--cut_ph1_angle")
        .help("Angle (tanθ) list of PH sum cut of PL1. Must be used with -ph1-th together.\n"
            "For example, \"-ph1-ang 0.1 0.2 0.5 -ph1-th 20 18 16\" defines cuts\n"
            "as [0.0,0.1): ≧20, [0.1,0.2): ≧18, [0.2,0.5): ≧16.")
        .default_value(std::vector<double>())
        .nargs(0, 100)
        .scan<'g', double>();
    parser.add_argument("-ph1-th", "--cut_ph1_threshold")
        .help("Threshold list of PH sum cut of PL1. Must be used with -ph1-ang together.\n"
            "For example, \"-ph1-ang 0.1 0.2 0.5 -ph1-th 20 18 16\" defines cuts\n"
            "as [0.0,0.1): ≧20, [0.1,0.2): ≧18, [0.2,0.5): ≧16.")
        .default_value(std::vector<int>())
        .nargs(0, 100)
        .scan<'d', int>();
    parser.add_argument("-vph0-ang", "--cut_vph0_angle")
        .help("Angle (tanθ) list of VPH sum cut of PL0. Must be used with -vph0-th together.\n"
            "For example, \"-vph0-ang 0.1 0.2 0.5 -vph0-th 60 40 20\" defines cuts\n"
            "as [0.0,0.1): ≧60, [0.1,0.2): ≧40, [0.2,0.5): ≧20.")
        .default_value(std::vector<double>())
        .nargs(0, 100)
        .scan<'g', double>();
    parser.add_argument("-vph0-th", "--cut_vph0_threshold")
        .help("Threshold list of VPH sum cut of PL0. Must be used with -vph0-ang together.\n"
            "For example, \"-vph0-ang 0.1 0.2 0.5 -vph0-th 60 40 20\" defines cuts\n"
            "as [0.0,0.1): ≧60, [0.1,0.2): ≧40, [0.2,0.5): ≧20.")
        .default_value(std::vector<int>())
        .nargs(0, 100)
        .scan<'d', int>();
    parser.add_argument("-vph1-ang", "--cut_vph1_angle")
        .help("Angle (tanθ) list of VPH sum cut of PL1. Must be used with -vph1-th together.\n"
            "For example, \"-vph1-ang 0.1 0.2 0.5 -vph1-th 60 40 20\" defines cuts\n"
            "as [0.0,0.1): ≧60, [0.1,0.2): ≧40, [0.2,0.5): ≧20.")
        .default_value(std::vector<double>())
        .nargs(0, 100)
        .scan<'g', double>();
    parser.add_argument("-vph1-th", "--cut_vph1_threshold")
        .help("Threshold list of VPH sum cut of PL1. Must be used with -vph1-ang together.\n"
            "For example, \"-vph1-ang 0.1 0.2 0.5 -vph1-th 60 40 20\" defines cuts\n"
            "as [0.0,0.1): ≧60, [0.1,0.2): ≧40, [0.2,0.5): ≧20.")
        .default_value(std::vector<int>())
        .nargs(0, 100)
        .scan<'d', int>();

    parser.add_group("More detailed options");
    parser.add_argument("-cor0", "--angle_correction0")
        .help("Base track angle correction factor of PL0.")
        .default_value(1.0)
        .scan<'g', double>();
    parser.add_argument("-cor1", "--angle_correction1")
        .help("Base track angle correction factor of PL1.")
        .default_value(1.0)
        .scan<'g', double>();
    parser.add_argument("--track_density_range")
        .help("Track density (/mm2) range for plotting. (default: auto)\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.0, 0.0}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("--angle_max")
        .help("Maximum angle range (tanθ) for plotting.")
        .default_value(5.0)
        .scan<'g', double>();
    parser.add_argument("--angle_resolution")
        .help("Resolution of angle hists (tanθ).")
        .default_value(0.05)
        .scan<'g', double>();
    parser.add_argument("--dxyz_cut_ph")
        .help("dx, dy, dz noise cut PH sum threshold (keep ≧ VAR)")
        .default_value(40)
        .scan<'i', int>();

    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << "\nError: " << err.what() << std::endl;
        std::cerr << parser;
        std::cerr << "\nError: " << err.what() << std::endl;
        std::exit(1);
    }
}

int ReadYaml(const std::string &param_file, std::unordered_map<std::string, bool> &on_plot,
    std::unordered_map<std::string, bool> &off_plot, argparse::ArgumentParser &parser, int &font_number,
    bool &hideGrid, std::string &palette_arg, int &NContours, bool &invertpalette, bool &negatepalette,
    std::vector<double> &cutX0, std::vector<double> &cutX1,std::vector<double> &cutY0,
    std::vector<double> &cutY1, std::vector<double> &cutAx0, std::vector<double> &cutAx1,
    std::vector<double> &cutAy0, std::vector<double> &cutAy1, std::vector<double> &cutTan0,
    std::vector<double> &cutTan1, std::vector<double> &cutPH0ang, std::vector<int> &cutPH0th,
    std::vector<double> &cutPH1ang, std::vector<int> &cutPH1th,
    std::vector<double> &cutVPH0ang, std::vector<int> &cutVPH0th,
    std::vector<double> &cutVPH1ang, std::vector<int> &cutVPH1th,
    double &angle_correction0, double &angle_correction1, std::vector<double> &TD_range,
    double &angle_max, double &angle_resolution, int &dxyz_cutPH, bool &failedFlag)
{
    failedFlag = true;
    {
        try {
            YAML::Node params = YAML::LoadFile(param_file);
            for (const auto &plot : params["on_plot"].as<std::vector<std::string>>(std::vector<std::string>())) {
                on_plot[plot] = true;
            }
            for (const auto &plot : params["off_plot"].as<std::vector<std::string>>(std::vector<std::string>())) {
                if (!on_plot.count(plot)) { // YAMLファイルのoff_plotからon_plotと被るものを除外(コマンドライン引数の優先)
                    off_plot[plot] = true;
                } else {
                    on_plot.erase(plot); // 被ったものはon_plotフラグをfalseにする
                }
            }
            if (!parser.is_used("-f")) {
                font_number = params["font_number"].as<int>(4);
            }
            hideGrid = params["hide_grid"].as<bool>(false);
            if (parser.is_used("-g")) {
                hideGrid = !hideGrid;
            }
            global_darkmode = params["dark_mode"].as<bool>(false);
            if (parser.is_used("-d")) {
                global_darkmode = !global_darkmode;
            }
            if (!parser.is_used("-p")) {
                palette_arg = params["palette"].as<std::string>("kBird");
            }
            if (!parser.is_used("-c")) {
                NContours = params["contours"].as<int>(256);
            }
            invertpalette = params["invert_palette"].as<bool>(false);
            if (parser.is_used("-i")) {
                invertpalette = !invertpalette;
            }
            negatepalette = params["negate_palette"].as<bool>(false);
            if (parser.is_used("-n")) {
                negatepalette = !negatepalette;
            }
            if (!parser.is_used("-x0")) {
                cutX0 = params["cut_x0"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-x1")) {
                cutX1 = params["cut_x1"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-y0")) {
                cutY0 = params["cut_y0"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-y1")) {
                cutY1 = params["cut_y1"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-ax0")) {
                cutAx0 = params["cut_ax0"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-ax1")) {
                cutAx1 = params["cut_ax1"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-ay0")) {
                cutAy0 = params["cut_ay0"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-ay1")) {
                cutAy1 = params["cut_ay1"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-tan0")) {
                cutTan0 = params["cut_tan0"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-tan1")) {
                cutTan1 = params["cut_tan1"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("-ph0-ang")) {
                cutPH0ang = params["cut_ph0_angle"].as<std::vector<double>>(std::vector<double>());
            }
            if (!parser.is_used("-ph0-th")) {
                cutPH0th = params["cut_ph0_threshold"].as<std::vector<int>>(std::vector<int>({0, 0}));
            }
            if (!parser.is_used("-ph1-ang")) {
                cutPH1ang = params["cut_ph1_angle"].as<std::vector<double>>(std::vector<double>());
            }
            if (!parser.is_used("-ph1-th")) {
                cutPH1th = params["cut_ph1_threshold"].as<std::vector<int>>(std::vector<int>({0, 0}));
            }
            if (!parser.is_used("-vph0-ang")) {
                cutVPH0ang = params["cut_vph0_angle"].as<std::vector<double>>(std::vector<double>());
            }
            if (!parser.is_used("-vph0-th")) {
                cutVPH0th = params["cut_vph0_threshold"].as<std::vector<int>>(std::vector<int>());
            }
            if (!parser.is_used("-vph1-ang")) {
                cutVPH1ang = params["cut_vph1_angle"].as<std::vector<double>>(std::vector<double>());
            }
            if (!parser.is_used("-vph1-th")) {
                cutVPH1th = params["cut_vph1_threshold"].as<std::vector<int>>(std::vector<int>());
            }
            if (!parser.is_used("--angle_correction0")) {
                angle_correction0 = params["angle_correction0"].as<double>(1.0);
            }
            if (!parser.is_used("--angle_correction1")) {
                angle_correction1 = params["angle_correction1"].as<double>(1.0);
            }
            if (!parser.is_used("--track_density_range")) {
                TD_range = params["track_density_range"].as<std::vector<double>>(std::vector<double>({0.0, 0.0}));
            }
            if (!parser.is_used("--angle_max")) {
                angle_max = params["angle_max"].as<double>(6.0);
            }
            if (!parser.is_used("--angle_resolution")) {
                angle_resolution = params["angle_resolution"].as<double>(0.1);
            }
            if (!parser.is_used("--dxyz_cut_ph")) {
                dxyz_cutPH = params["dxyz_cut_ph"].as<int>(40);
            }
        } catch (const YAML::Exception &e) {
            std::cerr << "\nError: Failed to load parameter file: " << param_file << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }
    failedFlag = false;
    return 0;
}

void ReadLinklet(std::string linkletfile, TTree* tree, TTree* subtree,
    uint32_t &pl0, uint32_t &pl1, size_t &entries, double &gap,
    const std::vector<double> &cutX0, const std::vector<double> &cutX1,
    const std::vector<double> &cutY0, const std::vector<double> &cutY1,
    const std::vector<double> &cutAx0, const std::vector<double> &cutAx1,
    const std::vector<double> &cutAy0, const std::vector<double> &cutAy1,
    const std::vector<double> &cutTan0, const std::vector<double> &cutTan1,
    const std::vector<double> &cutDar, const double cutDal,
    const std::vector<double> &cutDr, const double cutDl,
    const std::vector<double> &cutMd, const double cutOa,
    const std::vector<double> &cutPH0ang, const std::vector<int> &cutPH0th,
    const std::vector<double> &cutPH1ang, const std::vector<int> &cutPH1th,
    const std::vector<double> &cutVPH0ang, const std::vector<int> &cutVPH0th,
    const std::vector<double> &cutVPH1ang, const std::vector<int> &cutVPH1th,
    const double angle_correction0, const double angle_correction1, const int dxyz_cutPH)
{
    // TTreeのブランチを作成
    uint8_t ph0, ph1;
    uint32_t vph0, vph1;
    double x0, x1, y0, y1, ax0, ax1, ay0, ay1, tan0, tan1, ar, al, dx, dy, dz, dax, day, dr, dl, md, oa;
    // tree
    const std::array<std::pair<const char*, void*>, 2> uint8Branches = {{
        {"ph0", &ph0}, {"ph1", &ph1}
    }};
    for (const auto& branch : uint8Branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/B").c_str());
    }
    const std::array<std::pair<const char*, void*>, 4> uint32Branches = {{
        {"pl0", &pl0}, {"pl1", &pl1}, {"vph0", &vph0}, {"vph1", &vph1}
    }};
    for (const auto& branch : uint32Branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/I").c_str());
    }
    const std::array<std::pair<const char*, void*>, 21> doubleBranches = {{
        {"x0", &x0}, {"x1", &x1}, {"y0", &y0}, {"y1", &y1},
        {"ax0", &ax0}, {"ax1", &ax1}, {"ay0", &ay0}, {"ay1", &ay1},
        {"tan0", &tan0}, {"tan1", &tan1}, {"ar", &ar}, {"al", &al},
        {"dx", &dx}, {"dy", &dy}, {"dz", &dz}, {"dax", &dax}, {"day", &day},
        {"dr", &dr}, {"dl", &dl}, {"md", &md}, {"oa", &oa}
    }};
    for (const auto& branch : doubleBranches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }
    // subtree
    const std::array<std::pair<const char*, void*>, 5> doubleBranchesSub = {{
        {"x", &x0}, {"y", &y0}, {"dx", &dx}, {"dy", &dy}, {"dz", &dz}
    }};
    for (const auto& branch : doubleBranchesSub) {
        subtree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    std::ifstream ifs(linkletfile.c_str(), std::ios::in | std::ios::binary);
    if (!ifs.is_open()) {
        std::cerr << "\nError: Cannot open file '" << linkletfile << "'" << std::endl;
        std::exit(1);
    }
    entries = 0;
    double gap_sum = 0.0;

    const size_t sz_netscan = sizeof(vxx::linklet_t);
    const size_t sz_ninja = sizeof(NinjaFormat::output_format_link);
    const size_t buf_size = std::max(sz_netscan, sz_ninja);
    std::vector<char> buf(buf_size);

    // 最初の1レコードを読み、形式を自動判別してからループを回す
    if (!ifs.read(buf.data(), buf_size)) {
        std::cerr << "\nError: Cannot read file or file empty\n";
        std::exit(1);
    }
    size_t first_bytes = static_cast<size_t>(ifs.gcount());

    // NETSCAN形式らしさをチェック
    bool file_is_netscan = false;
    if (first_bytes >= sz_netscan) {
        vxx::linklet_t tmp;
        std::memcpy(&tmp, buf.data(), sz_netscan);
        if (NinjaFormat::PlausibleNetscan(tmp)) file_is_netscan = true;
    }

    // NINJA独自形式 (3D linlet) らしさをチェック
    bool file_is_ninja = false;
    if (!file_is_netscan && first_bytes >= sz_ninja) {
        NinjaFormat::output_format_link tmp;
        std::memcpy(&tmp, buf.data(), sz_ninja);
        if (NinjaFormat::PlausibleNinja(tmp)) file_is_ninja = true;
    }

    // 追加チェック: ファイルサイズが形式に合っているか
    if (!file_is_netscan && !file_is_ninja) {
        // ファイルサイズを取得して形式を判定
        try {
            auto total = static_cast<std::uintmax_t>(std::filesystem::file_size(linkletfile));
            if (total > 0) {
                if (total % sz_netscan == 0) file_is_netscan = true;
                else if (total % sz_ninja == 0) file_is_ninja = true;
                else file_is_netscan = true; // デフォルトでNETSCAN形式とみなす
            } else {
                std::cerr << "\nError: Cannot determine file format";
                std::exit(1);
            }
        } catch (const std::filesystem::filesystem_error &e) {
            std::cerr << "\nError: Cannot determine file format (filesystem error): " << e.what() << std::endl;
            std::exit(1);
        }
    }

    std::cout << "Detected file format: " << (
        file_is_netscan ? "NETSCAN (class linklet_t)" : "NINJA3D (class output_format_link)"
    ) << std::endl;

    auto process_linklet = [&](const vxx::linklet_t& l) {
        pl0 = l.pos[0] * 0.1;
        pl1 = l.pos[1] * 0.1;

        x0 = l.b[0].x;
        if (cutX0[0] < cutX0[1] && (x0 < cutX0[0] || x0 > cutX0[1])) return;
        x1 = l.b[1].x;
        if (cutX1[0] < cutX1[1] && (x1 < cutX1[0] || x1 > cutX1[1])) return;
        y0 = l.b[0].y;
        if (cutY0[0] < cutY0[1] && (y0 < cutY0[0] || y0 > cutY0[1])) return;
        y1 = l.b[1].y;
        if (cutY1[0] < cutY1[1] && (y1 < cutY1[0] || y1 > cutY1[1])) return;

        ax0 = l.b[0].ax * angle_correction0;
        if (cutAx0[0] < cutAx0[1] && (ax0 < cutAx0[0] || ax0 > cutAx0[1])) return;
        ax1 = l.b[1].ax * angle_correction1;
        if (cutAx1[0] < cutAx1[1] && (ax1 < cutAx1[0] || ax1 > cutAx1[1])) return;
        ay0 = l.b[0].ay * angle_correction0;
        if (cutAy0[0] < cutAy0[1] && (ay0 < cutAy0[0] || ay0 > cutAy0[1])) return;
        ay1 = l.b[1].ay * angle_correction1;
        if (cutAy1[0] < cutAy1[1] && (ay1 < cutAy1[0] || ay1 > cutAy1[1])) return;

        tan0 = std::sqrt(ax0 * ax0 + ay0 * ay0);
        if (cutTan0[0] < cutTan0[1] && (tan0 < cutTan0[0] || tan0 > cutTan0[1])) return;
        tan1 = std::sqrt(ax1 * ax1 + ay1 * ay1);
        if (cutTan1[0] < cutTan1[1] && (tan1 < cutTan1[0] || tan1 > cutTan1[1])) return;

        ph0 = static_cast<uint8_t>(l.b[0].m[0].ph * 0.0001 + l.b[0].m[1].ph * 0.0001);
        ph1 = static_cast<uint8_t>(l.b[1].m[0].ph * 0.0001 + l.b[1].m[1].ph * 0.0001);
        vph0 = static_cast<uint32_t>(l.b[0].m[0].ph % 10000 + l.b[0].m[1].ph % 10000);
        vph1 = static_cast<uint32_t>(l.b[1].m[0].ph % 10000 + l.b[1].m[1].ph % 10000);

        for (size_t i = 0; i < cutPH0ang.size() - 1; ++i) {
        if (tan0 >= cutPH0ang[i] && tan0 < cutPH0ang[i + 1] &&
            ph0 < cutPH0th[i]) return;
        }
        for (size_t i = 0; i < cutPH1ang.size() - 1; ++i) {
        if (tan1 >= cutPH1ang[i] && tan1 < cutPH1ang[i + 1] &&
            ph1 < cutPH1th[i]) return;
        }
        for (size_t i = 0; i < cutVPH0ang.size() - 1; ++i) {
        if (tan0 >= cutVPH0ang[i] && tan0 < cutVPH0ang[i + 1] &&
            vph0 < cutVPH0th[i]) return;
        }
        for (size_t i = 0; i < cutVPH1ang.size() - 1; ++i) {
        if (tan1 >= cutVPH1ang[i] && tan1 < cutVPH1ang[i + 1] &&
            vph1 < cutVPH1th[i]) return;
        }

        dax = ax1 - ax0;
        day = ay1 - ay0;

        ar = (dax * ax0 + day * ay0) / tan0;
        if (cutDar.size() >= 2 && (cutDar[0] != 0.0 || cutDar[1] != 0.0)) {
            double allowed_dar = cutDar[0] * tan0 + cutDar[1];
            if (std::abs(ar) > allowed_dar) return;
        }
        al = (dax * ay0 - day * ax0) / tan0;
        if (cutDal > 0.0) {
            if (std::abs(al) > cutDal) return;
        }

        dx = l.dx;
        dy = l.dy;
        dz = l.b[1].z - l.b[0].z;

        dr = CalcTrackPair::RadialDistance(l.b[0], l.b[1]);
        if (cutDr.size() >= 2 && (cutDr[0] != 0.0 || cutDr[1] != 0.0)) {
            double allowed_dr = cutDr[0] * tan0 + cutDr[1];
            if (std::abs(dr) > allowed_dr) return;
        }
        dl = CalcTrackPair::LateralDistance(l.b[0], l.b[1]);
        if (cutDl > 0.0) {
            if (std::abs(dl) > cutDl) return;
        }
        md = CalcTrackPair::MinimumDistance(l.b[0], l.b[1]);
        if (cutMd.size() >= 2 && (cutMd[0] != 0.0 || cutMd[1] != 0.0)) {
            double allowed_md = cutMd[0] * tan0 + cutMd[1];
            if (md > allowed_md) return;
        }
        oa = CalcTrackPair::OpeningAngle(l.b[0], l.b[1]);
        if (cutOa > 0.0) {
            if (oa > cutOa) return;
        }

        tree->Fill();

        if ((ph0 + ph1) >= dxyz_cutPH && tan0 > 1.0 && tan0 < 1.1) subtree->Fill();

        entries++;
        gap_sum += dz;
    };

    ifs.clear();
    ifs.seekg(0, std::ios::beg);

    vxx::linklet_t l;
    if (file_is_netscan) { // NETSCAN formatとして読み込み
        const size_t rec = sz_netscan;
        std::vector<char> recbuf(rec);
        while (ifs.read(recbuf.data(), rec)) {
            std::memcpy(&l, recbuf.data(), rec);
            process_linklet(l);
        }
    } else { // NINJA format (3D linklet) として読み込み
        const size_t rec = sz_ninja;
        std::vector<char> recbuf(rec);
        while (ifs.read(recbuf.data(), rec)) {
            NinjaFormat::output_format_link nl;
            std::memcpy(&nl, recbuf.data(), rec);
            l = NinjaFormat::ConvertNinjaToNetscan(nl); // linklet_tに変換
            process_linklet(l);
        }
    }

    if (!ifs.eof()) { // ファイル終端以外の理由でループが終了した場合
        std::cerr << "\nError: Failed to read file '" << linkletfile << "'" << std::endl;
        std::exit(1);
    }
    if (entries == 0) {
        std::cerr << "\nError: No valid entries found in file '" << linkletfile << "'" << std::endl;
        std::exit(1);
    }
    ifs.close();

    gap = gap_sum / static_cast<double>(entries);
}

void position(TCanvas *c1, TTree *tree, const size_t entries, uint32_t pl, uint8_t idx,
    const std::vector<double> &TD_range, const double *AreaParam) noexcept
{
	gStyle->SetTitleOffset(1.1, "xyz");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    uint32_t bin = static_cast<uint32_t>(AreaParam[0]);
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    TH2D* position_2D = new TH2D(
        "position_2D", Form("Position PL%03d;x [mm];y [mm];/mm^{2}", pl),
        bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001
    );
    if (TD_range[1] > 0.0) position_2D->GetZaxis()->SetRangeUser(TD_range[0], TD_range[1]);
    tree->Draw(Form("y%d*0.001:x%d*0.001 >> position_2D", idx, idx), "", "goff");

    TLegend* pos_lg = new TLegend(0.6, 0.9, 0.75, 1.0);
    pos_lg->SetFillStyle(0);
    pos_lg->SetBorderSize(0);
    pos_lg->SetTextSize(0.04);
    pos_lg->SetTextColor(global_darkmode ? 0 : 1);
    pos_lg->AddEntry(position_2D, Form("Entries %.0f", position_2D->GetEntries()), "");

    c1->cd(1);
    position_2D->Draw("colz");
    pos_lg->Draw();

    for (int pad = 2; pad <= 3; ++pad) {
        c1->cd(pad);
        position_2D->SetFillColor((pad == 2) ? 91 : 90);
        TH1D* proj = (pad == 2) ? position_2D->ProjectionY() : position_2D->ProjectionX();
        proj->Draw((pad == 2) ? "hbar" : "bar");
        proj->SetTitle("");
    }

    c1->cd(4);
    gStyle->SetTitleOffset(1.5, "y");
    TH1D* track_density = new TH1D(
        "track_density", ";Track Density [/mm^{2}];Frequency", 100000, 0, 100000
    );
    int Xbins = position_2D->GetNbinsX();
    int Ybins = position_2D->GetNbinsY();
    double min_density = 0.0;
    double max_density = 0.0;
    for (int xBin = 0; xBin < Xbins; ++xBin) {
        for (int yBin = 0; yBin < Ybins; ++yBin) {
            double density = position_2D->GetBinContent(xBin + 1, yBin + 1);
            if (density > 0.0) track_density->Fill(density);
            if (max_density < density) max_density = density;
        }
    }
    if (TD_range[1] != 0.0) {
        min_density = TD_range[0];
        max_density = TD_range[1];
    }
    // ビン幅と合っていないと範囲設定が正しくできないため調整する
    min_density = std::floor(min_density / 1.0) * 1.0; // 整数に切り下げ(1はビン幅)
    max_density = std::ceil(max_density / 1.0) * 1.0;  // 整数に切り上げ(1はビン幅)
    track_density->GetXaxis()->SetRangeUser(min_density, max_density);
    track_density->SetFillStyle(0);
    track_density->SetLineWidth(2);
    track_density->Draw();
    MyUtil::PaintBins(track_density, min_density, max_density); // 各ビンをカラーパレットの色で塗る

	int density_entries = track_density->GetEntries();
    double density_mean = track_density->GetMean();
    double density_stddev = track_density->GetStdDev();
    TLegend* density_lg = new TLegend(0.67, 0.7, 0.9, 0.9);
    density_lg->SetFillStyle(0);
    density_lg->SetBorderSize(0);
    density_lg->SetTextSize(0.04);
    density_lg->SetTextColor(global_darkmode ? 0 : 1);
    density_lg->AddEntry(track_density, "Track Density [/mm^{2}]", "");
    density_lg->AddEntry(track_density, Form("%d areas", density_entries), "");
    density_lg->AddEntry(track_density, Form("Mean   %.2f", density_mean), "");
    density_lg->AddEntry(track_density, Form("Std Dev   %.2f", density_stddev), "");
    density_lg->Draw();
}

void angle(TCanvas *c1, TTree *tree, uint32_t pl, uint8_t idx,
    const double angle_max, const double angle_resolution) noexcept
{
	gStyle->SetTitleOffset(1.1, "xyz");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }
    c1->GetPad(4)->SetBottomMargin(0.15);

    uint32_t angle_bin = 2 / angle_resolution * angle_max;

    TString angtitle = Form(
        "Angle PL%03d;tan#it{#theta}_{x};tan#it{#theta}_{y};/(%g rad)^{2}", pl, angle_resolution
    );
    TH2D* angle_2D = new TH2D(
        "angle_2D", angtitle, angle_bin, -angle_max, angle_max, angle_bin, -angle_max, angle_max
    );
    tree->Draw(Form("ay%d:ax%d >> angle_2D", idx, idx), "", "goff");

    TLegend* ang_lg = new TLegend(0.6, 0.9, 0.75, 1.0);
    ang_lg->SetFillStyle(0);
    ang_lg->SetBorderSize(0);
    ang_lg->SetTextSize(0.04);
    ang_lg->SetTextColor(global_darkmode ? 0 : 1);
    ang_lg->AddEntry(angle_2D, Form("Entries %.0f", angle_2D->GetEntries()), "");

    c1->cd(1);
    angle_2D->Draw("colz");
	ang_lg->Draw();

    for (int pad = 2; pad <= 3; ++pad) {
        c1->cd(pad);
        angle_2D->SetFillColor((pad == 2) ? 91 : 90);
        TH1D* proj = (pad == 2) ? angle_2D->ProjectionY() : angle_2D->ProjectionX();
        proj->Draw((pad == 2) ? "hbar" : "bar");
        proj->SetTitle("");
    }

    c1->cd(4);
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.5, "y");

    TString ang1Dtitle = Form(
        ";#sqrt{tan^{2}#it{#theta}_{x}^{ }#plus tan^{2}#it{#theta}_{y}} (PL%.3d);Frequency",
        pl
    );
	TH1D* angle_1D = new TH1D("angle_1D", ang1Dtitle, angle_bin, 0.0, angle_max);
    angle_1D->SetFillColorAlpha(92, 0.7);
	tree->Draw(Form("tan%d>>angle_1D", idx));
	TLegend* ang_lg2 = new TLegend(0.45, 0.9, 0.6, 1.0);
	ang_lg2->SetFillStyle(0);
	ang_lg2->SetBorderSize(0);
	ang_lg2->SetTextSize(0.04);
    ang_lg2->SetTextColor(global_darkmode ? 0 : 1);
	ang_lg2->AddEntry(angle_1D, Form("Entries %.0f", angle_1D->GetEntries()), "");
	ang_lg2->Draw();
}

void difference_xy(TCanvas *c1, TTree *tree, double gap, uint32_t pl0, uint32_t pl1,
    const double angle_max, const double angle_resolution) noexcept
{
    gStyle->SetTitleOffset(1.1,"x");
    gStyle->SetTitleOffset(1.3,"y");

	c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
        c1->GetPad(pad)->SetLeftMargin(0.15);
        c1->GetPad(pad)->SetBottomMargin(0.13);
    }

    gap = std::abs(gap);
    double ang_range = 0.8;
	double pos_range = gap * 0.375;
    uint32_t angle_bin = 2 / angle_resolution * angle_max;

    auto drawHistogram = [&](const char* draw_expr, const char* hist_name,
                             const TString& title, int pad) {
        double range = (pad <= 2) ? ang_range : pos_range;
        TH2D* hist = new TH2D(
            hist_name, title, angle_bin, -1 * angle_max, angle_max,
            200, -1 * range, range
        );
        c1->cd(pad);
        draw_expr = Form("%s >> %s", draw_expr, hist_name);
        tree->Draw(draw_expr, "", "colz");
    };

    drawHistogram("dax:ax0", "angdiff_x",
                  Form("Angle difference x;tan#it{#theta}_{x} (PL%.3d);\
                   #Deltatan#it{#theta}_{x} (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
                  1);

    drawHistogram("day:ay0", "angdiff_y",
                  Form("Angle difference y;tan#it{#theta}_{y} (PL%.3d);\
                   #Deltatan#it{#theta}_{y} (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
                  2);

    drawHistogram("dx:ax0", "posdiff_x",
                  Form("Position difference x;tan#it{#theta}_{x} (PL%.3d);\
                   #Deltax [#mum] (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
                  3);

    drawHistogram("dy:ay0", "posdiff_y",
                  Form("Position difference y;tan#it{#theta}_{y} (PL%.3d);\
                   #Deltay [#mum] (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
                  4);
}

void difference_rl(TCanvas *c1, TTree *tree, double gap, uint32_t pl0, uint32_t pl1,
    const double angle_max, const double angle_resolution) noexcept
{
    gStyle->SetTitleOffset(1.3,"x");
    gStyle->SetTitleOffset(1.3,"y");

	c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
        c1->GetPad(pad)->SetLeftMargin(0.15);
        c1->GetPad(pad)->SetBottomMargin(0.14);
    }

    gap = std::abs(gap);
    double ang_range_r = 0.8;
	double ang_range_l = 0.06;
	double pos_range_r = gap * 0.1;
	double pos_range_l = gap * 0.07;
    uint32_t angle_bin = 1 / angle_resolution * angle_max;

    auto drawHistogram = [&](const char* draw_expr, const char* hist_name, const TString& title,
                           int pad, int xbins, double xmin, double xmax,
                           int ybins, double ymin, double ymax) {
        TH2D* hist = new TH2D(hist_name, title, xbins, xmin, xmax, ybins, ymin, ymax);
        c1->cd(pad);
        draw_expr = Form("%s >> %s", draw_expr, hist_name);
        tree->Draw(draw_expr, "", "colz");
    };

    drawHistogram(
        "ar:tan0", "angdiff_r",
        Form("Angle difference radial;\
            #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}} (PL%.3d);\
            #Deltatan#it{#theta}_{radial} (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
        1, angle_bin, 0, angle_max, 200, -ang_range_r, ang_range_r
    );
    drawHistogram(
        "al:tan0", "angdiff_l",
        Form("Angle difference lateral;\
            #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}} (PL%.3d);\
            #Deltatan#it{#theta}_{lateral} (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
        2, angle_bin, 0, angle_max, 200, -ang_range_l, ang_range_l
    );
    drawHistogram(
        "dr:tan0", "posdiff_r",
        Form("Position difference radial;\
            #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}} (PL%.3d);\
            #Deltaradial [#mum] (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
        3, angle_bin, 0, angle_max, 200, -pos_range_r, pos_range_r
    );
    drawHistogram(
        "dl:tan0", "posdiff_l",
        Form("Position difference lateral;\
            #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}} (PL%.3d);\
            #Deltalateral [#mum] (PL%.3d#minus PL%.3d)", pl0, pl1, pl0),
        4, angle_bin, 0, angle_max, 200, -pos_range_l, pos_range_l
    );
}

void difference_xyrl(TCanvas *c1, TTree *tree, double gap, uint32_t pl0, uint32_t pl1) noexcept
{
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleOffset(1.4,"y");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.235);
        c1->GetPad(pad)->SetLeftMargin(0.23);
    }

    gap = std::abs(gap);
	double ang_max_xy = 0.02;
	double ang_range_r = 0.06;
	double ang_range_l = 0.015;
	double pos_max_xy = gap * 0.02;
	double pos_range_r = gap * 0.03;
	double pos_range_l = gap * 0.02;

    auto draw2DHist = [&](const char* draw_expr, const char* hist_name, const TString& title, int pad,
                          int xbins, double xmin, double xmax,
                          int ybins, double ymin, double ymax) {
        TH2D* hist = new TH2D(hist_name, title, xbins, xmin, xmax, ybins, ymin, ymax);
        c1->cd(pad);
        draw_expr = Form("%s >> %s", draw_expr, hist_name);
        tree->Draw(draw_expr, "", "colz");
    };

    draw2DHist(
        "day:dax",
        "angdiff_xy",
        Form("Angle difference xy;#Deltatan#it{#theta}_{x} (PL%.3d#minus PL%.3d);\
            #Deltatan#it{#theta}_{y} (PL%.3d#minus PL%.3d)", pl1, pl0, pl1, pl0),
        1, 100, -ang_max_xy, ang_max_xy, 100, -ang_max_xy, ang_max_xy
    );
    draw2DHist(
        "ar:al",
        "angdiff_rl",
        Form("Angle difference rl;#Deltatan#it{#theta}_{lateral} (PL%.3d#minus PL%.3d);\
            #Deltatan#it{#theta}_{radial} (PL%.3d#minus PL%.3d)", pl1, pl0, pl1, pl0),
        2, 100, -ang_range_l, ang_range_l, 100, -ang_range_r, ang_range_r
    );
    draw2DHist(
        "dy:dx",
        "posdiff_xy",
        Form("Position difference xy;#Deltax [#mum] (PL%.3d#minus PL%.3d);\
            #Deltay [#mum] (PL%.3d#minus PL%.3d)", pl1, pl0, pl1, pl0),
        3, 100, -pos_max_xy, pos_max_xy, 100, -pos_max_xy, pos_max_xy
    );
    draw2DHist(
        "dr:dl",
        "posdiff_rl",
        Form("Position difference rl;#Deltalateral [#mum] (PL%.3d#minus PL%.3d);\
            #Deltaradial [#mum] (PL%.3d#minus PL%.3d)", pl1, pl0, pl1, pl0),
        4, 100, -pos_range_l, pos_range_l, 100, -pos_range_r, pos_range_r
    );
}

void difference_1D(TCanvas *c1, TTree *tree, int type, const std::vector<double> &angle_list,
    int start_num, std::vector<double> &mean, std::vector<double> &sigma, std::vector<double> &sigma_err,
    bool &noEntries, double gap, uint32_t pl0, uint32_t pl1) noexcept
{
    double ang_low, ang_up, range, fit_range, fit_center;
	int Entries[6], histcolor;
	TString formula, name, title, draw;
    const char* sqrt_tan = "#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}";
	TCut angle_cut;
	noEntries = false;

	// formula
    switch (type) {
        case 0: formula = "dax"; break;
        case 1: formula = "dx"; break;
        case 2: formula = "day"; break;
        case 3: formula = "dy"; break;
        case 4: formula = "ar"; break;
        case 5: formula = "dr"; break;
        case 6: formula = "al"; break;
        case 7: formula = "dl"; break;
        default: formula = ""; break;
    }

	TH1D** hist = new TH1D*[6];
	TLegend *lg[6];

	gStyle->SetTitleOffset(1.1,"x");

	c1->Divide(3, 2);
    for (int pad = 1; pad <= 6; ++pad) {
         c1->GetPad(pad)->SetRightMargin(0.1);
       c1->GetPad(pad)->SetLeftMargin(0.11);
    }

    gap = std::abs(gap);

    // plot
    for (int i = 0; i < 6; ++i)
	{
        int n = start_num + i;
        if (n + 1 >= static_cast<int>(angle_list.size())) { noEntries = true; break; }
        ang_low = angle_list[n];
        ang_up = angle_list[n + 1];

		// histgram range, title, angle range
        switch (type) {
            case 0: // ax
                histcolor = 90;
                if (gap < 3000)
                    range = 0.012 * ang_up * ang_up + 0.03 * ang_up + 0.01;
                else
                    range = 0.03 * ang_up * ang_up + 0.012 * ang_up + 0.02;
                title = Form("Angle difference x (%.1f #leq |tan#it{#theta}_{x}| (PL%.3d) < %.1f);\
                    #Deltatan#it{#theta}_{x} (PL%.3d#minus PL%.3d);", ang_low, pl0, ang_up, pl1, pl0);
                angle_cut = Form("TMath::Abs(ax0)>=%.1f&&TMath::Abs(ax0)<%.1f", ang_low, ang_up);
                break;
            case 1: // x
                histcolor = 91;
                if (gap < 3000)
                    range = (0.0056 * ang_up * ang_up + 0.040 * ang_up + 0.007) * gap;
                else
                    range = (0.02 * ang_up * ang_up + 0.02 * ang_up + 0.05) * gap;
                title = Form("Position difference x (%.1f #leq |tan#it{#theta}_{x}| (PL%.3d) < %.1f);\
                    #Deltax [#mum] (PL%.3d#minus PL%.3d);", ang_low, pl0, ang_up, pl1, pl0);
                angle_cut = Form("TMath::Abs(ax0)>=%.1f&&TMath::Abs(ax0)<%.1f", ang_low, ang_up);
                break;
            case 2: // ay
                histcolor = 90;
                if (gap < 3000)
                    range = 0.012 * ang_up * ang_up + 0.03 * ang_up + 0.01;
                else
                    range = (0.012 * ang_up * ang_up + 0.03 * ang_up + 0.02) * gap * 0.0003;
                title = Form("Angle difference y (%.1f #leq |tan#it{#theta}_{y}| (PL%.3d) < %.1f);\
                    #Deltatan#it{#theta}_{y} (PL%.3d#minus PL%.3d);", ang_low, pl0, ang_up, pl1, pl0);
                angle_cut = Form("TMath::Abs(ay0)>=%.1f&&TMath::Abs(ay0)<%.1f", ang_low, ang_up);
                break;
            case 3: // y
                histcolor = 91;
                if (gap < 3000)
                    range = (0.0063 * ang_up * ang_up + 0.03 * ang_up + 0.018) * gap;
                else
                    range = (0.0109 * ang_up * ang_up + 0.0323 * ang_up + 0.02) * gap;
                title = Form("Position difference y (%.1f #leq |tan#it{#theta}_{y}| (PL%.3d) < %.1f);\
                    #Deltay [#mum] (PL%.3d#minus PL%.3d);", ang_low, pl0, ang_up, pl1, pl0);
                angle_cut = Form("TMath::Abs(ay0)>=%.1f&&TMath::Abs(ay0)<%.1f", ang_low, ang_up);
                break;
            case 4: // ar
                histcolor = 90;
                if (gap < 3000)
                    range = 0.012 * ang_up * ang_up + 0.03 * ang_up + 0.01;
                else
                    range = (range = 0.012 * ang_up * ang_up + 0.03 * ang_up + 0.01) * gap * 0.00035;
                title = Form(
                    "Angle difference r (%.1f #leq %s(PL%.3d) < %.1f);\
                    #Deltatan#it{#theta}_{radial} (PL%.3d#minus PL%.3d);",
                    ang_low, sqrt_tan, pl0, ang_up, pl1, pl0);
                angle_cut = Form("tan0>=%.1f&&tan0<%.1f", ang_low, ang_up);
                break;
            case 5: // r
                histcolor = 91;
                if (gap < 3000)
                    range = (0.007 * ang_up + 0.02) * gap;
                else
                    range = (0.007 * ang_up + 0.02) * gap;
                title = Form("Position difference r (%.1f #leq %s(PL%.3d) < %.1f);\
                    #Deltaradial [#mum] (PL%.3d#minus PL%.3d);",
                    ang_low, sqrt_tan, pl0, ang_up, pl1, pl0);
                angle_cut = Form("tan0>=%.1f&&tan0<%.1f", ang_low, ang_up);
                break;
            case 6: // al
                histcolor = 90;
                if (gap < 3000)
                    range = 0.001 * ang_up * ang_up + 0.0005 * ang_up + 0.015;
                else
                    range = 0.01 * ang_up + 0.02;
                title = Form("Angle difference l (%.1f #leq %s(PL%.3d) < %.1f);\
                    #Deltatan#it{#theta}_{lateral} (PL%.3d#minus PL%.3d);",
                    ang_low, sqrt_tan, pl0, ang_up, pl1, pl0);
                angle_cut = Form("tan0>=%.1f&&tan0<%.1f", ang_low, ang_up);
                break;
            case 7: // l
                histcolor = 91;
                if (gap < 3000)
                    range = (0.0007 * ang_up + 0.0107) * gap;
                else
                    if (ang_up < 0.25)
                        range = 100;
                    else
                        range = 200;
                title = Form("Position difference l (%.1f #leq %s(PL%.3d) < %.1f);\
                    #Deltalateral [#mum] (PL%.3d#minus PL%.3d);",
                    ang_low, sqrt_tan, pl0, ang_up, pl1, pl0);
                angle_cut = Form("tan0>=%.1f&&tan0<%.1f", ang_low, ang_up);
                break;
            default:
                break;
        }

		name = Form("hist[%d]", i);
		draw = formula + ">>" + name;

		hist[i] = new TH1D(name, title, 100, -1 * range, range);
		tree->Draw(draw, angle_cut, "goff");

		if (hist[i]->GetEntries() < 1) { noEntries = true; break; }

		c1->cd(i + 1);
		hist[i]->Draw();
        hist[i]->SetFillColorAlpha(histcolor, 0.7);
        hist[i]->SetLineColor(histcolor);

		if (hist[i]->GetEntries() < 100) continue;

		fit_range = 3 * hist[i]->GetStdDev();
		fit_center = hist[i]->GetMean();
        TF1* gaus = new TF1("gaus", "gaus", -1 * fit_range + fit_center, fit_range + fit_center);

		// first fitting range
        switch (type) {
            case 0: // ax
                fit_range *= (ang_low < 2.0) ? 1.5 : 2.0;
                break;
            case 1: // x
                fit_range *= (ang_low < 2.0) ? 1.5 : 2.0;
                break;
            case 2: // ay
                if (ang_low > 0.3) fit_range *= 1.5;
                break;
            case 3: // y
                if (ang_low > 0.9) fit_range *= 1.5;
                break;
            case 4: // ar
                fit_range *= 1.5;
                break;
            case 5: // r
                fit_range *= 1.5;
                break;
            case 6: // al
                if (ang_low < 0.3) fit_range *= 2.0;
                break;
            case 7: // l
                fit_range *= 1.5;
                break;
            default:
                break;
        }

		hist[i]->Fit("gaus", "q", "", -1 * fit_range + fit_center, fit_range + fit_center);

		// repeat fitting
		for (int j = 0; j < 5; ++j)
		{
            mean[n] = gaus->GetParameter(1);
            sigma[n] = gaus->GetParameter(2);

			fit_range = sigma[n];
			fit_center = mean[n];

            switch (type) {
                case 0: // ax
                    if (ang_low < 2.0) fit_range *= 1.5;
                    else fit_range *= 2.0;
                    break;
                case 1: // x
                    if (ang_low < 2.0) fit_range *= 1.5;
                    else fit_range *= 2.0;
                    break;
                case 2: // ay
                    if (ang_low > 0.3) fit_range *= 1.5;
                    break;
                case 3: // y
                    if (ang_low < 0.1 && hist[i]->GetStdDev() > 20) fit_range *= -0.15 * j + 1.0;
                    else if (ang_low < 0.6 && hist[i]->GetStdDev() > 25) fit_range *= 1.0;
                    else if (ang_low < 1.0 && hist[i]->GetStdDev() > 25) fit_range *= 1.2;
                    else fit_range *= 1.5;
                    break;
                case 4: // ar
                    fit_range *= 1.5;
                    break;
                case 5: // r
                    fit_range *= 1.5;
                    break;
                case 6: // al
                    if (ang_low < 0.3) fit_range *= 2.0;
                    else fit_range *= 1.2;
                    break;
                case 7: // l
                    fit_range *= 1.2;
                    break;
            }

			hist[i]->Fit("gaus", "q", "", -1 * fit_range + fit_center, fit_range + fit_center);
		}

		mean[n] = gaus->GetParameter(1);
		sigma[n] = gaus->GetParameter(2);
		sigma_err[n] = gaus->GetParError(2);

        double range_7sigma = 7 * sigma[n];
        range = std::min(range, range_7sigma);
        // hist生成時のrangeをはみ出ないよう1ビンだけ内側に範囲を狭める
        double binw = hist[i]->GetXaxis()->GetBinWidth(1);
        double display_range = range;
        if (display_range > binw) display_range -= binw;
        hist[i]->GetXaxis()->SetRangeUser(-1 * display_range, display_range);

        gaus->SetRange(hist[i]->GetXaxis()->GetXmin(), hist[i]->GetXaxis()->GetXmax());
        gaus->SetLineColor(1);
        gaus->Draw("same");

		Entries[i] = hist[i]->GetEntries();

		lg[i] = new TLegend(0.1, 0.7, 0.27, 0.86);

        lg[i]->SetFillStyle(0);
        lg[i]->SetBorderSize(0);
        lg[i]->SetTextSize(0.04);
        lg[i]->AddEntry(hist[i], Form("Entries   %d", Entries[i]), "");
        lg[i]->AddEntry(hist[i], Form("Mean    %.2g", mean[n]), "");
        lg[i]->AddEntry(hist[i], Form("Sigma   %.2g", sigma[n]), "");
        lg[i]->Draw();
	}
}

void deviation(TCanvas *c1, TTree *tree, const double angle_max, const std::vector<double>& angle_list,
    const std::vector<std::vector<double>>& sigma, const std::vector<std::vector<double>>& sigma_err,
    double gap, uint8_t pl0, bool xy, bool dividedSqrt2) noexcept
{
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleOffset(1.5,"y");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
        c1->GetPad(pad)->SetLeftMargin(0.15);
        c1->GetPad(pad)->SetBottomMargin(0.15);
    }

    double max_ax, max_ay, max_x, max_y, max_ar, max_al, max_r, max_l;
    gap = std::abs(gap);

    if (dividedSqrt2) {
        if (gap < 3000) {
            max_ax = max_ay = 0.08;
            max_x = max_y = 40;
            max_ar = 0.08;
            max_al = 0.004;
            max_r = 10;
            max_l = 3;
        } else {
            max_ax = max_ay = 0.16;
            max_x = max_y = 350;
            max_ar = 0.20;
            max_al = 0.004;
            max_r = 60;
            max_l = 20;
        }
    } else {
        if (gap < 3000) {
            max_ax = max_ay = 0.12;
            max_x = max_y = 50;
            max_ar = 0.12;
            max_al = 0.006;
            max_r = 20;
            max_l = 4;
        } else {
            max_ax = max_ay = 0.24;
            max_x = max_y = 500;
            max_ar = 0.24;
            max_al = 0.006;
            max_r = 90;
            max_l = 30;
        }
    }

    TGraphErrors* dev[4];

    const int n_ranges = static_cast<int>(angle_list.size()) - 1;
    std::vector<double> x(n_ranges), y(n_ranges), ex(n_ranges), ey(n_ranges);
    double bin_offset;

    for (int t = 0; t < 4; ++t)
    {
        for (int i = 0; i < n_ranges; ++i)
        {
            double low = angle_list[i];
            double high = angle_list[i + 1];
            double width = high - low;
            if (width <= 0.1001) bin_offset = width * 0.5; // ~0.05 for 0.1 bins
            else bin_offset = width * 0.5; // ~0.1 for 0.2 bins

            x[i] = (low + high) * 0.5;
            y[i] = sigma[xy ? t : t + 4][i];
            ex[i] = bin_offset;
            ey[i] = sigma_err[xy ? t : t + 4][i];

            if (dividedSqrt2)
            {
                y[i] *= 0.70710678118654752440084436210485;
                ey[i] *= 0.70710678118654752440084436210485;
            }
		}

        dev[t] = new TGraphErrors(n_ranges, x.data(), y.data(), ex.data(), ey.data());
    }

    auto drawDeviation = [&](int pad, const char* hist_name, const char* title,
        double max_yaxis, int graph_idx) {
        TH2D* hist = new TH2D(hist_name, title, 10, 0, angle_max, 10, 0, max_yaxis);
        c1->cd(pad);
        tree->Draw(Form("x0:x1 >> %s", hist_name), "", "axis");
        dev[graph_idx]->Draw("same E Z 0");
        dev[graph_idx]->SetLineColor(93);
        dev[graph_idx]->SetLineWidth(1);
    };

    const char* signdiv = dividedSqrt2 ? " /^{ }#sqrt{2}" : "";

    if (xy) {
        drawDeviation(1, "angdev_x",
            Form("Angle deviation^{ }#sigma%s ( x );|tan#it{#theta}_{x}| (PL%.3d);\
            Angle deviation^{ }#sigma (#Deltatan#it{#theta}_{x} )", signdiv, pl0), max_ax, 0);
        drawDeviation(2, "angdev_y",
            Form("Angle deviation^{ }#sigma%s ( y );|tan#it{#theta}_{y}| (PL%.3d);\
            Angle deviation^{ }#sigma (#Deltatan#it{#theta}_{y} )", signdiv, pl0), max_ay, 2);
        drawDeviation(3, "posdev_x",
            Form("Position deviation^{ }#sigma%s ( x );|tan#it{#theta}_{x}| (PL%.3d);\
            Position deviation^{ }#sigma (#Deltax [#mum] )", signdiv, pl0), max_x, 1);
        drawDeviation(4, "posdev_y",
            Form("Position deviation^{ }#sigma%s ( y );|tan#it{#theta}_{y}| (PL%.3d);\
            Position deviation^{ }#sigma (#Deltay [#mum] )", signdiv, pl0), max_y, 3);
    } else {
        drawDeviation(1, "angdev_r",
            Form("Angle deviation^{ }#sigma%s ( radial );%s (PL%.3d);\
            Angle deviation^{ }#sigma (#Deltatan#it{#theta}_{radial} )", signdiv,
            "#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}", pl0), max_ar, 0);
        drawDeviation(2, "angdev_l",
            Form("Angle deviation^{ }#sigma%s ( lateral );%s (PL%.3d);\
            Angle deviation^{ }#sigma (#Deltatan#it{#theta}_{lateral} )", signdiv,
            "#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}", pl0), max_al, 2);
        drawDeviation(3, "posdev_r",
            Form("Position deviation^{ }#sigma%s ( radial );%s (PL%.3d);\
            Position deviation^{ }#sigma (#Deltaradial [#mum] )", signdiv,
            "#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}", pl0), max_r, 1);
        drawDeviation(4, "posdev_l",
            Form("Position deviation^{ }#sigma%s ( lateral );%s (PL%.3d);\
            Position deviation^{ }#sigma (#Deltalateral [#mum] )", signdiv,
            "#sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}}", pl0), max_l, 3);
    }
}

void md_oa(TCanvas *c1, TTree *tree, uint8_t pl, const size_t entries,
    const double angle_max, const double angle_resolution) noexcept
{
	gStyle->SetTitleOffset(1.2,"x");
	gStyle->SetTitleOffset(1.3,"y");

	c1->Divide(2, 2);
    for (int pad = 1; pad <= 2; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
        c1->GetPad(pad)->SetLeftMargin(0.15);
        c1->GetPad(pad)->SetBottomMargin(0.13);
    }
    for (int pad = 3; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
        c1->GetPad(pad)->SetLeftMargin(0.15);
        c1->GetPad(pad)->SetBottomMargin(0.15);
    }

    uint32_t angle_bin = 1 / angle_resolution * angle_max;
    double md_range = 200;
    double oa_range = 0.12;

	TH1D* hist_md = new TH1D(
        "hist_md", "Minimum distance;Minimum distance [#mum];Frequency", 100, 0, md_range
    );
	TH1D* hist_oa = new TH1D(
        "hist_oa", "Opening angle;Opening angle [rad];Frequency", 100, 0, oa_range
    );

	c1->cd(1);
	tree->Draw("md >> hist_md");

	c1->cd(2);
	tree->Draw("oa >> hist_oa");

    TH2D* hist_md2 = new TH2D;
    TH2D* hist_oa2 = new TH2D;

    auto draw2DHistogram = [&](TH2D* hist, int pad, const char* yaxis, const char* hist_name,
                             const TString& title) {
        double range = (yaxis == "md") ? md_range : oa_range;
        hist = new TH2D(
            hist_name, title, angle_bin, 0, angle_max,
            100, 0, range
        );
        c1->cd(pad);
        const char* draw_expr = Form("%s:tan0 >> %s", yaxis, hist_name);
        tree->Draw(draw_expr, "", "colz");
    };

    draw2DHistogram(
        hist_md2, 3, "md", "hist_md2",
        Form("Minimum distance;\
        #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}} (PL%.3d);\
        Minimum distance [#mum]", pl)
    );
    draw2DHistogram(
        hist_oa2, 4, "oa", "hist_oa2",
        Form("Opening angle;\
        #sqrt{tan^{2}#it{#theta}_{x} #plus tan^{2}#it{#theta}_{y}} (PL%.3d);\
        Opening angle [rad]", pl)
    );

    TLatex show_entries;
    show_entries.SetTextAlign(22);
    show_entries.SetTextSize(0.045);
    show_entries.SetTextColor(global_darkmode ? 0 : 1);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->cd(pad);
        show_entries.DrawLatexNDC((pad <= 2) ? 0.75 : 0.8, 0.96, Form("Entries %d", entries));
    }
}

void dxdydz(TCanvas *c1, TTree *tree, const double *AreaParam) noexcept
{
    int bin = static_cast<int>(AreaParam[5]);
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];
    double pitch = AreaParam[6];
    double pitch_half = pitch * 0.5;

    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    TString dz_title = Form(
        "Gap at each %.1f^{ }#times %.1f mm^{2};x [mm];y [mm];"
        "Average of Gap [#mum] at each %.1f^{ }#times %.1f mm^{2}",
        pitch, pitch, pitch, pitch
    );
    TString dz_1D_title = Form(
        "Gap at each %.1f^{ }#times %.1f mm^{2}     ;"
        "Average of Gap [#mum] at each %.1f^{ }#times %.1f mm^{2};Frequency",
        pitch, pitch, pitch, pitch
    );
    TH2D* dz_2D = new TH2D("dz_2D", dz_title, bin, LowX*0.001, UpX *0.001, bin, LowY*0.001, UpY*0.001);
    TH1D* dz_1D = new TH1D("dz_1D", dz_1D_title, 20000000, -10000000, 10000000); // 100 m が上限
    TH1D* dz_temp = new TH1D("dz_temp", "dz_temp", 20000000, -10000000, 10000000); // 100 m が上限

    TString dx_title = Form(
        "#Deltax (1.0 < tan#it{#theta} < 1.1, at mid-plane);x [mm];y [mm];"
        "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2}",
        pitch, pitch
    );
    TString dx_1D_title = Form(
        "#Deltax (1.0 < tan#it{#theta} < 1.1, at mid-plane);"
        "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
        pitch, pitch
    );
    TH2D* dx_2D = new TH2D("dx_2D", dx_title, bin, LowX*0.001, UpX *0.001, bin, LowY*0.001, UpY*0.001);
    TH1D* dx_1D = new TH1D("dx_1D", dx_1D_title, 8000, -1000, 1000);
    TH1D* dx_temp = new TH1D("dx_temp", "dx_temp", 100, -5000, 5000);

    TString dy_title = Form(
        "#Deltay (1.0 < tan#it{#theta} < 1.1, at mid-plane)     ;x [mm];y [mm];"
        "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2}",
        pitch, pitch
    );
    TString dy_1D_title = Form(
        "#Deltay (1.0 < tan#it{#theta} < 1.1, at mid-plane)     ;"
        "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
        pitch, pitch
    );
    TH2D* dy_2D = new TH2D("dy_2D", dy_title, bin, LowX*0.001, UpX *0.001, bin, LowY*0.001, UpY*0.001);
    TH1D* dy_1D = new TH1D("dy_1D", dy_1D_title, 8000, -1000, 1000);
    TH1D* dy_temp = new TH1D("dy_temp", "dy_temp", 100, -5000, 5000);

    for (int ix = 0; ix <= bin; ++ix) {
        for (int iy = 0; iy <= bin; ++iy) {
            // 各ビンの領域を取得
            double xcenter = dz_2D->GetXaxis()->GetBinCenter(ix);
            double ycenter = dz_2D->GetYaxis()->GetBinCenter(iy);
            TCut area = Form(
                "(%f-x*0.001)*(%f-x*0.001)<%f*%f && (%f-y*0.001)*(%f-y*0.001)<%f*%f",
                xcenter, xcenter, pitch_half, pitch_half,
                ycenter, ycenter, pitch_half, pitch_half
            );

            // 各ビンの領域に対してdz, dx, dyのヒストグラムを作成
            tree->Draw("dz >> dz_temp", area, "goff");
            tree->Draw("dx >> dx_temp", area, "goff");
            tree->Draw("dy >> dy_temp", area, "goff");

            if (dz_temp->GetEntries() == 0) continue;

            // 各ヒストグラムの平均値を取得し、2Dヒストグラムと1Dヒストグラムに格納
            double area_gap = dz_temp->GetMean();
            dz_1D->Fill(area_gap);
            dz_2D->SetBinContent(ix, iy, area_gap);

            double dx_pitch = dx_temp->GetMean();
            dx_1D->Fill(dx_pitch);
            dx_2D->SetBinContent(ix, iy, dx_pitch);

            double dy_pitch = dy_temp->GetMean();
            dy_1D->Fill(dy_pitch);
            dy_2D->SetBinContent(ix, iy, dy_pitch);
        }
    }
    gDirectory->Delete("d*_temp");

    // ヒストグラムの範囲設定用のラムダ式
    auto setRange = [](TH1D* hist, TH2D* hist2D, double mean, double range) {
        hist->GetXaxis()->SetRangeUser(mean - range, mean + range);
        hist2D->GetZaxis()->SetRangeUser(mean - range, mean + range);
    };

    // 5σの範囲を設定
    double dz_5sigma = 5 * dz_1D->GetStdDev();
    double dz_1D_mean = dz_1D->GetMean();
    double dxdy_5sigma = std::max(dx_1D->GetStdDev(), dy_1D->GetStdDev()) * 5;
    setRange(dz_1D, dz_2D, dz_1D_mean, dz_5sigma);
    setRange(dx_1D, dx_2D, 0.0, dxdy_5sigma);
    setRange(dy_1D, dy_2D, 0.0, dxdy_5sigma);

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    c1->cd(1);
    dz_2D->Draw("colz1");

    c1->cd(2);
    dz_1D->SetFillStyle(0);
    dz_1D->SetLineWidth(2);
    dz_1D->Draw();
    MyUtil::PaintBins(dz_1D, dz_1D_mean - dz_5sigma, dz_1D_mean + dz_5sigma); // 各ビンをカラーパレットの色で塗る

    TLegend* dz_lg = new TLegend(0.68, 0.7, 0.9, 0.9);
    dz_lg->SetFillStyle(0);
    dz_lg->SetBorderSize(0);
    dz_lg->SetTextSize(0.04);
    dz_lg->SetTextColor(global_darkmode ? 0 : 1);
    dz_lg->AddEntry(dz_1D, Form("Areas      %.0f", dz_1D->GetEntries()), "");
    dz_lg->AddEntry(dz_1D, Form("Mean      %.2f [#mum]", dz_1D_mean), "");
    dz_lg->AddEntry(dz_1D, Form("Std Dev   %.2f [#mum]", dz_1D->GetStdDev()), "");
    dz_lg->Draw();

    c1->cd(3);
    dx_2D->Draw("colz1"); // colz1は0のビンを塗りつぶさない

    c1->cd(4);
    dy_2D->Draw("colz1"); // colz1は0のビンを塗りつぶさない
}

void dxdy(TCanvas *c1, TTree *tree, const double *AreaParam) noexcept
{
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    // gDirectoryからヒストグラムを取得。なければ新規作成
    TH2D* dx_2D = (TH2D*)gDirectory->Get("dx_2D");
    TH1D* dx_1D = (TH1D*)gDirectory->Get("dx_1D");
    TH2D* dy_2D = (TH2D*)gDirectory->Get("dy_2D");
    TH1D* dy_1D = (TH1D*)gDirectory->Get("dy_1D");
    if (!dx_2D || !dx_1D || !dy_2D || !dy_1D) {
        int bin = static_cast<int>(AreaParam[5]);
        double LowX = AreaParam[1];
        double UpX  = AreaParam[2];
        double LowY = AreaParam[3];
        double UpY  = AreaParam[4];
        double pitch = AreaParam[6];
        double pitch_half = pitch * 0.5;

        TString dx_title = Form(
            "#Deltax (1.0 < tan#it{#theta} < 1.1, at mid-plane);x [mm];y [mm];"
            "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2}",
            pitch, pitch
        );
        TString dx_1D_title = Form(
            "#Deltax (1.0 < tan#it{#theta} < 1.1, at mid-plane);"
            "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
            pitch, pitch
        );
        TH1D* dx_temp = new TH1D("dx_temp", "dx_temp", 100, -1000, 1000);
        dx_2D = new TH2D("dx_2D", dx_title, bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001);
        dx_1D = new TH1D("dx_1D", dx_1D_title, 500, -100, 100);

        TString dy_title = Form(
            "#Deltay (1.0 < tan#it{#theta} < 1.1, at mid-plane)     ;x [mm];y [mm];"
            "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2}",
            pitch, pitch
        );
        TString dy_1D_title = Form(
            "#Deltay (1.0 < tan#it{#theta} < 1.1, at mid-plane)     ;"
            "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
            pitch, pitch
        );
        dy_2D = new TH2D("dy_2D", dy_title, bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001);
        dy_1D = new TH1D("dy_1D", dy_1D_title, 500, -100, 100);
        TH1D* dy_temp = new TH1D("dy_temp", "dy_temp", 100, -1000, 1000);

        for (int ix = 0; ix <= bin; ++ix) {
            for (int iy = 0; iy <= bin; ++iy) {
                double xcenter = dx_2D->GetXaxis()->GetBinCenter(ix);
                double ycenter = dx_2D->GetYaxis()->GetBinCenter(iy);
                TCut area = Form(
                    "(%f-x*0.001)*(%f-x*0.001)<%f*%f && (%f-y*0.001)*(%f-y*0.001)<%f*%f",
                    xcenter, xcenter, pitch_half, pitch_half,
                    ycenter, ycenter, pitch_half, pitch_half
                );

                tree->Draw("dx >> dx_temp", area, "goff");
                tree->Draw("dy >> dy_temp", area, "goff");

                if (dx_temp->GetEntries() == 0) continue;

                double dx_pitch = dx_temp->GetMean();
                dx_1D->Fill(dx_pitch);
                dx_2D->SetBinContent(ix, iy, dx_pitch);

                double dy_pitch = dy_temp->GetMean();
                dy_1D->Fill(dy_pitch);
                dy_2D->SetBinContent(ix, iy, dy_pitch);
            }
        }
        gDirectory->Delete("d*_temp");

        // ヒストグラムの範囲設定用のラムダ式
        auto setRange = [](TH1D* hist, TH2D* hist2D, double mean, double range) {
            hist->GetXaxis()->SetRangeUser(mean - range, mean + range);
            hist2D->GetZaxis()->SetRangeUser(mean - range, mean + range);
        };

        double dxdy_5sigma = std::max(dx_1D->GetStdDev(), dy_1D->GetStdDev()) * 5;
        setRange(dx_1D, dx_2D, 0.0, dxdy_5sigma);
        setRange(dy_1D, dy_2D, 0.0, dxdy_5sigma);
    }

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    c1->cd(1);
    dx_2D->Draw("colz1"); // colz1は0のビンを塗りつぶさない

    c1->cd(2);
    dy_2D->Draw("colz1"); // colz1は0のビンを塗りつぶさない

    // ヒストグラムとMeanなどの情報を描画するためのラムダ式
    auto drawHistogramWithLegend = [](
        TCanvas* canvas, int pad, TH1D* hist, double hist_min, double hist_max,
        const char* legend_title, double legend_x1, double legend_y1, double legend_x2, double legend_y2
    ) {
        canvas->cd(pad);
        hist->SetFillStyle(0);
        hist->SetLineWidth(2);
        hist->Draw();
        MyUtil::PaintBins(hist, hist_min, hist_max); // 各ビンをカラーパレットの色で塗る

        TLegend* legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->SetTextColor(global_darkmode ? 0 : 1);
        legend->AddEntry(hist, Form("Areas      %.0f", hist->GetEntries()), "");
        legend->AddEntry(hist, Form("Mean      %.2f [#mum]", hist->GetMean()), "");
        legend->AddEntry(hist, Form("Std Dev   %.2f [#mum]", hist->GetStdDev()), "");
        legend->Draw();
    };

    // dxとdyの5σのうち、大きい方を基準にして範囲を設定
    double dxdy_5sigma = std::max(dx_1D->GetStdDev(), dy_1D->GetStdDev()) * 5;

    drawHistogramWithLegend(c1, 3, dx_1D, -dxdy_5sigma, dxdy_5sigma, "dx_1D", 0.74, 0.7, 0.96, 0.9);
    drawHistogramWithLegend(c1, 4, dy_1D, -dxdy_5sigma, dxdy_5sigma, "dy_1D", 0.68, 0.7, 0.9, 0.9);
}
