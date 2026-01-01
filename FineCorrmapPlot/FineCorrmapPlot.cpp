#include <iostream>
#include <fstream>
#include <sstream>
#include <csignal>
#include <queue>
#include <set>
#include <algorithm>
#include <cmath>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TError.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TEllipse.h>
#include <TGaxis.h>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <argparse/argparse.hpp>

#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>

void sig_map(TCanvas *c1, TTree *tree, const double *AreaParam_d, const int *AreaParam_i,
    const std::vector<std::pair<int, int>> &hole_area, uint32_t entries, double edgecut) noexcept;
void gap_rotz(TCanvas *c1, TTree *tree, const double *AreaParam_d,
    uint32_t entries, double nominal_gap) noexcept;
void shift_xy(TCanvas *c1, TTree *tree, const double *AreaParam_d, uint32_t entries) noexcept;
void shrink(TCanvas *c1, TTree *tree, const double *AreaParam_d, uint32_t entries) noexcept;

namespace {
    // Ctrl+Cで終了したときの処理用にグローバル変数を定義
    TCanvas *global_c1 = nullptr;
    std::string global_output = "";

    // その他のグローバル変数
    bool global_darkmode = false;
}

// Ctrl+Cで終了したときの処理
void handleSIGINT(int) {
    std::cout << "\n** *Catched Ctrl+C: ";

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

void parse_arguments(argparse::ArgumentParser& parser, int argc, char *argv[]) {
    parser.set_usage_max_line_width(80);
    parser.add_description("Tips: You can combine single-character arguments.\n"
        "      For example, \"-d -p kBird -i\" is equivalent to \"-dpi kBird\".");

    // 必須引数: corrmapファイル
    parser.add_argument("fine_corrmap")
        .help("Path to the fine corrmap file to be processed.")
        .required();
    // オプション引数: 出力ファイル名
    parser.add_argument("-o", "--output")
        .help("Output PDF file name. [default: local_corrmap.pdf]")
        .default_value(std::string());
    // オプション引数: その他の設定
    parser.add_argument("-e", "--edge_cut")
        .help("Edge cut value in hole map (mm).")
        .default_value(10.0)
        .scan<'g', double>();
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

    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << "\nError: " << err.what() << std::endl;
        std::cerr << parser;
        std::cerr << "\nError: " << err.what() << std::endl;
        std::exit(1);
    }
}

int main(int argc, char *argv[])
{
    // Ctrl+Cで終了したときの処理を設定
    std::signal(SIGINT, handleSIGINT);

    // 処理時間を計測
    TStopwatch sw;

    // argparseを使用して引数を解析
    std::cout << "\nInitializing..." << std::endl;
    argparse::ArgumentParser parser("CorrmapPlot.exe", "1.0.0");
    parse_arguments(parser, argc, argv);

    // 引数を取得
    const auto corrmapfile = parser.get<std::string>("fine_corrmap");
    const auto output_arg = parser.get<std::string>("--output");
    const auto edgecut = parser.get<double>("--edge_cut");
    auto font_number = parser.get<int>("--font_number");
    auto hideGrid = parser.get<bool>("--hide_grid");
    global_darkmode = parser.get<bool>("--dark_mode");
    auto palette_arg = parser.get<std::string>("--palette");
    auto NContours = parser.get<int>("--contours");
    auto invertpalette = parser.get<bool>("--invert_palette");
    auto negatepalette = parser.get<bool>("--negate_palette");

    // 出力ファイル名の設定
    const std::string output = (output_arg.empty())
        ? (corrmapfile + ".pdf")
        : ((output_arg.size() > 4 && output_arg.substr(output_arg.size() - 4) == ".pdf")
            ? output_arg
            : output_arg + ".pdf");

    // 引数を利用する変数を設定
    const int font_code = 10 * font_number + 2;

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
    gROOT->GetColor(96)->SetRGB(r1 * 0.9, g1 * 0.9, b1 * 0.9);
    gROOT->GetColor(97)->SetRGB(r3 * 0.9, g3 * 0.9, b3 * 0.9);
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
    gStyle->SetOptStat(0);               // 統計boxの表示をオフ
    gStyle->SetPadRightMargin(0.1);      // Pad右側のマージン
    gStyle->SetPadLeftMargin(0.1);       // Pad左側のマージン
    gStyle->SetPadTopMargin(0.1);        // Pad上側のマージン
    gStyle->SetPadBottomMargin(0.11);    // Pad下側のマージン
    gStyle->SetLabelOffset(0.008,"xyz"); // 軸ラベル(数値)と軸の距離
    gStyle->SetTitleOffset(1.1,"xyz");   // 軸titleと軸の距離
    gStyle->SetTitleY(0.985);            // タイトルのy方向の位置
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

    // エラーメッセージ未満のROOTのメッセージを非表示に設定
    gErrorIgnoreLevel = kError;

    // ファイル読み込み準備
    std::ifstream ifs(corrmapfile.c_str());
    std::string line;

    if (!ifs) {
        std::cerr << "\n Error: Cannot open file: " << corrmapfile << std::endl;
        return 1;
    }

    int id, sig, ix, iy;
    double x, y, z, dx, dy, dz, rot_x, rot_y, rot_z,
        shr_x, shr_y, shr_z, shr_xy, // shr_x=shr_yのためプロット時にはshr_xyに統一
        yx_shear, zx_shear, zy_shear,
        dx_, dy_;

    // ファイルの読み込み準備
    int ixmin = INT_MAX;
    int ixmax = -INT_MAX;
    int iymin = INT_MAX;
    int iymax = -INT_MAX;
    uint32_t entries = 0;
    std::map<int, std::vector<double>> ix_to_x_list;
    std::set<int> ix_set;
    std::set<int> iy_set;
    double rotz_sum = 0.0;

    // TTreeを作成
    TTree *tree = new TTree("tree", "");
    TTree *tree2 = new TTree("tree2", "");
    const std::array<std::pair<const char*, void*>, 3> int_branches = {{
        {"sig", &sig}, {"ix", &ix}, {"iy", &iy}
    }};
    const std::array<std::pair<const char*, void*>, 6> double_branches = {{
        {"dx_", &dx_}, {"dy_", &dy_}, {"dz", &dz},
        {"rot_z", &rot_z}, {"shr_xy", &shr_xy}, {"shr_z", &shr_z}
    }};
    for (const auto& branch : int_branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/I").c_str());
    }
    for (const auto& branch : double_branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }
    const std::array<std::pair<const char*, void*>, 2> int_branches2 = {{
        {"ix", &ix}, {"iy", &iy}
    }};
    const std::array<std::pair<const char*, void*>, 2> double_branches2 = {{
        {"dx", &dx}, {"dy", &dy}
    }};
    for (const auto& branch : int_branches2) {
        tree2->Branch(branch.first, branch.second, (std::string(branch.first) + "/I").c_str());
    }
    for (const auto& branch : double_branches2) {
        tree2->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    // ファイルの読み込み
    while (std::getline(ifs, line)) {
        std::istringstream iss(line);
        if (!(iss >> id >> ix >> iy >> sig >> x >> y >> z
                  >> rot_x >> rot_y >> rot_z >> shr_xy >> shr_y >> shr_z
                  >> yx_shear >> zx_shear >> zy_shear >> dx_ >> dy_ >> dz)) {
            std::cerr << "Warning : Malformed line skipped: " << line << std::endl;
            continue;
        }

        if (ix < ixmin) ixmin = ix;
        if (ix > ixmax) ixmax = ix;
        if (iy < iymin) iymin = iy;
        if (iy > iymax) iymax = iy;

        rotz_sum += rot_z;

        x *= 0.001; // mm単位に変換
        ix_to_x_list[ix].push_back(x);

        rot_z *= 1000; // mrad単位に変換

        ix_set.insert(ix);
        iy_set.insert(iy);

        tree->Fill();
        entries++;
    }
    ifs.close();

    if (entries == 0) {
        std::cerr << "\n Error : No parameter is found in " << corrmapfile << std::endl;
        return 1;
    }

    double nominal_gap = z;
    double rotz_mean = rotz_sum / entries;

    // 隣接する ix グループ間の距離 (= view_width) を計算
    std::map<int, double> ix_to_avg_x;
    double view_width = DBL_MAX;
    for (const auto& pair : ix_to_x_list) {
        double sum = 0.0;
        for (double val : pair.second) {
            sum += val;
        }
        ix_to_avg_x[pair.first] = sum / pair.second.size();
    }
    if (ix_to_avg_x.size() >= 2) {
        auto it = ix_to_avg_x.begin();
        auto prev_it = it;
        ++it;
        for (; it != ix_to_avg_x.end(); ++it, ++prev_it) {
            double diff = std::abs(it->second - prev_it->second);
            diff = std::round(diff / 0.1) * 0.1;  // 100 µm (0.1 mm) の位で丸める
            if (diff > 0.0 && diff < view_width) view_width = diff;
        }
    }

    // データが存在する領域をリストアップ
    std::vector<bool> has_data((ixmax - ixmin + 1) * (iymax - iymin + 1), false);
    for (int i = 0; i < entries; ++i) {
        tree->GetEntry(i);
        int idx = (ix - ixmin) * (iymax - iymin + 1) + (iy - iymin);
        has_data[idx] = true;
    }
    // ix_set と iy_set の全組み合わせからデータが存在しない領域を抽出
    std::vector<std::pair<int, int>> hole_area;
    for (int i : ix_set) {
        for (int j : iy_set) {
            int idx = (i - ixmin) * (iymax - iymin + 1) + (j - iymin);
            if (!has_data[idx]) {
                hole_area.emplace_back(i, j);
            }
        }
    }

    // dx, dyをグローバル座標系に変換
    for (int i = 0; i < entries; ++i) {
        tree->GetEntry(i);
        double theta = - rotz_mean;
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        dx =  cos_theta * dx_ - sin_theta * dy_;
        dy =  sin_theta * dx_ + cos_theta * dy_;
        tree2->Fill();
    }

    // 情報表示
    std::cout << "\n Corrmap File : " << corrmapfile << std::endl;
    std::cout << " # of Areas   : " << entries << std::endl;

    // プロット開始
    TDatime time_now;
    std::string Time_Now = fmt::format(
        "{:d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}",
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(),
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    sw.Start();
	std::cout << " Plot start   : " << Time_Now << std::endl;

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas *c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());
    global_c1 = c1;
    global_output = output;

	// プログレスバーの初期化
	int page = 0;
    constexpr int total = 4; // 合計ページ数
	MyUtil::ShowProgress(page, total);

    // データの座標の範囲を取得し、表示範囲とビンの数を決定する
    // フィルムの長辺の端から1cm外側までを最大表示範囲とし、縦横比を正しく保って表示する
    const int RangeIX = ixmax - ixmin; // データ領域
    const double RangeX = (RangeIX + 1) * view_width;
    const int RangeIY = iymax - iymin; // データ領域
    const double RangeY = (RangeIY + 1) * view_width;
    double LowX, UpX, LowY, UpY, bin; // 表示範囲とビンの数
    constexpr double margin = 10.0; // 10mmマージン
    int margin_num = static_cast<int>(margin / view_width);
    if (RangeIX >= RangeIY) {
        LowX = (ixmin - margin_num) * view_width;
        UpX = (ixmax + margin_num) * view_width;
        LowY = (iymin - (RangeIX - RangeIY + 2 * margin_num) * 0.5) * view_width;
        UpY = (iymax + (RangeIX - RangeIY + 2 * margin_num) * 0.5) * view_width;
        bin = RangeIX + 2 * margin_num;
    } else {
        LowX = (ixmin - (RangeIY - RangeIX + 2 * margin_num) * 0.5) * view_width;
        UpX = (ixmax + (RangeIY - RangeIX + 2 * margin_num) * 0.5) * view_width;
        LowY = (iymin - margin_num) * view_width;
        UpY = (iymax + margin_num) * view_width;
        bin = RangeIY + 2 * margin_num;
    }
    const double AreaParam_d[8] = {bin, LowX, UpX, LowY, UpY, RangeX, RangeY, view_width};
    const int AreaParam_i[3] = {ixmin, iymin, margin_num};

    sig_map(c1, tree, AreaParam_d, AreaParam_i, hole_area, entries, edgecut);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");
    gDirectory->Delete("*_2D");

    gap_rotz(c1, tree, AreaParam_d, entries, nominal_gap);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");
    gDirectory->Delete("*_2D");

    shift_xy(c1, tree2, AreaParam_d, entries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");
    gDirectory->Delete("*_2D");
    gDirectory->Delete("tree2");

    shrink(c1, tree, AreaParam_d, entries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");
    gDirectory->Delete("*_2D");
    gDirectory->Delete("tree");

    // PDFファイルを閉じる
    c1->Print((output + "]").c_str());
	if (page < total) page = total;
    MyUtil::ShowProgress(page, total);

    // プロット終了
    time_now.Set();
    Time_Now = fmt::format(
        "{:d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}",
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(),
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    double elapsed_time = sw.RealTime();
    double cpu_time = sw.CpuTime();
    std::cout << fmt::format(
        "\n Plot end     : {} - Elapsed {:.2f} [s] (CPU: {:.2f} [s])",
        Time_Now, elapsed_time, cpu_time
    ) << std::endl;
    std::cout << " Output       : " << output << std::endl;

    delete c1;
    gDirectory->Delete("tree");

    return 0;
}

void sig_map(TCanvas *c1, TTree *tree, const double *AreaParam_d, const int *AreaParam_i,
    const std::vector<std::pair<int, int>> &hole_area, uint32_t entries, double edgecut) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double bin  = AreaParam_d[0];
    double LowX = AreaParam_d[1];
    double UpX  = AreaParam_d[2];
    double LowY = AreaParam_d[3];
    double UpY  = AreaParam_d[4];
    double view_width = AreaParam_d[7];

    int ixmin = AreaParam_i[0];;
    int iymin = AreaParam_i[1];
    int margin_num = AreaParam_i[2];

    int ix, iy, sig;
    tree->SetBranchAddress("ix", &ix);
    tree->SetBranchAddress("iy", &iy);
    tree->SetBranchAddress("sig", &sig);

    // 2Dヒストグラム作成
    // alignment_map_interpolation.exeによって補完された領域には-1が入っているため、
    // 別途sig_interp_2Dを用意して描画する
    TH2D *sig_2D = new TH2D("sig_2D", ";x [mm];y [mm];Signal",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *sig_interp_2D = new TH2D("sig_interp_2D", ";x [mm];y [mm]",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *datamap_2D = new TH2D("datamap_2D", ";x [mm];y [mm]",
        bin, LowX, UpX, bin, LowY, UpY);

    for (int i = 0; i < entries; ++i) {
        tree->GetEntry(i);
        int binX = ix - ixmin + 1 + margin_num;
        int binY = iy - iymin + 1 + margin_num;
        if (sig > 0) {
            sig_2D->SetBinContent(binX, binY, sig);
            datamap_2D->SetBinContent(binX, binY, 1.0);
        }
        else if (sig == -1) sig_interp_2D->SetBinContent(binX, binY, 1.0);
    }

    // データが存在しない領域を描画するヒストグラム作成
    TH2D *hole_2D = new TH2D("hole_2D", ";x [mm];y [mm]",
        bin, LowX, UpX, bin, LowY, UpY);
    for (const auto& pair : hole_area) {
        int binX = pair.first - ixmin + 1 + margin_num;
        int binY = pair.second - iymin + 1 + margin_num;
        hole_2D->SetBinContent(binX, binY, 1.0);
    }

    // 描画
    c1->cd(1);
    sig_interp_2D->SetFillColor(93);
    sig_interp_2D->Draw("box"); // sig_2Dより先に描画しないとbox間に隙間ができる
    sig_2D->Draw("colz same0"); // sameだとsig_interp_2D基準のカラースケールになるためsame0を使用
    sig_2D->Draw("axis same");

    // タイトル
    TLatex title1;
    title1.SetTextAlign(22);
    title1.SetTextSize(0.06);
    title1.SetTextColor(global_darkmode ? 0 : 1);
    title1.DrawLatexNDC(0.5, 0.95, "Signal");

    // 1Dヒストグラム作成
    const double sig_max = tree->GetMaximum("sig");
    TH1D *sig_1D = new TH1D("sig_1D", ";Signal;Area", 100, 0.0, sig_max);

    tree->Draw("sig >> sig_1D", "", "goff");
    c1->cd(3);
    sig_1D->SetFillStyle(0);
    sig_1D->SetLineWidth(2);
    sig_1D->Draw();
    MyUtil::PaintBins(sig_1D, 0.0, sig_max); // 各ビンをカラーパレットの色で塗る

    TLatex title3;
    title3.SetTextAlign(22);
    title3.SetTextSize(0.06);
    title3.SetTextColor(global_darkmode ? 0 : 1);
    title3.DrawLatexNDC(0.5, 0.95, "Signal");

    TLegend *lg = new TLegend(0.73, 0.7, 0.95, 0.9);
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.04);
    lg->SetTextColor(global_darkmode ? 0 : 1);
    lg->AddEntry(sig_1D, fmt::format("{} areas", entries).c_str(), "");
    lg->AddEntry(sig_1D, fmt::format("Mean      {:.2f}", sig_1D->GetMean()).c_str(), "");
    lg->AddEntry(sig_1D, fmt::format("Std Dev   {:.2f}", sig_1D->GetStdDev()).c_str(), "");
    lg->Draw();

    c1->cd(2);
    sig_interp_2D->Draw("box");
    sig_2D->Draw("same0");
    sig_2D->Draw("axis same");

    c1->cd(4);
    if (sig_interp_2D->GetEntries() > 0) sig_interp_2D->Draw("box");
    datamap_2D->SetFillColor(93);
    datamap_2D->Draw("box same");
    datamap_2D->Draw("axis same");

    // edge cut 領域の描画
    double edge_LowX = LowX + margin_num * view_width + edgecut;
    double edge_UpX = UpX - margin_num * view_width - edgecut;
    double edge_LowY = LowY + margin_num * view_width + edgecut;
    double edge_UpY = UpY - margin_num * view_width - edgecut;

    TBox *edge_areas[4] = {
        new TBox(LowX, LowY, edge_LowX, edge_UpY),
        new TBox(edge_UpX, edge_LowY, UpX, UpY),
        new TBox(LowX, edge_UpY, edge_UpX, UpY),
        new TBox(edge_LowX, LowY, UpX, edge_LowY)
    };
    for (int i = 0; i < 4; ++i) {
        edge_areas[i]->SetFillStyle(1001);
        edge_areas[i]->SetFillColorAlpha(94, 0.7);
        edge_areas[i]->SetLineStyle(0);
    }

    for (int i = 1; i <= 2; ++i) {
        c1->cd(i * 2);
        for (int j = 0; j < 4; ++j) {
            edge_areas[j]->Draw("same");
        }
    }

    // 穴検出&クラスタリング
    // hole_areaからビン座標に変換して穴セルを抽出
    std::vector<std::pair<double, double>> hole_cells;

    for (const auto& pair : hole_area) {
        int binX = pair.first - ixmin + 1 + margin_num;
        int binY = pair.second - iymin + 1 + margin_num;
        double xc = hole_2D->GetXaxis()->GetBinCenter(binX);
        double yc = hole_2D->GetYaxis()->GetBinCenter(binY);
        hole_cells.push_back({xc, yc});
    }

    int holecells = hole_cells.size();

    std::vector<bool> visited(holecells, false);
    std::vector<std::vector<int>> clusters;
    double threshold_sq = view_width * 2.0;
    threshold_sq *= threshold_sq;  // 距離比較を高速化するため二乗値を事前計算

    // Breadth-First Search (BFS)を用いて各holecellをクラスタリング
    // 2つの点のx, y方向の距離が view_width * 2.0 より近ければクラスタリングする
    // edgecutで指定した範囲外の点はクラスタリングしない

    for (int i = 0; i < holecells; ++i) {
        if (visited[i]) continue;

        double x0 = hole_cells[i].first;
        double y0 = hole_cells[i].second;
        if (x0 < edge_LowX || x0 > edge_UpX || y0 < edge_LowY || y0 > edge_UpY) continue;

        visited[i] = true;
        std::vector<int> cluster;
        std::queue<int> q;
        q.push(i);

        while (!q.empty()) {
            int j = q.front();
            q.pop();
            cluster.push_back(j);

            double x1 = hole_cells[j].first;
            double y1 = hole_cells[j].second;

            for (int k = i + 1; k < holecells; ++k) {
                if (visited[k]) continue;

                double x2 = hole_cells[k].first;
                double y2 = hole_cells[k].second;

                if (x2 < edge_LowX || x2 > edge_UpX || y2 < edge_LowY || y2 > edge_UpY) {
                    continue;
                }

                // 2つの点のx、y方向の距離が view_width * 2.0 より近ければクラスタリング
                double dx = std::abs(x1 - x2);
                double dy = std::abs(y1 - y2);
                if (dx < view_width * 2.0 && dy < view_width * 2.0) {
                    visited[k] = true;
                    q.push(k);
                }
            }
        }
        clusters.push_back(cluster);
    }

    // クラスターごとに楕円を描画
    for (size_t i = 0; i < clusters.size(); ++i) {
        const std::vector<int>& cluster = clusters[i];

        // クラスター内の座標の中央値を計算
        std::vector<double> x_vals, y_vals;
        for (int idx : cluster) {
            x_vals.push_back(hole_cells[idx].first);
            y_vals.push_back(hole_cells[idx].second);
        }

        std::sort(x_vals.begin(), x_vals.end());
        std::sort(y_vals.begin(), y_vals.end());

        double centerX, centerY;
        size_t size = x_vals.size();
        if (size % 2 == 1) {
            centerX = x_vals[size / 2];
            centerY = y_vals[size / 2];
        } else {
            centerX = (x_vals[size / 2 - 1] + x_vals[size / 2]) / 2.0;
            centerY = (y_vals[size / 2 - 1] + y_vals[size / 2]) / 2.0;
        }

        // 中心からの最大距離を計算
        double dist_max = 0;
        if (cluster.size() == 1) {
            dist_max = 0.5 * view_width + 5.0;
        } else {
            for (int idx : cluster) {
                double x = hole_cells[idx].first;
                double y = hole_cells[idx].second;
                double d = std::hypot(x - centerX, y - centerY) + 0.5 * view_width + 5.0;
                if (d > dist_max) dist_max = d;
            }
        }

        // 楕円を描画
        TEllipse *ellipse = new TEllipse(centerX, centerY, dist_max, dist_max);
        ellipse->SetLineColor(global_darkmode ? 0 : 1);
        ellipse->SetLineWidth(1);
        ellipse->SetFillStyle(0);

        c1->cd(2);
        ellipse->Draw("same");
        c1->cd(4);
        ellipse->Draw("same");
    }

    TLegend *lg2 = new TLegend(0.665, 0.7, 0.93, 0.98);
    lg2->SetFillStyle(0);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.04);
    lg2->SetTextColor(global_darkmode ? 0 : 1);
    lg2->AddEntry(datamap_2D, fmt::format("Edge cut: {:.1f} mm", edgecut).c_str(), "");

    c1->cd(2);
    lg2->Draw();

    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    std::string title2_text = "Data Hole Map";
    title2->DrawLatexNDC(0.43, 0.95, title2_text.c_str());

    c1->cd(4);
    lg2->Draw();

    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    std::string title4_text = fmt::format("Number of Holes: {}", clusters.size());
    title4->DrawLatexNDC(0.43, 0.95, title4_text.c_str());
}

void gap_rotz(TCanvas *c1, TTree *tree, const double *AreaParam,
    uint32_t entries, double nominal_gap) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double bin  = AreaParam[0];
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];
    double view_width = AreaParam[7];

    // |nominal_gap|が1000 µmを超える場合mmスケールで表示
    bool mmScale = false;
    if (std::abs(nominal_gap) > 1000.0) {
        mmScale = true;
        nominal_gap /= 1000.0;
    }

    // ブランチアドレスを設定
    tree->ResetBranchAddresses();
    int ix, iy;
    double dz, rot_z;
    tree->SetBranchAddress("ix", &ix);
    tree->SetBranchAddress("iy", &iy);
    tree->SetBranchAddress("dz", &dz);
    tree->SetBranchAddress("rot_z", &rot_z);

    // 2Dヒストグラム作成
    TH2D *gap_2D = new TH2D("gap_2D",
        mmScale ? ";x [mm];y [mm];Gap [mm]" : ";x [mm];y [mm];Gap [#mum]",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *rotz_2D = new TH2D("rotz_2D", ";x [mm];y [mm];Rotation z-axis [mrad]",
        bin, LowX, UpX, bin, LowY, UpY);

    for (uint32_t i = 0; i < entries; ++i) {
        tree->GetEntry(i);
        double x = ix * view_width;
        double y = iy * view_width;
        int bx = gap_2D->GetXaxis()->FindBin(x);
        int by = gap_2D->GetYaxis()->FindBin(y);
        if (bx > 0 && bx <= bin && by > 0 && by <= bin) {
            gap_2D->SetBinContent(
                bx, by, mmScale ? dz / 1000.0 + nominal_gap : dz + nominal_gap);
            rotz_2D->SetBinContent(bx, by, rot_z);
        }
    }

    // 1binが最小単位を下回らないように調整
    const double gap_unit = mmScale ? 1e-5 : 0.01;
    constexpr double rotz_unit = 1e-4;

    TH1D *gap_temp = new TH1D("gap_temp", "", 2000000, -10000.0, 10000.0);
    tree->Draw(
        fmt::format(mmScale ? "dz / 1000.0  + {} >> gap_temp"
            : "dz + {} >> gap_temp", nominal_gap).c_str(),
        "", "goff");
    double gap_mean = gap_temp->GetMean();
    double gap_stddev = gap_temp->GetStdDev();
    double gap_low = std::max(
        gap_mean - 5 * gap_stddev,
        mmScale ? tree->GetMinimum("dz") / 1000.0 - gap_stddev + nominal_gap
                : tree->GetMinimum("dz") - gap_stddev + nominal_gap
    );
    double gap_up  = std::min(
        gap_mean + 5 * gap_stddev,
        mmScale ? tree->GetMaximum("dz") / 1000.0 + gap_stddev + nominal_gap
                : tree->GetMaximum("dz") + gap_stddev + nominal_gap
    );
    double gap_range = gap_up - gap_low;
    int gap_bin = (gap_range < gap_unit * 100)
        ? std::floor(gap_range / gap_unit)
        : 100;
    gDirectory->Delete("gap_temp");

    TH1D *rotz_temp = new TH1D("rotz_temp", "", 200000000, -10000.0, 10000.0);
    tree->Draw("rot_z >> rotz_temp", "", "goff");
    double rotz_mean = rotz_temp->GetMean();
    double rotz_stddev = rotz_temp->GetStdDev();
    double rotz_low = std::max(rotz_mean - 5 * rotz_stddev, tree->GetMinimum("rot_z") - rotz_stddev);
    double rotz_up  = std::min(rotz_mean + 5 * rotz_stddev, tree->GetMaximum("rot_z") + rotz_stddev);
    double rotz_range = rotz_up - rotz_low;
    int rotz_bin = (rotz_range < rotz_unit * 100)
        ? std::floor(rotz_range / rotz_unit)
        : 100;
        gDirectory->Delete("rotz_temp");

    gap_2D->GetZaxis()->SetRangeUser(gap_low, gap_up);
    rotz_2D->GetZaxis()->SetRangeUser(rotz_low, rotz_up);

    // 1Dヒストグラム作成
    TH1D *gap_1D = new TH1D("gap_1D",
        mmScale ? ";Gap [mm];Area" : ";Gap [#mum];Area", gap_bin, gap_low, gap_up);
    TH1D *rotz_1D = new TH1D("rotz_1D", ";Rotation z-axis [mrad];Area", rotz_bin, rotz_low, rotz_up);

    tree->Draw(fmt::format(mmScale ? "dz / 1000.0  + {}>> gap_1D"
        : "dz + {} >> gap_1D", nominal_gap).c_str(), "", "goff");
    tree->Draw("rot_z >> rotz_1D", "", "goff");

    // pad1: dz 2D分布
    c1->cd(1);
    gap_2D->Draw("colz1");
    TLatex *title1 = new TLatex();
    title1->SetTextAlign(22);
    title1->SetTextSize(0.06);
    title1->SetTextColor(global_darkmode ? 0 : 1);
    title1->DrawLatexNDC(0.5, 0.95,
        fmt::format(mmScale ? "Gap (Nominal: {:.3f} mm)"
            : "Gap (Nominal: {:.0f}^{{ }}#mum)", nominal_gap).c_str());

    // pad2: rot_z 2D分布
    c1->cd(2);
    rotz_2D->Draw("colz1");
    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    title2->DrawLatexNDC(0.43, 0.95, "Rotation z-axis");

    // pad3: dz 1D
    c1->cd(3);
    gap_1D->SetFillStyle(0);
    gap_1D->SetLineWidth(2);
    gap_1D->Draw();
    MyUtil::PaintBins(gap_1D, gap_low, gap_up);

    TLatex *title3 = new TLatex();
    title3->SetTextAlign(22);
    title3->SetTextSize(0.06);
    title3->SetTextColor(global_darkmode ? 0 : 1);
    title3->DrawLatexNDC(0.5, 0.95,
        fmt::format(mmScale ? "Gap (Nominal: {:.3f} mm)"
            : "Gap (Nominal: {:.0f}^{{ }}#mum)", nominal_gap).c_str());

    TLegend *lg1 = new TLegend(0.73, 0.7, 0.95, 0.9);
    lg1->SetFillStyle(0);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(0.04);
    lg1->SetTextColor(global_darkmode ? 0 : 1);
    lg1->AddEntry(gap_1D, fmt::format("{} areas", entries).c_str(), "");
    lg1->AddEntry(gap_1D, fmt::format("Mean      {:.2f}", gap_1D->GetMean()).c_str(), "");
    lg1->AddEntry(gap_1D, fmt::format("Std Dev   {:.2f}", gap_1D->GetStdDev()).c_str(), "");
    lg1->Draw();

    // pad4: rot_z 1D
    c1->cd(4);
    rotz_1D->SetFillStyle(0);
    rotz_1D->SetLineWidth(2);
    rotz_1D->Draw();
    MyUtil::PaintBins(rotz_1D, rotz_low, rotz_up);

    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    title4->DrawLatexNDC(0.43, 0.95, "Rotation z-axis");

    TLegend *lg2 = new TLegend(0.665, 0.7, 0.885, 0.9);
    lg2->SetFillStyle(0);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.04);
    lg2->SetTextColor(global_darkmode ? 0 : 1);
    lg2->AddEntry(rotz_1D, fmt::format("{} areas", entries).c_str(), "");
    lg2->AddEntry(rotz_1D, fmt::format("Mean      {:.2f}", rotz_1D->GetMean()).c_str(), "");
    lg2->AddEntry(rotz_1D, fmt::format("Std Dev   {:.2f}", rotz_1D->GetStdDev()).c_str(), "");
    lg2->Draw();
}

void shift_xy(TCanvas *c1, TTree *tree, const double *AreaParam_d, uint32_t entries) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    TGaxis::SetMaxDigits(3);

    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double bin  = AreaParam_d[0];
    double LowX = AreaParam_d[1];
    double UpX  = AreaParam_d[2];
    double LowY = AreaParam_d[3];
    double UpY  = AreaParam_d[4];
    double view_width = AreaParam_d[7];

    // ブランチアドレスを設定
    int ix, iy;
    double dx, dy;
    tree->SetBranchAddress("ix", &ix);
    tree->SetBranchAddress("iy", &iy);
    tree->SetBranchAddress("dx", &dx);
    tree->SetBranchAddress("dy", &dy);

    TH1D *dx_temp = new TH1D("dx_temp", "", 20000000, -100000.0, 100000.0);
    TH1D *dy_temp = new TH1D("dy_temp", "", 20000000, -100000.0, 100000.0);
    tree->Draw("dx >> dx_temp", "", "goff");
    tree->Draw("dy >> dy_temp", "", "goff");
    double dx_stddev = dx_temp->GetStdDev();
    double dy_stddev = dy_temp->GetStdDev();
    double dx_mean = dx_temp->GetMean();
    double dy_mean = dy_temp->GetMean();
    gDirectory->Delete("dx_temp");
    gDirectory->Delete("dy_temp");

    double dx_low = dx_mean - 5 * dx_stddev;
    double dx_up  = dx_mean + 5 * dx_stddev;
    double dy_low = dy_mean - 5 * dy_stddev;
    double dy_up  = dy_mean + 5 * dy_stddev;

    double mmScale = false;
    if (dx_low < -1000.0 || dx_up > 1000.0 || dy_low < -1000.0 || dy_up > 1000.0) {
        mmScale = true;
        dx_low /= 1000.0;
        dx_up  /= 1000.0;
        dy_low /= 1000.0;
        dy_up  /= 1000.0;
    }

    // 1binが最小単位を下回らないように調整
    const double dxy_unit = mmScale ? 1e-5 : 0.01;

    double dx_range = dx_up - dx_low;
    int dx_bin = (dx_range < dxy_unit * 100) ? std::floor(dx_range / dxy_unit) : 100;
    double dy_range = dy_up - dy_low;
    int dy_bin = (dy_range < dxy_unit * 100) ? std::floor(dy_range / dxy_unit) : 100;

    // 2Dヒストグラム作成
    TH2D *dx_2D = new TH2D("dx_2D",
        mmScale ? ";x [mm];y [mm];x shift [mm]" : ";x [mm];y [mm];x shift [#mum]",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *dy_2D = new TH2D("dy_2D",
        mmScale ? ";x [mm];y [mm];y shift [mm]" : ";x [mm];y [mm];y shift [#mum]",
        bin, LowX, UpX, bin, LowY, UpY);

    for (uint32_t i = 0; i < entries; ++i) {
        tree->GetEntry(i);
        double x = ix * view_width;
        double y = iy * view_width;
        int bx = dx_2D->GetXaxis()->FindBin(x);
        int by = dx_2D->GetYaxis()->FindBin(y);
        if (bx > 0 && bx <= bin && by > 0 && by <= bin) {
            dx_2D->SetBinContent(bx, by, mmScale ? dx / 1000.0 : dx);
            dy_2D->SetBinContent(bx, by, mmScale ? dy / 1000.0 : dy);
        }
    }

    dx_2D->GetZaxis()->SetRangeUser(dx_low, dx_up);
    dy_2D->GetZaxis()->SetRangeUser(dy_low, dy_up);

    // 1Dヒストグラム作成
    TH1D *dx_1D = new TH1D("dx_1D",
        mmScale ? ";x shift [mm];Area" : ";x shift [#mum];Area",
        dx_bin, dx_low, dx_up);
    TH1D *dy_1D = new TH1D("dy_1D",
        mmScale ? ";y shift [mm];Area" : ";y shift [#mum];Area",
        dy_bin, dy_low, dy_up);

    tree->Draw(mmScale ? "dx / 1000.0 >> dx_1D" : "dx >> dx_1D", "", "goff");
    tree->Draw(mmScale ? "dy / 1000.0 >> dy_1D" : "dy >> dy_1D", "", "goff");

    // pad1: dx 2D分布
    c1->cd(1);
    dx_2D->Draw("colz1");
    TLatex *title1 = new TLatex();
    title1->SetTextAlign(22);
    title1->SetTextSize(0.06);
    title1->SetTextColor(global_darkmode ? 0 : 1);
    title1->DrawLatexNDC(0.5, 0.95, "x shift");

    // pad2: dy 2D分布
    c1->cd(2);
    dy_2D->Draw("colz1");
    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    title2->DrawLatexNDC(0.43, 0.95, "y shift");

    // pad3: dx 1D
    c1->cd(3);
    dx_1D->SetFillStyle(0);
    dx_1D->SetLineWidth(2);
    dx_1D->Draw();
    MyUtil::PaintBins(dx_1D, dx_low, dx_up);

    TLatex *title3 = new TLatex();
    title3->SetTextAlign(22);
    title3->SetTextSize(0.06);
    title3->SetTextColor(global_darkmode ? 0 : 1);
    title3->DrawLatexNDC(0.5, 0.95, "x shift");

    TLegend *lg1 = new TLegend(0.73, 0.7, 0.95, 0.9);
    lg1->SetFillStyle(0);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(0.04);
    lg1->SetTextColor(global_darkmode ? 0 : 1);
    lg1->AddEntry(dx_1D, fmt::format("{} areas", entries).c_str(), "");
    lg1->AddEntry(dx_1D, fmt::format("Mean      {:.2f}", dx_1D->GetMean()).c_str(), "");
    lg1->AddEntry(dx_1D, fmt::format("Std Dev   {:.2f}", dx_1D->GetStdDev()).c_str(), "");
    lg1->Draw();

    // pad4: dy 1D
    c1->cd(4);
    dy_1D->SetFillStyle(0);
    dy_1D->SetLineWidth(2);
    dy_1D->Draw();
    MyUtil::PaintBins(dy_1D, dy_low, dy_up);

    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    title4->DrawLatexNDC(0.43, 0.95, "y shift");

    TLegend *lg2 = new TLegend(0.665, 0.7, 0.885, 0.9);
    lg2->SetFillStyle(0);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.04);
    lg2->SetTextColor(global_darkmode ? 0 : 1);
    lg2->AddEntry(dy_1D, fmt::format("{} areas", entries).c_str(), "");
    lg2->AddEntry(dy_1D, fmt::format("Mean      {:.2f}", dy_1D->GetMean()).c_str(), "");
    lg2->AddEntry(dy_1D, fmt::format("Std Dev   {:.2f}", dy_1D->GetStdDev()).c_str(), "");
    lg2->Draw();
}

void shrink(TCanvas *c1, TTree *tree, const double *AreaParam_d, uint32_t entries) noexcept
{
        c1->Divide(2, 2);
        for (int pad = 1; pad <= 4; ++pad) {
            c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
            c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
        }

        TGaxis::SetMaxDigits(3);

        gStyle->SetTitleOffset(1.4, "x");
        gStyle->SetTitleOffset(1.4, "y");
        gStyle->SetTitleOffset(1.2, "z");

        double bin  = AreaParam_d[0];
        double LowX = AreaParam_d[1];
        double UpX  = AreaParam_d[2];
        double LowY = AreaParam_d[3];
        double UpY  = AreaParam_d[4];
        double view_width = AreaParam_d[7];

        // ブランチアドレスを設定
        tree->ResetBranchAddresses();
        int ix, iy;
        double shr_xy, shr_z;
        tree->SetBranchAddress("ix", &ix);
        tree->SetBranchAddress("iy", &iy);
        tree->SetBranchAddress("shr_xy", &shr_xy);
        tree->SetBranchAddress("shr_z", &shr_z);

        TH1D *shr_xy_temp = new TH1D("shr_xy_temp", "", 20000, -100.0, 100.0);
        TH1D *shr_z_temp = new TH1D("shr_z_temp", "", 20000, -100.0, 100.0);
        tree->Draw("shr_xy >> shr_xy_temp", "", "goff");
        tree->Draw("shr_z >> shr_z_temp", "", "goff");
        double shr_xy_stddev = shr_xy_temp->GetStdDev();
        double shr_z_stddev = shr_z_temp->GetStdDev();
        double shr_xy_mean = shr_xy_temp->GetMean();
        double shr_z_mean = shr_z_temp->GetMean();
        gDirectory->Delete("shr_xy_temp");
        gDirectory->Delete("shr_z_temp");

        double shr_xy_low = shr_xy_mean - 5 * shr_xy_stddev;
        double shr_xy_up  = shr_xy_mean + 5 * shr_xy_stddev;
        double shr_z_low = shr_z_mean - 5 * shr_z_stddev;
        double shr_z_up  = shr_z_mean + 5 * shr_z_stddev;

        // 1binが最小単位を下回らないように調整
        constexpr double shr_unit = 0.000001;

        double shr_xy_range = shr_xy_up - shr_xy_low;
        int shr_xy_bin = (shr_xy_range < shr_unit * 100) ? std::floor(shr_xy_range / shr_unit) : 100;
        double shr_z_range = shr_z_up - shr_z_low;
        int shr_z_bin = (shr_z_range < shr_unit * 100) ? std::floor(shr_z_range / shr_unit) : 100;

        // 2Dヒストグラム作成
        TH2D *shr_xy_2D = new TH2D("shr_xy_2D", ";x [mm];y [mm];Shrink xy",
            bin, LowX, UpX, bin, LowY, UpY);
        TH2D *shr_z_2D = new TH2D("shr_z_2D", ";x [mm];y [mm];Shrink z",
            bin, LowX, UpX, bin, LowY, UpY);

        for (uint32_t i = 0; i < entries; ++i) {
            tree->GetEntry(i);
            double x = ix * view_width;
            double y = iy * view_width;
            int bx = shr_xy_2D->GetXaxis()->FindBin(x);
            int by = shr_xy_2D->GetYaxis()->FindBin(y);
            if (bx > 0 && bx <= bin && by > 0 && by <= bin) {
                shr_xy_2D->SetBinContent(bx, by, shr_xy);
                shr_z_2D->SetBinContent(bx, by, shr_z);
            }
        }

        shr_xy_2D->GetZaxis()->SetRangeUser(shr_xy_low, shr_xy_up);
        shr_z_2D->GetZaxis()->SetRangeUser(shr_z_low, shr_z_up);

        // 1Dヒストグラム作成
        TH1D *shr_xy_1D = new TH1D("shr_xy_1D", ";Shrink xy;Area",
            shr_xy_bin, shr_xy_low, shr_xy_up);
        TH1D *shr_z_1D = new TH1D("shr_z_1D", ";Shrink z;Area",
            shr_z_bin, shr_z_low, shr_z_up);

        tree->Draw("shr_xy >> shr_xy_1D", "", "goff");
        tree->Draw("shr_z >> shr_z_1D", "", "goff");

        // pad1: shr_xy 2D分布
        c1->cd(1);
        shr_xy_2D->Draw("colz1");
        TLatex *title1 = new TLatex();
        title1->SetTextAlign(22);
        title1->SetTextSize(0.06);
        title1->SetTextColor(global_darkmode ? 0 : 1);
        title1->DrawLatexNDC(0.5, 0.95, "Shrink xy");

        // pad2: shr_z 2D分布
        c1->cd(2);
        shr_z_2D->Draw("colz1");
        TLatex *title2 = new TLatex();
        title2->SetTextAlign(22);
        title2->SetTextSize(0.06);
        title2->SetTextColor(global_darkmode ? 0 : 1);
        title2->DrawLatexNDC(0.43, 0.95, "Shrink z");

        // pad3: shr_xy 1D
        c1->cd(3);
        shr_xy_1D->SetFillStyle(0);
        shr_xy_1D->SetLineWidth(2);
        shr_xy_1D->Draw();
        MyUtil::PaintBins(shr_xy_1D, shr_xy_low, shr_xy_up);

        TLatex *title3 = new TLatex();
        title3->SetTextAlign(22);
        title3->SetTextSize(0.06);
        title3->SetTextColor(global_darkmode ? 0 : 1);
        title3->DrawLatexNDC(0.5, 0.95, "Shrink xy");

        TLegend *lg1 = new TLegend(0.73, 0.7, 0.95, 0.9);
        lg1->SetFillStyle(0);
        lg1->SetBorderSize(0);
        lg1->SetTextSize(0.04);
        lg1->SetTextColor(global_darkmode ? 0 : 1);
        lg1->AddEntry(shr_xy_1D, fmt::format("{} areas", entries).c_str(), "");
        lg1->AddEntry(shr_xy_1D, fmt::format("Mean      {:.4f}", shr_xy_1D->GetMean()).c_str(), "");
        lg1->AddEntry(shr_xy_1D, fmt::format("Std Dev   {:.4f}", shr_xy_1D->GetStdDev()).c_str(), "");
        lg1->Draw();

        // pad4: shr_z 1D
        c1->cd(4);
        shr_z_1D->SetFillStyle(0);
        shr_z_1D->SetLineWidth(2);
        shr_z_1D->Draw();
        MyUtil::PaintBins(shr_z_1D, shr_z_low, shr_z_up);

        TLatex *title4 = new TLatex();
        title4->SetTextAlign(22);
        title4->SetTextSize(0.06);
        title4->SetTextColor(global_darkmode ? 0 : 1);
        title4->DrawLatexNDC(0.43, 0.95, "Shrink z");

        TLegend *lg2 = new TLegend(0.665, 0.7, 0.885, 0.9);
        lg2->SetFillStyle(0);
        lg2->SetBorderSize(0);
        lg2->SetTextSize(0.04);
        lg2->SetTextColor(global_darkmode ? 0 : 1);
        lg2->AddEntry(shr_z_1D, fmt::format("{} areas", entries).c_str(), "");
        lg2->AddEntry(shr_z_1D, fmt::format("Mean      {:.4f}", shr_z_1D->GetMean()).c_str(), "");
        lg2->AddEntry(shr_z_1D, fmt::format("Std Dev   {:.4f}", shr_z_1D->GetStdDev()).c_str(), "");
        lg2->Draw();
    }
