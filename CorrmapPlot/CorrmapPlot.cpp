#include <iostream>
#include <fstream>
#include <sstream>
#include <csignal>
#include <queue>
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
#include <TGraph2D.h>
#include <TGraph.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TColor.h>
#include <TArrow.h>
#include <TPaveText.h>
#include <TCut.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TGaxis.h>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <argparse/argparse.hpp>

#include <MyHeader/Affine.hpp>
#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>

void sig_map(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept;
void bg_hole(TCanvas *c1, TTree *tree, const double *AreaParam, const uint32_t pl[2],
    uint32_t entries, double area_l[5], double edgecut) noexcept;
void dz_rms(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept;
void shift_xy(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries, const Affine& afp_g,
    double resolution, double arr, double ref, TTree *tree2, TGraph *dxdy_graph[4]) noexcept;
void shift_xys(TCanvas *c1, TTree *tree2, const double *AreaParam, uint32_t entries,
    double resolution, TGraph *dxdy_graph[4]) noexcept;
void rms_xy(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept;
void shift_axay(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries, const Affine& afa_g,
    double resolution, double arr, double ref) noexcept;
void shift_axays(TCanvas *c1, TTree *tree, const double *AreaParam, uint32_t entries) noexcept;
void rms_axay(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept;

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
    parser.add_argument("local_corrmap")
        .help("Path to the local corrmap file (output of ali-l.exe) to be processed.")
        .required();
    parser.add_argument("global_corrmap")
        .help("Path to the global corrmap file (output of ali-g.exe) to be processed.")
        .required();
    // オプション引数: 出力ファイル名
    parser.add_argument("-o", "--output")
        .help("Output PDF file name. [default: local_corrmap.pdf]")
        .default_value(std::string());
    // オプション引数: その他の設定
    parser.add_argument("-res", "--resolution")
        .help("Arrow plot resolution (0.0 - 1.0).")
        .default_value(1.0)
        .scan<'g', double>();
    parser.add_argument("-posref", "--position_reference")
        .help("Reference arrow length of position difference (mm).")
        .default_value(0.1)
        .scan<'g', double>();
    parser.add_argument("-angref", "--angle_reference")
        .help("Reference arrow length of angle difference (mrad).")
        .default_value(10.0)
        .scan<'g', double>();
    parser.add_argument("-arr", "--arrow_length")
        .help("Scale factor for arrow lengths.")
        .default_value(1.0)
        .scan<'g', double>();
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
    argparse::ArgumentParser parser("DcPlot.exe", "1.0.0");
    parse_arguments(parser, argc, argv);

    // 引数を取得
    const auto localfile = parser.get<std::string>("local_corrmap");
    const auto globalfile = parser.get<std::string>("global_corrmap");
    const auto output_arg = parser.get<std::string>("--output");
    const auto resolution = parser.get<double>("--resolution");
    const double pos_reference = parser.get<double>("--position_reference");
    const double ang_reference = parser.get<double>("--angle_reference");
    const auto arr = parser.get<double>("--arrow_length");
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
        ? (localfile + ".pdf")
        : ((output_arg.size() > 4 && output_arg.substr(output_arg.size() - 4) == ".pdf")
            ? output_arg
            : output_arg + ".pdf");

    // 引数を利用する変数を設定
    const int font_code = 10 * font_number + 2;
    const double extz = arr * 250000.0;

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
    std::ifstream ifs_g(globalfile.c_str());
    std::ifstream ifs_l(localfile.c_str());
    std::string line;

    if (!ifs_g) {
        std::cerr << "\n Error: Cannot open file: " << globalfile << std::endl;
        return 1;
    }
    if (!ifs_l) {
        std::cerr << "\n Error: Cannot open file: " << localfile << std::endl;
        return 1;
    }

    uint32_t id, pos0, pos1;
    double x, y, xmin, xmax, ymin, ymax,
        a_pos, b_pos, c_pos, d_pos, p_pos, q_pos,
        a_ang, b_ang, c_ang, d_ang, p_ang, q_ang,
        dz, sig, bg, signif, rms_x, rms_y, rms_ax, rms_ay,
        dmy0, dmy1;

    // globalファイルの読み込み
    Affine afp_g, afa_g;
    int line_count = 0;
    while (std::getline(ifs_g, line)) {
        std::istringstream iss(line);
        if (!(iss >> id >> pos0 >> pos1 >> xmin >> xmax >> ymin >> ymax
                  >> a_pos >> b_pos >> c_pos >> d_pos >> p_pos >> q_pos
                  >> a_ang >> b_ang >> c_ang >> d_ang >> p_ang >> q_ang
                  >> dz >> sig >> bg >> signif >> rms_x >> rms_y
                  >> dmy0 >> dmy1 >> rms_ax >> rms_ay)) {
            std::cerr << "Warning : Malformed line skipped: " << line << std::endl;
            continue;
        }

        // mm単位に変換
        p_pos /= 1000.0;
        q_pos /= 1000.0;

        // mrad単位に変換
        p_ang *= 1000.0;
        q_ang *= 1000.0;

        afp_g = Affine(a_pos, b_pos, c_pos, d_pos, p_pos, q_pos);
        afa_g = Affine(a_ang, b_ang, c_ang, d_ang, p_ang, q_ang);

        line_count++;
    }
    ifs_g.close();
    if (line_count != 1) {
        std::cerr <<
            "\n Attention : Multiple global parameters exist. Last line is selected."
            << std::endl;
    }
    if (line_count == 0) {
        std::cerr << "\n Error : No parameter is found in " << globalfile << std::endl;
        return 1;
    }

    uint32_t g_pos0 = pos0;
    uint32_t g_pos1 = pos1;
    uint32_t pl[2];
    pl[0] = g_pos0 / 10;
    pl[1] = g_pos1 / 10;

    // localファイルの読み込み準備
    double xmin_l = DBL_MAX;
    double xmax_l = -DBL_MAX;
    double ymin_l = DBL_MAX;
    double ymax_l = -DBL_MAX;
    double view_width = DBL_MAX;

    // TTreeを作成
    TTree *tree = new TTree("tree", "");
    const std::array<std::pair<const char*, void*>, 22> double_branches = {{
        {"x", &x}, {"y", &y}, {"a_pos", &a_pos}, {"b_pos", &b_pos}, {"c_pos", &c_pos},
        {"d_pos", &d_pos}, {"p_pos", &p_pos}, {"q_pos", &q_pos}, {"a_ang", &a_ang},
        {"b_ang", &b_ang}, {"c_ang", &c_ang}, {"d_ang", &d_ang}, {"p_ang", &p_ang},
        {"q_ang", &q_ang}, {"dz", &dz}, {"sig", &sig}, {"bg", &bg}, {"signif", &signif},
        {"rms_x", &rms_x}, {"rms_y", &rms_y}, {"rms_ax", &rms_ax}, {"rms_ay", &rms_ay}
    }};
    for (const auto& branch : double_branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    // localファイルの読み込み
    line_count = 0;
    while (std::getline(ifs_l, line)) {
        std::istringstream iss(line);
        if (!(iss >> id >> pos0 >> pos1 >> xmin >> xmax >> ymin >> ymax
                  >> a_pos >> b_pos >> c_pos >> d_pos >> p_pos >> q_pos
                  >> a_ang >> b_ang >> c_ang >> d_ang >> p_ang >> q_ang
                  >> dz >> sig >> bg >> signif >> rms_x >> rms_y
                  >> dmy0 >> dmy1 >> rms_ax >> rms_ay)) {
            std::cerr << "Warning : Malformed line skipped: " << line << std::endl;
            continue;
        }

        if (pos0 != g_pos0 || pos1 != g_pos1) {
            std::cerr << "Warning : Position IDs do not match global parameters. Line skipped: "
                      << line << std::endl;
            continue;
        }

        if (xmin < xmin_l) xmin_l = xmin;
        if (xmax > xmax_l) xmax_l = xmax;
        if (ymin < ymin_l) ymin_l = ymin;
        if (ymax > ymax_l) ymax_l = ymax;

        double width = xmax - xmin;
        if (width < view_width) view_width = width;

        // mm単位に変換
        p_pos /= 1000.0;
        q_pos /= 1000.0;
        x = (xmax + xmin) / 2000.0;
        y = (ymax + ymin) / 2000.0;

        // mrad単位に変換
        p_ang *= 1000.0;
        q_ang *= 1000.0;

        tree->Fill();
        line_count++;
    }
    ifs_l.close();

    double area_l[5] = {xmin_l, xmax_l, ymin_l, ymax_l, view_width};

    uint32_t lentries = line_count;
    if (lentries == 0) {
        std::cerr << "\n Error : No parameter is found in " << localfile << std::endl;
        return 1;
    }

    // 情報表示
    std::cout << "\n Local File  : " << localfile << std::endl;
    std::cout << " Global File : " << globalfile << std::endl;
    std::cout << " Arrow       : " << fmt::format("{:.1f}", arr) << std::endl;
    std::cout << " Resolution  : " << fmt::format("{:.1f}", resolution) << std::endl;
    std::cout << " # of Areas  : " << lentries << std::endl;

    // プロット開始
    TDatime time_now;
    std::string Time_Now = fmt::format(
        "{:d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}",
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(),
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    sw.Start();
	std::cout << " Plot start  : " << Time_Now << std::endl;

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas *c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());
    global_c1 = c1;
    global_output = output;

	// プログレスバーの初期化
	int page = 0;
    constexpr int total = 9; // 合計ページ数
	MyUtil::ShowProgress(page, total);

    // データの座標の範囲を取得し、表示範囲とビンの数を決定する
    // フィルムの長辺の端から1cm外側までを最大表示範囲とし、縦横比を正しく保って表示する
    const double MinX = tree->GetMinimum("x");
    const double MaxX = tree->GetMaximum("x");
    const double MinY = tree->GetMinimum("y");
    const double MaxY = tree->GetMaximum("y");
    const double RangeX = MaxX - MinX; // データ領域
    const double RangeY = MaxY - MinY; // データ領域
    double LowX, UpX, LowY, UpY, bin; // 表示範囲とビンの数
    if (RangeX >= RangeY) {
        LowX = MinX - 10;
        UpX = MaxX + 10;
        LowY = MinY - (RangeX - RangeY + 20) * 0.5;
        UpY = MaxY + (RangeX - RangeY + 20) * 0.5;
        bin = (RangeX + 20) * 1000 / view_width;
    } else {
        LowX = MinX - (RangeY - RangeX + 20) * 0.5;
        UpX = MaxX + (RangeY - RangeX + 20) * 0.5;
        LowY = MinY - 10;
        UpY = MaxY + 10;
        bin = (RangeY + 20) * 1000 / view_width;
    }
    const double AreaParam[7] = {bin, LowX, UpX, LowY, UpY, RangeX, RangeY};

    sig_map(c1, tree, AreaParam, pl, lentries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");

    bg_hole(c1, tree, AreaParam, pl, lentries, area_l, edgecut);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("sig*");
    gDirectory->Delete("bg*");
    gDirectory->Delete("hole*");

    dz_rms(c1, tree, AreaParam, pl, lentries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("dz*");

    TTree *tree2 = new TTree("tree2", "");
    double dx, dy;
    tree2->Branch("dx", &dx, "dx/D");
    tree2->Branch("dy", &dy, "dy/D");

    TGraph *dxdy_graph[4] = {
        new TGraph(),
        new TGraph(),
        new TGraph(),
        new TGraph()
    };

    shift_xy(c1, tree, AreaParam, pl, lentries, afp_g, resolution, arr,
        pos_reference, tree2, dxdy_graph);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("shift*");

    shift_xys(c1, tree2, AreaParam, lentries, resolution, dxdy_graph);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    delete tree2;
    for (int i = 0; i < 4; ++i) {
        delete dxdy_graph[i];
    }
    gDirectory->Delete("shift*");

    rms_xy(c1, tree, AreaParam, pl, lentries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("rms*");

    shift_axay(c1, tree, AreaParam, pl, lentries, afa_g, resolution, arr, ang_reference);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("shift*");

    shift_axays(c1, tree, AreaParam, lentries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("shift*");

    rms_axay(c1, tree, AreaParam, pl, lentries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("rms*");

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
        "\n Plot end    : {} - Elapsed {:.2f} [s] (CPU: {:.2f} [s])",
        Time_Now, elapsed_time, cpu_time
    ) << std::endl;
    std::cout << " Output      : " << output << std::endl;

    delete c1;
    gDirectory->Delete("tree");

    return 0;
}

void sig_map(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept
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

    double x, y, sig, signif;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("sig", &sig);
    tree->SetBranchAddress("signif", &signif);

    // 2Dヒストグラム作成
    TH2D *sig_hist[2];
    for (int i = 0; i < 2; ++i) {
        std::string title_text = (i == 0)
            ? "Signal"
            : "Significance";
        sig_hist[i] = new TH2D(
            fmt::format("sig_hist{}", i+1).c_str(),
            fmt::format(";x [mm];y [mm];{}", title_text).c_str(),
            bin, LowX, UpX,
            bin, LowY, UpY
        );

        // 各ビンの合計値とカウントを保存して平均化
        std::vector<double> sum(bin * bin, 0.0);
        std::vector<int> count(bin * bin, 0);

        for (int j = 0; j < entries; ++j) {
            tree->GetEntry(j);
            int binX = sig_hist[i]->GetXaxis()->FindBin(x);
            int binY = sig_hist[i]->GetYaxis()->FindBin(y);
            if (binX >= 1 && binX <= bin && binY >= 1 && binY <= bin) {
                int idx = (binY - 1) * bin + (binX - 1);
                sum[idx] += (i == 0) ? sig : signif;
                count[idx] += 1;
            }
        }

        // 平均値をヒストグラムにセット
        for (int bx = 1; bx <= bin; ++bx) {
            for (int by = 1; by <= bin; ++by) {
                int idx = (by - 1) * bin + (bx - 1);
                if (count[idx] > 0) {
                    sig_hist[i]->SetBinContent(bx, by, sum[idx] / count[idx]);
                } else {
                    sig_hist[i]->SetBinContent(bx, by, 0);
                }
            }
        }

        // 描画
        c1->cd(i + 1);

        sig_hist[i]->Draw("colz");

        // タイトル
        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatexNDC(
            (i == 0) ? 0.5 : 0.43,
            0.95,
            fmt::format("{} PL{:03d}-PL{:03d}", title_text, pl[0], pl[1]).c_str());
    }

    // 1Dヒストグラム作成
    const double sig_max = tree->GetMaximum("sig");
    const double signif_max = tree->GetMaximum("signif");
    TH1D *sig_1D = new TH1D("sig_1D", ";Signal;Area", 100, 0.0, sig_max);
    TH1D *signif_1D = new TH1D("signif_1D", ";Significance;Area", 100, 0.0, signif_max);

    for (int i = 0; i < 2; ++i) {
        std::string draw_variable = (i == 0) ? "sig" : "signif";
        double max_value = (i == 0) ? sig_max : signif_max;
        TH1D *hist = (i == 0) ? sig_1D : signif_1D;
        std::string title_1D = (i == 0) ? "Signal" : "Significance";
        tree->Draw(fmt::format("{} >> {}", draw_variable, hist->GetName()).c_str(), "", "goff");
        c1->cd(i + 3);
        hist->SetFillStyle(0);
        hist->SetLineWidth(2);
        hist->Draw();
        MyUtil::PaintBins(hist, 0.0, max_value); // 各ビンをカラーパレットの色で塗る

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatexNDC((i == 0) ? 0.5 : 0.43, 0.95, title_1D.c_str());

        TLegend *lg = new TLegend(
            (i == 0) ? 0.73 : 0.665, 0.7,
            (i == 0) ? 0.95 : 0.885, 0.9
        );
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.04);
        lg->SetTextColor(global_darkmode ? 0 : 1);
        lg->AddEntry(hist, fmt::format("{} areas", entries).c_str(), "");
        lg->AddEntry(hist, fmt::format("Mean      {:.2f}", hist->GetMean()).c_str(), "");
        lg->AddEntry(hist, fmt::format("Std Dev   {:.2f}", hist->GetStdDev()).c_str(), "");
        lg->Draw();
    }
}

void bg_hole(TCanvas *c1, TTree *tree, const double *AreaParam, const uint32_t pl[2],
    uint32_t entries, double area_l[5], double edgecut) noexcept
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

    double xmin_l = area_l[0];
    double xmax_l = area_l[1];
    double ymin_l = area_l[2];
    double ymax_l = area_l[3];
    double view_width = area_l[4];

    double x, y, bg, sig;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("bg", &bg);
    tree->SetBranchAddress("sig", &sig);

    // bg の2Dヒストグラム作成
    TH2D *bg_hist = new TH2D(
        "bg_hist", ";x [mm];y [mm];Background",
        bin, LowX, UpX,
        bin, LowY, UpY
    );

    // 各ビンの合計値とカウントを保存して平均化
    std::vector<double> sum(bin * bin, 0.0);
    std::vector<int> count(bin * bin, 0);

    for (int j = 0; j < entries; ++j) {
        tree->GetEntry(j);
        int binX = bg_hist->GetXaxis()->FindBin(x);
        int binY = bg_hist->GetYaxis()->FindBin(y);
        if (binX >= 1 && binX <= bin && binY >= 1 && binY <= bin) {
            int idx = (binY - 1) * bin + (binX - 1);
            sum[idx] += bg;
            count[idx] += 1;
        }
    }

    // 平均値をヒストグラムにセット
    for (int bx = 1; bx <= bin; ++bx) {
        for (int by = 1; by <= bin; ++by) {
            int idx = (by - 1) * bin + (bx - 1);
            if (count[idx] > 0) {
                bg_hist->SetBinContent(bx, by, sum[idx] / count[idx]);
            } else {
                bg_hist->SetBinContent(bx, by, 0);
            }
        }
    }

    // 描画
    c1->cd(2);
    bg_hist->Draw("colz1");

    TLatex *title = new TLatex();
    title->SetTextAlign(22);
    title->SetTextSize(0.06);
    title->SetTextColor(global_darkmode ? 0 : 1);
    title->DrawLatexNDC(0.43, 0.95, fmt::format("Background PL{:03d}-PL{:03d}", pl[0], pl[1]).c_str());

    // bg の1Dヒストグラム作成
    const double bg_max = tree->GetMaximum("bg");
    const double bg_min = tree->GetMinimum("bg");
    TH1D *bg_1D = new TH1D("bg_1D", ";Background;Area", 20, bg_min, bg_max);
    tree->Draw(fmt::format("bg >> bg_1D").c_str(), "", "goff");

    c1->cd(4);
    bg_1D->SetFillStyle(0);
    bg_1D->SetLineWidth(2);
    bg_1D->Draw();
    MyUtil::PaintBins(bg_1D, bg_min, bg_max); // 各ビンをカラーパレットの色で塗る

    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    title2->DrawLatexNDC(0.43, 0.95, "Background");

    TLegend *lg = new TLegend(0.665, 0.7, 0.885, 0.9);
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.04);
    lg->SetTextColor(global_darkmode ? 0 : 1);
    lg->AddEntry(bg_1D, fmt::format("{} areas", entries).c_str(), "");
    lg->AddEntry(bg_1D, fmt::format("Mean      {:.2f}", bg_1D->GetMean()).c_str(), "");
    lg->AddEntry(bg_1D, fmt::format("Std Dev   {:.2f}", bg_1D->GetStdDev()).c_str(), "");
    lg->Draw();

    // gDirectoryからsig_hist1を取得
    TH2D *sig_hist1 = (TH2D*)gDirectory->Get("sig_hist1");
    c1->cd(1);
    sig_hist1->Draw("col");

    TLatex *title3 = new TLatex();
    title3->SetTextAlign(22);
    title3->SetTextSize(0.06);
    title3->SetTextColor(global_darkmode ? 0 : 1);
    title3->DrawLatexNDC(0.5, 0.95, "Data Hole Map");

    // hole の2Dヒストグラム作成
    TH2D *hole_hist = new TH2D(
        "hole_hist", ";x [mm];y [mm]",
        bin, LowX, UpX,
        bin, LowY, UpY
    );

    std::vector<int> count_data(bin * bin, 0);

    for (int i = 0; i < entries; ++i) {
        tree->GetEntry(i);
        int binX = hole_hist->GetXaxis()->FindBin(x);
        int binY = hole_hist->GetYaxis()->FindBin(y);
        if (binX >= 1 && binX <= bin && binY >= 1 && binY <= bin) {
            int idx = (binY - 1) * bin + (binX - 1);
            count_data[idx] += 1;
        }
    }

    // 平均値をヒストグラムにセット
    for (int bx = 1; bx <= bin; ++bx) {
        for (int by = 1; by <= bin; ++by) {
            int idx = (by - 1) * bin + (bx - 1);
            if (count_data[idx] > 0) {
                hole_hist->SetBinContent(bx, by, 1);
            } else {
                hole_hist->SetBinContent(bx, by, 0);
            }
        }
    }

    c1->cd(3);
    hole_hist->SetFillColor(93);
    hole_hist->Draw("box");

    // edge cut 領域の描画
    double edge_LowX = xmin_l * 0.001 + edgecut;
    double edge_UpX = xmax_l * 0.001 - edgecut;
    double edge_LowY = ymin_l * 0.001 + edgecut;
    double edge_UpY = ymax_l * 0.001 - edgecut;

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

    for (int i = 0; i < 2; ++i) {
        c1->cd(i * 2 + 1);
        for (int j = 0; j < 4; ++j) {
            edge_areas[j]->Draw("same");
        }
    }

    // 穴検出&クラスタリング
    // TH2Dから穴（count_data == 0）の座標を抽出
    struct HoleCell {
        double x, y;
        int binX, binY;
    };
    std::vector<HoleCell> hole_cells;

    for (int bx = 1; bx <= bin; ++bx) {
        for (int by = 1; by <= bin; ++by) {
            int idx = (by - 1) * bin + (bx - 1);
            if (count_data[idx] == 0) {
                // ビンの中心座標を取得
                double xc = hole_hist->GetXaxis()->GetBinCenter(bx);
                double yc = hole_hist->GetYaxis()->GetBinCenter(by);
                hole_cells.push_back({xc, yc, bx, by});
            }
        }
    }

    int holecells = hole_cells.size();
    if (holecells == 0) return;

    std::vector<bool> visited(holecells, false);
    std::vector<std::vector<int>> clusters;

    // Breadth-First Search (BFS)を用いて各holecellをクラスタリング
    // 2つの点のx, y方向の距離が view_width * 2.0 より近ければクラスタリングする
    // edgecutで指定した範囲外の点はクラスタリングしない

    view_width /= 1000.0; // mm単位に変換

    for (int i = 0; i < holecells; ++i) {
        if (visited[i]) continue;
        visited[i] = true;

        double x0 = hole_cells[i].x;
        double y0 = hole_cells[i].y;
        if (x0 < edge_LowX || x0 > edge_UpX || y0 < edge_LowY || y0 > edge_UpY) continue;

        std::vector<int> cluster;
        std::queue<int> q;
        q.push(i);

        while (!q.empty()) {
            int j = q.front();
            q.pop();
            cluster.push_back(j);

            for (int k = 0; k < holecells; ++k) {
                if (visited[k]) continue;

                double x1 = hole_cells[j].x;
                double y1 = hole_cells[j].y;
                double x2 = hole_cells[k].x;
                double y2 = hole_cells[k].y;

                if (x1 < edge_LowX || x1 > edge_UpX || y1 < edge_LowY || y1 > edge_UpY ||
                    x2 < edge_LowX || x2 > edge_UpX || y2 < edge_LowY || y2 > edge_UpY) {
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
        for (size_t j = 0; j < cluster.size(); ++j) {
            x_vals.push_back(hole_cells[cluster[j]].x);
            y_vals.push_back(hole_cells[cluster[j]].y);
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
            for (size_t j = 0; j < cluster.size(); ++j) {
                double x = hole_cells[cluster[j]].x;
                double y = hole_cells[cluster[j]].y;
                double d = std::hypot(x - centerX, y - centerY) + 0.5 * view_width + 5.0;
                if (d > dist_max) dist_max = d;
            }
        }

        // 楕円を描画
        TEllipse *ellipse = new TEllipse(centerX, centerY, dist_max, dist_max);
        ellipse->SetLineColor(global_darkmode ? 0 : 1);
        ellipse->SetLineWidth(1);
        ellipse->SetFillStyle(0);

        c1->cd(1);
        ellipse->Draw("same");
        c1->cd(3);
        ellipse->Draw("same");
    }

    TLegend *lg2 = new TLegend(0.76, 0.7, 0.93, 0.98);
    lg2->SetFillStyle(0);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.04);
    lg2->SetTextColor(global_darkmode ? 0 : 1);
    lg2->AddEntry(bg_1D, fmt::format("Edge cut: {0:.1f} mm", edgecut).c_str(), "");
    c1->cd(1);
    lg2->Draw();
    c1->cd(3);
    lg2->Draw();

    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    std::string title4_text = fmt::format("Number of Holes: {}", clusters.size());
    title4->DrawLatexNDC(0.5, 0.95, title4_text.c_str());
}

void dz_rms(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    int bin  = static_cast<int>(AreaParam[0]);
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    double x, y, dz, rms_x, rms_y;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("dz", &dz);
    tree->SetBranchAddress("rms_x", &rms_x);
    tree->SetBranchAddress("rms_y", &rms_y);

    // 3つの2Dヒストグラムを作成
    TH2D *dz_2D = new TH2D("dz_2D", ";x [mm];y [mm];dz",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *rmsx_2D = new TH2D("rmsx_2D", ";x [mm];y [mm];RMS x",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *rmsy_2D = new TH2D("rmsy_2D", ";x [mm];y [mm];RMS y",
        bin, LowX, UpX, bin, LowY, UpY);

    // 各ビンの合計値とカウントを保存して平均化
    std::vector<double> sum_dz(bin * bin, 0.0);
    std::vector<double> sum_rmsx(bin * bin, 0.0);
    std::vector<double> sum_rmsy(bin * bin, 0.0);
    std::vector<int> count(bin * bin, 0);

    for (int j = 0; j < entries; ++j) {
        tree->GetEntry(j);
        int binX = dz_2D->GetXaxis()->FindBin(x);
        int binY = dz_2D->GetYaxis()->FindBin(y);
        if (binX >= 1 && binX <= bin && binY >= 1 && binY <= bin) {
            int idx = (binY - 1) * bin + (binX - 1);
            sum_dz[idx] += dz;
            sum_rmsx[idx] += rms_x;
            sum_rmsy[idx] += rms_y;
            count[idx] += 1;
        }
    }

    // 平均値をヒストグラムにセット
    for (int bx = 1; bx <= bin; ++bx) {
        for (int by = 1; by <= bin; ++by) {
            int idx = (by - 1) * bin + (bx - 1);
            if (count[idx] > 0) {
                dz_2D->SetBinContent(bx, by, sum_dz[idx] / count[idx]);
                if (sum_rmsx[idx] < 0 || sum_rmsy[idx] < 0) continue;
                rmsx_2D->SetBinContent(bx, by, sum_rmsx[idx] / count[idx]);
                rmsy_2D->SetBinContent(bx, by, sum_rmsy[idx] / count[idx]);
            }
        }
    }

    // dz の1Dヒストグラム作成
    TH1D *dz_temp = new TH1D("dz_temp", "", 200000, -100000, 100000);
    tree->Draw("dz >> dz_temp", "", "goff");
    double dz_mean = dz_temp->GetMean();
    double dz_3sigma = 3 * dz_temp->GetStdDev();
    double range_low = dz_mean - dz_3sigma;
    double range_up  = dz_mean + dz_3sigma;
    gDirectory->Delete("dz_temp");

    TH1D *dz_1D = new TH1D("dz_1D", ";dz;Area", 50, range_low, range_up);
    tree->Draw("dz >> dz_1D", "", "goff");
    dz_2D->GetZaxis()->SetRangeUser(range_low, range_up);

    // pad1: dz 2D分布
    c1->cd(1);
    dz_2D->Draw("colz1");
    TLatex *title1 = new TLatex();
    title1->SetTextAlign(22);
    title1->SetTextSize(0.06);
    title1->SetTextColor(global_darkmode ? 0 : 1);
    std::string title1_text = fmt::format("dz PL{:03d}-PL{:03d}", pl[0], pl[1]);
    title1->DrawLatexNDC(0.5, 0.95, title1_text.c_str());

    // pad2: dz 1D分布
    c1->cd(2);
    dz_1D->SetFillStyle(0);
    dz_1D->SetLineWidth(2);
    dz_1D->Draw();
    MyUtil::PaintBins(dz_1D, range_low, range_up);

    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    title2->DrawLatexNDC(0.43, 0.95, "dz");

    TLegend *lg = new TLegend(0.665, 0.7, 0.885, 0.9);
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.04);
    lg->SetTextColor(global_darkmode ? 0 : 1);

    std::string areas_text = fmt::format("{} areas", entries);
    std::string mean_text = fmt::format("Mean      {:.2f}", dz_1D->GetMean());
    std::string std_text = fmt::format("Std Dev   {:.2f}", dz_1D->GetStdDev());

    lg->AddEntry(dz_1D, areas_text.c_str(), "");
    lg->AddEntry(dz_1D, mean_text.c_str(), "");
    lg->AddEntry(dz_1D, std_text.c_str(), "");
    lg->Draw();

    double rms_max = std::max(tree->GetMaximum("rms_x"), tree->GetMaximum("rms_y")) + 0.2;
    rmsx_2D->GetZaxis()->SetRangeUser(0.0, rms_max);
    rmsy_2D->GetZaxis()->SetRangeUser(0.0, rms_max);

    // pad3: rms_x 2D分布
    c1->cd(3);
    rmsx_2D->Draw("colz");
    TLatex *title3 = new TLatex();
    title3->SetTextAlign(22);
    title3->SetTextSize(0.06);
    title3->SetTextColor(global_darkmode ? 0 : 1);
    std::string title3_text = fmt::format("RMS x PL{:03d}-PL{:03d}", pl[0], pl[1]);
    title3->DrawLatexNDC(0.5, 0.95, title3_text.c_str());

    // pad4: rms_y 2D分布
    c1->cd(4);
    rmsy_2D->Draw("colz");
    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    std::string title4_text = fmt::format("RMS y PL{:03d}-PL{:03d}", pl[0], pl[1]);
    title4->DrawLatexNDC(0.43, 0.95, title4_text.c_str());
}

void shift_xy(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries, const Affine& afp_g,
    double resolution, double arr, double ref, TTree *tree2, TGraph *dxdy_graph[4]) noexcept
{
    double LowX = AreaParam[1];
    double UpX = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY = AreaParam[4];
    double RangeX = AreaParam[5];
    double RangeY = AreaParam[6];

    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(1.2, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    // resolutionの値が0や>1でも必ず1以上のステップ数になるように保護
    double res = (std::isfinite(resolution) && resolution > 0.0) ? resolution : 1.0;
    if (res > 1.0) res = 1.0;
    uint32_t loop_step = static_cast<uint32_t>(std::max(1.0, static_cast<double>(std::lround(1.0 / res))));
    uint32_t loop_count = (entries + loop_step - 1) / loop_step;

    if (!afp_g.isInvertible()) {
        std::cerr << "Warning: global affine is non-invertible; shift_xy skipped." << std::endl;
        return;
    }
    const Affine afp_g_inv = afp_g.inverse();

    std::string shiftxy_title;
    if (resolution == 1.0) {
        shiftxy_title = fmt::format("Position Shift PL{:03d}-PL{:03d};x [mm];y [mm]", pl[0], pl[1]);
    } else {
        shiftxy_title = fmt::format("Position Shift PL{:03d}-PL{:03d} ({} / {} areas);x [mm];y [mm]",
                                    pl[0], pl[1], loop_count, entries);
    }
    TH2D *shiftxy_frame = new TH2D("shiftxy_frame", shiftxy_title.c_str(),
                                   100, LowX, UpX, 100, LowY, UpY);
    shiftxy_frame->Draw();

    // extrapolated z
    double extz = 1000. * 0.3 * arr;

    // 既存のブランチアドレスをリセットし、正しいポインタに付け替える
    tree->ResetBranchAddresses();

    double x, y, a_pos, b_pos, c_pos, d_pos, p_pos, q_pos;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("a_pos", &a_pos);
    tree->SetBranchAddress("b_pos", &b_pos);
    tree->SetBranchAddress("c_pos", &c_pos);
    tree->SetBranchAddress("d_pos", &d_pos);
    tree->SetBranchAddress("p_pos", &p_pos);
    tree->SetBranchAddress("q_pos", &q_pos);

    double dx = 0.0, dy = 0.0;
    tree2->SetBranchAddress("dx", &dx);
    tree2->SetBranchAddress("dy", &dy);

    TGraph *gdx_x = dxdy_graph[0];
    TGraph *gdx_y = dxdy_graph[1];
    TGraph *gdy_x = dxdy_graph[2];
    TGraph *gdy_y = dxdy_graph[3];

    for (int i = 0; i < entries; i += loop_step) {
        tree->GetEntry(i);

        // calculate displacement due to local correction
        Affine afp(a_pos, b_pos, c_pos, d_pos, p_pos, q_pos);
        Affine afp_l = afp * afp_g_inv;
        afp_l.transform(x, y, &dx, &dy);
        dx -= x;
        dy -= y;
        dx *= 1000.0; // mm -> um
        dy *= 1000.0; // mm -> um

        tree2->Fill();

        // SetPoint(i, xc, dx) だと (0, 0) に点が打たれてしまうため GetN() を使っている
        gdx_x->SetPoint(gdx_x->GetN(), x, dx);
        gdx_y->SetPoint(gdx_y->GetN(), y, dx);
        gdy_x->SetPoint(gdy_x->GetN(), x, dy);
        gdy_y->SetPoint(gdy_y->GetN(), y, dy);

        TArrow *arrow = new TArrow(x, y, x + dx * 0.001 * extz, y + dy * 0.001 * extz, 0.01, ">");
        arrow->SetFillColor(93);
        arrow->SetLineColor(93);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(1);
        arrow->SetLineStyle(1);
        arrow->Draw();
    }

    // reference
    double ref_length = ref * extz;
    std::string ref_text = fmt::format("{} mm", ref);
    TPaveText *ref_pt = new TPaveText(
        UpX + RangeX * 0.08, LowY + RangeY * 0.07,
        UpX + RangeX * 0.28, LowY + RangeY * 0.17,
        "br"
    );
    ref_pt->SetTextColor(global_darkmode ? 0 : 1);
    ref_pt->SetFillColor(global_darkmode ? 1 : 0);
    ref_pt->SetLineColor(global_darkmode ? 1 : 0);
    ref_pt->SetShadowColor(global_darkmode ? 1 : 0);
    TText *text_xy = ref_pt->AddText(ref_text.c_str());
    ref_pt->Draw();

    TArrow *arrowX = new TArrow(
        UpX + RangeX * 0.05, LowY + RangeY * 0.05,
        UpX + RangeX * 0.05 + ref_length, LowY + RangeY * 0.05, 0.01,
        ">"
    );
    arrowX->SetFillColor(global_darkmode ? 0 : 1);
    arrowX->SetLineColor(global_darkmode ? 0 : 1);
    arrowX->SetFillStyle(1001);
    arrowX->SetLineWidth(2);
    arrowX->SetLineStyle(1);
    arrowX->Draw();

    TArrow *arrowY = new TArrow(
        UpX + RangeX * 0.05, LowY + RangeY * 0.05,
        UpX + RangeX * 0.05, LowY + RangeY * 0.05 + ref_length, 0.01,
        ">"
    );
    arrowY->SetFillColor(global_darkmode ? 0 : 1);
    arrowY->SetLineColor(global_darkmode ? 0 : 1);
    arrowY->SetFillStyle(1001);
    arrowY->SetLineWidth(2);
    arrowY->SetLineStyle(1);
    arrowY->Draw();
}

void shift_xys(TCanvas *c1, TTree *tree2, const double *AreaParam, uint32_t entries,
    double resolution, TGraph *dxdy_graph[4]) noexcept
{
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.5, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    c1->Clear();
    c1->Divide(3, 2);

    // 既存のブランチアドレスをリセットし、正しいポインタに付け替える
    tree2->ResetBranchAddresses();

    double dx, dy;
    tree2->SetBranchAddress("dx", &dx);
    tree2->SetBranchAddress("dy", &dy);

    // ループ数表示のための計算
    double res = (std::isfinite(resolution) && resolution > 0.0) ? resolution : 1.0;
    if (res > 1.0) res = 1.0;
    uint32_t loop_step = static_cast<uint32_t>(std::max(1.0, static_cast<double>(std::lround(1.0 / res))));
    uint32_t loop_count = (entries + loop_step - 1) / loop_step;

    // マーカー設定
    constexpr int markerstyle = 20;
    constexpr double markersize = 0.2;

    TH1D *shift_x_temp = new TH1D("shift_x_temp", "", 40000, -1000.0, 1000.0);
    tree2->Draw("dx>>shift_x_temp", "", "goff");
    double sx_mean = shift_x_temp->GetMean();
    double sx_stddev = shift_x_temp->GetStdDev();
    double sx_range_low = sx_mean - 5 * sx_stddev;
    double sx_range_up  = sx_mean + 5 * sx_stddev;
    gDirectory->Delete("shift_x_temp");
    TH1D *shift_y_temp = new TH1D("shift_y_temp", "", 40000, -1000.0, 1000.0);
    tree2->Draw("dy>>shift_y_temp", "", "goff");
    double sy_mean = shift_y_temp->GetMean();
    double sy_stddev = shift_y_temp->GetStdDev();
    double sy_range_low = sy_mean - 5 * sy_stddev;
    double sy_range_up  = sy_mean + 5 * sy_stddev;
    gDirectory->Delete("shift_y_temp");
    double range_low = std::min(sx_range_low, sy_range_low);
    double range_up  = std::max(sx_range_up, sy_range_up);

    // pad1: x shift ヒスト
    c1->cd(1);
    std::string shiftx_title;
    if (resolution == 1.0) {
        shiftx_title = "x shift;x shift [mm];";
    } else {
        shiftx_title = fmt::format("x shift ({} / {} areas);x shift [mm];", loop_count, entries);
    }

    TH1D *shift_x = new TH1D("shift_x", shiftx_title.c_str(), 100, range_low, range_up);
    tree2->Draw("dx>>shift_x", "", "goff");
    shift_x->SetFillColor(90);
    shift_x->Draw();

    TLegend *sx_lg = new TLegend(0.55, 0.75, 0.75, 0.85);
    sx_lg->SetFillStyle(0);
    sx_lg->SetBorderSize(0);
    sx_lg->SetTextSize(0.04);
    sx_lg->AddEntry(shift_x, Form("Mean   %.2f", sx_mean), "");
    sx_lg->AddEntry(shift_x, Form("Std Dev   %.2f", sx_stddev), "");
    sx_lg->Draw();

    // 元配列
    TGraph *gdx_x = dxdy_graph[0];
    TGraph *gdx_y = dxdy_graph[1];
    TGraph *gdy_x = dxdy_graph[2];
    TGraph *gdy_y = dxdy_graph[3];

    // pad2: x shift : x 散布図（mm）
    c1->cd(2);
    std::string dx_x_title;
    if (resolution == 1.0) {
        dx_x_title = "x shift : x;x [mm];x shift [mm]";
    } else {
        dx_x_title = fmt::format("x shift : x ({} / {} areas);x [mm];x shift [mm]", loop_count, entries);
    }
    gdx_x->SetTitle(dx_x_title.c_str());
    gdx_x->SetMarkerStyle(markerstyle);
    gdx_x->SetMarkerSize(markersize);
    gdx_x->SetMarkerColor(90);
    gdx_x->GetXaxis()->SetLimits(LowX, UpX);
    gdx_x->GetYaxis()->SetRangeUser(range_low, range_up);
    gdx_x->Draw("AP");

    // pad3: x shift : y 散布図（mm）
    c1->cd(3);
    std::string dx_y_title;
    if (resolution == 1.0) {
        dx_y_title = "x shift : y;y [mm];x shift [mm]";
    } else {
        dx_y_title = fmt::format("x shift : y ({} / {} areas);y [mm];x shift [mm]", loop_count, entries);
    }
    gdx_y->SetTitle(dx_y_title.c_str());
    gdx_y->SetMarkerStyle(markerstyle);
    gdx_y->SetMarkerSize(markersize);
    gdx_y->SetMarkerColor(90);
    gdx_y->GetXaxis()->SetLimits(LowY, UpY);
    gdx_y->GetYaxis()->SetRangeUser(range_low, range_up);
    gdx_y->Draw("AP");

    // pad4: y shift ヒスト
    c1->cd(4);
    std::string shifty_title;
    if (resolution == 1.0) {
        shifty_title = "y shift;y shift [mm];";
    } else {
        shifty_title = fmt::format("y shift ({} / {} areas);y shift [mm];", loop_count, entries);
    }

    TH1D *shift_y = new TH1D("shift_y", shifty_title.c_str(), 100, range_low, range_up);
    tree2->Draw("dy>>shift_y", "", "goff");
    shift_y->SetFillColor(91);
    shift_y->Draw();

    TLegend *sy_lg = new TLegend(0.55, 0.75, 0.75, 0.85);
    sy_lg->SetFillStyle(0);
    sy_lg->SetBorderSize(0);
    sy_lg->SetTextSize(0.04);
    sy_lg->AddEntry(shift_y, Form("Mean   %.2f", sy_mean), "");
    sy_lg->AddEntry(shift_y, Form("Std Dev   %.2f", sy_stddev), "");
    sy_lg->Draw();

    // pad5: y shift : x 散布図（mm）
    c1->cd(5);
    std::string dy_x_title;
    if (resolution == 1.0) {
        dy_x_title = "y shift : x;x [mm];y shift [mm]";
    } else {
        dy_x_title = fmt::format("y shift : x ({} / {} areas);x [mm];y shift [mm]", loop_count, entries);
    }
    gdy_x->SetTitle(dy_x_title.c_str());
    gdy_x->SetMarkerStyle(markerstyle);
    gdy_x->SetMarkerSize(markersize);
    gdy_x->SetMarkerColor(91);
    gdy_x->GetXaxis()->SetLimits(LowX, UpX);
    gdy_x->GetYaxis()->SetRangeUser(range_low, range_up);
    gdy_x->Draw("AP");

    // pad6: y shift : y 散布図（mm）
    c1->cd(6);
    std::string dy_y_title;
    if (resolution == 1.0) {
        dy_y_title = "y shift : y;y [mm];y shift [mm]";
    } else {
        dy_y_title = fmt::format("y shift : y ({} / {} areas);y [mm];y shift [mm]", loop_count, entries);
    }
    gdy_y->SetTitle(dy_y_title.c_str());
    gdy_y->SetMarkerStyle(markerstyle);
    gdy_y->SetMarkerSize(markersize);
    gdy_y->SetMarkerColor(91);
    gdy_y->GetXaxis()->SetLimits(LowY, UpY);
    gdy_y->GetYaxis()->SetRangeUser(range_low, range_up);
    gdy_y->Draw("AP");
}

void rms_xy(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    // dz_rmsで作成された2DヒストグラムをgDirectoryから取得
    TH2D *rmsx_2D = (TH2D*)gDirectory->Get("rmsx_2D");
    TH2D *rmsy_2D = (TH2D*)gDirectory->Get("rmsy_2D");

    // pad1: RMS x 2D分布
    c1->cd(1);
    rmsx_2D->Draw("colz1");
    TLatex *title1 = new TLatex();
    title1->SetTextAlign(22);
    title1->SetTextSize(0.06);
    title1->SetTextColor(global_darkmode ? 0 : 1);
    title1->DrawLatexNDC(0.5, 0.95, fmt::format("RMS x PL{:03d}-PL{:03d}", pl[0], pl[1]).c_str());

    // pad2: RMS y 2D分布
    c1->cd(2);
    rmsy_2D->Draw("colz");
    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    title2->DrawLatexNDC(0.5, 0.95, fmt::format("RMS y PL{:03d}-PL{:03d}", pl[0], pl[1]).c_str());

    // 1Dヒストグラム作成
    TH1D *rmsx_1D = new TH1D("rmsx_1D", ";RMS x;Area", 10000, 0.0, 1000.0);
    TH1D *rmsy_1D = new TH1D("rmsy_1D", ";RMS y;Area", 10000, 0.0, 1000.0);

    tree->Draw("rms_x >> rmsx_1D", "rms_x > 0", "goff");
    tree->Draw("rms_y >> rmsy_1D", "rms_y > 0", "goff");

    double rms_max = std::max(tree->GetMaximum("rms_x"), tree->GetMaximum("rms_y")) + 0.2;

    // pad3: RMS x 1D
    c1->cd(3);
    rmsx_1D->GetXaxis()->SetRangeUser(0.0, rms_max);
    rmsx_1D->SetFillStyle(0);
    rmsx_1D->SetLineWidth(2);
    rmsx_1D->Draw();
    MyUtil::PaintBins(rmsx_1D, 0.0, rms_max);
    TLatex *title3 = new TLatex();
    title3->SetTextAlign(22);
    title3->SetTextSize(0.06);
    title3->SetTextColor(global_darkmode ? 0 : 1);
    title3->DrawLatexNDC(0.5, 0.95, "RMS x");

    TLegend *lg1 = new TLegend(0.73, 0.7, 0.95, 0.9);
    lg1->SetFillStyle(0);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(0.04);
    lg1->SetTextColor(global_darkmode ? 0 : 1);
    lg1->AddEntry(rmsx_1D, fmt::format("{} areas", entries).c_str(), "");
    lg1->AddEntry(rmsx_1D, fmt::format("Mean      {:.2f}", rmsx_1D->GetMean()).c_str(), "");
    lg1->AddEntry(rmsx_1D, fmt::format("Std Dev   {:.2f}", rmsx_1D->GetStdDev()).c_str(), "");
    lg1->Draw();

    // pad4: RMS y 1D
    c1->cd(4);
    rmsy_1D->GetXaxis()->SetRangeUser(0.0, rms_max);
    rmsy_1D->SetFillStyle(0);
    rmsy_1D->SetLineWidth(2);
    rmsy_1D->Draw();
    MyUtil::PaintBins(rmsy_1D, 0.0, rms_max);
    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    title4->DrawLatexNDC(0.43, 0.95, "RMS y");

    TLegend *lg2 = new TLegend(0.665, 0.7, 0.885, 0.9);
    lg2->SetFillStyle(0);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.04);
    lg2->SetTextColor(global_darkmode ? 0 : 1);
    lg2->AddEntry(rmsy_1D, fmt::format("{} areas", entries).c_str(), "");
    lg2->AddEntry(rmsy_1D, fmt::format("Mean      {:.2f}", rmsy_1D->GetMean()).c_str(), "");
    lg2->AddEntry(rmsy_1D, fmt::format("Std Dev   {:.2f}", rmsy_1D->GetStdDev()).c_str(), "");
    lg2->Draw();
}

void shift_axay(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries, const Affine& afa_g,
    double resolution, double arr, double ref) noexcept
{
    double LowX = AreaParam[1];
    double UpX = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY = AreaParam[4];
    double RangeX = AreaParam[5];
    double RangeY = AreaParam[6];

    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(1.2, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    // resolutionの値が0や>1でも必ず1以上のステップ数になるように保護
    double res = (std::isfinite(resolution) && resolution > 0.0) ? resolution : 1.0;
    if (res > 1.0) res = 1.0;
    uint32_t loop_step = static_cast<uint32_t>(std::max(1.0, static_cast<double>(std::lround(1.0 / res))));
    uint32_t loop_count = (entries + loop_step - 1) / loop_step;

    if (!afa_g.isInvertible()) {
        std::cerr << "Warning: global angular affine is non-invertible; shift_axay skipped." << std::endl;
        return;
    }
    const Affine afa_g_inv = afa_g.inverse();

    std::string shiftaxay_title;
    if (resolution == 1.0) {
        shiftaxay_title = fmt::format("Angle Shift PL{:03d}-PL{:03d};x [mm];y [mm]", pl[0], pl[1]);
    } else {
        shiftaxay_title = fmt::format("Angle Shift PL{:03d}-PL{:03d} ({} / {} areas);x [mm];y [mm]",
                                      pl[0], pl[1], loop_count, entries);
    }
    TH2D *shiftaxay_frame = new TH2D("shiftaxay_frame", shiftaxay_title.c_str(),
                                    100, LowX, UpX, 100, LowY, UpY);
    shiftaxay_frame->Draw();

    double extz = 1.5 * arr;

    // 既存のブランチアドレスをリセットし、正しいポインタに付け替える
    tree->ResetBranchAddresses();

    double x, y, a_ang, b_ang, c_ang, d_ang, p_ang, q_ang;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("a_ang", &a_ang);
    tree->SetBranchAddress("b_ang", &b_ang);
    tree->SetBranchAddress("c_ang", &c_ang);
    tree->SetBranchAddress("d_ang", &d_ang);
    tree->SetBranchAddress("p_ang", &p_ang);
    tree->SetBranchAddress("q_ang", &q_ang);

    for (uint32_t i = 0; i < entries; i += loop_step) {
        tree->GetEntry(i);

        // calculate angular displacement due to local correction
        Affine afa(a_ang, b_ang, c_ang, d_ang, p_ang, q_ang);
        Affine afa_l = afa * afa_g_inv;
        double dx = 0.0, dy = 0.0;
        afa_l.transform(0.0, 0.0, &dx, &dy);

        TArrow *arrow = new TArrow(
            x, y,
            x + dx * extz,
            y + dy * extz,
            0.01,
            ">"
        );
        arrow->SetFillColor(93);
        arrow->SetLineColor(93);
        arrow->SetFillStyle(1001);
        arrow->SetLineWidth(1);
        arrow->SetLineStyle(1);
        arrow->Draw();
    }

    // reference (ref is expected in rad; display in mrad)
    double ref_length = ref * extz;
    std::string ref_text = fmt::format("{:.3g} mrad", ref);
    TPaveText *ref_pt = new TPaveText(
        UpX + RangeX * 0.08, LowY + RangeY * 0.07,
        UpX + RangeX * 0.28, LowY + RangeY * 0.17,
        "br"
    );
    ref_pt->SetTextColor(global_darkmode ? 0 : 1);
    ref_pt->SetFillColor(global_darkmode ? 1 : 0);
    ref_pt->SetLineColor(global_darkmode ? 1 : 0);
    ref_pt->SetShadowColor(global_darkmode ? 1 : 0);
    TText *text_axay = ref_pt->AddText(ref_text.c_str());
    (void)text_axay;
    ref_pt->Draw();

    TArrow *arrowX = new TArrow(
        UpX + RangeX * 0.05, LowY + RangeY * 0.05,
        UpX + RangeX * 0.05 + ref_length, LowY + RangeY * 0.05, 0.01,
        ">"
    );
    arrowX->SetFillColor(global_darkmode ? 0 : 1);
    arrowX->SetLineColor(global_darkmode ? 0 : 1);
    arrowX->SetFillStyle(1001);
    arrowX->SetLineWidth(2);
    arrowX->SetLineStyle(1);
    arrowX->Draw();

    TArrow *arrowY = new TArrow(
        UpX + RangeX * 0.05, LowY + RangeY * 0.05,
        UpX + RangeX * 0.05, LowY + RangeY * 0.05 + ref_length, 0.01,
        ">"
    );
    arrowY->SetFillColor(global_darkmode ? 0 : 1);
    arrowY->SetLineColor(global_darkmode ? 0 : 1);
    arrowY->SetFillStyle(1001);
    arrowY->SetLineWidth(2);
    arrowY->SetLineStyle(1);
    arrowY->Draw();
}

void shift_axays(TCanvas *c1, TTree *tree, const double *AreaParam, uint32_t entries) noexcept
{
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.5, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    c1->Clear();
    c1->Divide(3, 2);

    // 既存のブランチアドレスをリセットし、正しいポインタに付け替える
    tree->ResetBranchAddresses();

    double x, y, p_ang, q_ang;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("p_ang", &p_ang);
    tree->SetBranchAddress("q_ang", &q_ang);

    // マーカー設定
    constexpr int markerstyle = 20;
    constexpr double markersize = 0.2;

    TH1D *shift_ax_temp = new TH1D("shift_ax_temp", "", 2000, -1000.0, 1000.0);
    tree->Draw("p_ang >> shift_ax_temp", "", "goff");
    double sax_mean = shift_ax_temp->GetMean();
    double sax_stddev = shift_ax_temp->GetStdDev();
    double sax_range_low = sax_mean - 5 * sax_stddev;
    double sax_range_up  = sax_mean + 5 * sax_stddev;
    gDirectory->Delete("shift_ax_temp");
    TH1D *shift_ay_temp = new TH1D("shift_ay_temp", "", 2000, -1000.0, 1000.0);
    tree->Draw("q_ang >> shift_ay_temp", "", "goff");
    double say_mean = shift_ay_temp->GetMean();
    double say_stddev = shift_ay_temp->GetStdDev();
    double say_range_low = say_mean - 5 * say_stddev;
    double say_range_up  = say_mean + 5 * say_stddev;
    gDirectory->Delete("shift_ay_temp");
    double range_low = std::min(sax_range_low, say_range_low);
    double range_up  = std::max(sax_range_up, say_range_up);

    // pad1: ax shift ヒストグラム
    c1->cd(1);
    std::string shiftax_title = "#it{#theta}_{x} shift;#it{#theta}_{x} shift [mrad];";
    TH1D *shift_ax = new TH1D("shift_ax", shiftax_title.c_str(), 100, range_low, range_up);
    tree->Draw("p_ang >> shift_ax", "", "goff");
    shift_ax->SetFillColor(90);
    shift_ax->Draw();

    TLegend *sax_lg = new TLegend(0.5, 0.75, 0.75, 0.85);
    sax_lg->SetFillStyle(0);
    sax_lg->SetBorderSize(0);
    sax_lg->SetTextSize(0.04);
    sax_lg->AddEntry(shift_ax, fmt::format("Mean   {:.2g}", sax_mean).c_str(), "");
    sax_lg->AddEntry(shift_ax, fmt::format("Std Dev   {:.2g}", sax_stddev).c_str(), "");
    sax_lg->Draw();

    // TGraph生成: p_ang のデータセット
    TGraph *gdax_x = new TGraph();
    TGraph *gdax_y = new TGraph();
    TGraph *gday_x = new TGraph();
    TGraph *gday_y = new TGraph();

    for (uint32_t i = 0; i < entries; ++i) {
        tree->GetEntry(i);

        // SetPoint(i, xc, p_ang) だと (0, 0) に点が打たれてしまうため GetN() を使っている
        gdax_x->SetPoint(gdax_x->GetN(), x, p_ang);
        gdax_y->SetPoint(gdax_y->GetN(), y, p_ang);
        gday_x->SetPoint(gday_x->GetN(), x, q_ang);
        gday_y->SetPoint(gday_y->GetN(), y, q_ang);
    }

    // pad2: ax shift : x 散布図
    c1->cd(2);
    std::string dax_x_title = "#it{#theta}_{x} shift : x;x [mm];#it{#theta}_{x} shift [mrad]";
    gdax_x->SetTitle(dax_x_title.c_str());
    gdax_x->SetMarkerStyle(markerstyle);
    gdax_x->SetMarkerSize(markersize);
    gdax_x->SetMarkerColor(90);
    gdax_x->GetXaxis()->SetLimits(LowX, UpX);
    gdax_x->GetYaxis()->SetRangeUser(range_low, range_up);
    gdax_x->Draw("AP");

    // pad3: ax shift : y 散布図
    c1->cd(3);
    std::string dax_y_title = "#it{#theta}_{x} shift : y;y [mm];#it{#theta}_{x} shift [mrad]";
    gdax_y->SetTitle(dax_y_title.c_str());
    gdax_y->SetMarkerStyle(markerstyle);
    gdax_y->SetMarkerSize(markersize);
    gdax_y->SetMarkerColor(90);
    gdax_y->GetXaxis()->SetLimits(LowY, UpY);
    gdax_y->GetYaxis()->SetRangeUser(range_low, range_up);
    gdax_y->Draw("AP");

    // pad4: ay shift ヒストグラム
    c1->cd(4);
    std::string shiftay_title = "#it{#theta}_{y} shift;#it{#theta}_{y} shift [mrad];";
    TH1D *shift_ay = new TH1D("shift_ay", shiftay_title.c_str(), 100, range_low, range_up);
    tree->Draw("q_ang >> shift_ay", "", "goff");
    shift_ay->SetFillColor(91);
    shift_ay->Draw();

    TLegend *say_lg = new TLegend(0.5, 0.75, 0.75, 0.85);
    say_lg->SetFillStyle(0);
    say_lg->SetBorderSize(0);
    say_lg->SetTextSize(0.04);
    say_lg->AddEntry(shift_ay, fmt::format("Mean   {:.2g}", say_mean).c_str(), "");
    say_lg->AddEntry(shift_ay, fmt::format("Std Dev   {:.2g}", say_stddev).c_str(), "");
    say_lg->Draw();

    // pad5: ay shift : x 散布図
    c1->cd(5);
    std::string day_x_title = "#it{#theta}_{y} shift : x;x [mm];#it{#theta}_{y} shift [mrad]";
    gday_x->SetTitle(day_x_title.c_str());
    gday_x->SetMarkerStyle(markerstyle);
    gday_x->SetMarkerSize(markersize);
    gday_x->SetMarkerColor(91);
    gday_x->GetXaxis()->SetLimits(LowX, UpX);
    gday_x->GetYaxis()->SetRangeUser(range_low, range_up);
    gday_x->Draw("AP");

    // pad6: ay shift : y 散布図
    c1->cd(6);
    std::string day_y_title = "#it{#theta}_{y} shift : y;y [mm];#it{#theta}_{y} shift [mrad]";
    gday_y->SetTitle(day_y_title.c_str());
    gday_y->SetMarkerStyle(markerstyle);
    gday_y->SetMarkerSize(markersize);
    gday_y->SetMarkerColor(91);
    gday_y->GetXaxis()->SetLimits(LowY, UpY);
    gday_y->GetYaxis()->SetRangeUser(range_low, range_up);
    gday_y->Draw("AP");
}

void rms_axay(TCanvas *c1, TTree *tree, const double *AreaParam,
    const uint32_t pl[2], uint32_t entries) noexcept
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

    double bin = AreaParam[0];
    double LowX = AreaParam[1];
    double UpX = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY = AreaParam[4];

    // 2D ヒストグラムの作成
    TH2D *rmsax_2D = new TH2D("rmsax_2D", ";x [mm];y [mm];RMS#it{#theta}_{x} [mrad]",
        bin, LowX, UpX, bin, LowY, UpY);
    TH2D *rmsay_2D = new TH2D("rmsay_2D", ";x [mm];y [mm];RMS#it{#theta}_{y} [mrad]",
        bin, LowX, UpX, bin, LowY, UpY);

    // 既存のブランチアドレスをリセット
    tree->ResetBranchAddresses();
    double x, y, rms_ax, rms_ay;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("rms_ax", &rms_ax);
    tree->SetBranchAddress("rms_ay", &rms_ay);

    // データを2Dヒストグラムにセット
    for (uint32_t j = 0; j < entries; ++j) {
        tree->GetEntry(j);
        int bx = rmsax_2D->GetXaxis()->FindBin(x);
        int by = rmsax_2D->GetYaxis()->FindBin(y);
        if (bx > 0 && bx <= bin && by > 0 && by <= bin) {
            rmsax_2D->SetBinContent(bx, by, rms_ax);
            rmsay_2D->SetBinContent(bx, by, rms_ay);
        }
    }

    TH1D *rmsax_temp = new TH1D("rmsax_temp", "", 10000000, 0.0, 1000.0);
    TH1D *rmsay_temp = new TH1D("rmsay_temp", "", 10000000, 0.0, 1000.0);
    tree->Draw("rms_ax >> rmsax_temp", "rms_ax > 0", "goff");
    tree->Draw("rms_ay >> rmsay_temp", "rms_ay > 0", "goff");
    double rms_stddev = std::max(rmsax_temp->GetStdDev(), rmsay_temp->GetStdDev());
    double rms_max = std::max(tree->GetMaximum("rms_ax"), tree->GetMaximum("rms_ay")) + rms_stddev;
    // rmsの最小単位が0.0001のため、1binが0.0001未満にならないように調整
    int rms_bin = (rms_max < 0.01) ? std::floor(rms_max / 0.0001) : 100;
    gDirectory->Delete("*_temp");

    rmsax_2D->GetZaxis()->SetRangeUser(0.0, rms_max);
    rmsay_2D->GetZaxis()->SetRangeUser(0.0, rms_max);

    // 1Dヒストグラム作成
    TH1D *rmsax_1D = new TH1D("rmsax_1D", ";RMS#it{#theta}_{x} [mrad];Area", rms_bin, 0.0, rms_max);
    TH1D *rmsay_1D = new TH1D("rmsay_1D", ";RMS#it{#theta}_{y} [mrad];Area", rms_bin, 0.0, rms_max);

    tree->Draw("rms_ax >> rmsax_1D", "rms_ax > 0", "goff");
    tree->Draw("rms_ay >> rmsay_1D", "rms_ay > 0", "goff");

    // pad1: RMS ax 2D分布
    c1->cd(1);
    rmsax_2D->Draw("colz");
    TLatex *title1 = new TLatex();
    title1->SetTextAlign(22);
    title1->SetTextSize(0.06);
    title1->SetTextColor(global_darkmode ? 0 : 1);
    title1->DrawLatexNDC(0.5, 0.95,
        fmt::format("RMS#it{{#theta}}_{{x}} PL{:03d}-PL{:03d}", pl[0], pl[1]).c_str());

    // pad2: RMS ay 2D分布
    c1->cd(2);
    rmsay_2D->Draw("colz");
    TLatex *title2 = new TLatex();
    title2->SetTextAlign(22);
    title2->SetTextSize(0.06);
    title2->SetTextColor(global_darkmode ? 0 : 1);
    title2->DrawLatexNDC(0.43, 0.95,
        fmt::format("RMS#it{{#theta}}_{{y}} PL{:03d}-PL{:03d}", pl[0], pl[1]).c_str());

    // pad3: RMS ax 1D
    c1->cd(3);
    rmsax_1D->GetXaxis()->SetRangeUser(0.0, rms_max);
    rmsax_1D->SetFillStyle(0);
    rmsax_1D->SetLineWidth(2);
    rmsax_1D->Draw();
    MyUtil::PaintBins(rmsax_1D, 0.0, rms_max);

    TLatex *title3 = new TLatex();
    title3->SetTextAlign(22);
    title3->SetTextSize(0.06);
    title3->SetTextColor(global_darkmode ? 0 : 1);
    title3->DrawLatexNDC(0.5, 0.95, "RMS#it{#theta}_{x}");

    TLegend *lg1 = new TLegend(0.73, 0.7, 0.95, 0.9);
    lg1->SetFillStyle(0);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(0.04);

    lg1->SetTextColor(global_darkmode ? 0 : 1);
    lg1->AddEntry(rmsax_1D, fmt::format("{} areas", entries).c_str(), "");
    lg1->AddEntry(rmsax_1D, fmt::format("Mean      {:.2g}", rmsax_1D->GetMean()).c_str(), "");
    lg1->AddEntry(rmsax_1D, fmt::format("Std Dev   {:.2g}", rmsax_1D->GetStdDev()).c_str(), "");
    lg1->Draw();

    // pad4: RMS ay 1D
    c1->cd(4);
    rmsay_1D->GetXaxis()->SetRangeUser(0.0, rms_max);
    rmsay_1D->SetFillStyle(0);
    rmsay_1D->SetLineWidth(2);
    rmsay_1D->Draw();
    MyUtil::PaintBins(rmsay_1D, 0.0, rms_max);

    TLatex *title4 = new TLatex();
    title4->SetTextAlign(22);
    title4->SetTextSize(0.06);
    title4->SetTextColor(global_darkmode ? 0 : 1);
    title4->DrawLatexNDC(0.43, 0.95, "RMS#it{#theta}_{y}");

    TLegend *lg2 = new TLegend(0.665, 0.7, 0.885, 0.9);
    lg2->SetFillStyle(0);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.04);
    lg2->SetTextColor(global_darkmode ? 0 : 1);
    lg2->AddEntry(rmsay_1D, fmt::format("{} areas", entries).c_str(), "");
    lg2->AddEntry(rmsay_1D, fmt::format("Mean      {:.2g}", rmsay_1D->GetMean()).c_str(), "");
    lg2->AddEntry(rmsay_1D, fmt::format("Std Dev   {:.2g}", rmsay_1D->GetStdDev()).c_str(), "");
    lg2->Draw();
}