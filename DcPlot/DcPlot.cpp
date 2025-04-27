#include <iostream>
#include <fstream>
#include <sstream>
#include <csignal>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TError.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TColor.h>
#include <TArrow.h>
#include <TPaveText.h>
#include <TCut.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include <argparse/argparse.hpp>

#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>

void dist_arrow(TCanvas *c1, TTree *tree[2], const double *AreaParam, const double resolution,
    const double reference, int32_t *entries, const double extz, const uint8_t face) noexcept;
void dist_map(TCanvas *c1, TTree *tree[2], const double *AreaParam, const double dist_max,
    int32_t *entries, const double markersize) noexcept;
void shr_map(TCanvas *c1, TTree *tree[2], const double *AreaParam, const std::vector<double> shr_range,
    int32_t *entries, const double markersize) noexcept;
void sig_bg_map(TCanvas *c1, TTree *tree, const double *AreaParam, double sn_max[2],
    int32_t *entries, const double markersize) noexcept;
void dist_shr_map(TCanvas *c1, TTree *tree[2], const double dist_max,
    const std::vector<double> shr_range, const double sig_max, int32_t *entries) noexcept;

namespace {
    // Ctrl+Cで終了したときの処理用にグローバル変数を定義
    TCanvas *global_c1 = nullptr;
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

void parse_arguments(argparse::ArgumentParser& parser, int argc, char *argv[]) {
    parser.set_usage_max_line_width(80);
    parser.add_description("Tips: You can combine single-character arguments.\n"
        "      For example, \"-d -p kBird -i\" is equivalent to \"-dpi kBird\".");

    // 必須引数: bvxxファイル
    parser.add_argument("input_dc")
        .help("Path to the bvxx file to be processed.")
        .required();
    // オプション引数: 出力ファイル名
    parser.add_argument("-o", "--output")
        .help("Output PDF file name. [default: input_dc.pdf]")
        .default_value(std::string());
    // オプション引数: その他の設定
    parser.add_argument("-res", "--resolution")
        .help("Arrow plot resolution (0.0 - 1.0).")
        .default_value(1.0)
        .scan<'g', double>();
    parser.add_argument("-ref", "--reference")
        .help("Reference arrow length (mrad).")
        .default_value(50.0)
        .scan<'g', double>();
    parser.add_argument("-arr", "--arrow_length")
        .help("Plot arrow length.")
        .default_value(1.0)
        .scan<'g', double>();
    parser.add_argument("-dist", "--distortion_max")
        .help("Maximum distortion range.")
        .default_value(0.2)
        .scan<'g', double>();
    parser.add_argument("-shr", "--shrink_range")
        .help("Shrink range.\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({0.7, 1.8}))
        .nargs(2)
        .scan<'g', double>();
    parser.add_argument("-sig", "--signal_max")
        .help("Maximum signal range (default: maximum signal value).")
        .scan<'g', double>();
    parser.add_argument("-bg", "--background_max")
        .help("Maximum background range (default: maximum background value).")
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
    const auto dcfile = parser.get<std::string>("input_dc");
    const auto output_arg = parser.get<std::string>("--output");
    const auto resolution = parser.get<double>("--resolution");
    const auto reference = parser.get<double>("--reference");
    const auto arr = parser.get<double>("--arrow_length");
    auto dist_max = parser.get<double>("--distortion_max");
    auto shr_range = parser.get<std::vector<double>>("--shrink_range");
    auto font_number = parser.get<int>("--font_number");
    auto hideGrid = parser.get<bool>("--hide_grid");
    global_darkmode = parser.get<bool>("--dark_mode");
    auto palette_arg = parser.get<std::string>("--palette");
    auto NContours = parser.get<int>("--contours");
    auto invertpalette = parser.get<bool>("--invert_palette");
    auto negatepalette = parser.get<bool>("--negate_palette");

    // 出力ファイル名の設定
    const std::string output = (output_arg.empty())
        ? (dcfile + ".pdf")
        : ((output_arg.size() > 4 && output_arg.substr(output_arg.size() - 4) == ".pdf")
            ? output_arg
            : output_arg + ".pdf");

    // 引数を利用する変数を設定
    const int font_code = 10 * font_number + 2;
    const double extz = arr * 250000.0;
    double sn_max[2];
    sn_max[0] = parser.is_used("-sig") ? parser.get<double>("--signal_max") : 0.0;
    sn_max[1] = parser.is_used("-bg") ? parser.get<double>("--background_max") : 0.0;

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

    // TTreeの作成
    TTree *tree[2];
    tree[0] = new TTree("tree[0]", ""); // face1
    tree[1] = new TTree("tree[1]", ""); // face2

    // TTreeのブランチを作成
    uint32_t id, pos, pos_;
    double x, y, xmin, xmax, ymin, ymax, shr, ddz, shr_, dax, day, dz, sig, bg, signif, dist;
    double dmy0, dmy1, dmy2, dmy3, dmy4, dmy5, dmy6;
    const std::array<std::pair<const char*, void*>, 9> branches = {{
        {"x", &x}, {"y", &y}, {"shr", &shr}, {"dax", &dax}, {"day", &day}, {"dz", &dz},
        {"sig", &sig}, {"bg", &bg}, {"dist", &dist}
    }};
    for (const auto& branch : branches) {
        tree[0]->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
        tree[1]->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    // ファイルの読み込み
    std::ifstream ifs(dcfile.c_str());
    if (!ifs) {
        std::cerr << "Error: Cannot open file: " << dcfile << std::endl;
        return 1;
    }
    int32_t entries[2] = {0, 0};
    double view = DBL_MAX;
    std::string line;
    while (std::getline(ifs, line)) {
        std::istringstream iss(line);
        if (!(iss >> id >> pos >> pos_ >> xmin >> xmax >> ymin >> ymax
                  >> dmy0 >> dmy1 >> dmy2 >> dmy3 >> dmy4 >> dmy5
                  >> shr >> ddz >> dmy6 >> shr_ >> dax >> day >> dz >> sig >> bg >> signif)) {
            std::cerr << "Warning: Malformed line skipped: " << line << std::endl;
            continue;
        }

        x = (xmax + xmin) / 2000.0;
        y = (ymax + ymin) / 2000.0;
        dist = sqrt(dax * dax + day * day);
        view = std::min(view, xmax - xmin);

        uint8_t num = (pos & 1) ^ 1;
        entries[num]++;
        tree[num]->Fill();
    }
    ifs.close();

    // 情報表示
    std::cout << "\n File       : " << dcfile << std::endl;
    std::cout << " Arrow      : " << Form("%.1f", arr) << std::endl;
    std::cout << " Resolution : " << Form("%.1f", resolution) << std::endl;
    std::cout << " View       : " << view << std::endl;
    std::cout << " # of Areas : " << entries[0] << " (face1) & " << entries[1] << " (face2)" << std::endl;

    // プロット開始
    TDatime time_now;
    std::string Time_Now = Form(
        "%d-%02d-%02d %02d:%02d:%02d",
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(),
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    sw.Start();
	std::cout << " Plot start : " << Time_Now << std::endl;

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas *c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());
    global_c1 = c1;
    global_output = output;

	// プログレスバーの初期化
	int page = 0;
    constexpr int total = 10; // 合計ページ数
	MyUtil::ShowProgress(page, total);

    // データの座標の範囲を取得し、表示範囲とビンの数を決定する
    // フィルムの長辺の端から1cm外側までを最大表示範囲とし、縦横比を正しく保って表示する
    const double MinX = tree[0]->GetMinimum("x");
    const double MaxX = tree[0]->GetMaximum("x");
    const double MinY = tree[0]->GetMinimum("y");
    const double MaxY = tree[0]->GetMaximum("y");
    const double RangeX = MaxX - MinX; // データ領域
    const double RangeY = MaxY - MinY; // データ領域
    double LowX, UpX, LowY, UpY, bin; // 表示範囲とビンの数
    double markersize; // マーカーサイズ
    if (RangeX >= RangeY) {
        LowX = MinX - 10;
        UpX = MaxX + 10;
        LowY = MinY - (RangeX - RangeY + 20) * 0.5;
        UpY = MaxY + (RangeX - RangeY + 20) * 0.5;
        bin = (RangeX + 20) * 1000 / view;
        if (RangeX < 100) {
            markersize = view * 0.00023;
        } else if (RangeX < 150) {
            markersize = view * 0.00011;
        } else if (RangeX < 300) {
            markersize = view * 0.00006;
        }
    } else {
        LowX = MinX - (RangeY - RangeX + 20) * 0.5;
        UpX = MaxX + (RangeY - RangeX + 20) * 0.5;
        LowY = MinY - 10;
        UpY = MaxY + 10;
        bin = (RangeY + 20) * 1000 / view;
        if (RangeY < 100) {
            markersize = view * 0.00023;
        } else if (RangeY < 150) {
            markersize = view * 0.00011;
        } else if (RangeY < 300) {
            markersize = view * 0.00006;
        }
    }
    const double AreaParam[7] = {bin, LowX, UpX, LowY, UpY, RangeX, RangeY};

    dist_arrow(c1, tree, AreaParam, resolution, reference, entries, extz, 1);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);

    dist_arrow(c1, tree, AreaParam, resolution, reference, entries, extz, 2);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);

    dist_map(c1, tree, AreaParam, dist_max, entries, markersize);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");

    shr_map(c1, tree, AreaParam, shr_range, entries, markersize);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");

    sig_bg_map(c1, tree[0], AreaParam, sn_max, entries, markersize);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("*_1D");

    dist_shr_map(c1, tree, dist_max, shr_range, sn_max[0], entries);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, total);
    gDirectory->Delete("dist*");
    gDirectory->Delete("shr*");

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
    double elapsed_time = sw.RealTime();
    double cpu_time = sw.CpuTime();
    std::cout << Form(
        "\n Plot end   : %s - Elapsed %.2f [s] (CPU: %.2f [s])",
        Time_Now.data(), elapsed_time, cpu_time
    ) << std::endl;
    std::cout << " Output     : " << output << std::endl;

    delete c1;
    gDirectory->Delete("tree*");

    return 0;
}

void dist_arrow(TCanvas *c1, TTree *tree[2], const double *AreaParam, const double resolution,
    const double reference, int32_t *entries, const double extz, const uint8_t face) noexcept
{
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];
    double RangeX = AreaParam[5];
    double RangeY = AreaParam[6];
    uint8_t num = face - 1;
    const uint32_t loop_step = static_cast<int>(1.0 / resolution);
    const uint32_t loop_count = (entries[num] + loop_step - 1) / loop_step;

    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(1.6, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    std::string title;
    if (resolution == 1.0) {
        title = Form("Distortion %d;x [mm];y [mm]", face);
    } else {
        title = Form("Distortion %d (%d / %d areas);x [mm];y [mm]", face, loop_count, entries[num]);
    }

    TH1 *frame = c1->DrawFrame(LowX, LowY, UpX, UpY, title.data());

    double x, y, dax, day;
    tree[num]->SetBranchAddress("x", &x);
    tree[num]->SetBranchAddress("y", &y);
    tree[num]->SetBranchAddress("dax", &dax);
    tree[num]->SetBranchAddress("day", &day);

    for (int i = 0; i < entries[num]; i += loop_step) {
        tree[num]->GetEntry(i);
        TArrow *arrow = new TArrow(x, y, x + dax * 0.001 * extz, y + day * 0.001 * extz, 0.01, ">");
        arrow->SetLineColor(95 + face);
        arrow->SetFillColor(95 + face);
        arrow->SetLineWidth(1);
        arrow->Draw();
    }

    const double ref_length = 0.000001 * reference * extz;
    TPaveText *ref_pt = new TPaveText(
        UpX + RangeX * 0.08, LowY + RangeY * 0.07,
        UpX + RangeX * 0.28, LowY + RangeY * 0.17,
        "br"
    );
    ref_pt->SetTextColor(global_darkmode ? 0 : 1);
    ref_pt->SetFillColor(global_darkmode ? 1 : 0);
    ref_pt->SetLineColor(global_darkmode ? 1 : 0);
    ref_pt->SetShadowColor(global_darkmode ? 1 : 0);
    const char *shiftxy_text = Form("%.0f mrad", reference);
    TText *text_xy = ref_pt->AddText(shiftxy_text);
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

void dist_map(TCanvas *c1, TTree *tree[2], const double *AreaParam, const double dist_max,
    int32_t *entries, const double markersize) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    TGraph2D *gframe = new TGraph2D;
    gframe->SetPoint(0, LowX, LowY, 0);
    gframe->SetPoint(1, UpX, UpY, 0);
    gframe->SetMarkerColor(global_darkmode ? 1 : 0);
    gframe->GetHistogram()->GetXaxis()->SetTickLength(-0.03);
    gframe->GetHistogram()->GetYaxis()->SetTickLength(-0.03);
    gframe->GetXaxis()->SetLabelOffset(0.05);
    gframe->GetYaxis()->SetLabelOffset(0.1);
    gframe->GetHistogram()->SetMinimum(0.0);
    gframe->GetHistogram()->SetMaximum(dist_max);

    for (int i = 0; i < 2; ++i) {
        TGraph2D *dist_graph = new TGraph2D;
        double x, y, dist;
        tree[i]->SetBranchAddress("x", &x);
        tree[i]->SetBranchAddress("y", &y);
        tree[i]->SetBranchAddress("dist", &dist);

        for (int j = 0; j < entries[i]; ++j) {
            tree[i]->GetEntry(j);
            dist_graph->SetPoint(dist_graph->GetN(), x, y, dist);
        }

        c1->cd(i + 1);
        gPad->SetTheta(90);
        gPad->SetPhi(1.e-14);
        gPad->SetGrid(0, 0);
        gframe->SetTitle(";x [mm];y [mm];");
        gframe->Draw("p");
        dist_graph->SetMarkerStyle(21);
        dist_graph->SetMarkerSize(markersize);
        dist_graph->Draw("pcolzsame");

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatex(0.0, 0.65, Form("|Distortion %d|", i + 1));

        TLine *l0, *l1, *l2, *l3;
        l0 = new TLine(-0.574, -0.574, -0.574, 0.574);
        l1 = new TLine(0.574, -0.574, 0.574, 0.574);
        l2 = new TLine(-0.574, -0.574, 0.574, -0.574);
        l3 = new TLine(-0.574, 0.574, 0.574, 0.574);
        for (TLine *line : {l0, l1, l2, l3}) {
            if (line) {
                line->SetLineWidth(1);
                line->SetLineColor(global_darkmode ? 0 : 1);
                line->Draw("same");
            }
        }
    }

    TH1D *dist1_1D = new TH1D(
        "dist1_1D",
        ";|Distortion 1| = #sqrt{#Deltatan^{2}#it{#theta}_{x1}#plus#Deltatan^{2}#it{#theta}_{y1}}"
        ";Area",
        100, 0.0, dist_max
    );
    TH1D *dist2_1D = new TH1D(
        "dist2_1D",
        ";|Distortion 2| = #sqrt{#Deltatan^{2}#it{#theta}_{x2}#plus#Deltatan^{2}#it{#theta}_{y2}}"
        ";Area",
        100, 0.0, dist_max
    );

    for (int i = 0; i < 2; ++i) {
        TH1D *hist = (i == 0) ? dist1_1D : dist2_1D;
        tree[i]->Draw(Form("dist >> dist%d_1D", i + 1), "", "goff");
        c1->cd(i + 3);
        hist->SetFillStyle(0);
        hist->SetLineWidth(2);
        hist->Draw();
        MyUtil::PaintBins(hist, 0.0, dist_max); // 各ビンをカラーパレットの色で塗る

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatexNDC((i == 0) ? 0.5 : 0.43, 0.95, Form("|Distortion %d|", i + 1));

        TLegend *lg = new TLegend(
            (i == 0) ? 0.73 : 0.665,
            0.7,
            (i == 0) ? 0.95 : 0.885,
            0.9
        );
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.04);
        lg->SetTextColor(global_darkmode ? 0 : 1);
        lg->AddEntry(hist, Form("%d areas", entries[i]), "");
        lg->AddEntry(hist, Form("Mean      %.5f", hist->GetMean()), "");
        lg->AddEntry(hist, Form("Std Dev   %.5f", hist->GetStdDev()), "");
        lg->Draw();
    }
}

void shr_map(TCanvas *c1, TTree *tree[2], const double *AreaParam, const std::vector<double> shr_range,
    int32_t *entries, const double markersize) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    TGraph2D *gframe = new TGraph2D;
    gframe->SetPoint(0, LowX, LowY, 0);
    gframe->SetPoint(1, UpX, UpY, 0);
    gframe->SetMarkerColor(global_darkmode ? 1 : 0);
    gframe->GetHistogram()->GetXaxis()->SetTickLength(-0.03);
    gframe->GetHistogram()->GetYaxis()->SetTickLength(-0.03);
    gframe->GetXaxis()->SetLabelOffset(0.05);
    gframe->GetYaxis()->SetLabelOffset(0.1);
    gframe->GetHistogram()->SetMinimum(shr_range[0]);
    gframe->GetHistogram()->SetMaximum(shr_range[1]);

    for (int i = 0; i < 2; ++i) {
        TGraph2D *shr_graph = new TGraph2D;
        double x, y, shr, dist;
        tree[i]->SetBranchAddress("x", &x);
        tree[i]->SetBranchAddress("y", &y);
        tree[i]->SetBranchAddress("shr", &shr);

        // distブランチは使わないが、設定しないとなぜかshrにdistが入る
        tree[i]->SetBranchAddress("dist", &dist);

        for (int j = 0; j < entries[i]; ++j) {
            tree[i]->GetEntry(j);
            shr_graph->SetPoint(shr_graph->GetN(), x, y, shr);
        }

        c1->cd(i + 1);
        gPad->SetTheta(90);
        gPad->SetPhi(1.e-14);
        gPad->SetGrid(0, 0);
        gframe->SetTitle(";x [mm];y [mm];");
        gframe->Draw("p");
        shr_graph->SetMarkerStyle(21);
        shr_graph->SetMarkerSize(markersize);
        shr_graph->Draw("pcolzsame");

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatex(0.0, 0.65, Form("Shrink %d", i + 1));

        TLine *l0, *l1, *l2, *l3;
        l0 = new TLine(-0.574, -0.574, -0.574, 0.574);
        l1 = new TLine(0.574, -0.574, 0.574, 0.574);
        l2 = new TLine(-0.574, -0.574, 0.574, -0.574);
        l3 = new TLine(-0.574, 0.574, 0.574, 0.574);
        for (TLine *line : {l0, l1, l2, l3}) {
            if (line) {
                line->SetLineWidth(1);
                line->SetLineColor(global_darkmode ? 0 : 1);
                line->Draw("same");
            }
        }
    }

    TH1D *shr1_1D = new TH1D("shr1_1D", ";Shrink 1;Area", 100, shr_range[0], shr_range[1]);
    TH1D *shr2_1D = new TH1D("shr2_1D", ";Shrink 2;Area", 100, shr_range[0], shr_range[1]);

    for (int i = 0; i < 2; ++i) {
        TH1D *hist = (i == 0) ? shr1_1D : shr2_1D;
        tree[i]->Draw(Form("shr >> shr%d_1D", i + 1), "", "goff");
        c1->cd(i + 3);
        hist->SetFillStyle(0);
        hist->SetLineWidth(2);
        hist->Draw();
        MyUtil::PaintBins(hist, shr_range[0], shr_range[1]); // 各ビンをカラーパレットの色で塗る

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatexNDC((i == 0) ? 0.5 : 0.43, 0.95, Form("Shrink %d", i + 1));

        TLegend *lg = new TLegend(
            (i == 0) ? 0.73 : 0.665,
            0.7,
            (i == 0) ? 0.95 : 0.885,
            0.9
        );
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.04);
        lg->SetTextColor(global_darkmode ? 0 : 1);
        lg->AddEntry(hist, Form("%d areas", entries[i]), "");
        lg->AddEntry(hist, Form("Mean      %.5f", hist->GetMean()), "");
        lg->AddEntry(hist, Form("Std Dev   %.5f", hist->GetStdDev()), "");
        lg->Draw();
    }
}

void sig_bg_map(TCanvas *c1, TTree *tree, const double *AreaParam, double sn_max[2],
    int32_t *entries, const double markersize) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    double max_array[2];
    max_array[0] = (sn_max[0] > 0.0) ? sn_max[0] : tree->GetMaximum("sig");
    sn_max[0] = max_array[0]; // dist_shr_mapに使う
    max_array[1] = (sn_max[1] > 0.0) ? sn_max[1] : tree->GetMaximum("bg");

    for (int i = 0; i < 2; ++i) {
        TGraph2D *gframe = new TGraph2D;
        gframe->SetPoint(0, LowX, LowY, 0);
        gframe->SetPoint(1, UpX, UpY, 0);
        gframe->SetMarkerColor(global_darkmode ? 1 : 0);
        gframe->GetHistogram()->GetXaxis()->SetTickLength(-0.03);
        gframe->GetHistogram()->GetYaxis()->SetTickLength(-0.03);
        gframe->GetXaxis()->SetLabelOffset(0.05);
        gframe->GetYaxis()->SetLabelOffset(0.1);
        gframe->GetHistogram()->SetMinimum(0.0);
        gframe->GetHistogram()->SetMaximum(max_array[i]);

        TGraph2D *graph = new TGraph2D;
        double x, y, sig, bg, dist, shr;
        tree->SetBranchAddress("x", &x);
        tree->SetBranchAddress("y", &y);
        tree->SetBranchAddress("sig", &sig);
        tree->SetBranchAddress("bg", &bg);

        // distブランチとshrブランチは使わないが、設定しないと出力がおかしくなる
        tree->SetBranchAddress("dist", &dist);
        tree->SetBranchAddress("shr", &shr);

        for (int j = 0; j < entries[i]; ++j) {
            tree->GetEntry(j);
            graph->SetPoint(
                graph->GetN(), x, y,
                (i == 0) ? sig : bg
            );
        }

        c1->cd(i + 1);
        gPad->SetTheta(90);
        gPad->SetPhi(1.e-14);
        gPad->SetGrid(0, 0);
        gframe->SetTitle(";x [mm];y [mm];");
        gframe->Draw("p");
        graph->SetMarkerStyle(21);
        graph->SetMarkerSize(markersize);
        graph->Draw("pcolzsame");

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatex(0.0, 0.65, (i == 0) ? "Signal" : "Background");

        TLine *l0, *l1, *l2, *l3;
        l0 = new TLine(-0.574, -0.574, -0.574, 0.574);
        l1 = new TLine(0.574, -0.574, 0.574, 0.574);
        l2 = new TLine(-0.574, -0.574, 0.574, -0.574);
        l3 = new TLine(-0.574, 0.574, 0.574, 0.574);
        for (TLine *line : {l0, l1, l2, l3}) {
            if (line) {
                line->SetLineWidth(1);
                line->SetLineColor(global_darkmode ? 0 : 1);
                line->Draw("same");
            }
        }
    }

    TH1D *sig_1D = new TH1D("sig_1D", ";Signal;Area", 80, 0.0, max_array[0]);
    TH1D *bg_1D = new TH1D("bg_1D", ";Background;Area", 20, 0.0, max_array[1]);

    for (int i = 0; i < 2; ++i) {
        TH1D *hist = (i == 0) ? sig_1D : bg_1D;
        tree->Draw(
            (i == 0) ? "sig >> sig_1D" : "bg >> bg_1D",
            "", "goff"
        );
        c1->cd(i + 3);
        hist->SetFillStyle(0);
        hist->SetLineWidth(2);
        hist->Draw();
        MyUtil::PaintBins(hist, 0.0, max_array[i]); // 各ビンをカラーパレットの色で塗る

        TLatex title;
        title.SetTextAlign(22);
        title.SetTextSize(0.06);
        title.SetTextColor(global_darkmode ? 0 : 1);
        title.DrawLatexNDC((i == 0) ? 0.5 : 0.43, 0.95, (i == 0) ? "Signal" : "Background");

        TLegend *lg = new TLegend(
            (i == 0) ? 0.73 : 0.665,
            0.7,
            (i == 0) ? 0.95 : 0.885,
            0.9
        );
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.04);
        lg->SetTextColor(global_darkmode ? 0 : 1);
        lg->AddEntry(hist, Form("%d areas", entries[i]), "");
        lg->AddEntry(hist, Form("Mean      %.2f", hist->GetMean()), "");
        lg->AddEntry(hist, Form("Std Dev   %.2f", hist->GetStdDev()), "");
        lg->Draw();
    }
}

void dist_shr_map(TCanvas *c1, TTree *tree[2], const double dist_max,
    const std::vector<double> shr_range, const double sig_max, int32_t *entries) noexcept
{
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    c1->Divide(2, 2);

    auto createAndDrawHist = [&](int pad, const char* name, const char* main_title,
                                 const char* title, const char* drawExpr,
                                 double xMin, double xMax, double yMin, double yMax, TTree* tree) {
        c1->cd(pad);
        gPad->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        gPad->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);

        TH2D* hist = new TH2D(name, title, 100, xMin, xMax, 100, yMin, yMax);
        tree->Draw(Form("%s >> %s", drawExpr, name), "", "goff");
        hist->Draw("colz");

        TLatex hist_title;
        hist_title.SetTextAlign(22);
        hist_title.SetTextSize(0.06);
        hist_title.SetTextColor(global_darkmode ? 0 : 1);
        hist_title.DrawLatexNDC((pad % 2 == 0) ? 0.43 : 0.5, 0.95, main_title);
    };

    createAndDrawHist(1, "dist1", "|Distortion 1| vs Signal", ";Signal;|Distortion 1|;Area",
                      "dist:sig", 0.0, sig_max, 0.0, dist_max, tree[0]);
    createAndDrawHist(2, "dist2", "|Distortion 2| vs Signal", ";Signal;|Distortion 2|;Area",
                      "dist:sig", 0.0, sig_max, 0.0, dist_max, tree[1]);
    createAndDrawHist(3, "shr1", "Shrink 1 vs Signal", ";Signal;Shrink 1;Area",
                      "shr:sig", 0.0, sig_max, shr_range[0], shr_range[1], tree[0]);
    createAndDrawHist(4, "shr2", "Shrink 2 vs Signal", ";Signal;Shrink 2;Area",
                      "shr:sig", 0.0, sig_max, shr_range[0], shr_range[1], tree[1]);
}
