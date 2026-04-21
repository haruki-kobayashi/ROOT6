/*
    predictionファイル（CSV, TSV, またはバイナリ形式）から読み込んだ
    prediction/reconstructionデータに対して、プレートごとの検出効率を
    複数形式でプロットしてPDF出力するプログラム。

    処理フロー:
    1. コマンドライン引数を解析し、表示範囲・カット条件・描画設定を取得
    2. 入力フォーマットを判定（未指定時は拡張子から自動判定）
    3. predictionファイルを読み込み（PredRecord.hpp の FileReader を使用）
    4. 各プレートについて each_plate_efficiency を実行
       - 位置2D効率、角度2D効率、角度1D効率（Wilson区間エラーバー付き）を描画
       - PH/VPH分布を描画
       - 角度1D効率データ（EfficiencyRecord）を収集
    5. 全プレート比較として all_plate_map を描画
       - 角度可変ビン × Plate番号の TH2Poly 効率マップを描画
    6. 全プレート重ね描きとして all_plate_overlay を描画
       - 収集済み EfficiencyRecord を使って1D効率曲線を重ね描き
    7. すべてのページをPDFとして出力し、必要に応じてテキストも出力
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <csignal>
#include <cmath>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TError.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>

#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <argparse/argparse.hpp>

#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>
#include <MyHeader/PredRecord.hpp>
#include <MyHeader/ConfidenceInterval.hpp>

// ========================================
// Efficiency計算用構造体
// ========================================

struct EfficiencyRecord {
    int pl;
    double angle_min, angle_max;
    int pred, found;
    double eff, eff_err_low, eff_err_high;
};


// ========================================
// グローバル変数
// ========================================

namespace {
    // Ctrl+Cで終了したときの処理用にグローバル変数を定義
    TCanvas *global_c1 = nullptr;
    std::string global_output = "";

    // その他のグローバル変数
    bool global_darkmode = false;

    constexpr double kPi = 3.14159265358979323846;

    inline double ToDisplayAngle(double tan_theta, bool use_degree)
    {
        return use_degree ? std::atan(tan_theta) * 180.0 / kPi : tan_theta;
    }

    inline std::vector<double> ToDisplayAngleBins(const std::vector<double> &angle_bins, bool use_degree)
    {
        std::vector<double> out;
        out.reserve(angle_bins.size());
        for (double v : angle_bins) {
            out.push_back(ToDisplayAngle(v, use_degree));
        }
        return out;
    }

    inline std::string AngleAxisTitle(bool use_degree)
    {
        return use_degree ? "#it{#theta} [#circ]" : "tan#it{#theta}";
    }

    inline std::string AngleAxisTitleXY(bool use_degree, const char *axis)
    {
        return use_degree
            ? fmt::format("#it{{#theta}}_{{{}}} [#circ]", axis)
            : fmt::format("tan#it{{#theta}}_{{{}}}", axis);
    }
}

// ========================================
// Ctrl+Cで終了したときの処理
// ========================================

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

// ========================================
// コマンドライン引数解析関数
// ========================================

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[]) {
    parser.set_usage_max_line_width(80);
    parser.add_description("Tips: You can combine single-character arguments.\n"
        "      For example, \"-d -p kBird -i\" is equivalent to \"-dpi kBird\".");

    // 必須引数: predファイル
    parser.add_argument("input_pred")
        .help("Path to the pred file (CSV, TSV, or binary) to be processed.")
        .required();

    parser.add_group("Optional arguments");

    // オプション引数: 入力フォーマット
    parser.add_argument("-fmt", "--format")
        .help("Input file format: csv, tsv, or bin. [default: auto-detect from extension]")
        .default_value(std::string());

    // オプション引数: 出力ファイル名
    parser.add_argument("-o", "--output")
        .help("Output PDF file name. [default: <input_pred>.pdf]")
        .default_value(std::string());

    // オプション引数: 角度ビン
    parser.add_argument("-a", "--angle_bins")
        .help("Angle bins to divide, separated by spaces.\n"
            "e.g., '--angle_bins 0.1 0.3 0.5' creates bins [0,0.1), [0.1,0.3), [0.3,0.5).\n"
            "[default: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0\n"
            "1.2 1.4 1.6 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0]")
        .default_value(std::vector<double>({
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
            1.2, 1.4, 1.6, 1.8, 2.0,
            2.5, 3.0, 3.5, 4.0, 4.5, 5.0
        }))
        .nargs(1, 10000)
        .scan<'g', double>();

    // オプション引数: エラーバーの長さ
        parser.add_argument("--error_sigma")
        .help("Error sigma for Wilson score interval. [default: 1.0]")
        .default_value(1.0)
        .scan<'g', double>();

    parser.add_argument("-%", "--percent")
        .help("Display efficiency values as percentages (0-100%) "
            "instead of decimals. [default: false]")
        .flag();

    parser.add_argument("-deg", "--angle_degree")
        .help("Plot angle axes in degrees instead of tan#theta. [default: false]")
        .flag();

    parser.add_argument("-txt", "-text")
        .help("Output text file name for efficiency values. [default: none]")
        .default_value(std::string());

    // プロット範囲
    parser.add_group("Plot range options");

    parser.add_argument("-efflow", "-effmin", "--efficiency_lower_limit")
        .help("Lower limit of efficiency display range.")
        .default_value(0.8)
        .scan<'g', double>();

    parser.add_argument("-angmax", "--angle_max")
        .help("Maximum angle to plot.")
        .default_value(5.0)
        .scan<'g', double>();

    parser.add_argument("--pl_range")
        .help("Plate number range to plot.")
        .default_value(std::vector<int>({0, INT_MAX}))
        .nargs(2)
        .scan<'i', int>();

    parser.add_argument("--ph_range")
        .help("Pulse Height range to plot.")
        .default_value(std::vector<double>({9.5, 32.5}))
        .nargs(2)
        .scan<'g', double>();

    parser.add_argument("--vph_range")
        .help("Volume Pulse Height range to plot.")
        .default_value(std::vector<double>({0.0, 250.0}))
        .nargs(2)
        .scan<'g', double>();

    // カット機能
    parser.add_group("Cut options");

    parser.add_argument("-x", "--cut_x")
        .help("Cut X range (mm).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({-DBL_MAX, DBL_MAX}))
        .nargs(2)
        .scan<'g', double>();

    parser.add_argument("-y", "--cut_y")
        .help("Cut Y range (mm).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({-DBL_MAX, DBL_MAX}))
        .nargs(2)
        .scan<'g', double>();

    parser.add_argument("-ax", "--cut_ax")
        .help("Cut ax range (tan_theta).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({-DBL_MAX, DBL_MAX}))
        .nargs(2)
        .scan<'g', double>();

    parser.add_argument("-ay", "--cut_ay")
        .help("Cut ay range (tan_theta).\n"
            "First = minimum, Second = maximum.")
        .default_value(std::vector<double>({-DBL_MAX, DBL_MAX}))
        .nargs(2)
        .scan<'g', double>();

    parser.add_group("Other options");

    // オプション引数: フォント設定
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
        .default_value(20)
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
        std::cerr << parser;
        std::cerr << "\nError: " << err.what() << std::endl;
        std::exit(1);
    }
}

// ========================================
// プロット関数の前方宣言
// ========================================

std::vector<EfficiencyRecord> each_plate_efficiency(
    TCanvas *c1, const std::vector<pred::Record> &records, int pl,
    const std::vector<double> &angle_bins, double eff_lower_lim,
    double angle_max, double error_sigma, bool use_percent, bool use_degree,
    const std::vector<double> &ph_range, const std::vector<double> &vph_range,
    const std::vector<double> &cut_x, const std::vector<double> &cut_y,
    const std::vector<double> &cut_ax, const std::vector<double> &cut_ay
) noexcept;

void all_plate_map(
    TCanvas *c1, const std::vector<int> &plates,
    const std::vector<EfficiencyRecord> &all_efficiency_data,
    const std::vector<double> &angle_bins, double eff_lower_lim, bool use_percent, bool use_degree
) noexcept;

void all_plate_overlay(
    TCanvas *c1, const std::vector<int> &plates,
    const std::vector<EfficiencyRecord> &all_efficiency_data,
    const std::vector<double> &angle_bins, double eff_lower_lim, bool use_percent, bool use_degree
) noexcept;

// ========================================
// メイン関数
// ========================================

int main(int argc, char* argv[])
{
    // Ctrl+Cで終了したときの処理を設定
    std::signal(SIGINT, handleSIGINT);

    // 処理時間を計測
    TStopwatch sw;

    // argparseを使用して引数を解析
    std::cout << "\n Initializing..." << std::endl;
    argparse::ArgumentParser parser("PredEfficiencyPlot.exe", "1.0.0");
    parse_arguments(parser, argc, argv);

    // 引数を取得
    const auto predfile = parser.get<std::string>("input_pred");
    auto input_format = parser.get<std::string>("--format");
    const auto output_arg = parser.get<std::string>("--output");
    auto angle_bins = parser.get<std::vector<double>>("--angle_bins");
    const auto error_sigma = parser.get<double>("--error_sigma");
    const auto eff_lower_lim = parser.get<double>("--efficiency_lower_limit");
    auto angle_max = parser.get<double>("--angle_max");
    const auto pl_range = parser.get<std::vector<int>>("--pl_range");
    const auto ph_range = parser.get<std::vector<double>>("--ph_range");
    const auto vph_range = parser.get<std::vector<double>>("--vph_range");
    auto cutX = parser.get<std::vector<double>>("--cut_x");
    auto cutY = parser.get<std::vector<double>>("--cut_y");
    auto cutAx = parser.get<std::vector<double>>("--cut_ax");
    auto cutAy = parser.get<std::vector<double>>("--cut_ay");
    auto font_number = parser.get<int>("--font_number");
    auto hideGrid = parser.get<bool>("--hide_grid");
    auto use_percent = parser.get<bool>("--percent");
    auto use_degree = parser.get<bool>("--angle_degree");
    const auto text_output = parser.get<std::string>("-txt");
    global_darkmode = parser.get<bool>("--dark_mode");
    auto palette_arg = parser.get<std::string>("--palette");
    auto NContours = parser.get<int>("--contours");
    auto invertpalette = parser.get<bool>("--invert_palette");
    auto negatepalette = parser.get<bool>("--negate_palette");

    // angle_maxの絶対値を取り、角度カット範囲に適用
    angle_max = std::abs(angle_max);
    cutAx[0] = std::max(cutAx[0], -angle_max);
    cutAx[1] = std::min(cutAx[1], angle_max);
    cutAy[0] = std::max(cutAy[0], -angle_max);
    cutAy[1] = std::min(cutAy[1], angle_max);

    // angle_binsからangle_max以上のビンを削除
    while (!angle_bins.empty() && angle_bins.back() > angle_max) {
        angle_bins.pop_back();
    }

    // ビンの境界が単調増加であることを確認
    for (size_t i = 1; i < angle_bins.size(); ++i) {
        if (angle_bins[i] <= angle_bins[i-1]) {
            std::cerr << fmt::format(
                "Error: Angle bins must be in strictly increasing order.\n"
                "       bins[{}] = {:.2f} <= bins[{}] = {:.2f}\n",
                i, angle_bins[i], i-1, angle_bins[i-1]
            );
            std::exit(1);
        }
    }

    // 入力フォーマットの自動判定
    if (input_format.empty()) {
        std::cout << " Input format not specified. Auto-detecting from file extension.\n";
        std::string ext = "";
        if (predfile.size() > 4) {
            ext = predfile.substr(predfile.size() - 4);
        }
        if (ext == ".csv") {
            input_format = "csv";
        } else if (ext == ".tsv") {
            input_format = "tsv";
        } else if (ext == ".bin") {
            input_format = "bin";
        } else {
            std::cerr << "Error: Unable to auto-detect input format from file extension.\n";
            std::exit(1);
        }
    }

    // フォーマット文字列を小文字に統一
    for (auto& c : input_format) {
        c = std::tolower(static_cast<unsigned char>(c));
    }

    // フォーマットの妥当性チェック
    if (input_format != "csv" && input_format != "," &&
        input_format != "tsv" && input_format != "\\t" &&
        input_format != "bin" && input_format != "binary") {
        std::cerr << fmt::format(
            "Error: Invalid input format ({}) specified. Use 'csv', 'tsv', or 'bin'.\n",
            input_format
        );
        std::exit(1);
    }

    // 出力ファイル名の設定
    const std::string output = (output_arg.empty())
        ? (predfile + ".pdf")
        : ((output_arg.size() > 4 && output_arg.substr(output_arg.size() - 4) == ".pdf")
            ? output_arg
            : output_arg + ".pdf");

    global_output = output;

    // 引数を利用する変数を設定
    const int font_code = 10 * font_number + 2;

    // カラーパレットの設定
    std::variant<int, std::string> Palette;
    try {
        Palette = std::stoi(palette_arg);
    } catch (const std::invalid_argument&) {
        Palette = palette_arg;
    }
    std::visit([NContours](auto&& arg) {
        MyPalette::SetPalette(arg, NContours);
    }, Palette);
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
    gROOT->GetColor(94)->SetRGB(200./255., 200./255., 200./255.);
    gROOT->GetColor(95)->SetRGB(60./255., 60./255., 60./255.);
    gROOT->GetColor(96)->SetRGB(r1 * 0.9, g1 * 0.9, b1 * 0.9);
    gROOT->GetColor(97)->SetRGB(r3 * 0.9, g3 * 0.9, b3 * 0.9);

    // カラーの設定
    gStyle->SetHistFillColor(0);                           // ヒストグラム
    gStyle->SetHistLineColor(93);                          // ヒストグラムの枠
    gStyle->SetFuncColor(93);                              // グラフ
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
    gStyle->SetOptStat(0);                // 統計boxの表示をオフ
    gStyle->SetPadRightMargin(0.19);      // Pad右側のマージン
    gStyle->SetPadLeftMargin(0.12);       // Pad左側のマージン
    gStyle->SetPadTopMargin(0.1);         // Pad上側のマージン
    gStyle->SetPadBottomMargin(0.11);     // Pad下側のマージン
    gStyle->SetLabelOffset(0.008, "xyz"); // 軸ラベル（数値）と軸の距離
    gStyle->SetLabelSize(0.05, "xyz");    // 軸ラベルのサイズ
    gStyle->SetTitleSize(0.05, "xyz");    // 軸タイトルのサイズ
    gStyle->SetTitleOffset(0.9, "x");     // x軸タイトルのオフセット
    gStyle->SetTitleOffset(1.0, "y");     // y軸タイトルのオフセット
    gStyle->SetTitleOffset(1.25, "z");    // z軸タイトルのオフセット
    gStyle->SetTitleY(0.985);             // タイトルのY位置
    gStyle->SetTitleFontSize(0.07);       // タイトルのフォントサイズ

    // グリッドの設定
    gStyle->SetPadGridX(!hideGrid); // グリッドの表示
    gStyle->SetPadGridY(!hideGrid); // グリッドの表示
    gStyle->SetPadTickX(1);         // 上側x軸の目盛り表示
    gStyle->SetPadTickY(1);         // 右側y軸の目盛り表示

    // フォント設定
    gStyle->SetStatFont(font_code);
    gStyle->SetLabelFont(font_code, "xyz");
    gStyle->SetTitleFont(font_code, "xyz");
    gStyle->SetTitleFont(font_code, "");
    gStyle->SetTextFont(font_code);
    gStyle->SetLegendFont(font_code);

    // エラーメッセージ未満のROOTのメッセージを非表示に設定
    gErrorIgnoreLevel = kError;

    // predictionファイルの読み込み
    std::cout << " Reading prediction file...";
    pred::FileReader reader(input_format);
    std::vector<pred::Record> records;
    try {
        records = reader.ReadAll(predfile);
    } catch (const std::exception& e) {
        std::cerr << "\nError: Failed to read prediction file: " << e.what() << std::endl;
        std::exit(1);
    }
    if (records.empty()) {
        std::cerr << "\nError: No records found in prediction file.\n";
        std::exit(1);
    }

    std::cout << fmt::format(" Loaded {} records.\n", records.size());

    // プレートの一覧を取得
    std::set<int> plate_set;
    for (const auto& record : records) {
        if (record.pred.pl < pl_range[0] || record.pred.pl > pl_range[1]) {
            continue; // 範囲外のプレートはスキップ
        }
        plate_set.insert(record.pred.pl);
    }
    std::vector<int> plates(plate_set.begin(), plate_set.end());
    std::sort(plates.begin(), plates.end());

    std::cout << fmt::format(" Found {} plates\n", plates.size());

    // Canvas作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas *c1 = new TCanvas("c1");
    global_c1 = c1;
    c1->Print((output + "[").c_str());

    TDatime starttime;
    int year = starttime.GetYear();
    int month = starttime.GetMonth();
    int day = starttime.GetDay();
    int hour = starttime.GetHour();
    int minute = starttime.GetMinute();
    int second = starttime.GetSecond();

    std::string StartTime = fmt::format("{:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}",
        year, month, day, hour, minute, second);
    std::cout << " Plot start : " << StartTime << std::endl;

	// プログレスバーの初期化
	int page = 0;
    const int total = static_cast<int>(plates.size()) + 2;
	MyUtil::ShowProgress(page, total);

    // 各プレートごとにプロット
    std::vector<EfficiencyRecord> all_efficiency_data;
    for (int pl : plates) {
        auto eff_data = each_plate_efficiency(c1, records, pl, angle_bins, eff_lower_lim,
            angle_max, error_sigma, use_percent, use_degree, ph_range, vph_range, cutX, cutY, cutAx, cutAy);
        all_efficiency_data.insert(all_efficiency_data.end(), eff_data.begin(), eff_data.end());
        c1->Print(output.c_str());
        c1->Clear();
        MyUtil::ShowProgress(page, total);
    }

    // 全プレート比較 (TH2Poly, 可変角度ビン)
    all_plate_map(c1, plates, all_efficiency_data, angle_bins, eff_lower_lim, use_percent, use_degree);
    c1->Print(output.c_str());
    c1->Clear();
    MyUtil::ShowProgress(page, total);

    // 全プレートの1D角度分布を重ね描き
    all_plate_overlay(c1, plates, all_efficiency_data, angle_bins, eff_lower_lim, use_percent, use_degree);
    c1->Print(output.c_str());
    c1->Clear();
    MyUtil::ShowProgress(page, total);

    // PDFファイルを閉じる
    c1->Print((output + "]").c_str());
	if (page < total) page = total;
    MyUtil::ShowProgress(page, total);

    TDatime endtime;
    year = endtime.GetYear();
    month = endtime.GetMonth();
    day = endtime.GetDay();
    hour = endtime.GetHour();
    minute = endtime.GetMinute();
    second = endtime.GetSecond();

    std::string EndTime = fmt::format("{:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}",
        year, month, day, hour, minute, second);

    double elapsedtime = sw.CpuTime();
    std::cout << "\n Plot end   : " << EndTime << " - Elapsed " << elapsedtime << " [s] (CPU)" << std::endl;
    std::cout << " Output     : " << output << std::endl;

    // テキストファイル出力
    if (!text_output.empty()) {
        std::ofstream ofs(text_output);
        if (!ofs.is_open()) {
            std::cerr << "Error: Cannot open text output file: " << text_output << std::endl;
            std::exit(1);
        }
        ofs << "pl angle_min angle_max prediction found efficiency error_low error_high\n";

        // プレートとangle_minでソート
        std::sort(all_efficiency_data.begin(), all_efficiency_data.end(),
            [](const EfficiencyRecord& a, const EfficiencyRecord& b) {
                if (a.pl != b.pl) return a.pl < b.pl;
                return a.angle_min < b.angle_min;
            });

        for (const auto& er : all_efficiency_data) {
            ofs << fmt::format("{:>4}", er.pl) << " "
                << fmt::format("{:>6.2f}", er.angle_min) << " "
                << fmt::format("{:>6.2f}", er.angle_max) << " "
                << fmt::format("{:>10}", er.pred) << " "
                << fmt::format("{:>10}", er.found) << " "
                << fmt::format("{:>9.6f}", er.eff) << " "
                << fmt::format("{:>9.6f}", er.eff_err_low) << " "
                << fmt::format("{:>9.6f}", er.eff_err_high) << "\n";
        }
        ofs.close();
        std::cout << " Text file  : " << text_output << std::endl;
    }

    return 0;
}

// ========================================
// プロット関数の実装
// ========================================

std::vector<EfficiencyRecord> each_plate_efficiency(
    TCanvas *c1, const std::vector<pred::Record> &records, int pl,
    const std::vector<double> &angle_bins, double eff_lower_lim,
    double angle_max, double error_sigma, bool use_percent, bool use_degree,
    const std::vector<double> &ph_range, const std::vector<double> &vph_range,
    const std::vector<double> &cut_x, const std::vector<double> &cut_y,
    const std::vector<double> &cut_ax, const std::vector<double> &cut_ay
) noexcept
{
    // プレート内のレコードをフィルタリング
    std::vector<const pred::Record*> filtered_records;
    double x_min = DBL_MAX, x_max = -DBL_MAX;
    double y_min = DBL_MAX, y_max = -DBL_MAX;

    for (const auto& record : records) {
        if (record.pred.pl != pl) continue;

        // カット条件の適用
        if (record.pred.x * 0.001 < cut_x[0] || record.pred.x * 0.001 > cut_x[1]) continue;
        if (record.pred.y * 0.001 < cut_y[0] || record.pred.y * 0.001 > cut_y[1]) continue;
        if (record.pred.ax < cut_ax[0] || record.pred.ax > cut_ax[1]) continue;
        if (record.pred.ay < cut_ay[0] || record.pred.ay > cut_ay[1]) continue;

        filtered_records.push_back(&record);
        x_min = std::min(x_min, record.pred.x * 0.001);
        x_max = std::max(x_max, record.pred.x * 0.001);
        y_min = std::min(y_min, record.pred.y * 0.001);
        y_max = std::max(y_max, record.pred.y * 0.001);
    }

    if (filtered_records.empty()) {
        // データなし - 空のプレートホルダーを表示
        std::string hist_title = fmt::format("Plate {:03d} has no data", pl);
        TLatex *text = new TLatex(0.5, 0.5, hist_title.c_str());
        text->SetNDC();
        text->SetTextAlign(22);
        text->Draw();
        return {};
    }

    // Canvasを分割
    c1->Divide(2, 1);
    c1->GetPad(1)->Divide(1, 2);
    c1->GetPad(2)->Divide(1, 3);
    for (int i = 1; i <= 3; ++i) {
        c1->GetPad(2)->GetPad(i)->SetBottomMargin(0.16);
    }

    const std::vector<double> angle_bins_disp = ToDisplayAngleBins(angle_bins, use_degree);
    const double angle_max_disp = ToDisplayAngle(angle_max, use_degree);
    const std::string angle_axis = AngleAxisTitle(use_degree);
    const std::string angle_axis_x = AngleAxisTitleXY(use_degree, "x");
    const std::string angle_axis_y = AngleAxisTitleXY(use_degree, "y");

    // プロット範囲とビン幅の決定
    // フィルムの長辺の端から1cm外側までを最大表示範囲とし、縦横比を正しく保って表示
    // ビン幅は1mm
    double LowX, UpX, LowY, UpY;
    int bin;
    constexpr double bin_width = 1.0;
    const double RangeX = x_max - x_min;
    const double RangeY = y_max - y_min;
    if (RangeX >= RangeY) {
        LowX = x_min - 10.0;
        UpX = x_max + 10.0;
        LowY = y_min - (RangeX - RangeY + 20.0) * 0.5;
        UpY = y_max + (RangeX - RangeY + 20.0) * 0.5;
        bin = static_cast<int>((RangeX + 20.0) * bin_width);
    } else {
        LowX = x_min - (RangeY - RangeX + 20.0) * 0.5;
        UpX = x_max + (RangeY - RangeX + 20.0) * 0.5;
        LowY = y_min - 10.0;
        UpY = y_max + 10.0;
        bin = static_cast<int>((RangeY + 20.0) * bin_width);
    }

    // ========================================
    // 1. efficiency 2D位置分布
    // ========================================
    c1->GetPad(1)->cd(1);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(0.2);

    gStyle->SetTitleOffset(0.8, "x");
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetTitleOffset(1.1, "z");
    gStyle->SetTitleY(1.025);
    gStyle->SetTitleSize(0.06, "xyz");
    gStyle->SetTitleFontSize(0.11);

    TH2D *pos0 = new TH2D("pos0", "all_prediction", bin, LowX, UpX, bin, LowY, UpY);
    TH2D *pos1 = new TH2D("pos1", "found_prediction", bin, LowX, UpX, bin, LowY, UpY);

    for (const auto* rec : filtered_records) {
        pos0->Fill(rec->pred.x * 0.001, rec->pred.y * 0.001);
        if (pred::found(*rec)) {
            pos1->Fill(rec->pred.x * 0.001, rec->pred.y * 0.001);
        }
    }

    std::string z_title = use_percent ? "Efficiency (%)" : "Efficiency";
    std::string postitle = fmt::format("PL{:03d};x [mm];y [mm];{}", pl, z_title);
    TH2D *pos = new TH2D("pos", postitle.c_str(), bin, LowX, UpX, bin, LowY, UpY);
    pos->Divide(pos1, pos0);
    if (use_percent) pos->Scale(100.0);
    pos->Draw("colz");
    double eff_lower = use_percent ? eff_lower_lim * 100.0 : eff_lower_lim;
    double eff_upper = use_percent ? 100.0 : 1.0;
    pos->GetZaxis()->SetRangeUser(eff_lower, eff_upper);

    // ========================================
    // 2. efficiency 2D角度分布
    // ========================================
    c1->GetPad(1)->cd(2);
    gPad->SetRightMargin(0.25);
    gPad->SetLeftMargin(0.2);

    TH2D *ang0 = new TH2D("ang0", "all_prediction",
        100, -angle_max_disp, angle_max_disp, 100, -angle_max_disp, angle_max_disp);
    TH2D *ang1 = new TH2D("ang1", "found_prediction",
        100, -angle_max_disp, angle_max_disp, 100, -angle_max_disp, angle_max_disp);

    for (const auto* rec : filtered_records) {
        ang0->Fill(ToDisplayAngle(rec->pred.ax, use_degree), ToDisplayAngle(rec->pred.ay, use_degree));
        if (pred::found(*rec)) {
            ang1->Fill(ToDisplayAngle(rec->pred.ax, use_degree), ToDisplayAngle(rec->pred.ay, use_degree));
        }
    }

    std::string ang_title = fmt::format(";{};{};{}", angle_axis_x, angle_axis_y, z_title);
    TH2D *ang = new TH2D("ang", ang_title.c_str(),
        100, -angle_max_disp, angle_max_disp, 100, -angle_max_disp, angle_max_disp);
    ang->Divide(ang1, ang0);
    if (use_percent) ang->Scale(100.0);
    ang->Draw("colz");
    ang->GetXaxis()->SetNdivisions(212);
    ang->GetYaxis()->SetNdivisions(212);
    ang->SetTitleOffset(0.8, "y");
    ang->GetZaxis()->SetRangeUser(eff_lower, eff_upper);

    // ========================================
    // 3. efficiency 1D角度分布
    // ========================================
    c1->GetPad(2)->cd(1);

    const int bin_num = static_cast<int>(angle_bins_disp.size());
    std::vector<double> xbins_vec = {0.0};
    xbins_vec.insert(xbins_vec.end(), angle_bins_disp.begin(), angle_bins_disp.end());

    double* xbins = xbins_vec.data();

    TH1D *eff0 = new TH1D("eff0", "all_prediction", bin_num, xbins);
    TH1D *eff1 = new TH1D("eff1", "found_prediction", bin_num, xbins);

    for (const auto* rec : filtered_records) {
        double tan_theta = ToDisplayAngle(std::hypot(rec->pred.ax, rec->pred.ay), use_degree);
        eff0->Fill(tan_theta);
        if (pred::found(*rec)) {
            eff1->Fill(tan_theta);
        }
    }

    TH1D *eff = new TH1D("eff", fmt::format(";{};{}", angle_axis, z_title).c_str(),
                         bin_num, xbins);
    eff->Divide(eff1, eff0);
    if (use_percent) eff->Scale(100.0);
    eff->Draw();
    eff->SetLineWidth(0);
    eff->GetYaxis()->SetNdivisions(505);
    eff->SetTitleOffset(0.8, "x");
    eff->SetTitleOffset(0.7, "y");
    eff->SetTitleSize(0.09, "x");
    eff->SetTitleSize(0.09, "y");
    eff->SetLabelSize(0.08, "xy");
    eff->GetYaxis()->SetRangeUser(eff_lower, eff_upper);

    // エラーバーの追加
    TGraphAsymmErrors *eff_err = new TGraphAsymmErrors();
    for (int i = 1; i <= bin_num; ++i) {
        int n_all = static_cast<int>(eff0->GetBinContent(i));
        int n_found = static_cast<int>(eff1->GetBinContent(i));
        if (n_all > 0) {
            double efficiency = static_cast<double>(n_found) / n_all;
            auto [ci_low, ci_high] = Confidence::wilson_interval(n_all, n_found, error_sigma);
            double err_low = efficiency - ci_low;
            double err_high = ci_high - efficiency;
            if (use_percent) {
                efficiency *= 100.0;
                err_low *= 100.0;
                err_high *= 100.0;
            }

            double bin_center = eff->GetBinCenter(i);
            double bin_half_width = eff->GetBinWidth(i) / 2.0;
            eff_err->SetPoint(i - 1, bin_center, efficiency);
            eff_err->SetPointError(i - 1, bin_half_width, bin_half_width, err_low, err_high);
        }
    }
    eff_err->Draw("same E Z 0");
    eff_err->SetLineWidth(2);
    eff_err->SetLineColor(93);

    TLegend *description = new TLegend(0.15, 0.94, 0.55, 0.97);
    description->SetFillStyle(0);
    description->SetBorderSize(0);
    description->SetTextSize(0.08);
    description->AddEntry(
        eff,
        fmt::format("Error bars: Wilson score interval ({}#sigma)", error_sigma).c_str(),
        "");
    description->Draw();

    // efficiencyデータの収集
    std::vector<EfficiencyRecord> efficiency_data;
    for (int i = 1; i <= bin_num; ++i) {
        int n_all = static_cast<int>(eff0->GetBinContent(i));
        int n_found = static_cast<int>(eff1->GetBinContent(i));
        if (n_all > 0) {
            double efficiency = static_cast<double>(n_found) / n_all;
            auto [ci_low, ci_high] = Confidence::wilson_interval(n_all, n_found, error_sigma);
            double err_low = efficiency - ci_low;
            double err_high = ci_high - efficiency;

            EfficiencyRecord er;
            er.pl = pl;
            er.angle_min = (i == 1) ? 0.0 : angle_bins_disp[i - 2];
            er.angle_max = angle_bins_disp[i - 1];
            er.pred = n_all;
            er.found = n_found;
            er.eff = efficiency;
            er.eff_err_low = err_low;
            er.eff_err_high = err_high;
            efficiency_data.push_back(er);
        }
    }

    // ========================================
    // 4. 検出された飛跡のPH角度分布
    // ========================================
    c1->GetPad(2)->cd(2);

    const int ph_bins = static_cast<int>(ph_range[1] - ph_range[0]);
    TH2D *ph = new TH2D(
        "ph", fmt::format(";{};PH (detected tracks)", angle_axis).c_str(),
        100, 0.0, angle_max_disp, ph_bins, ph_range[0], ph_range[1]
    );

    for (const auto* rec : filtered_records) {
        for (const auto& reco : rec->recos) {
            ph->Fill(ToDisplayAngle(std::hypot(reco.ax, reco.ay), use_degree), reco.ph * 0.0001);
        }
    }
    ph->Draw("colz");
    ph->SetTitleOffset(0.8, "x");
    ph->SetTitleOffset(0.6, "y");
    ph->SetTitleOffset(1.0, "z");
    ph->SetTitleSize(0.09, "x");
    ph->SetTitleSize(0.09, "y");
    ph->SetTitleSize(0.065, "z");
    ph->SetLabelSize(0.08, "xyz");

    // ========================================
    // 5. 検出された飛跡のVPH角度分布
    // ========================================
    c1->GetPad(2)->cd(3);

    const int vph_bins = 50;
    TH2D *vph = new TH2D(
        "vph", fmt::format(";{};VPH (detected tracks)", angle_axis).c_str(),
        100, 0.0, angle_max_disp, vph_bins, vph_range[0], vph_range[1]
    );

    for (const auto* rec : filtered_records) {
        for (const auto& reco : rec->recos) {
            vph->Fill(ToDisplayAngle(std::hypot(reco.ax, reco.ay), use_degree), reco.ph % 10000);
        }
    }
    vph->Draw("colz");
    vph->SetTitleOffset(0.8, "x");
    vph->SetTitleOffset(0.6, "y");
    vph->SetTitleOffset(1.0, "z");
    vph->SetTitleSize(0.09, "x");
    vph->SetTitleSize(0.09, "y");
    vph->SetTitleSize(0.065, "z");
    vph->SetLabelSize(0.08, "xyz");

    return efficiency_data;
}

void all_plate_map(
    TCanvas *c1, const std::vector<int> &plates,
    const std::vector<EfficiencyRecord> &all_efficiency_data,
    const std::vector<double> &angle_bins, double eff_lower_lim, bool use_percent, bool use_degree
) noexcept
{
    if (plates.empty() || angle_bins.empty()) {
        c1->cd();
        TLatex *text = new TLatex(0.5, 0.5, "No all-plate data to plot");
        text->SetNDC();
        text->SetTextAlign(22);
        text->Draw();
        return;
    }

    c1->cd();

    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetTitleOffset(0.9, "x");
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetTitleOffset(1.25, "z");
    gStyle->SetTitleY(0.985);
    gStyle->SetTitleFontSize(0.07);

    const int min_pl = plates.front();
    const int max_pl = plates.back();
    const int pl_num = max_pl - min_pl + 1;
    const auto plate_to_display_y = [min_pl, max_pl](int pl) {
        return static_cast<double>(min_pl + max_pl - pl);
    };
    int division = pl_num + 1;
    if (pl_num > 20) division = 220;

    const double z_lower = use_percent ? eff_lower_lim * 100.0 : eff_lower_lim;
    const double z_upper = use_percent ? 100.0 : 1.0;
    const std::vector<double> angle_bins_disp = ToDisplayAngleBins(angle_bins, use_degree);
    const std::string angle_axis = AngleAxisTitle(use_degree);

    std::string z_title = use_percent ? "Efficiency (%)" : "Efficiency";
    std::string all_title = fmt::format(
        "All plate efficiency map;tan#it{{#theta}};Plate number;{}", z_title
    );

    TH2Poly *allPL = new TH2Poly();
    allPL->SetName("allPL");
    allPL->SetTitle(all_title.c_str());
    allPL->SetTitle(fmt::format("All plate efficiency map;{};Plate number;{}", angle_axis, z_title).c_str());

    std::map<std::pair<int, int>, double> eff_map;
    for (const auto& er : all_efficiency_data) {
        int angle_index = -1;
        for (int i = 0; i < static_cast<int>(angle_bins_disp.size()); ++i) {
            if (std::abs(er.angle_max - angle_bins_disp[i]) < 1e-9) {
                angle_index = i;
                break;
            }
        }
        if (angle_index < 0) continue;

        double value = use_percent ? er.eff * 100.0 : er.eff;
        eff_map[{er.pl, angle_index}] = value;
    }

    double x_low = 0.0;
    for (int ai = 0; ai < static_cast<int>(angle_bins_disp.size()); ++ai) {
        double x_up = angle_bins_disp[ai];
        for (int pl : plates) {
            const double y_center = plate_to_display_y(pl);
            const int bin_idx = allPL->AddBin(x_low, y_center - 0.5, x_up, y_center + 0.5);
            auto it = eff_map.find({pl, ai});
            if (it != eff_map.end()) {
                allPL->SetBinContent(bin_idx, it->second);
            }
        }
        x_low = x_up;
    }

    allPL->SetMinimum(z_lower);
    allPL->SetMaximum(z_upper);

    // 描画枠の設定：TH2Polyのみでmin_pl-0.5からmax_pl+0.5の範囲設定が難しいため
    TH1F *frame = gPad->DrawFrame(
        0.0, min_pl - 0.5, angle_bins_disp.back(), max_pl + 0.5,
        all_title.c_str()
    );
    frame->SetLabelSize(0.04, "x");
    frame->SetLabelSize(0.045, "z");

    // Y軸の主目盛りをビン境界(半整数)に合わせて、グリッド位置を固定
    frame->GetYaxis()->SetNdivisions(pl_num, false);

    // frame側のYラベルは非表示（境界グリッド専用の軸にする）
    frame->GetYaxis()->SetLabelSize(0.0);

    allPL->Draw("colz same");
    gPad->RedrawAxis("g"); // グリッドを前面に出す

    // 整数ラベル用の補助Y軸（ビン中心: min_pl..max_pl）
    const double x_axis = 0.0;
    TGaxis *yaxis_center = new TGaxis(
        x_axis, max_pl,
        x_axis, min_pl,
        min_pl, max_pl,
        std::max(pl_num - 1, 1),
        "S"
    );
    yaxis_center->SetNdivisions(division);
    yaxis_center->SetLabelSize(0.045);
    yaxis_center->SetLabelOffset(-0.008);
    yaxis_center->SetLabelFont(gStyle->GetLabelFont("y"));
    yaxis_center->SetLabelColor(gStyle->GetLabelColor("y"));
    yaxis_center->SetLineColor(gStyle->GetAxisColor("y"));
    yaxis_center->SetTickSize(0.0);
    yaxis_center->Draw();

    frame->Draw("same axis"); // 軸を前面に出す
}

void all_plate_overlay(
    TCanvas *c1, const std::vector<int> &plates,
    const std::vector<EfficiencyRecord> &all_efficiency_data,
    const std::vector<double> &angle_bins, double eff_lower_lim, bool use_percent, bool use_degree
) noexcept
{
    if (plates.empty() || angle_bins.empty()) {
        c1->cd();
        TLatex *text = new TLatex(0.5, 0.5, "No all-plate data to plot");
        text->SetNDC();
        text->SetTextAlign(22);
        text->Draw();
        return;
    }

    c1->cd();

    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetTitleOffset(0.9, "x");
    gStyle->SetTitleOffset(1.0, "y");
    gStyle->SetTitleY(0.985);
    gStyle->SetTitleFontSize(0.07);

    const int min_pl = plates.front();
    const int max_pl = plates.back();
    gStyle->SetNumberContours(max_pl - min_pl + 1);

    const double y_lower = use_percent ? eff_lower_lim * 100.0 : eff_lower_lim;
    const double y_upper = use_percent ? 100.0 : 1.0;
    const std::string y_title = use_percent ? "Efficiency (%)" : "Efficiency";
    const std::vector<double> angle_bins_disp = ToDisplayAngleBins(angle_bins, use_degree);
    const std::string angle_axis = AngleAxisTitle(use_degree);

    std::vector<double> xbins_vec = {0.0};
    xbins_vec.insert(xbins_vec.end(), angle_bins_disp.begin(), angle_bins_disp.end());
    const double *xbins = xbins_vec.data();
    const int bin_num = static_cast<int>(angle_bins_disp.size());

    TH2D *axis_hist = new TH2D(
        "all_plate_overlay_axis",
        fmt::format("All plate efficiency;{};{};Plate number", angle_axis, y_title).c_str(),
        bin_num, xbins, 1, y_lower, y_upper
    );
    axis_hist->Fill(0.0, 0.0); // ダミーデータで軸を描画
    axis_hist->SetLabelSize(0.04, "x");
    axis_hist->SetLabelSize(0.045, "y");
    axis_hist->GetYaxis()->SetRangeUser(y_lower, y_upper);
    axis_hist->GetZaxis()->SetRangeUser(min_pl - 0.5, max_pl + 0.5);
    axis_hist->Draw("colz");

    const int n_colors = std::max(gStyle->GetNumberOfColors(), 1);
    std::map<int, std::map<int, const EfficiencyRecord*>> plate_bins;
    for (const auto& er : all_efficiency_data) {
        int angle_index = -1;
        for (int i = 0; i < static_cast<int>(angle_bins_disp.size()); ++i) {
            if (std::abs(er.angle_max - angle_bins_disp[i]) < 1e-9) {
                angle_index = i;
                break;
            }
        }
        if (angle_index < 0) continue;
        plate_bins[er.pl][angle_index] = &er;
    }

    for (int pl : plates) {
        TH1D *hist = new TH1D(Form("allPL_overlay_eff_%03d", pl), "", bin_num, xbins);

        const int color_index = (max_pl > min_pl)
            ? static_cast<int>(std::lround((n_colors - 1.0) * (pl - min_pl) / static_cast<double>(max_pl - min_pl)))
            : 0;
        const int color = gStyle->GetColorPalette(std::clamp(color_index, 0, n_colors - 1));

        TGraphAsymmErrors *eff_err = new TGraphAsymmErrors();
        auto it = plate_bins.find(pl);
        if (it != plate_bins.end()) {
            for (const auto& [angle_index, er] : it->second) {
                const int bin_index = angle_index + 1;
                const double efficiency = use_percent ? er->eff * 100.0 : er->eff;
                const double err_low = use_percent ? er->eff_err_low * 100.0 : er->eff_err_low;
                const double err_high = use_percent ? er->eff_err_high * 100.0 : er->eff_err_high;

                hist->SetBinContent(bin_index, efficiency);

                const double bin_center = hist->GetBinCenter(bin_index);
                const double bin_half_width = hist->GetBinWidth(bin_index) / 2.0;
                eff_err->SetPoint(bin_index - 1, bin_center, efficiency);
                eff_err->SetPointError(bin_index - 1, bin_half_width, bin_half_width, err_low, err_high);
            }
        }

        eff_err->SetLineWidth(2);
        eff_err->SetLineColorAlpha(color, 0.6);
        eff_err->Draw("same E Z 0");
    }
}
