#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <set>
#include <variant>
#include <array>
#include <csignal>

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
#include <TGraph.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TCut.h>
#include <TLegend.h>
#include <TColor.h>
#include <TLeaf.h>

#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>
#include <ROOT6/SensorArray.hpp>
#include <ROOT6/StderrSuppressor.hpp>
#include "main.hpp"

// Ctrl+Cで終了したときの処理用にグローバル変数を定義
TCanvas* global_c1 = nullptr;
std::string global_output = "";

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
    // Ctrl+Cで終了したときの処理
    std::signal(SIGINT, handleSIGINT);

    // 引数の数を確認
    if (argc < 2) {
        std::cerr << "\n Usage: " << argv[0] << " bvxx_file [opsion]\n" << std::endl;
        return 1;
    }

    // コマンドライン引数を解析
    const std::string bvxxfile = argv[1];

    // 一通り完成したらargparseで受け取れるように変更する
    std::variant<int, std::string, MyPalette::Palette, EColorPalette> Palette;
    if (argc > 2) {
        try {
            Palette = std::stoi(argv[2]); // argv[2]を整数として解釈
        } catch (const std::invalid_argument&) {
            Palette = argv[2]; // 認識できない場合は文字列として扱う
        }
    } else {
        Palette = kBird; // デフォルトのパレットはkBird
    }
    const std::string outputfile = (argc > 3) ? argv[3] : bvxxfile;
    const std::string output = (outputfile.size() > 4 &&
                                outputfile.substr(outputfile.size() - 4) == ".pdf")
                                ? outputfile
                                : outputfile + ".pdf";
    int pl = -1;
    const double TD_range[2] = {0.0, 0.0};
    const double angle_max = 6.0;
    const double angle_resolution = 0.1;
    const double da_cut_slope = 0.08;
    const double da_cut_intercept = 0.02;
    const uint8_t da_cutPH = 9;
    const double dlat_range = 0.05;
    const double drad_range = 1.0;
    const uint32_t vph_range = 100;
    const int32_t ranking_range_min = 30;
    const bool invertpalette = false;
    const int font_number = 4;
    const bool showGrid = true;
    int vph_standard = -1;
    const RankingParams ranking_params_array[][3] = {
        {
            {0.0, 0.1, 120, 0.10, 0.05}, 
            {0.1, 0.2,  70, 0.10, 0.05}, 
            {0.2, 0.3,  45, 0.10, 0.05}
        },
        {
            {0.3, 0.4, 25, 0.14, 0.05}, 
            {0.4, 0.5, 20, 0.15, 0.05}, 
            {0.5, 0.6, 15, 0.15, 0.05}
        },
        {
            {0.6, 0.7, 10, 0.20, 0.05}, 
            {0.7, 0.8,  0, 0.20, 0.05}, 
            {0.8, 0.9,  0, 0.20, 0.05}
        },
        {
            {0.9, 1.0, 0, 0.25, 0.05}, 
            {1.0, 1.1, 0, 0.30, 0.05}, 
            {1.1, 1.3, 0, 0.35, 0.05}
        },
        {
            {1.3, 1.5, 0, 0.35, 0.05}, 
            {1.5, 2.0, 0, 0.50, 0.05}, 
            {2.0, 2.5, 0, 0.60, 0.05}
        },
        {
            {2.5, 3.0, 0, 0.70, 0.05}, 
            {3.0, 4.0, 0, 0.90, 0.05}, 
            {4.0, 5.0, 0, 1.00, 0.05}
        }
    };

    const int font_code = 10 * font_number + 2;
    const TString da_cutX = Form(
        "(%f * (ax < 0 ? -ax : ax) + %f)",
        da_cut_slope, da_cut_intercept
    );
    const TString da_cutY = Form(
        "(%f * (ay < 0 ? -ay : ay) + %f)",
        da_cut_slope, da_cut_intercept
    );

    // 処理時間を計測
    TStopwatch sw;

    // エラーメッセージ未満のROOTのメッセージを非表示に設定
    gErrorIgnoreLevel = kError;

    // TTreeの作成
    std::cout << "\nReading bvxx file: ";
    TTree* tree = new TTree("tree", "");
    TTree* subtree = new TTree("subtree", "");

    // TTreeのブランチを作成
    uint8_t ph1, ph2;
    uint32_t ShotID1, ViewID1, ImagerID1, ShotID2, ViewID2, ImagerID2, vph1, vph2;
    double x, y, ax, ay, ax1, ay1, ax2, ay2, dax1, day1, dax2, day2, dx, dy, dz, tan, lin, linl;
    // tree
    const std::array<std::pair<const char*, void*>, 2> uint8Branches = {{
        {"ph1", &ph1}, {"ph2", &ph2}
    }};
    for (const auto& branch : uint8Branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/B").c_str());
    }
    const std::array<std::pair<const char*, void*>, 4> uint32Branches = {{
        {"ImagerID1", &ImagerID1}, {"ImagerID2", &ImagerID2}, 
        {"vph1", &vph1}, {"vph2", &vph2}
    }};
    for (const auto& branch : uint32Branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/I").c_str());
    }
    const std::array<std::pair<const char*, void*>, 15> doubleBranches = {{
        {"x", &x}, {"y", &y}, {"ax", &ax}, {"ay", &ay}, 
        {"ax1", &ax1}, {"ay1", &ay1}, {"ax2", &ax2}, {"ay2", &ay2}, 
        {"dax1", &dax1}, {"day1", &day1}, {"dax2", &dax2}, {"day2", &day2}, 
        {"tan", &tan}, {"lin", &lin}, {"linl", &linl}
    }};
    for (const auto& branch : doubleBranches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }
    // subtree
    const std::array<std::pair<const char*, void*>, 5> doubleBranchesSub = {{
        {"x", &x}, {"y", &y}, {"dx", &dx}, {"dy", &dy}, {"dz", &dz}
    }};
    for (const auto& branch : doubleBranchesSub) {
        subtree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    // bvxxファイルの読み込み
    vxx::BvxxReader br;
    std::set<uint32_t> uniqueViewID2, uniqueViewID1;
    constexpr uint8_t NumberOfImager = 72;
    uint32_t NoTcount[2][72] = {0};
    uint32_t NoTcount_same[72] = {0};
    double sensor_dr_sum[72] = {0.0};
    double sensor_dz_sum[72] = {0.0};

    StderrSuppressor sup; // 標準エラー出力を抑制するためのオブジェクト
    bool readSuccess = false;
    if (pl == -1) { // PL番号が指定されていない場合
        // PL番号を0から999まで試す。正しいPL番号が見つかるまでのエラー出力をStderrSuppressorで抑制する
        sup.suppress(); // 標準エラー出力の抑制開始
        for (int pli = 0; pli < 1000; ++pli) {
            if (br.Begin(bvxxfile, pli, 0)) {
                sup.restore(); // 見つかったら標準エラー出力の抑制を解除
                std::cout << "PL" << Form("%03d", pli) << " Loading..." << std::endl;
                pl = pli;

                vxx::HashEntry h;
                vxx::base_track_t b;

                while (br.NextHashEntry(h)) {
                    while (br.NextBaseTrack(b)) {
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
                        dx = (ax - 0.5 * (ax1 + ax2)) * dz; // microtrack1, 2をベース中央に内挿したときの位置ずれ
                        dy = (ay - 0.5 * (ay1 + ay2)) * dz; // microtrack1, 2をベース中央に内挿したときの位置ずれ
                        tan = sqrt(ax * ax + ay * ay);
                        lin = sqrt(dax1 * dax1 + day1 * day1 + dax2 * dax2 + day2 * day2); // 飛跡の直線性
                        linl = sqrt(
                            ((ax * ay1 - ay * ax1) / tan) * ((ax * ay1 - ay * ax1) / tan) +
                            ((ax * ay2 - ay * ax2) / tan) * ((ax * ay2 - ay * ax2) / tan)
                        ); // 飛跡のlateral方向の直線性

                        uniqueViewID2.insert(ViewID2); // 視野数のカウント
                        uniqueViewID1.insert(ViewID1); // 視野数のカウント
                        NoTcount[0][ImagerID2]++; // Imagerごとの飛跡本数カウント。NoTcount[i]とLayeriが対応
                        NoTcount[1][ImagerID1]++; // Imagerごとの飛跡本数カウント。NoTcount[i]とLayeriが対応

                        // 両面のmicro trackが同じImagerで検出された場合
                        if (ImagerID1 == ImagerID2) {
                            NoTcount_same[ImagerID1]++; // Imagerごとの飛跡本数カウント
                            sensor_dr_sum[ImagerID2] += sqrt(dx*dx + dy*dy); // Imagerごとの平面方向の位置ずれ
                            sensor_dz_sum[ImagerID2] += dz; // Imagerごとのz方向の位置ずれ
                        }

                        tree->Fill();

                        if ((ph1+ph2)>24 && tan>1.0 && tan<1.1) subtree->Fill();
                    }
                }
                br.End();
                readSuccess = true;
                break; // 成功したらループを抜ける
            }
        }
        if (!readSuccess) {
            std::cerr << "Error: Failed to load bvxx file. Please try using the -pl option." << std::endl;
            return 1;
        }
    } else { // PL番号が指定されている場合
        if (br.Begin(bvxxfile, pl, 0)) {
            std::cout << "PL" << Form("%03d", pl) << " Loading..." << std::endl;

            vxx::HashEntry h;
            vxx::base_track_t b;

            while (br.NextHashEntry(h)) {
                while (br.NextBaseTrack(b)) {
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
                    dx = (ax - 0.5 * (ax1 + ax2)) * dz; // microtrack1, 2をベース中央に内挿したときの位置ずれ
                    dy = (ay - 0.5 * (ay1 + ay2)) * dz; // microtrack1, 2をベース中央に内挿したときの位置ずれ
                    tan = sqrt(ax * ax + ay * ay);
                    lin = sqrt(dax1 * dax1 + day1 * day1 + dax2 * dax2 + day2 * day2); // 飛跡の直線性
                    linl = sqrt(
                        ((ax * ay1 - ay * ax1) / tan) * ((ax * ay1 - ay * ax1) / tan) +
                        ((ax * ay2 - ay * ax2) / tan) * ((ax * ay2 - ay * ax2) / tan)
                    ); // 飛跡のlateral方向の直線性

                    uniqueViewID2.insert(ViewID2); // 視野数のカウント
                    uniqueViewID1.insert(ViewID1); // 視野数のカウント
                    NoTcount[0][ImagerID2]++; // Imagerごとの飛跡本数カウント。NoTcount[i]とLayeriが対応
                    NoTcount[1][ImagerID1]++; // Imagerごとの飛跡本数カウント。NoTcount[i]とLayeriが対応

                    // 両面のmicro trackが同じImagerで検出された場合
                    if (ImagerID1 == ImagerID2) {
                        NoTcount_same[ImagerID1]++; // Imagerごとの飛跡本数カウント
                        sensor_dr_sum[ImagerID2] += sqrt(dx*dx + dy*dy); // Imagerごとの平面方向の位置ずれ
                        sensor_dz_sum[ImagerID2] += dz; // Imagerごとのz方向の位置ずれ
                    }

                    tree->Fill();

                    if ((ph1+ph2)>24 && tan>1.0 && tan<1.1) subtree->Fill();
                }
            }
            br.End();
        } else {
            std::cerr << "Error: Failed to load bvxx file. Please check if the correct -pl is specified." << std::endl;
            return 1;
        }
    } // bvxxファイルの読み込み終了

    int fieldsOfView[2] = {uniqueViewID2.size(), uniqueViewID1.size()}; // 視野数。fieldsOfView[i]とLayeriが対応

    double elapsed_time = sw.RealTime();
    double cpu_time = sw.CpuTime();
    std::cout << TString::Format(
        "Track data successfully loaded. - Elapsed %.2f [s] (CPU: %.2f [s])", 
        elapsed_time, cpu_time
    ) << std::endl;

    // 情報表示
    const size_t entries = tree->GetEntriesFast();
    std::cout << "\n Input file : " << bvxxfile << std::endl;
    std::cout << " # of BT    : " << entries << " tracks" << std::endl;

    // プロット開始
    TDatime time_now;
    TString Time_Now = Form(
        "%d-%02d-%02d %02d:%02d:%02d", 
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(), 
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    sw.Start();
	std::cout << " Plot start : " << Time_Now << std::endl;

	// プログレスバーの初期化
	int page = 0;
	const int total = 40; // 合計ページ数
	MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    // カラーパレットの設定
    MyPalette::SetPalette(Palette);
    if (invertpalette) MyPalette::InvertPalette();
    // カラー定義
    Float_t r1, g1, b1, r2, g2, b2, r3, g3, b3;
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.15))->GetRGB(r1, g1, b1);
    gROOT->GetColor(90)->SetRGB(r1, g1, b1);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.85))->GetRGB(r3, g3, b3);
    gROOT->GetColor(91)->SetRGB(r3, g3, b3);
    gROOT->GetColor(gStyle->GetColorPalette(256 * 0.5))->GetRGB(r2, g2, b2);
    gROOT->GetColor(92)->SetRGB(r2, g2, b2);
    gROOT->GetColor(93)->SetRGB(r2 * 0.6, g2 * 0.6, b2 * 0.6);
    gStyle->SetHistFillColor(92);
    gStyle->SetHistLineColor(93);
    gStyle->SetFuncColor(93);
    gROOT->GetColor(94)->SetRGB( 200./255., 200./255., 203./255.); // light gray
    gStyle->SetGridColor(94);

    // スタイルの設定
    gStyle->SetPadRightMargin(0.1);      // Pad右側のマージン
    gStyle->SetPadLeftMargin(0.1);       // Pad左側のマージン
    gStyle->SetPadTopMargin(0.1);        // Pad上側のマージン
    gStyle->SetPadBottomMargin(0.11);    // Pad下側のマージン
    gStyle->SetLabelOffset(0.008,"xyz"); // 軸ラベル（数値）と軸の距離
    gStyle->SetTitleOffset(1.1,"xyz");   // 軸titleと軸の距離
    gStyle->SetTitleY(0.985);            // タイトルのy方向の位置
    gStyle->SetHistFillStyle(1011); // ヒストグラム内部のパターン 1011=塗りつぶし
    gStyle->SetHistLineStyle(1);    // ヒストグラムの線種 1=直線
    gStyle->SetHistLineWidth(1);    // ヒストグラムの線幅 pixel
    gStyle->SetPadGridX(showGrid); // グリッドの表示
    gStyle->SetPadGridY(showGrid); // グリッドの表示
    gStyle->SetPadTickX(1); // 上側x軸の目盛り表示
    gStyle->SetPadTickY(1); // 右側y軸の目盛り表示
    gStyle->SetStatFont(font_code);         // 統計box内のフォント
    gStyle->SetLabelFont(font_code, "xyz"); // 軸ラベルのフォント
    gStyle->SetTitleFont(font_code, "xyz"); // 軸titleのフォント
    gStyle->SetTitleFont(font_code, "");    // titleのフォント
    gStyle->SetTextFont(font_code);         // textのフォント
    gStyle->SetLegendFont(font_code);       // 凡例のフォント

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());
    global_c1 = c1;
    global_output = output;

    // データの座標の範囲を取得し、表示範囲とビンの数を決定する
    // フィルムの長辺の端から1cm外側までを最大表示範囲とし、縦横比を正しく保って表示する
    // 位置分布等のビン幅は1mm。dx, dy, dzのプロットは計算時間短縮のためフィルムサイズによって変える
    const int MinX = tree->GetMinimum("x");
    const int MaxX = tree->GetMaximum("x");
    const int MinY = tree->GetMinimum("y");
    const int MaxY = tree->GetMaximum("y");
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

    position(c1, tree, pl, AreaParam, TD_range);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    position_projection(c1, tree, entries, pl, TD_range, AreaParam);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("pos*");
    gDirectory->Delete("track_density");
	MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    angle(c1, tree, pl, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    angle_projection(c1, tree, pl, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("ang*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    d_angle(c1, tree);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
	MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    d_angle_Ncut(c1, tree, da_cutX, da_cutY, da_cutPH);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
	MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    d_angle_rl(c1, tree, angle_max, dlat_range, drad_range, 1);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    d_angle_rl(c1, tree, angle_max, dlat_range, drad_range, 2);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*a*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    phvph_2D(c1, tree, vph_range, angle_max, angle_resolution);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*ph*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    for (int i = 0; i < 10; ++i) // 0.0-0.1 ~ 0.9-1.0
    {
        phvph_1D(c1, tree, vph_range, i, 0.1);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*ph*");
        MyUtil::ShowProgress(page, static_cast<double>(page) / total);
    }

    const int phvph_loop = static_cast<int>((angle_max - 0.1) * 2) + 2;
    for (int i = 2; i < phvph_loop; ++i) // 1.0-1.1 ~
    {
        phvph_1D(c1, tree, vph_range, i, 0.5);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("*ph*");
        MyUtil::ShowProgress(page, static_cast<double>(page) / total);
    }

    // ランキングプロットにおけるy軸(VPH)の描画範囲の基準値vph_standardを決定
    if (vph_standard < 0) { // vph_standardが指定されていない場合
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
                vph_sigma = vph_temp->GetStdDev();
                range = vph_sigma * 3;
                if (vph_entries < 1000) cut_type = 1;
            }

            range_max = range_min + range;

            if (cut_type == 0) {
                cutvph_temp = Form(
                "(vph1+vph2)>%d && (vph1+vph2)<%d && tan>1.5 && tan<1.51 && lin>0.03 && lin<0.05"
                "&& dax1<0.02 && dax2<0.02 && day1<0.02 && day2<0.02", 
                range_min, range_max
                );
            } else {
                cutvph_temp = Form(
                "(vph1+vph2)>%d && (vph1+vph2)<%d && tan>1.5 && tan<1.51 && lin>0.08 && lin<0.10"
                "&& dax1<0.02 && dax2<0.02 && day1<0.02 && day2<0.02", 
                range_min, range_max
                );
            }
            tree->Draw("(vph1+vph2)>>vph_temp", cutvph_temp, "goff");
            vph_mean = vph_temp->GetMean();
            if (range_min < 15) break;
        } while (vph_mean < range_min + 0.5 * vph_sigma);

        TF1* gaus = new TF1("gaus", "gaus", 0.0, 1000.0);
        vph_temp->Fit("gaus", "q 0", "", range_min, vph_mean + vph_sigma);
        vph_standard = gaus->GetParameter(1) + 5;
        if (vph_standard < ranking_range_min) vph_standard = ranking_range_min;
        gDirectory->Delete("vph_temp");
        delete gaus;
    }

    for (const auto& params : ranking_params_array) {
        ranking(c1, tree, params, vph_standard);
        c1->Print(output.c_str()); c1->Clear();
        gDirectory->Delete("rank*");
        MyUtil::ShowProgress(page, static_cast<double>(page) / total);
    }

    dxdydz(c1, subtree, AreaParam);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("*_temp");
    gDirectory->Delete("dz*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    dxdy(c1, subtree, AreaParam);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("d*");
	MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    gDirectory->Delete("subtree");

    sensor_not(c1, fieldsOfView, NoTcount);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("not*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    sensor_drdz(c1, NoTcount_same, sensor_dr_sum, sensor_dz_sum);
    c1->Print(output.c_str()); c1->Clear();
    gDirectory->Delete("dz*");
    MyUtil::ShowProgress(page, static_cast<double>(page) / total);

    gDirectory->Delete("sersor_array");

    // PDFファイルを閉じる
    c1->Print((output + "]").c_str());
	if (page < total) page = total;
    MyUtil::ShowProgress(page, 1.0);

    // プロット終了
    time_now.Set();
    Time_Now = Form(
        "%d-%02d-%02d %02d:%02d:%02d", 
        time_now.GetYear(), time_now.GetMonth(), time_now.GetDay(), 
        time_now.GetHour(), time_now.GetMinute(), time_now.GetSecond()
    );
    elapsed_time = sw.RealTime();
    cpu_time = sw.CpuTime();
    std::cout << TString::Format(
        "\n Plot end   : %s - Elapsed %.2f [s] (CPU: %.2f [s])", 
        Time_Now.Data(), elapsed_time, cpu_time
    ) << std::endl;
    std::cout << " Output     : " << output << std::endl;

    delete c1;
    gDirectory->Delete("tree");

    return 0;
}

void position(TCanvas *c1, TTree *tree, int pl, const double *AreaParam, const double *TD_range) noexcept
{
    uint32_t bin = static_cast<uint32_t>(AreaParam[0]);
    double LowX = AreaParam[1];
    double UpX  = AreaParam[2];
    double LowY = AreaParam[3];
    double UpY  = AreaParam[4];

    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.2, "y");
    gStyle->SetTitleOffset(1.6, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    TH2D* position_2D = new TH2D(
        "position_2D", Form("Position PL%03d;x [mm];y [mm];/mm^{2}", pl),
        bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001
    );
    if (TD_range[1] > 0.0) position_2D->GetZaxis()->SetRangeUser(TD_range[0], TD_range[1]);
    tree->Draw("y*0.001:x*0.001 >> position_2D", "", "colz");

	TLegend* pos_lg = new TLegend(0.6, 0.9, 0.75, 1.0);
    pos_lg->SetName("pos_lg");
    gDirectory->Add(pos_lg);
	pos_lg->SetFillStyle(0);
	pos_lg->SetBorderSize(0);
	pos_lg->SetTextSize(0.04);
	pos_lg->AddEntry(position_2D, Form("Entries %.0f", position_2D->GetEntries()), "");
	pos_lg->Draw();
}

void position_projection(
    TCanvas *c1, TTree *tree, const size_t entries, int pl, const double *TD_range, const double *AreaParam
) noexcept
{
    gStyle->SetOptStat("");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    // gDirectoryからposition_2D, pos_lgを取得。なければ新規作成
    TH2D* position_2D = (TH2D*)gDirectory->Get("position_2D");
    TLegend* pos_lg = (TLegend*)gDirectory->Get("pos_lg");
    if (!position_2D || !pos_lg) {
        uint32_t bin = static_cast<uint32_t>(AreaParam[0]);
        double LowX = AreaParam[1];
        double UpX  = AreaParam[2];
        double LowY = AreaParam[3];
        double UpY  = AreaParam[4];

        position_2D = new TH2D(
            "position_2D", Form("Position PL%03d;x [mm];y [mm];/mm^{2}", pl),
            bin, LowX*0.001, UpX*0.001, bin, LowY*0.001, UpY*0.001
        );
        if (TD_range[1] > 0.0) position_2D->GetZaxis()->SetRangeUser(TD_range[0], TD_range[1]);
        tree->Draw("y*0.001:x*0.001 >> position_2D", "", "goff");

        pos_lg = new TLegend(0.6, 0.9, 0.75, 1.0);
        pos_lg->SetFillStyle(0);
        pos_lg->SetBorderSize(0);
        pos_lg->SetTextSize(0.04);
        pos_lg->AddEntry(position_2D, Form("Entries %.0f", position_2D->GetEntries()), "");
    }

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
        "track_density", ";Track Density [/mm^{2}];Frequency", 10000, 0, 100000
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
    density_lg->AddEntry(track_density, "Track Density [/mm^{2}]", "");
    density_lg->AddEntry(track_density, Form("%d areas", density_entries), "");
    density_lg->AddEntry(track_density, Form("Mean   %.2f", density_mean), "");
    density_lg->AddEntry(track_density, Form("Std Dev   %.2f", density_stddev), "");
    density_lg->Draw();
}

void angle(TCanvas *c1, TTree *tree, int pl, const double angle_max, const double angle_resolution) noexcept
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(0.9, "y");
    gStyle->SetTitleOffset(1.6, "z");
    c1->SetRightMargin(0.235);
    c1->SetLeftMargin(0.23);

    uint32_t angle_bin = 2 / angle_resolution * angle_max;

    TString angtitle = Form(
        "Angle PL%03d;tan#it{#theta}_{x};tan#it{#theta}_{y};/(%g rad)^{2}", pl, angle_resolution
    );
    TH2D* angle_2D = new TH2D(
        "angle_2D", angtitle, angle_bin, -angle_max, angle_max, angle_bin, -angle_max, angle_max
    );
    tree->Draw("ay:ax >> angle_2D", "", "colz");

	TLegend* ang_lg = new TLegend(0.6, 0.9, 0.75, 1.0);
    ang_lg->SetName("ang_lg");
    gDirectory->Add(ang_lg);
	ang_lg->SetFillStyle(0);
	ang_lg->SetBorderSize(0);
	ang_lg->SetTextSize(0.04);
	ang_lg->AddEntry(angle_2D, Form("Entries %.0f", angle_2D->GetEntries()), "");
	ang_lg->Draw();
}

void angle_projection(TCanvas *c1, TTree *tree, int pl, const double angle_max, const double angle_resolution) noexcept
{
    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    // gDirectoryからangle_2D, ang_lgを取得。なければ新規作成
	TH2D* angle_2D = (TH2D*)gDirectory->Get("angle_2D");
    TLegend* ang_lg = (TLegend*)gDirectory->Get("ang_lg");
    if (!angle_2D || !ang_lg) {
        uint32_t angle_bin = 2 / angle_resolution * angle_max;

        TString angtitle = Form(
            "Angle PL%03d;tan#it{#theta}_{x};tan#it{#theta}_{y};/(%g rad)^{2}", pl, angle_resolution
        );
        TH2D* angle_2D = new TH2D(
            "angle_2D", angtitle, angle_bin, -angle_max, angle_max, angle_bin, -angle_max, angle_max
        );
        tree->Draw("ay:ax >> angle_2D", "", "goff");

        ang_lg = new TLegend(0.6, 0.9, 0.75, 1.0);
        ang_lg->SetFillStyle(0);
        ang_lg->SetBorderSize(0);
        ang_lg->SetTextSize(0.04);
        ang_lg->AddEntry(angle_2D, Form("Entries %.0f", angle_2D->GetEntries()), "");
    }

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
    gStyle->SetOptStat("");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.5, "y");

    uint32_t angle_bin = angle_max / angle_resolution;

    TString ang1Dtitle = ";#sqrt{tan^{2}#it{#theta}_{x}^{ }#plus tan^{2}#it{#theta}_{y}};Frequency";
	TH1D* angle_1D = new TH1D("angle_1D", ang1Dtitle, angle_bin, 0.0, angle_max);
    angle_1D->SetFillColorAlpha(92, 0.7);
	tree->Draw("tan>>angle_1D");
	TLegend* ang_lg2 = new TLegend(0.45, 0.9, 0.6, 1.0);
	ang_lg2->SetFillStyle(0);
	ang_lg2->SetBorderSize(0);
	ang_lg2->SetTextSize(0.04);
	ang_lg2->AddEntry(angle_1D, Form("Entries %.0f", angle_1D->GetEntries()), "");
	ang_lg2->Draw();
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

    // ヒストグラムを作成して描画するためのラムダ式
    auto createAndDrawHistogram = [&](int pad, const char* name, const char* title, const char* drawExpr) {
        c1->cd(pad);
        TH2D* hist = new TH2D(name, title, 100, -2.0, 2.0, 100, -0.1, 0.1);
        tree->Draw((std::string(drawExpr) + " >> " + name).c_str(), "", "colz");
        hist->Draw("colz");
    };

    createAndDrawHistogram(1, 
        "axdax1", 
        "tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x1} : tan#it{#theta}_{x};"
        "tan#it{#theta}_{x};tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x1}",
        "dax1:ax"
    );
    createAndDrawHistogram(2, 
        "ayday1", 
        "tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y1} : tan#it{#theta}_{y};"
        "tan#it{#theta}_{y};tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y1}", 
        "day1:ay"
    );
    createAndDrawHistogram(3, 
        "axdax2", 
        "tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x2} : tan#it{#theta}_{x};"
        "tan#it{#theta}_{x};tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x2}", 
        "dax2:ax"
    );
    createAndDrawHistogram(4, 
        "ayday2", 
        "tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y2} : tan#it{#theta}_{y};"
        "tan#it{#theta}_{y};tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y2}", 
        "day2:ay"
    );
}

void d_angle_Ncut(
    TCanvas *c1, TTree *tree, const TString da_cutX, const TString da_cutY, const uint8_t da_cutPH
) noexcept
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

    // ヒストグラムを作成して描画するためのラムダ式
    auto createAndDrawHistogram = [&](
        int pad, const char* name, const char* title, const char* drawExpr, const TCut& cut
    ) {
        c1->cd(pad);
        TH2D* hist = new TH2D(name, title, 100, -2.0, 2.0, 100, -0.1, 0.1);
        tree->Draw((std::string(drawExpr) + " >> " + name).c_str(), cut, "colz");
        hist->Draw("colz");
    };

    TCut cut_temp = Form(
        "dax2*dax2<%s*%s&&day2*day2<%s*%s&&ph1>%d&&ph2>%d", 
        da_cutX.Data(), da_cutX.Data(), 
        da_cutY.Data(), da_cutY.Data(), 
        da_cutPH, da_cutPH
    );
    createAndDrawHistogram(
        1, 
        "axdax1", 
        "tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x1} : tan#it{#theta}_{x} (Noise cut);"
        "tan#it{#theta}_{x};tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x1}", 
        "dax1:ax", 
        cut_temp
    );
    createAndDrawHistogram(
        2, 
        "ayday1", 
        "tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y1} : tan#it{#theta}_{y} (Noise cut);"
        "tan#it{#theta}_{y};tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y1}", 
        "day1:ay", 
        cut_temp
    );

    cut_temp = Form(
        "dax1*dax1<%s*%s&&day1*day1<%s*%s&&ph1>%d&&ph2>%d", 
        da_cutX.Data(), da_cutX.Data(), 
        da_cutY.Data(), da_cutY.Data(), 
        da_cutPH, da_cutPH
    );
    createAndDrawHistogram(
        3, 
        "axdax2", 
        "tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x2} : tan#it{#theta}_{x} (Noise cut);"
        "tan#it{#theta}_{x};tan#it{#theta}_{x}^{ }#minus tan#it{#theta}_{x2}", 
        "dax2:ax", 
        cut_temp
    );
    createAndDrawHistogram(
        4, 
        "ayday2", 
        "tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y2} : tan#it{#theta}_{y} (Noise cut);"
        "tan#it{#theta}_{y};tan#it{#theta}_{y}^{ }#minus tan#it{#theta}_{y2}", 
        "day2:ay", 
        cut_temp
    );
}

void d_angle_rl(
    TCanvas *c1, TTree *tree, const double angle_max, const double dlat_range, const double drad_range, const int face
) noexcept
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

    // ヒストグラムを作成して描画するためのラムダ式
    auto createAndDrawHistogram = [&](
        int pad, const char* name, const char* title, const char* drawExpr, double xMax, double yMin, double yMax
    ) {
        c1->cd(pad);
        TH2D* hist = new TH2D(name, title, 200, 0.0, xMax, 200, yMin, yMax);
        tree->Draw(Form(drawExpr, suffix.Data(), suffix.Data(), suffix.Data()), "", "colz");
    };

    createAndDrawHistogram(
        1, 
        "lat", 
        Form(
            "#Deltalateral %s;tan#it{#theta};"
            "#frac{tan#it{#theta}_{x%s}^{ }#times tan#it{#theta}_{y}^{ }"
            "#plus tan#it{#theta}_{y%s}^{ }#times tan#it{#theta}_{x}}"
            "{#sqrt{tan^{2}#it{#theta}_{x}^{ }#plus tan^{2}#it{#theta}_{y}}}",
            suffix.Data(),
            suffix.Data(),
            suffix.Data()
        ),
        "(ax%s*ay-ay%s*ax)/tan:tan>>lat", 
        angle_max, 
        -dlat_range, 
        dlat_range
    );
    createAndDrawHistogram(
        2, 
        "rad", 
        Form(
            "#Deltaradial %s;tan#it{#theta};"
            "#frac{tan#it{#theta}_{x%s}^{ }#times tan#it{#theta}_{x}^{ }"
            "#plus tan#it{#theta}_{y%s}^{ }#times tan#it{#theta}_{y}}"
            "{#sqrt{tan^{2}#it{#theta}_{x}^{ }#plus tan^{2}#it{#theta}_{y}}}^{ }"
            "#minus #sqrt{tan^{2}#it{#theta}_{x}^{ }#plus tan^{2}#it{#theta}_{y}}",
            suffix.Data(),
            suffix.Data(),
            suffix.Data()
        ),
        "(ax%s*ax+ay%s*ay)/tan-tan:tan>>rad", 
        angle_max, 
        -drad_range, 
        drad_range
    );
}

void phvph_2D(
    TCanvas *c1, TTree *tree, const uint32_t vph_range, const double angle_max, const double angle_resolution
) noexcept
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.5, "y");

    c1->Divide(3, 2);
    for (int pad = 1; pad <= 6; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
    }

    uint32_t phvph_bin = static_cast<uint32_t>(angle_max / angle_resolution);

    // ヒストグラムを作成して描画するためのラムダ式
    auto createAndDraw2DHistogram = [&](
        int pad, const char* name, const char* title, const char* drawExpr,
        int yBins, double yMin, double yMax, int NDiv
    ) {
        c1->cd(pad);
        TH2D* hist = new TH2D(name, title, phvph_bin, 0.0, angle_max, yBins, yMin, yMax);
        tree->Draw(Form("%s:tan >> %s", drawExpr, name), "", "colz");
        hist->GetYaxis()->SetNdivisions(NDiv);
    };

    createAndDraw2DHistogram(
        1, "phs0", "PHsum;tan#it{#theta};PHsum", "(ph1+ph2)", 27, 5.5, 32.5, 14
    );
    createAndDraw2DHistogram(
        2, "phs1", "PH1;tan#it{#theta};PH1", "ph1", 11, 5.5, 16.5, 11
    );
    createAndDraw2DHistogram(
        3, "phs2", "PH2;tan#it{#theta};PH2", "ph2", 11, 5.5, 16.5, 11
    );
    createAndDraw2DHistogram(
        4, "vphs0", "VPHsum;tan#it{#theta};VPHsum", "(vph1+vph2)", 2 * vph_range, 0.5, 2 * vph_range + 0.5, 10
    );
    createAndDraw2DHistogram(
        5, "vphs1", "VPH1;tan#it{#theta};VPH1", "vph1", vph_range, 0.5, vph_range + 0.5, 10
    );
    createAndDraw2DHistogram(
        6, "vphs2", "VPH2;tan#it{#theta};VPH2", "vph2", vph_range, 0.5, vph_range + 0.5, 10
    );
}

void phvph_1D(TCanvas *c1, TTree *tree, const uint32_t vph_range, const uint8_t i, const double interval) noexcept
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");

    double tan_low = i * interval;
    double tan_up = i * interval + 0.1;
    TString range = Form("tan>=%.1f&&tan<%.1f", tan_low, tan_up);

    c1->Divide(3, 2);

    // ヒストグラムを作成して描画するためのラムダ式
    auto createAndDrawHistogram = [&](
        int pad, const char* name, const char* title, const char* drawExpr,
        int bins, double xMin, double xMax, const double legendCoords[4]
    ) {
        c1->cd(pad);
        TString histTitle = Form("%s (%.1f^{ }#leq tan#it{#theta} < %.1f);%s;", title, tan_low, tan_up, title);
        TH1D* hist = new TH1D(name, histTitle, bins, xMin, xMax);
        hist->SetFillColorAlpha(92, 0.7);
        tree->Draw(Form("%s >> %s", drawExpr, name), range);

        TLegend* lg = new TLegend(legendCoords[0], legendCoords[1], legendCoords[2], legendCoords[3]);
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.04);
        lg->AddEntry(hist, Form("Entries    %.0f", hist->GetEntries()), "");
        lg->AddEntry(hist, Form("Mean      %.1f", hist->GetMean()), "");
        lg->Draw();
    };

    // EntriesとMeanを表示する座標を指定
    const double legendCoords1[4] = {0.12, 0.76, 0.22, 0.86};
    const double legendCoords2[4] = {0.55, 0.76, 0.65, 0.86};

    // PHsum, PH1, PH2
    createAndDrawHistogram(1, "phsum", "PHsum", "ph1+ph2", 32, 0.5, 32.5, legendCoords1);
    createAndDrawHistogram(2, "ph1", "PH1", "ph1", 16, 0.5, 16.5, legendCoords1);
    createAndDrawHistogram(3, "ph2", "PH2", "ph2", 16, 0.5, 16.5, legendCoords1);

    // VPHsum, VPH1, VPH2
    createAndDrawHistogram(4, "vphsum", "VPHsum", "vph1+vph2", 2 * vph_range, 0, 2 * vph_range, legendCoords2);
    createAndDrawHistogram(5, "vph1", "VPH1", "vph1", vph_range, 0, vph_range, legendCoords2);
    createAndDrawHistogram(6, "vph2", "VPH2", "vph2", vph_range, 0, vph_range, legendCoords2);
}

void ranking(TCanvas *c1, TTree *tree, const RankingParams (&params)[3], uint32_t vph_standard) noexcept
{
    gStyle->SetOptStat("e");
    gStyle->SetStatFormat("6.2f");
    gStyle->SetStatY(0.15);
    gStyle->SetStatX(0.85);
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.15);
    gStyle->SetTitleOffset(1.4, "x");
    gStyle->SetTitleOffset(1.5, "y");

    c1->Divide(3, 2);
    for (int pad = 1; pad <= 6; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.15);
    }

    // 角度範囲が同じ2つのランキングプロットを描画するためのラムダ式
    auto createAndDrawRank = [&](int pad, const RankingParams& param) {
        if (param.tan_low >= param.tan_up) return;

        TCut range = Form("tan>=%.1f&&tan<%.1f", param.tan_low, param.tan_up);

        // x, y
        c1->cd(pad);
        TString rank_title = Form(
            "2D TRank (%.1f^{ }#leq tan#it{#theta} < %.1f);"
            "#sqrt{(dtan#it{#theta}_{x1})^{2}^{ }#plus (dtan#it{#theta}_{y1})^{2}^{ }"
            "#plus (dtan#it{#theta}_{x2})^{2}^{ }#plus (dtan#it{#theta}_{y2})^{2}};"
            "VPHsum",
            param.tan_low,
            param.tan_up
        );
        TH2D* rank = new TH2D(
            "rank",
            rank_title,
            50,
            0.0,
            param.xy_lin_max,
            static_cast<int>(vph_standard + param.vph_max_plus),
            0.0,
            static_cast<double>(vph_standard + param.vph_max_plus)
        );
        tree->Draw("(vph1+vph2):lin >> rank", range, "colz");

        // lateral
        c1->cd(pad + 3);
        TString rankl_title = Form(
            "2D TRank Lateral (%.1f^{ }#leq tan#it{#theta} < %.1f);"
            "Angular difference (Lateral);VPHsum",
            param.tan_low,
            param.tan_up
        );
        TH2D* rankl = new TH2D(
            "rankl",
            rankl_title,
            50,
            0.0,
            param.lat_lin_max,
            static_cast<int>(vph_standard + param.vph_max_plus),
            0.0,
            static_cast<double>(vph_standard + param.vph_max_plus)
        );
        rankl->GetXaxis()->SetNdivisions(505);
        tree->Draw("(vph1+vph2):linl >> rankl", range, "colz");
    };

    for (int i = 0; i < 3; ++i) {
        createAndDrawRank(i + 1, params[i]);
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

    gStyle->SetOptStat("");
    gStyle->SetStatFormat(".4g");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "z");

    TString dz_title = Form(
        "#Deltaz (1.0 < tan#it{#theta} < 1.1);x [mm];y [mm];"
        "Average of#Deltaz [#mum] at each %.1f^{ }#times %.1f mm^{2}",
        pitch,
        pitch
    );
    TString dz_1D_title = Form(
        "#Deltaz at each %.1f^{ }#times %.1f mm^{2}     ;"
        "Average of#Deltaz [#mum] at each %.1f^{ }#times %.1f mm^{2};Frequency",
        pitch,
        pitch,
        pitch,
        pitch
    );
    TH2D* dz_2D = new TH2D("dz_2D", dz_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
    TH1D* dz_1D = new TH1D("dz_1D", dz_1D_title, 5000, 0, 5000);
    TH1D* dz_temp = new TH1D("dz_temp", "dz_temp", 10000, 0, 10000);

    TString dx_title = Form(
        "#Deltax (1.0 < tan#it{#theta} < 1.1, interpolated);x [mm];y [mm];"
        "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2}",
        pitch,
        pitch
    );
    TString dx_1D_title = Form(
        "#Deltax (1.0 < tan#it{#theta} < 1.1, interpolated);"
        "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
        pitch,
        pitch
    );
    TH2D* dx_2D = new TH2D("dx_2D", dx_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
    TH1D* dx_1D = new TH1D("dx_1D", dx_1D_title, 500, -100, 100);
    TH1D* dx_temp = new TH1D("dx_temp", "dx_temp", 100, -1000, 1000);

    TString dy_title = Form(
        "#Deltay (1.0 < tan#it{#theta} < 1.1, interpolated)     ;x [mm];y [mm];"
        "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2}",
        pitch,
        pitch
    );
    TString dy_1D_title = Form(
        "#Deltay (1.0 < tan#it{#theta} < 1.1, interpolated)     ;"
        "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
        pitch,
        pitch
    );
    TH2D* dy_2D = new TH2D("dy_2D", dy_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
    TH1D* dy_1D = new TH1D("dy_1D", dy_1D_title, 500, -100, 100);
    TH1D* dy_temp = new TH1D("dy_temp", "dy_temp", 100, -1000, 1000);

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
            double thickness = dz_temp->GetMean();
            dz_1D->Fill(thickness);
            dz_2D->SetBinContent(ix, iy, thickness);

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
    double dxdy_5sigma = 5 * std::max(dx_1D->GetStdDev(), dy_1D->GetStdDev());
    setRange(dz_1D, dz_2D, dz_1D_mean, dz_5sigma);
    setRange(dx_1D, dx_2D, 0.0, dxdy_5sigma);
    setRange(dy_1D, dy_2D, 0.0, dxdy_5sigma);

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 4; ++pad) {
        c1->GetPad(pad)->SetRightMargin((pad % 2 == 0) ? 0.3 : 0.235);
        c1->GetPad(pad)->SetLeftMargin((pad % 2 == 0) ? 0.165 : 0.23);
    }

    c1->cd(1);
    dz_2D->Draw("colz");

    c1->cd(2);
    dz_1D->SetFillStyle(0);
    dz_1D->SetLineWidth(2);
    dz_1D->Draw();
    MyUtil::PaintBins(dz_1D, dz_1D_mean - dz_5sigma, dz_1D_mean + dz_5sigma); // 各ビンをカラーパレットの色で塗る

    TLegend* dz_lg = new TLegend(0.68, 0.7, 0.9, 0.9);
    dz_lg->SetFillStyle(0);
    dz_lg->SetBorderSize(0);
    dz_lg->SetTextSize(0.04);
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
    gStyle->SetOptStat("");
    gStyle->SetStatFormat(".4g");
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
            "#Deltax (1.0 < tan#it{#theta} < 1.1, interpolated);x [mm];y [mm];"
            "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2}",
            pitch,
            pitch
        );
        TString dx_1D_title = Form(
            "#Deltax (1.0 < tan#it{#theta} < 1.1, interpolated);"
            "Average of#Deltax [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
            pitch,
            pitch
        );
        TH1D* dx_temp = new TH1D("dx_temp", "dx_temp", 100, -1000, 1000);
        dx_2D = new TH2D("dx_2D", dx_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
        dx_1D = new TH1D("dx_1D", dx_1D_title, 500, -100, 100);

        TString dy_title = Form(
            "#Deltay (1.0 < tan#it{#theta} < 1.1, interpolated)     ;x [mm];y [mm];"
            "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2}",
            pitch,
            pitch
        );
        TString dy_1D_title = Form(
            "#Deltay (1.0 < tan#it{#theta} < 1.1, interpolated)     ;"
            "Average of#Deltay [#mum] in each %.1f^{ }#times %.1f mm^{2};Frequency",
            pitch,
            pitch
        );
        dy_2D = new TH2D("dy_2D", dy_title, bin, LowX * 0.001, UpX * 0.001, bin, LowY * 0.001, UpY * 0.001);
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

        double dxdy_5sigma = 5 * std::max(dx_1D->GetStdDev(), dy_1D->GetStdDev());
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
        int hist_bin_min = hist->FindFixBin(hist_min);
        int hist_bin_max = hist->FindFixBin(hist_max);
        for (int i = hist_bin_min; i <= hist_bin_max; ++i) {
            if (hist->GetBinContent(i) == 0) continue;
            int ci = gStyle->GetColorPalette(256 * (hist->GetBin(i) - hist_bin_min) / (hist_bin_max - hist_bin_min));
            MyUtil::PaintBin(hist, i, ci);
        }

        TLegend* legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->AddEntry(hist, Form("Areas      %.0f", hist->GetEntries()), "");
        legend->AddEntry(hist, Form("Mean      %.2f [#mum]", hist->GetMean()), "");
        legend->AddEntry(hist, Form("Std Dev   %.2f [#mum]", hist->GetStdDev()), "");
        legend->Draw();
    };

    // dxとdyの5σのうち、大きい方を基準にして範囲を設定
    double dxdy_5sigma = 5 * std::max(dx_1D->GetStdDev(), dy_1D->GetStdDev());

    drawHistogramWithLegend(c1, 3, dx_1D, -dxdy_5sigma, dxdy_5sigma, "dx_1D", 0.74, 0.7, 0.96, 0.9);
    drawHistogramWithLegend(c1, 4, dy_1D, -dxdy_5sigma, dxdy_5sigma, "dy_1D", 0.68, 0.7, 0.9, 0.9);
}

void sensor_not(TCanvas *c1, const int fieldsOfView[2], const uint32_t NoTcount[2][72]) noexcept
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 2; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.12);
        c1->GetPad(pad)->SetLeftMargin(0.1);
        c1->GetPad(pad)->SetTopMargin(0.12);
        c1->GetPad(pad)->SetBottomMargin(0.08);
    }
    // 非対称なパッドを作成するためのラムダ式
    auto createPad = [](
        const char* name, const char* title, double x1, double y1, double x2, double y2,
        double rightMargin, double leftMargin
    ) {
        TPad* pad = new TPad(name, title, x1, y1, x2, y2);
        pad->SetRightMargin(rightMargin);
        pad->SetLeftMargin(leftMargin);
        pad->Draw();
        return pad;
    };
    TPad* pad31 = createPad("pad31", "pad31", 0.01, 0.0, 0.4, 0.5, 0.01, 0.12);
    TPad* pad32 = createPad("pad32", "pad32", 0.4, 0.0, 0.5, 0.5, 0.3, 0.05);
    TPad* pad41 = createPad("pad41", "pad41", 0.51, 0.0, 0.9, 0.5, 0.01, 0.12);
    TPad* pad42 = createPad("pad42", "pad42", 0.9, 0.0, 1.0, 0.5, 0.3, 0.05);

    TH2D* not0_2D = new TH2D(
        "not0_2D", "Number of Micro Tracks per View for Each Imager (Layer0)", 8, -0.5, 7.5, 9, -0.5, 8.5
    );
    TH1D* not0_1D = new TH1D("not0_1D", "not0_1D", 4000000, -10000000.0, 10000000.0);
    TGraph* not0_graph = new TGraph();
    TH2D* not1_2D = new TH2D(
        "not1_2D", "Number of Micro Tracks per View for Each Imager (Layer1)", 8, -0.5, 7.5, 9, -0.5, 8.5
    );
    TH1D* not1_1D = new TH1D("not1_1D", "not1_1D", 4000000, -10000000.0, 10000000.0);
    TGraph* not1_graph = new TGraph();

    TH2I* sensor_array = new TH2I("sensor_array", "Sensor Array", 8, -0.5, 7.5, 9, -0.5, 8.5);
    sensor_array->SetMarkerSize(1.8);
    sensor_array->SetMarkerColor(93);
    sensor_array->SetBarOffset(-0.25);

    for (int ix = 0; ix < 8; ++ix) {
        for (int iy = 0; iy < 9; ++iy) {
            int sensor_num = ix * 9 + iy;
            auto index = SensorArray::GetIndex(sensor_num);

            // ImagerID表示用の2Dヒストグラムを作成
            sensor_array->SetBinContent(index.x + 1, index.y + 1, sensor_num);

            // 各Imagerの平均飛跡本数を取得
            double not0_value = NoTcount[0][sensor_num] / static_cast<double>(fieldsOfView[0]);
            double not1_value = NoTcount[1][sensor_num] / static_cast<double>(fieldsOfView[1]);

            if (not0_value != 0.0) {
                not0_1D->Fill(not0_value);
                not0_2D->SetBinContent(index.x + 1, index.y + 1, not0_value);
            }
            not0_graph->SetPoint(not0_graph->GetN(), sensor_num, not0_value);

            if (not1_value != 0.0) {
                not1_1D->Fill(not1_value);
                not1_2D->SetBinContent(index.x + 1, index.y + 1, not1_value);
            }
            not1_graph->SetPoint(not1_graph->GetN(), sensor_num, not1_value);
        }
    }

    double not0_mean = not0_1D->GetMean();
    double not0_5sigma = 5 * not0_1D->GetStdDev();
    double not0_range[2] = {not0_mean - not0_5sigma, not0_mean + not0_5sigma};
    double not1_mean = not1_1D->GetMean();
    double not1_5sigma = 5 * not1_1D->GetStdDev();
    double not1_range[2] = {not1_mean - not1_5sigma, not1_mean + not1_5sigma};
    double not_range[2] = {std::min(not0_range[0], not1_range[0]), std::max(not0_range[1], not1_range[1])};

    // 2Dヒストグラムを作成して描画するためのラムダ式
    auto configure2DHistogram = [&](TH2D* hist, int pad) {
        c1->cd(pad);
        gPad->SetGrid(0, 0);
        gPad->SetTicks(0, 0);
        hist->SetMarkerSize(2.0);
        hist->SetBarOffset(0.2);
        hist->GetXaxis()->SetNdivisions(8);
        hist->GetYaxis()->SetNdivisions(9);
        hist->GetXaxis()->SetTickLength(-0.015);
        hist->GetYaxis()->SetTickLength(-0.01);
        hist->GetXaxis()->SetLabelOffset(0.02);
        hist->GetYaxis()->SetLabelOffset(0.02);
        hist->GetZaxis()->SetRangeUser(not_range[0], not_range[1]);
        hist->Draw("colz1 text"); // colz1は0のビンを塗りつぶさない
        sensor_array->Draw("same text");
    };

    // グラフを作成して描画するためのラムダ式
    auto configureGraph = [&](TGraph* graph, TPad* pad, const char* title) {
        pad->cd();
        graph->GetYaxis()->SetTitleOffset(1.9);
        graph->GetXaxis()->SetLimits(-3.5, 74.5);
        graph->GetXaxis()->SetNdivisions(16);
        graph->GetYaxis()->SetRangeUser(not_range[0], not_range[1]);
        graph->SetTitle(title);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.5);
        graph->SetMarkerColor(93);
        graph->SetLineColor(93);
        graph->Draw("a p l");
    };

    // 1Dヒストグラムを作成して描画するためのラムダ式
    auto configure1DHistogram = [&](TH1D* hist, TPad* pad) {
        pad->cd();
        hist->GetXaxis()->SetRangeUser(not_range[0], not_range[1]);
        hist->GetYaxis()->SetNdivisions(5);
        hist->SetTitle("");
        hist->SetLabelColor(0, "xy");
        hist->SetFillColor(93);
        hist->Draw("hbar");
    };

    configure2DHistogram(not0_2D, 1);
    configure2DHistogram(not1_2D, 2);
    configureGraph(not0_graph, pad31, ";Imager ID;Number of Micro Tracks per View");
    configureGraph(not1_graph, pad41, ";Imager ID;Number of Micro Tracks per View");
    configure1DHistogram(not0_1D, pad32);
    configure1DHistogram(not1_1D, pad42);
}

void sensor_drdz(
    TCanvas *c1, const uint32_t NoTcount_same[72], const double dr_sum[72], const double dz_sum[72]
) noexcept
{
    gStyle->SetOptStat("");
    gStyle->SetTitleOffset(1.1, "x");

    c1->Divide(2, 2);
    for (int pad = 1; pad <= 2; ++pad) {
        c1->GetPad(pad)->SetRightMargin(0.12);
        c1->GetPad(pad)->SetLeftMargin(0.1);
        c1->GetPad(pad)->SetTopMargin(0.12);
        c1->GetPad(pad)->SetBottomMargin(0.08);
    }
    // 非対称なパッドを作成するためのラムダ式
    auto createPad = [](
        const char* name, const char* title, double x1, double y1, double x2, double y2,
        double rightMargin, double leftMargin
    ) {
        TPad* pad = new TPad(name, title, x1, y1, x2, y2);
        pad->SetRightMargin(rightMargin);
        pad->SetLeftMargin(leftMargin);
        pad->Draw();
        return pad;
    };
    TPad* pad31 = createPad("pad31", "pad31", 0.01, 0.0, 0.4, 0.5, 0.01, 0.12);
    TPad* pad32 = createPad("pad32", "pad32", 0.4, 0.0, 0.5, 0.5, 0.3, 0.05);
    TPad* pad41 = createPad("pad41", "pad41", 0.51, 0.0, 0.9, 0.5, 0.01, 0.12);
    TPad* pad42 = createPad("pad42", "pad42", 0.9, 0.0, 1.0, 0.5, 0.3, 0.05);

    TH2D* dr_2D = new TH2D(
        "dr_2D",
        "Average of^{ }#sqrt{#Deltax^{2}^{ }#plus #Deltay^{2}} (Same Imager on Both Layers)",
        8, -0.5, 7.5, 9, -0.5, 8.5
    );
    TH1D* dr_1D = new TH1D("dr_1D", "dr_1D", 400000, -10000.0, 10000.0);
    TGraph* dr_graph = new TGraph();
    TH2D* dz_2D = new TH2D(
        "dz_2D",
        "Average of^{ }#Deltaz (Same Imager on Both Layers)",
        8, -0.5, 7.5, 9, -0.5, 8.5
    );
    TH1D* dz_1D = new TH1D("dz_1D", "dz_1D", 200000, -10000.0, 10000.0);
    TGraph* dz_graph = new TGraph();

    // gDirectoryからsensor_arrayを取得。なければ新規作成
    bool got_sensor_array = true;
    TH2I* sensor_array = (TH2I*)gDirectory->Get("sensor_array");
    if (!sensor_array) {
        got_sensor_array = false;
        sensor_array = new TH2I("sensor_array", "Sensor Array", 8, -0.5, 7.5, 9, -0.5, 8.5);
        sensor_array->SetMarkerSize(1.8);
        sensor_array->SetMarkerColor(93);
        sensor_array->SetBarOffset(-0.25);
    }

    for (int ix = 0; ix < 8; ++ix) {
        for (int iy = 0; iy < 9; ++iy) {
            int sensor_num = ix * 9 + iy;
            auto index = SensorArray::GetIndex(sensor_num);

            if (!got_sensor_array) sensor_array->SetBinContent(index.x + 1, index.y + 1, sensor_num);

            // 各Imagerの平均dr, dzを取得
            double dr_value = 0.0, dz_value = 0.0;
            if (NoTcount_same[sensor_num] != 0) { // 0除算を避ける
                dr_value = dr_sum[sensor_num] / static_cast<double>(NoTcount_same[sensor_num]);
                dz_value = dz_sum[sensor_num] / static_cast<double>(NoTcount_same[sensor_num]);
            }

            if (dr_value != 0.0) {
                dr_1D->Fill(dr_value);
                dr_2D->SetBinContent(index.x + 1, index.y + 1, dr_value);
            }
            dr_graph->SetPoint(dr_graph->GetN(), sensor_num, dr_value);

            if (dz_value != 0.0) {
                dz_1D->Fill(dz_value);
                dz_2D->SetBinContent(index.x + 1, index.y + 1, dz_value);
            }
            dz_graph->SetPoint(dz_graph->GetN(), sensor_num, dz_value);
        }
    }

    double dr_mean = dr_1D->GetMean();
    double dr_5sigma = 5 * dr_1D->GetStdDev();
    double dr_range[2] = {dr_mean - dr_5sigma, dr_mean + dr_5sigma};
    double dz_mean = dz_1D->GetMean();
    double dz_5sigma = 5 * dz_1D->GetStdDev();
    double dz_range[2] = {dz_mean - dz_5sigma, dz_mean + dz_5sigma};

    // 2Dヒストグラムを作成して描画するためのラムダ式
    auto configure2DHistogram = [&](TH2D* hist, int pad, const double range[2]) {
        c1->cd(pad);
        gPad->SetGrid(0, 0);
        gPad->SetTicks(0, 0);
        hist->SetMarkerSize(2.0);
        hist->SetBarOffset(0.2);
        hist->GetXaxis()->SetNdivisions(8);
        hist->GetYaxis()->SetNdivisions(9);
        hist->GetXaxis()->SetTickLength(-0.015);
        hist->GetYaxis()->SetTickLength(-0.01);
        hist->GetXaxis()->SetLabelOffset(0.02);
        hist->GetYaxis()->SetLabelOffset(0.02);
        hist->GetZaxis()->SetRangeUser(range[0], range[1]);
        hist->Draw("colz1 text"); // colz1は0のビンを塗りつぶさない
        sensor_array->Draw("same text");
    };

    // グラフを作成して描画するためのラムダ式
    auto configureGraph = [&](TGraph* graph, TPad* pad, const char* title, const double offset, const double range[2]) {
        pad->cd();
        graph->GetYaxis()->SetTitleOffset(offset);
        graph->GetXaxis()->SetLimits(-3.5, 74.5);
        graph->GetXaxis()->SetNdivisions(16);
        graph->GetYaxis()->SetRangeUser(range[0], range[1]);
        graph->SetTitle(title);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(0.5);
        graph->SetMarkerColor(93);
        graph->SetLineColor(93);
        graph->Draw("a p l");
    };

    // 1Dヒストグラムを作成して描画するためのラムダ式
    auto configure1DHistogram = [&](TH1D* hist, TPad* pad, const double range[2]) {
        pad->cd();
        hist->GetXaxis()->SetRangeUser(range[0], range[1]);
        hist->GetYaxis()->SetNdivisions(5);
        hist->SetTitle("");
        hist->SetLabelColor(0, "xy");
        hist->SetFillColor(93);
        hist->Draw("hbar");
    };

    configure2DHistogram(dr_2D, 1, dr_range);
    configure2DHistogram(dz_2D, 2, dz_range);
    configureGraph(dr_graph, pad31, ";Imager ID;Average of^{ }#sqrt{#Deltax^{2}^{ }#plus #Deltay^{2}}", 1.7, dr_range);
    configureGraph(dz_graph, pad41, ";Imager ID;Average of^{ }#Deltaz", 1.9, dz_range);
    configure1DHistogram(dr_1D, pad32, dr_range);
    configure1DHistogram(dz_1D, pad42, dz_range);
}
