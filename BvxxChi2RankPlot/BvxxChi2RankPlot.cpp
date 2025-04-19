#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TString.h>
#include <TColor.h>
#include <TLine.h>
#include <TStopwatch.h>
#include <TError.h>

#include <VxxReader/VxxReader.h>
#include <VxxReader/KeywordArgs.h>
#include <VxxReader/Template.h>
#include <VxxReader/netscan_data_types_ui.h>
#include <argparse/argparse.hpp>

#include <ROOT6/MyUtil.hpp>
#include <ROOT6/MyPalette.hpp>
#include <ROOT6/StderrSuppressor.hpp>

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[]);
TTree* read_bvxx(const std::string& bvxxfile, int pl, int area);
TH2D* PlotRanking_2d(TCanvas *c1, TTree *tree, double ang_min, double ang_max,
    double thickness_2, double dx_2, double dz_2,
    const std::string& param_file, int i, const std::string& title
);

int main(int argc, char* argv[])
{
    // 処理時間を計測
    TStopwatch sw;

    // argparseを使用して引数を解析
    std::cout << "\nInitializing..." << std::endl;
    argparse::ArgumentParser parser("BvxxChi2RankPlot.exe", "1.0.0");
    parse_arguments(parser, argc, argv);

    // 引数を取得
    const auto bvxxfile = parser.get<std::string>("input_bvxx");
    const auto angle_acc_file = parser.get<std::string>("angle_acc_file");
    const auto param_file = parser.get<std::string>("param_file");
    const auto output_arg = parser.get<std::string>("--output");
    auto title = parser.get<std::string>("--title");
    const auto thickness = parser.get<double>("--thickness");
    auto font_number = parser.get<int>("--font_number");
    auto hideGrid = parser.get<bool>("--hide_grid");
    auto darkmode = parser.get<bool>("--dark_mode");
    auto palette_arg = parser.get<std::string>("--palette");
    auto NContours = parser.get<int>("--contours");
    auto invertpalette = parser.get<bool>("--invert_palette");
    auto negatepalette = parser.get<bool>("--negate_palette");

    // 出力PDFファイル名の設定
    const std::string output = (output_arg.empty())
        ? (bvxxfile + "_Chi2Rank" + title + ".pdf")
        : ((output_arg.size() > 4 && output_arg.substr(output_arg.size() - 4) == ".pdf")
            ? output_arg
            : output_arg + ".pdf");

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
    // カラーの設定
    gROOT->GetColor(94)->SetRGB( 200./255., 200./255., 200./255.);
    gROOT->GetColor(95)->SetRGB(  60./255.,  60./255.,  60./255.);
    gStyle->SetGridColor(darkmode ? 95 : 94);       // グリッド
    gStyle->SetCanvasColor(darkmode ? 1 : 0);       // キャンバス(全体の背景)
    gStyle->SetStatColor(darkmode ? 1 : 0);         // 統計box
    gStyle->SetAxisColor(darkmode ? 0 : 1, "xyz");  // 軸
    gStyle->SetLabelColor(darkmode ? 0 : 1, "xyz"); // 軸ラベル(数値)
    gStyle->SetTitleColor(darkmode ? 0 : 1, "xyz"); // 軸タイトル
    gStyle->SetTitleTextColor(darkmode ? 0 : 1);    // メインのタイトル
    gStyle->SetFrameLineColor(darkmode ? 0 : 1);    // 描画エリアの枠
    gStyle->SetStatTextColor(darkmode ? 0 : 1);     // 統計box内のtext
    gStyle->SetLineColor(darkmode ? 0 : 1);         // 統計boxの枠など

    // スタイルの設定
	gStyle->SetOptStat("e");             // 統計boxの表示
    gStyle->SetStatY(0.975);              // 統計boxのy方向の位置
    gStyle->SetPadRightMargin(0.1);      // Pad右側のマージン
    gStyle->SetPadLeftMargin(0.1);       // Pad左側のマージン
    gStyle->SetPadTopMargin(0.1);        // Pad上側のマージン
    gStyle->SetPadBottomMargin(0.11);    // Pad下側のマージン
    gStyle->SetLabelOffset(0.008,"xyz"); // 軸ラベル(数値)と軸の距離
    gStyle->SetTitleOffset(1.2, "x");    // x軸titleと軸の距離
    gStyle->SetTitleOffset(1.5, "y");    // y軸titleと軸の距離
    gStyle->SetTitleY(0.985);            // タイトルのy方向の位置
    gStyle->SetPadGridX(!hideGrid); // グリッドの表示
    gStyle->SetPadGridY(!hideGrid); // グリッドの表示
    gStyle->SetPadTickX(1);         // 上側x軸の目盛り表示
    gStyle->SetPadTickY(1);         // 右側y軸の目盛り表示
    const int font_code = 10 * font_number + 2;
    gStyle->SetStatFont(font_code);         // 統計box内のフォント
    gStyle->SetLabelFont(font_code, "xyz"); // 軸ラベルのフォント
    gStyle->SetTitleFont(font_code, "xyz"); // 軸titleのフォント
    gStyle->SetTitleFont(font_code, "");    // titleのフォント

    // エラーメッセージ未満のROOTのメッセージを非表示に設定
    gErrorIgnoreLevel = kError;

    // angle_acc_fileの読み込み
    std::ifstream ifs(angle_acc_file.c_str());
	int pl, area;
	double dx, dz;
	ifs >> pl >> area >> dx >> dz;
    ifs.close();

    // 2乗の形でしか使わないためあらかじめ計算しておく
    double dx_2 = dx * dx;
    double dz_2 = dz * dz;

    // TTreeの作成
    std::cout << "Starting to read bvxx file..." << std::endl;
    TTree* tree = read_bvxx(bvxxfile, pl, area);

    double elapsed_time = sw.RealTime();
    double cpu_time = sw.CpuTime();
    std::cout << Form(
        "Angle accuracy data & track data successfully loaded. - Elapsed %.2f [s] (CPU: %.2f [s])", 
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

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());

	// プログレスバーの初期化
	int page = 0;
    const int total = 4;// 合計ページ数
	MyUtil::ShowProgress(page, total);

	TH2D** hist = new TH2D*[6];
    // 2乗の形でしか使わないためあらかじめ計算しておく
    double thickness_2 = thickness * thickness;

    constexpr std::array<std::pair<double, double>, 24> angle_ranges = {{
        {0.0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6},
        {0.6, 0.8}, {0.8, 1.0}, {1.0, 1.2}, {1.2, 1.4}, {1.4, 1.6}, {1.6, 1.8},
        {1.8, 2.1}, {2.1, 2.4}, {2.4, 2.7}, {2.7, 3.0}, {3.0, 3.3}, {3.3, 3.6},
        {3.6, 4.0}, {4.0, 4.4}, {4.4, 4.8}, {4.8, 5.2}, {5.2, 5.6}, {5.6, 6.0}
    }};

    if (title != "") title = title + " / ";
    for (size_t i = 0; i < angle_ranges.size(); i += 6) {
        c1->Divide(3, 2);
        for (size_t j = 0; j < 6 && i + j < angle_ranges.size(); ++j) {
            hist[j] = PlotRanking_2d(
                c1, tree, angle_ranges[i + j].first, angle_ranges[i + j].second,
                thickness_2, dx_2, dz_2, param_file, j, title
            );
        }
        c1->Print(output.c_str()); c1->Clear();
        MyUtil::ShowProgress(page, total);
        for (size_t j = 0; j < 6 && i + j < angle_ranges.size(); ++j) {
            delete hist[j];
        }
    }

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
        "\n Plot end   : %s - Elapsed %.2f [s] (CPU: %.2f [s])", 
        Time_Now.Data(), elapsed_time, cpu_time
    ) << std::endl;
    std::cout << " Output     : " << output << std::endl;

    delete c1;
    gDirectory->Delete("tree");

    return 0;
}

TH2D* PlotRanking_2d(TCanvas *c1, TTree *tree, double ang_min, double ang_max,
    double thickness_2, double dx_2, double dz_2,
    const std::string& param_file, int i, const std::string& title
) {
    c1->cd(i + 1);
    gPad->SetRightMargin(0.15);
    TString hist_title = Form(
        "%s%.1f^{ }#leq tan#it{#theta} < %.1f;#it{#chi}^{2};VPH",
        title.c_str(), ang_min, ang_max
    );

    TH2D* hist;
	if (0.0 <= (ang_min + ang_max) / 2 && (ang_min + ang_max) / 2 < 0.1) {
		hist = new TH2D(Form("hist_%d", i), hist_title, 100, 0, 100, 100, 0, 200);
	} else if (0.1 <= (ang_min + ang_max) / 2 && (ang_min + ang_max) / 2 < 0.3) {
		hist = new TH2D(Form("hist_%d", i), hist_title, 60, 0, 60, 150, 0, 150);
	} else {
		hist = new TH2D(Form("hist_%d", i), hist_title, 60, 0, 60, 100, 0, 100);
	}

    int vph;
	double chi, tan, dr1_2, dr2_2, dl1_2, dl2_2, sigma_r_2, sigma_l_2;
    tree->SetBranchAddress("tan", &tan);
    tree->SetBranchAddress("dr1_2", &dr1_2);
    tree->SetBranchAddress("dr2_2", &dr2_2);
    tree->SetBranchAddress("dl1_2", &dl1_2);
    tree->SetBranchAddress("dl2_2", &dl2_2);
    tree->SetBranchAddress("vph", &vph);

    for (int i = 0; i < tree->GetEntries(); ++i) {

		tree->GetEntry(i);
		if (tan < ang_min || tan > ang_max) continue;

		sigma_r_2 = 2 * (dx_2 + tan * tan * dz_2) / thickness_2;
		sigma_l_2 = 2 * dx_2 / thickness_2;
		chi = (dl1_2 / sigma_l_2) + (dl2_2 / sigma_l_2) + (dr1_2 / sigma_r_2) + (dr2_2 / sigma_r_2);

        hist->Fill(chi, vph);
	}

	double x0, y0, x1, y1, x2, y2, angle_min, angle_max, num, effi;
	std::ifstream ifs_p(param_file.c_str());

	TLine* l0 = nullptr;
	TLine* l1 = nullptr;
	TLine* l2 = nullptr;
	TLine* l3 = nullptr;

	while (ifs_p >> angle_min >> angle_max >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> num >> effi) {
		if (ang_min <= (angle_min + angle_max) / 2 && (angle_min + angle_max) / 2 < ang_max) {
			l0 = new TLine(0, y0, x0, y0);
			l1 = new TLine(x0, y0, x1, y1);
			l2 = new TLine(x1, y1, x2, y2);
			if (0.0 <= (angle_min + angle_max) / 2 && (angle_min + angle_max) / 2 <= 0.1) {
                l3 = new TLine(x2, y2, 100, y2);
            } else {
                l3 = new TLine(x2, y2, 60, y2);
            }
			break;
		}
	}

	hist->Draw("colz");

    for (TLine* line : {l0, l1, l2, l3}) {
        if (line) {
            line->SetLineWidth(2);
            line->Draw("same");
        }
    }

	return hist;
}

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[])
{
    parser.set_usage_max_line_width(80);

    // 必須引数: bvxxファイル
    parser.add_argument("input_bvxx")
        .help("Path to the bvxx file to be processed.")
        .required();
    // 必須引数: 角度精度ファイル
    parser.add_argument("angle_acc_file")
        .help("Path to the angle accuracy file\n"
            "generated by \"BvxxAngleAccuracy.exe\" or \"Plot_angle_acc_ver2.C\".")
        .required();
    // 必須引数: パラメータファイル
    parser.add_argument("param_file")
        .help("Path to the parameter file\n"
            "generated by \"Basetrack_Rankingcut_chi2.exe\".")
        .required();
    // オプション引数: 出力PDFファイル名
    parser.add_argument("-o", "--output")
        .help("Output PDF file name. [default: [input_bvxx]_Chi2Rank[--title].pdf]")
        .default_value(std::string());
    // オプション引数: タイトル
    parser.add_argument("-t", "--title")
        .help("Title of histograms such as \"thick 0\".")
        .default_value(std::string());
    // オプション引数: 乳剤厚
    parser.add_argument("-thick", "--thickness")
        .help("Thickness of one emulsion layer.")
        .default_value(60.0)
        .scan<'g', double>();
    // オプション引数: その他の設定
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

TTree* read_bvxx(const std::string& bvxxfile, int pl, int area)
{
    TTree* tree = new TTree("tree", "");

    // ブランチを作成
    int vph;
    double tan, dr1_2, dl1_2, dr2_2, dl2_2;
    tree->Branch("vph", &vph, "vph/I");
    const std::array<std::pair<const char*, void*>, 5> branches = {{
        {"tan", &tan}, {"dr1_2", &dr1_2}, {"dl1_2", &dl1_2}, {"dr2_2", &dr2_2}, {"dl2_2", &dl2_2}
    }};
    for (const auto& branch : branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    // bvxxファイルの読み込み
    vxx::BvxxReader br;

    if (br.Begin(bvxxfile, pl, area)) {
        vxx::HashEntry h;
        vxx::base_track_t b;

        while (br.NextHashEntry(h)) {
            while (br.NextBaseTrack(b)) {
                double ax = b.ax;
                double ay = b.ay;
                double ax1 = b.m[0].ax;
                double ay1 = b.m[0].ay;
                double ax2 = b.m[1].ax;
                double ay2 = b.m[1].ay;
                int ph1 = b.m[0].ph;
                int ph2 = b.m[1].ph;

                vph = (ph1 + ph2) % 10000;

                tan = sqrt(ax * ax + ay * ay);

                double dr1 = ((ax1 - ax) * ax + (ay1 - ay) * ay) / tan;
                double dl1 = ((ax1 - ax) * ay - (ay1 - ay) * ax) / tan;
                double dr2 = ((ax2 - ax) * ax + (ay2 - ay) * ay) / tan;
                double dl2 = ((ax2 - ax) * ay - (ay2 - ay) * ax) / tan;

                // 2乗の形でしか使わないため、あらかじめ計算して2乗をtreeに保存しておく
                dr1_2 = dr1 * dr1;
                dl1_2 = dl1 * dl1;
                dr2_2 = dr2 * dr2;
                dl2_2 = dl2 * dl2;

                tree->Fill();
            }
        }
        br.End();
    } else {
        std::cerr <<
        "Error: Failed to load bvxx file due to incompatible PL or Area number."
        << std::endl;
        std::exit(1);
    } // bvxxファイルの読み込み終了

    return tree;
}
