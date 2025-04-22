#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TString.h>
#include <TColor.h>
#include <TLegend.h>

#include <VxxReader/VxxReader.h>
#include <VxxReader/KeywordArgs.h>
#include <VxxReader/Template.h>
#include <VxxReader/netscan_data_types_ui.h>
#include <argparse/argparse.hpp>

#include <ROOT6/StderrSuppressor.hpp>

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[]);
TTree* read_bvxx(const std::string& bvxxfile, int &pl);
void diff_hist(
    TCanvas *c1, std::string output, TTree *tree, const int Num, double *angle,
    double *sigma, const char* rl, int layer
);
void GetThresholdMinMax(TH1F *hist, double threshold, double &valMin, double &valMax);

int main(int argc, char* argv[])
{
    // argparseを使用して引数を解析
    argparse::ArgumentParser parser("BvxxAngleAccuracy.exe", "1.0.0");
    parse_arguments(parser, argc, argv);

    // 引数を取得
    const auto bvxxfile = parser.get<std::string>("input_bvxx");
    const auto output_txt = parser.get<std::string>("output_txt");
    const auto output_pdf = parser.get<std::string>("--output_pdf");
    auto pl = parser.get<int>("--pl_number");
    const auto area = parser.get<int>("--area_number");
    const auto thickness = parser.get<double>("--thickness");

    // 出力PDFファイル名の設定
    const std::string output = (output_pdf.empty())
        ? (bvxxfile + "_angle_accuracy.pdf")
        : ((output_pdf.size() > 4 && output_pdf.substr(output_pdf.size() - 4) == ".pdf")
            ? output_pdf
            : output_pdf + ".pdf");

    int fontid = 42;                     // 42=Helvetica, 132=Times
    gStyle->SetStatFont(fontid);         // 統計box内
    gStyle->SetLabelFont(fontid, "xyz"); // 軸ラベル（数値）
    gStyle->SetTitleFont(fontid, "xyz"); // 軸title
    gStyle->SetTitleFont(fontid, "");    // title
    gStyle->SetLegendFont(fontid);       // 凡例

    gROOT->GetColor(90)->SetRGB( 255./255.,  75./255.,   0./255.); // red
    gROOT->GetColor(91)->SetRGB(   0./255.,  90./255., 255./255.); // blue
    gROOT->GetColor(92)->SetRGB( 200./255., 200./255., 203./255.); // light gray
    gStyle->SetHistLineWidth(1);
    gStyle->SetLineColor(1);
    gStyle->SetTextColor(1);
	gStyle->SetFuncColor(1);
	gStyle->SetMarkerSize(0.7);
    gStyle->SetGridColor(92);

    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.11);
    gStyle->SetTitleOffset(1.1,"x");
    gStyle->SetTitleOffset(1.5,"y");
    gStyle->SetTitleY(0.985);

    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

	gStyle->SetOptStat("e");
	gStyle->SetOptFit(1102);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);

    // エラーメッセージ未満のROOTのメッセージを非表示に設定
    gErrorIgnoreLevel = kError;

    // TTreeの作成
    std::cout << "Starting to read bvxx file..." << std::endl;
    TTree* tree = read_bvxx(bvxxfile, pl);

    // キャンバスとPDFファイルの作成
    gStyle->SetPaperSize(TStyle::kA4);
    TCanvas* c1 = new TCanvas("c1");
    c1->Print((output + "[").c_str());

	constexpr int Num = 26;
	double angle[Num], sigma[4][Num];

    diff_hist(c1, output, tree, Num, angle, sigma[0], "Lateral", 1);
    diff_hist(c1, output, tree, Num, angle, sigma[1], "Lateral", 2);
    diff_hist(c1, output, tree, Num, angle, sigma[2], "Radial", 1);
    diff_hist(c1, output, tree, Num, angle, sigma[3], "Radial", 2);

    gDirectory->Delete("hist*");
    gDirectory->Delete("tree");

    TGraph **graph = new TGraph*[4];
    for (int i = 0; i < 4; ++i) {
        graph[i] = new TGraph;
        for (int j = 0; j < Num; ++j) {
            if (sigma[i][j] > 0) {
                graph[i]->SetPoint(graph[i]->GetN(), angle[j], sigma[i][j]);
            }
        }
        graph[i]->SetMarkerColor(i % 2 == 0 ? 90 : 91);
        graph[i]->SetMarkerStyle(i % 2 == 0 ? 20 : 25);
    }

    TH1* frame = c1->DrawFrame(0.0, 0.0, 5.0, 0.015, "Lateral angle accuracy;tan#it{#theta};#it{#sigma}_{lateral}");
	graph[0]->Draw("sameP");
	graph[1]->Draw("sameP");
    TLegend* lat_lg = new TLegend(0.17, 0.19, 0.37, 0.34);
    lat_lg->SetFillStyle(0);
    lat_lg->SetBorderSize(0);
    lat_lg->SetTextSize(0.04);
    lat_lg->AddEntry(graph[0], "Layer 1", "p");
    lat_lg->AddEntry(graph[1], "Layer 2", "p");
    lat_lg->Draw();
    c1->Print(output.c_str()); c1->Clear();

	frame = c1->DrawFrame(0.0, 0.0, 5.0, 0.2, "Radial angle accuracy;tan#it{#theta};#it{#sigma}_{radial}");
	graph[2]->Draw("sameP");
	graph[3]->Draw("sameP");
    TLegend* rad_lg = new TLegend(0.69, 0.19, 0.89, 0.34);
    rad_lg->SetFillStyle(0);
    rad_lg->SetBorderSize(0);
    rad_lg->SetTextSize(0.04);
    rad_lg->AddEntry(graph[2], "Layer 1", "p");
    rad_lg->AddEntry(graph[3], "Layer 2", "p");
    rad_lg->Draw();
    c1->Print(output.c_str()); c1->Clear();

    TGraph **graph2 = new TGraph*[2];
    for (int i = 0; i < 2; ++i) {
        graph2[i] = new TGraph;
        for (int j = 0; j < Num; ++j) {
            for (int k = 0; k < 2; ++k) {
                if (sigma[i * 2 + k][j] > 0) {
                    graph2[i]->SetPoint(graph2[i]->GetN(), angle[j], sigma[i * 2 + k][j]);
                }
            }
        }
        graph2[i]->SetMarkerColor(0);
        graph2[i]->SetMarkerStyle(1);
    }

	gStyle->SetOptFit(112);
    gStyle->SetStatFormat(".4f");
    gStyle->SetStatY(0.26);

	TF1* flat = new TF1("flat", Form("sqrt(2)/%f*[0]", thickness), 0.0, 5.0);
	flat->SetParameter(0, 1);
	graph2[0]->Fit(flat, "q0", "", 0.0, 5.0);
	double dx = flat->GetParameter(0);
    frame = c1->DrawFrame(0.0, 0.0, 5.0, 0.015, "Lateral angle accuracy;tan#it{#theta};#it{#sigma}_{Lateral}");
	graph2[0]->Draw("sameP");
	graph[0]->Draw("sameP");
	graph[1]->Draw("sameP");
    flat->Draw("same");
    TLegend* lat_lg2 = new TLegend(0.3, 0.25, 0.7, 0.45);
    lat_lg2->SetFillStyle(0);
    lat_lg2->SetBorderSize(0);
    lat_lg2->SetTextSize(0.05);
    lat_lg2->AddEntry(graph2[0], Form("#it{#sigma}_{Lateral} =^{ }#frac{#sqrt{2}}{%.0f} p_{0}", thickness), "");
    lat_lg2->Draw();
    lat_lg->Draw();
    c1->Print(output.c_str()); c1->Clear();

    gStyle->SetStatX(0.39);
    gStyle->SetStatY(0.9);

	TF1* frad = new TF1("frad", Form("sqrt(2*([0]^2+[1]^2*x^2))/%f", thickness), 0.0, 5.0);
    frad->FixParameter(0, dx);
	frad->SetParameter(1, 4);
	graph2[1]->Fit(frad, "q0", "", 0.0, 5.0);
	double dz= frad->GetParameter(1);
	frame = c1->DrawFrame(0.0, 0.0, 5.0, 0.2, "Radial angle accuracy;tan#it{#theta};#it{#sigma}_{Radial}");
	graph2[1]->Draw("sameP");
	graph[2]->Draw("sameP");
	graph[3]->Draw("sameP");
    frad->Draw("same");
    TLegend* rad_lg2 = new TLegend(0.12, 0.37, 0.32, 0.77);
    rad_lg2->SetFillStyle(0);
    rad_lg2->SetBorderSize(0);
    rad_lg2->SetTextSize(0.05);
    rad_lg2->AddEntry(
        graph2[0],
        Form(
            "#it{#sigma}_{Radial} =^{ }#frac{#sqrt{2}}{%.0f} #sqrt{p_{0}^{2}^{ }+ p_{1}^{2}^{ }tan^{2}#it{#theta}}",
            thickness
        ),
        ""
    );
    rad_lg2->Draw();
    rad_lg->Draw();
    c1->Print(output.c_str()); c1->Clear();

    for (int i = 0; i < 4; ++i) {
        delete graph[i];
    }
    delete[] graph;
    delete lat_lg;
    delete rad_lg;
    for (int i = 0; i < 2; ++i) {
        delete graph2[i];
    }
    delete[] graph2;
    gDirectory->Delete("graph*");
    delete flat;
    delete frad;
    delete lat_lg2;
    delete rad_lg2;

    c1->Print((output + "]").c_str());
    delete c1;

	std::ofstream ofs(output_txt);
	ofs << Form("%d %d %6.4f %6.4f", pl, area, dx, dz) << std::endl;
    ofs.close();

    std::cout << "Output.pdf: " << output << std::endl;
    std::cout << "Output.txt: " << output_txt << std::endl;
    return 0;
}

void parse_arguments(argparse::ArgumentParser& parser, int argc, char* argv[])
{
    parser.set_usage_max_line_width(80);

    // 必須引数: bvxxファイル
    parser.add_argument("input_bvxx")
        .help("Path to the bvxx file to be processed.")
        .required();
    // 必須引数: 出力テキストファイル名
    parser.add_argument("output_txt")
        .help("Path to the output text file.")
        .required();
    // オプション引数: 出力PDFファイル名
    parser.add_argument("-pdf", "--output_pdf")
        .help("Output PDF file name. [default: [input_bvxx]_angle_accuracy.pdf]")
        .default_value(std::string());
    // オプション引数: PL番号
    parser.add_argument("-pl", "--pl_number")
        .help("Plate number to be used. (default: auto)")
        .default_value(-1)
        .scan<'i', int>();
    // オプション引数: area番号
    parser.add_argument("-area", "--area_number")
        .help("Area number to be used.")
        .default_value(0)
        .scan<'i', int>();
    // オプション引数: 乳剤厚
    parser.add_argument("-thick", "--thickness")
        .help("Thickness of one emulsion layer.")
        .default_value(60.0)
        .scan<'g', double>();

    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << "\nError: " << err.what() << std::endl;
        std::cerr << parser;
        std::cerr << "\nError: " << err.what() << std::endl;
        std::exit(1);
    }
}

TTree* read_bvxx(const std::string& bvxxfile, int &pl)
{
    TTree* tree = new TTree("tree", "");

    // ブランチを作成
    double tan, dr1, dl1, dr2, dl2;
    const std::array<std::pair<const char*, void*>, 5> branches = {{
        {"tan", &tan}, {"dr1", &dr1}, {"dl1", &dl1}, {"dr2", &dr2}, {"dl2", &dl2}
    }};
    for (const auto& branch : branches) {
        tree->Branch(branch.first, branch.second, (std::string(branch.first) + "/D").c_str());
    }

    // bvxxファイルの読み込み
    vxx::BvxxReader br;

    // bvxxファイルを読み込むためのラムダ式
    auto readBaseTrack = [&](const vxx::base_track_t& b) {
        double ax = b.ax;
        double ay = b.ay;
        double ax1 = b.m[0].ax;
        double ay1 = b.m[0].ay;
        double ax2 = b.m[1].ax;
        double ay2 = b.m[1].ay;

        tan = sqrt(ax * ax + ay * ay);
        dr1 = ((ax1 - ax) * ax + (ay1 - ay) * ay) / tan;
        dl1 = ((ax1 - ax) * ay - (ay1 - ay) * ax) / tan;
        dr2 = ((ax2 - ax) * ax + (ay2 - ay) * ay) / tan;
        dl2 = ((ax2 - ax) * ay - (ay2 - ay) * ax) / tan;

        tree->Fill();
    };

    StderrSuppressor sup; // 標準エラー出力を抑制するためのオブジェクト
    bool readSuccess = false;
    if (pl == -1) { // PL番号が指定されていない場合
        // PL番号を0から999まで試す。正しいPL番号が見つかるまでのエラー出力をStderrSuppressorで抑制する
        sup.suppress(); // 標準エラー出力の抑制開始
        for (int pli = 0; pli < 1000; ++pli) {
            if (br.Begin(bvxxfile, pli, 0)) {
                sup.restore(); // 正しいPL番号が見つかったら標準エラー出力の抑制を解除
                pl = pli;

                vxx::HashEntry h;
                vxx::base_track_t b;

                while (br.NextHashEntry(h)) {
                    while (br.NextBaseTrack(b)) {
                        readBaseTrack(b);
                    }
                }
                br.End();
                readSuccess = true;
                break; // 成功したらループを抜ける
            }
        }
        if (!readSuccess) {
            std::cerr << "Error: Failed to load bvxx file. Please try using the -pl option." << std::endl;
            std::exit(1);
        }
    } else { // PL番号が指定されている場合
        if (br.Begin(bvxxfile, pl, 0)) {
            vxx::HashEntry h;
            vxx::base_track_t b;

            while (br.NextHashEntry(h)) {
                while (br.NextBaseTrack(b)) {
                    readBaseTrack(b);
                }
            }
            br.End();
        } else {
            std::cerr << "Error: Failed to load bvxx file. Please check if the correct -pl is specified." << std::endl;
            std::exit(1);
        }
    } // bvxxファイルの読み込み終了

    return tree;
}

void diff_hist(
    TCanvas *c1, std::string output, TTree *tree, const int Num, double *angle,
    double *sigma, const char* rl, int layer
) {
    double angmin, angmax, range;
    const char* branch_name = (rl == "Lateral") ? "dl" : "dr";
    TH1F** hist = new TH1F*[Num];
    TF1** gaus = new TF1*[Num];

    for (int i = 0; i < Num; ++i) {
        if (i < 14) {
            // tanθ 1.4 までは 0.1 刻み
            angmin = i * 0.1;
            angmax = angmin + 0.1;
        } else {
            // tanθ 1.4 からは 0.3 刻み
            angmin = 1.4 + (i - 14) * 0.3;
            angmax = angmin + 0.3;
        }
        angle[i] = (angmin + angmax) / 2;

        if (rl == "Lateral") {
            range = 0.05;
        } else if (angle[i] < 1.0) {
            range = angle[i] * 0.2 + 0.1;
        } else {
            range = angle[i] * 0.2;
        }

        TString name = Form("hist_%s_%d_%d", rl, layer, i);
        TString title = Form(
            "%s angle diff (%.1f^{ }#leq |tan#it{#theta}| < %.1f);%s angle difference [Layer%d];",
            rl, angmin, angmax, rl, layer
        );

        hist[i] = new TH1F(name, title, 100, -range, range);
        hist[i]->SetFillColorAlpha(layer == 1 ? 90 : 91, 0.4);
        hist[i]->SetLineColor(layer == 1 ? 90 : 91);

        tree->Draw(
            Form("%s%d >> %s", branch_name, layer, name.Data()),
            Form("%.1f <= tan && tan < %.1f && tan > 0.01", angmin, angmax),
            "goff"
        );

        if (hist[i]->GetEntries() >= 1000) {
            double min, max;
            double binmax = hist[i]->GetMaximum();
            GetThresholdMinMax(hist[i], 0.2 * binmax, min, max);
            gaus[i] = new TF1("gaus", "gaus", min, max);
            hist[i]->Fit(gaus[i], "q0", "", min, max);
            sigma[i] = gaus[i]->GetParameter(2);
        } else {
            sigma[i] = -1.0;
        }
    }

    for (int i = 0; i < Num; ++i) {
        if (i % 6 == 0) {
            c1->Divide(3, 2);
        }

        c1->cd(i % 6 + 1);
        hist[i]->Draw();
        gaus[i]->Draw("same");

        if (i % 6 == 5 || i == Num - 1) {
            c1->Print(output.c_str()); c1->Clear();
        }
    }
    for (int i = 0; i < Num; ++i) {
        if (hist[i]) delete hist[i];
        if (gaus[i]) delete gaus[i];
    }
    delete[] hist;
    delete[] gaus;
}

// thresholdを超える最小のビンのx座標とthresholdを下回らない最大のビンのx座標を取得する関数
void GetThresholdMinMax(TH1F *hist, double threshold, double &valMin, double &valMax)
{
    int Nbin = hist->GetNbinsX();
	int binMin = 0, binMax = 0;
	for (int i = 0; i < Nbin; i++) {
		if (hist->GetBinContent(i + 1) > threshold) {
			binMin = i;
            break;
		}
        if (i == Nbin - 1) binMin = 1;
	}
	for (int i = Nbin; i > 0; i--) {
		if (hist->GetBinContent(i) > threshold) {
			binMax = i;
			break;
		}
        if (i == 1) binMax = Nbin - 1;
	}
	valMin = hist->GetXaxis()->GetBinCenter(binMin);
	valMax = hist->GetXaxis()->GetBinCenter(binMax);
}
