#include <TROOT.h>
#include <TStyle.h>

void logon()
{
    // font
    int fontid = 42;                         // 42=Helvetica, 132=Times
    gStyle->SetStatFont(fontid);               // 統計box内
    gStyle->SetLabelFont(fontid, "xyz");       // 軸ラベル（数値）
    gStyle->SetTitleFont(fontid, "xyz");       // 軸title
    gStyle->SetTitleFont(fontid, "");          // title
    gStyle->SetTextFont(fontid);               // text
    gStyle->SetLegendFont(fontid);             // 凡例

    // redefine some colors
    gROOT->GetColor(0)->SetRGB( 255./255., 255./255., 255./255.);  // white
    gROOT->GetColor(1)->SetRGB(   0./255.,   0./255.,   0./255.);  // black
    gROOT->GetColor(89)->SetRGB( 255./255.,  75./255.,   0./255.); // red
    gROOT->GetColor(90)->SetRGB(   0./255.,  90./255., 255./255.); // blue
    gROOT->GetColor(91)->SetRGB( 246./255., 170./255.,   0./255.); // orange
    gROOT->GetColor(92)->SetRGB( 191./255., 228./255., 255./255.); // light blue
    gROOT->GetColor(93)->SetRGB( 255./255., 202./255., 128./255.); // light orange
    gROOT->GetColor(94)->SetRGB( 200./255., 200./255., 203./255.); // light gray

    // main color
    int colorid = 1;
    gStyle->SetAxisColor(colorid, "xyz");        // 軸
    gStyle->SetLabelColor(colorid, "xyz");       // 軸ラベル（数値）
    gStyle->SetTitleColor(colorid, "xyz");       // 軸title
    gStyle->SetTitleTextColor(colorid);          // title
    gStyle->SetFrameLineColor(colorid);          // 描画エリアの枠
    gStyle->SetStatTextColor(colorid);           // 統計box内のtext
    gStyle->SetLineColor(colorid);               // 統計boxの枠など
    gStyle->SetTextColor(colorid);               // text
	gStyle->SetFuncColor(colorid);				 // fitting function

    // grid
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetGridColor(94);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // margin and offset
    gStyle->SetPadRightMargin(0.1);      // Pad右側のマージン
    gStyle->SetPadLeftMargin(0.1);       // Pad左側のマージン
    gStyle->SetPadTopMargin(0.1);        // Pad上側のマージン
    gStyle->SetPadBottomMargin(0.11);    // Pad下側のマージン
    gStyle->SetLabelOffset(0.008,"xyz"); // 軸ラベル（数値）と軸の距離
    gStyle->SetTitleOffset(1.1,"xyz");   // 軸titleと軸の距離
    gStyle->SetTitleY(0.985);

    // histogram
    gStyle->SetHistFillColor(92);   // ヒストグラム内部の色
    gStyle->SetHistLineColor(1);    // ヒストグラムの線の色
    gStyle->SetHistFillStyle(1011); // ヒストグラム内部のパターン 1011=塗りつぶし
    gStyle->SetHistLineStyle(1);    // ヒストグラムの線種 1=直線
    gStyle->SetHistLineWidth(1);    // ヒストグラムの線幅 pixel

    // split z color
    gStyle->SetNumberContours(256); // z軸を256段階に色分け
}
