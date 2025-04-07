#include <cmath>
#include <sstream>
#include <VxxReader.h>
#include <KeywordArgs.h>
#include <Template.h>
#include <netscan_data_types_ui.h>
#include <TFile.h>
#include <TTree.h>

int main(int argc, char* argv[]) {
    // Check if the correct number of arguments is provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " basetrack.vxx pl [output.root]" << std::endl;
        return 1;
    }

    // Parse command line arguments
    std::string bvxxfile = argv[1];
    if (bvxxfile.substr(bvxxfile.find_last_of('.') + 1) != "vxx") {
        std::cerr << "Error: Input file must have a .vxx extension." << std::endl;
        return 1;
    }
    int pl_i = std::stoi(argv[2]);
    std::string outputfile = (argc == 4) ? argv[3] : bvxxfile.substr(0, bvxxfile.find_last_of('.')) + ".root";
    if (outputfile.substr(outputfile.find_last_of('.') + 1) != "root") {
        std::cerr << "Error: Output file must have a .root extension." << std::endl;
        return 1;
    }

    // Create the ROOT file and TTree
    TFile* fout = new TFile(outputfile.c_str(), "RECREATE");
    TTree* nt = new TTree("nt", "");

    int64_t rawid, rawid1, rawid2;
    int pl, isg1, ph1, pos1, col1, row1, zone1, isg2, ph2, pos2, col2, row2, zone2;
    double ax, ay, x, y;
    double ax1, ay1, z1;
    double ax2, ay2, z2;

    // Create branches for the tree
    const std::vector<std::pair<std::string, void*>> intBranches = {
        {"pl", &pl},
        {"rawid", &rawid},
        {"rawid1", &rawid1},
        {"rawid2", &rawid2},
        {"isg1", &isg1},
        {"ph1", &ph1},
        {"pos1", &pos1},
        {"col1", &col1},
        {"row1", &row1},
        {"zone1", &zone1},
        {"isg2", &isg2},
        {"ph2", &ph2},
        {"pos2", &pos2},
        {"col2", &col2},
        {"row2", &row2},
        {"zone2", &zone2}
    };
    const std::vector<std::pair<std::string, void*>> doubleBranches = {
        {"ax", &ax},
        {"ay", &ay},
        {"x", &x},
        {"y", &y},
        {"ax1", &ax1},
        {"ay1", &ay1},
        {"z1", &z1},
        {"ax2", &ax2},
        {"ay2", &ay2},
        {"z2", &z2}
    };
    for (const auto& branch : intBranches) {
        nt->Branch(branch.first.c_str(), branch.second, (branch.first + "/I").c_str());
    }
    for (const auto& branch : doubleBranches) {
        nt->Branch(branch.first.c_str(), branch.second, (branch.first + "/D").c_str());
    }

    // Process each base track and fill the tree
    vxx::BvxxReader br;
    int cnt = 0;
    std::ostringstream oss;
    if (br.Begin(bvxxfile, pl_i, 0))
    {
        vxx::HashEntry h;
        vxx::base_track_t b;

        while (br.NextHashEntry(h))
        {
            while (br.NextBaseTrack(b))
            {
                rawid = (int64_t)b.rawid;
                pl = b.pl;
                ax = b.ax;
                ay = b.ay;
                x = b.x;
                y = b.y;

                ph1 = b.m[0].ph;
                ax1 = b.m[0].ax;
                ay1 = b.m[0].ay;
                z1 = b.m[0].z;
                pos1 = b.m[0].pos;
                col1 = b.m[0].col;
                row1 = b.m[0].row;
                zone1 = b.m[0].zone;
                rawid1 = b.m[0].rawid;
                isg1 = b.m[0].isg;

                ph2 = b.m[1].ph;
                ax2 = b.m[1].ax;
                ay2 = b.m[1].ay;
                z2 = b.m[1].z;
                pos2 = b.m[1].pos;
                col2 = b.m[1].col;
                row2 = b.m[1].row;
                zone2 = b.m[1].zone;
                rawid2 = b.m[1].rawid;
                isg2 = b.m[1].isg;

                // Fill the tree with the data
                nt->Fill();

                ++cnt;
                if (cnt % 10000 == 0) {
                    oss << " Filled " << cnt << " tracks\r";
                    std::cout << oss.str() << std::flush;
                    oss.str(""); // Clear the string stream
                    oss.clear();
                }
            }
        }
        br.End();
    }

    // Save the tree to the file and close the file
    nt->AutoSave();
    fout->Close();

    std::cout << " Filled " << cnt << " tracks" << std::endl;
    std::cout << " Finished." << std::endl;

    return 0;
}