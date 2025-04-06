#include <cmath>
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

    // Read the BVXX file
    std::cout << "Reading BVXX file... " << std::endl;
    vxx::BvxxReader br;
    std::vector<vxx::base_track_t> btvec = br.ReadAll(bvxxfile, pl_i, 0);
    if (btvec.empty()) {
        std::cerr << "Error: Failed to read any data from the input file." << std::endl;
        return 1;
    }
    size_t entries = btvec.size();
    std::cout << "Read " << entries << " tracks." << std::endl;

    // Create the ROOT file and TTree
    TFile* fout = new TFile(outputfile.c_str(), "RECREATE");
    TTree* nt = new TTree("nt", "");

    int RawId, RawId1, RawId2, isg, isg1, isg2;
    double pl, ax, ay, x, y;
    double ph1, ax1, ay1, x1, y1, z1, pos1, col1, row1, zone1;
    double ph2, ax2, ay2, x2, y2, z2, pos2, col2, row2, zone2;

    // Create branches for the tree
    const std::vector<std::pair<std::string, void*>> intBranches = {
        {"RawId", &RawId},
        {"RawId1", &RawId1},
        {"RawId2", &RawId2},
        {"isg", &isg},
        {"isg1", &isg1},
        {"isg2", &isg2}
    };
    const std::vector<std::pair<std::string, void*>> doubleBranches = {
        {"pl", &pl},
        {"ax", &ax},
        {"ay", &ay},
        {"x", &x},
        {"y", &y},
        {"ph1", &ph1},
        {"ax1", &ax1},
        {"ay1", &ay1},
        {"x1", &x1},
        {"y1", &y1},
        {"z1", &z1},
        {"pos1", &pos1},
        {"col1", &col1},
        {"row1", &row1},
        {"zone1", &zone1},
        {"ph2", &ph2},
        {"ax2", &ax2},
        {"ay2", &ay2},
        {"x2", &x2},
        {"y2", &y2},
        {"z2", &z2},
        {"pos2", &pos2},
        {"col2", &col2},
        {"row2", &row2},
        {"zone2", &zone2}
    };
    for (const auto& branch : intBranches) {
        nt->Branch(branch.first.c_str(), branch.second, (branch.first + "/I").c_str());
    }
    for (const auto& branch : doubleBranches) {
        nt->Branch(branch.first.c_str(), branch.second, (branch.first + "/D").c_str());
    }

    // Fill the tree with data
    int cnt = 0;
    for (const auto& bt : btvec) {
        RawId = (int)bt.rawid;
        pl = (double)bt.pl;
        isg = (int)bt.isg;
        ax = (double)bt.ax;
        ay = (double)bt.ay;
        x = (double)bt.x;
        y = (double)bt.y;

        ph1 = (double)bt.m[0].ph;
        ax1 = (double)bt.m[0].ax;
        ay1 = (double)bt.m[0].ay;
        x1 = 0.0;
        y1 = 0.0;
        z1 = (double)bt.m[0].z;
        pos1 = (double)bt.m[0].pos;
        col1 = (double)bt.m[0].col;
        row1 = (double)bt.m[0].row;
        zone1 = (double)bt.m[0].zone;
        RawId1 = (int)bt.m[0].rawid;
        isg1 = (int)bt.m[0].isg;

        ph2 = (double)bt.m[1].ph;
        ax2 = (double)bt.m[1].ax;
        ay2 = (double)bt.m[1].ay;
        x2 = 0.0;
        y2 = 0.0;
        z2 = (double)bt.m[1].z;
        pos2 = (double)bt.m[1].pos;
        col2 = (double)bt.m[1].col;
        row2 = (double)bt.m[1].row;
        zone2 = (double)bt.m[1].zone;
        RawId2 = (int)bt.m[1].rawid;
        isg2 = (int)bt.m[1].isg;

        // Fill the tree with the data
        nt->Fill();

        if (cnt % 10000 == 0) {
            std::cout << "Filled " << cnt << " / " << entries << "\r" << std::flush;
        }
        ++cnt;
    }

    // Write the tree to the file and close it
    nt->AutoSave();
    fout->Close();

    std::cout << "Filled " << cnt << " / " << entries << std::endl;
    std::cout << "Finished." << std::endl;

    return 0;
}