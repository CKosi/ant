#include "calibration/DataManager.h"
#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "TH2.h"
#include "TH3.h"
#include "root-addons/analysis_codes/TwoPi0_MCSmearing_tools.h"
#include "TROOT.h"
#include "TRint.h"
#include "tree/TCalibrationData.h"
#include "calibration/modules/detail/TH2Storage.h"
#include "base/TH_ext.h"
#include <iostream>
#include <cstring>
#include "TDirectory.h"
#include "base/std_ext/string.h"

using namespace std;
using namespace ant;
using namespace  ant::calibration::gui;

TH2* GetHist(WrapTFileInput& file, const string& histname) {
    TH2* h = nullptr;
    file.GetObject(std_ext::formatter() <<  histname, h);
    if(!h) {
        LOG(ERROR) << histname << " not found in " << file.FileNames();
        exit(EXIT_FAILURE);
    }
    return h;
}

inline string ToUpper(const string& s) {
    string out = s;
    transform(s.begin(), s.end(),out.begin(), ::toupper);
    return out;
}

inline string ToLower(const string& s) {
    string out = s;
    transform(s.begin(), s.end(),out.begin(), ::tolower);
    return out;
}

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-mcsmearing", ' ', "0.1");
    auto cmd_verbose   = cmd.add<TCLAP::ValueArg<int>>   ("v","verbose", "Verbosity level (0..9)", false, 0,"level");
    auto cmd_batchmode = cmd.add<TCLAP::SwitchArg>       ("b","batch",   "Run in batch mode (no GUI, autosave)",false);
    auto cmd_setupname = cmd.add<TCLAP::ValueArg<string>>("s","setup",   "Setup name",       true, "", "setup");
    auto cmd_detector  = cmd.add<TCLAP::ValueArg<string>>("" ,"detector","Detector Name",    true, "", "detector");
    auto cmd_file      = cmd.add<TCLAP::ValueArg<string>>("" ,"file",    "Input file",       true, "", "file");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(!cmd_detector->isSet()) {
        LOG(ERROR) << "Detector not set";
        return EXIT_FAILURE;
    }

    const auto det = ToUpper(cmd_detector->getValue());

    // create TRint app early in order to have valid gStyle pointer...
    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    if(cmd_batchmode->isSet()) {
        fake_argv[fake_argc++] = strdup("-q");
    }
    auto app = new TRint("Ant-mcsmearing",&fake_argc,fake_argv,nullptr,0,true);


    const auto setup_name = cmd_setupname->getValue();
    auto setup = ExpConfig::Setup::Get(setup_name);
    if(setup == nullptr) {
        LOG(ERROR) << "Did not find setup instance for name " << setup_name;
        return 1;
    }

    auto manager = setup->GetCalibrationDataManager();
    manager->SetOverrideToDefault(true);

    const auto id = TID(0,0,{TID::Flags_t::MC});

    const string histname = std_ext::formatter() << "MCClusterECorr/h_EtrueErec_" << det;

    WrapTFileInput infile(cmd_file->getValue());
    const auto ecorr = GetHist(infile, histname);

    const string calName = std_ext::formatter() << det << "_ClusterECorr";

    TCalibrationData cdata(calName, id, id);
    calibration::detail::TH2Storage::Encode(ecorr, cdata);

    manager->Add(cdata,  Calibration::AddMode_t::AsDefault);


    app->Run(kTRUE);
    ExpConfig::Setup::Cleanup();
    setup = nullptr;
    manager = nullptr;

    return EXIT_SUCCESS;
}
