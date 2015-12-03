#include "Manager.h"

#include "ManagerWindow.h"

#include "calibration/gui/CalCanvas.h"

#include "tree/TDataRecord.h"
#include "tree/TAntHeader.h"

#include "base/interval.h"
#include "base/std_ext/misc.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TH2D.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

Manager::Manager(const std::vector<std::string>& inputfiles, unsigned avglength):
    buffer(avglength),
    state()
{
    BuildInputFiles(inputfiles);
}

Manager::~Manager()
{

}

void Manager::InitGUI(ManagerWindow* window_) {
    window = window_;
    module->InitGUI(window);

    if(buffer.GetSumLength()==0)
        window->SetProgressMax(1, nChannels-1);
    else
        window->SetProgressMax(input_files.size(), nChannels-1);
}


bool Manager::input_file_t::operator <(const Manager::input_file_t& o) const {
    return range.Start() < o.range.Start();
}

void Manager::BuildInputFiles(const vector<string>& filenames)
{
    if(filenames.empty())
        return;

    for(const auto& filename : filenames) {

        try {

            WrapTFileInput file;
            file.OpenFile(filename);

            TAntHeader* header = nullptr;
            file.GetObject("AntHeader", header);

            if(!header) {
                LOG(WARNING) << "No TAntHeader found in " << filename;
                continue;
            }

            if(SetupName.empty()) {
                SetupName = header->SetupName;
            }
            else if(SetupName != header->SetupName) {
                LOG(WARNING) << "Previously found setup name '" << SetupName
                             << "' does not match '" << header->SetupName << "' of file "
                             << filename;
                continue;
            }

            auto range = interval<TID>(header->FirstID, header->LastID);
            if(!range.IsSane()) {
                LOG(WARNING) << "Range " << range << " not sane in " << filename;
                continue;
            }
            input_files.emplace_back(filename, range);

        } catch (const std::runtime_error& e) {
            LOG(WARNING) << "Can't open " << filename << " " << e.what();
        }
    }

    input_files.sort();
    LOG(INFO) << "Loaded " << input_files.size()
              << " files from " << filenames.size() << " provided filenames ("
              << 100*(double)input_files.size()/filenames.size() << " %)";
    if(!input_files.empty()) {
        LOG(INFO) << "Files TID range from "
                  << input_files.front().range.Start()
                  << " to "
                  << input_files.back().range.Stop();
    }
}



void Manager::FillBufferFromFiles()
{
    while(buffer.Empty() && state.it_file != input_files.end()) {
        const input_file_t& file_input = *state.it_file;
        try
        {
            WrapTFileInput file;
            file.OpenFile(file_input.filename);

            auto hist = file.GetSharedHist(module->GetHistogramName());

            if(!hist) {
                LOG(WARNING) << "Histogram " << module->GetHistogramName() << " not found in " << file_input.filename;
            } else {
                LOG(INFO) << "Buffer filled with " << file_input.filename;
                buffer.Push(hist, file_input.range);
            }

        }
        catch (const std::runtime_error& e) {
            LOG(WARNING) << "Can't open " << file_input.filename << ": " << e.what();
        }
        state.it_file++;
    }

    if(state.it_file == input_files.end()) {
        VLOG(7) << "Reached end of files, processing remaining buffer";
        buffer.Finish();
    }
}

bool Manager::DoInit(int gotoSlice_)
{
    // use the sign of gotoSlice_ to detect if gotoSlice
    // is active at all
    unsigned gotoSlice = 0;
    if(gotoSlice_ >= 0) {
        gotoSlice = gotoSlice_;
        state.oneslice = true;
    }

    if(state.oneslice && gotoSlice>=input_files.size()) {
        LOG(ERROR) << "Requested slice " << gotoSlice << " not smaller than number of input files " << input_files.size();
        return false;
    }

    state.channel = -1;
    state.slice = 0;
    state.it_file = input_files.begin();

    nChannels = module->GetNumberOfChannels();
    if(nChannels==0) {
        LOG(ERROR) << "Module reports zero channels, nothing to do then.";
        return false;
    }

    FillBufferFromFiles();
    if(buffer.Empty()) {
        LOG(ERROR) << "Could not initially fill the buffer from given files";
        return false;
    }

    while(state.slice < gotoSlice) {
        state.slice++;
        buffer.GotoNextID();
        FillBufferFromFiles();
        if(buffer.Empty()) {
            LOG(ERROR) << "Could not move buffer to slice " << gotoSlice << ", stopped at " << state.slice;
            return false;
        }
    }

    module->StartSlice(buffer.CurrentID());
    return true;
}


Manager::RunReturn_t Manager::Run()
{
    // this statement is executed once the class goes out-of-scope
    std_ext::execute_on_destroy setProgress([this] () {
        window->SetProgress(state.slice, state.channel);
    });


    if(!state.breakpoint_finish && state.channel < nChannels && state.channel >= 0) {
        bool noskip = true;
        if(!state.breakpoint_fit) {

            const auto ret = module->DoFit(buffer.CurrentSum(), state.channel);
            noskip = ret != CalibModule_traits::DoFitReturn_t::Skip;

            if(ret == CalibModule_traits::DoFitReturn_t::Display
               || (!window->Mode.autoContinue && noskip)
               ) {
                VLOG(7) << "Displaying Fit...";
                module->DisplayFit();
                state.breakpoint_fit = true;
                return RunReturn_t::Wait;
            }
            else if(window->Mode.showEachFit && noskip) {
                module->DisplayFit();
            }
        }
        if(noskip)
            module->StoreFit(state.channel);
        state.breakpoint_fit = false;
    }

    if(state.breakpoint_finish
       || (state.channel >= nChannels && window->Mode.gotoNextSlice))
    {
        if(!state.breakpoint_finish) {
            VLOG(7) << "Finish module";
            state.breakpoint_finish = module->FinishSlice();
            if(!window->Mode.autoFinish && state.breakpoint_finish) {
                VLOG(7) << "Displaying finished range...";
                window->SetFinishMode(true);
                return RunReturn_t::Wait;
            }
        }
        module->StoreFinishSlice(buffer.CurrentID());
        window->SetFinishMode(false);

        if(state.oneslice) {
            LOG(INFO) << "Finished processing this one slice";
            return RunReturn_t::Exit;
        }

        state.breakpoint_finish = false;
        state.channel = 0;
        state.slice++;
        buffer.GotoNextID();

        // try refilling the worklist
        FillBufferFromFiles();
        if(buffer.Empty()) {
            /// \todo give module a chance to do something again here???
            LOG(INFO) << "Finished processing whole buffer";
            return RunReturn_t::Exit;
        }

        module->StartSlice(buffer.CurrentID());

    }
    else
    {
        // check if some specific channel is requested
        // then set the channel appropiately
        if(window->Mode.requestChannel<0) {
            state.channel += window->Mode.channelStep;
        }
        else {
            state.channel = window->Mode.requestChannel;
            window->Mode.requestChannel = -1; // request is handled...
        }

        // check if we run out of the range
        // stop running if we do so
        if(state.channel<0) {
            state.channel = 0;
            if(window->Mode.autoContinue)
                return RunReturn_t::Wait;
        }
        else if(state.channel>=nChannels && !window->Mode.gotoNextSlice) {
            state.channel = module->GetNumberOfChannels()-1;
            if(window->Mode.autoContinue)
                return RunReturn_t::Wait;
        }
    }

    // continue running by default
    return RunReturn_t::Continue;
}


