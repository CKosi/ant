#include "plot/root_draw.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObject.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "THStack.h"
#include <TVirtualPad.h>
#include <algorithm>

using namespace ant;
using namespace std;

const padoption::map_options_t padoption::map_options =
{
    {padoption_t::Legend, [] (TVirtualPad* p) {p->BuildLegend();} },
    {padoption_t::LogX, [] (TVirtualPad* p) {p->SetLogx();} },
    {padoption_t::LogY, [] (TVirtualPad* p) {p->SetLogy();} },
    {padoption_t::LogZ, [] (TVirtualPad* p) {p->SetLogz();} },
};


unsigned int ant::canvas::num = 0;

const ant::endcanvas endc;

TCanvas *ant::canvas::create(const string& title)
{
    stringstream s;
    s << "_canvas_" << setfill('0') << setw(3) << num++;
    name = s.str();

    TCanvas* c = new TCanvas(name.c_str(), title.c_str());
    return c;
}


TCanvas *ant::canvas::find()
{
    TObject* o = gROOT->FindObjectAny(name.c_str());
    TCanvas* c = dynamic_cast<TCanvas*>(o);
    if(c)
        return c;
    else
        return create();
}

canvas::canvas(const string &title) :
  name(), objs(), current_drawoption(), current_padoptions()
{
    create(title);
}

canvas::~canvas()
{

}

void canvas::cd()
{
    TCanvas* c = find();
    if(c) {
        c->cd();
    }
}

canvas &canvas::operator<<(root_drawable_traits &drawable)
{
    std::unique_ptr<root_drawable_traits> c(new drawable_container<root_drawable_traits*>(&drawable)) ;

    objs.emplace_back( move(c), current_drawoption, current_padoptions);

    return *this;
}

canvas &canvas::operator<<(TObject *hist)
{
    std::unique_ptr<root_drawable_traits> c(new drawable_container<TObject*>(hist)) ;

    objs.emplace_back( move(c), current_drawoption, current_padoptions);

    return *this;
}

canvas &canvas::operator<<(const endcanvas&)
{
    if(objs.empty()) {
        return *this;
    }

    TCanvas* c = find();

    if(c) {

        const int cols = ceil(sqrt(objs.size()));
        const int rows = ceil((double)objs.size()/(double)cols);

        c->Divide(cols,rows);
        int pad=1;
        for(const auto& o : objs) {
            TVirtualPad* vpad = c->cd(pad++);
            // draw the object
            get<0>(o)->Draw(get<1>(o).c_str());
            // set pad options
            for(const auto& o_ : get<2>(o)) {
                const auto& it = padoption::map_options.find(o_);
                // silently ignore not implemented pad options
                if(it == padoption::map_options.end())
                    continue;
                it->second(vpad);
            }
        }
    }
    return *this;
}

canvas &canvas::operator<<(const drawoption &c)
{
    current_drawoption = c.Option();
    return *this;
}

canvas &canvas::operator<<(const padoption::set &c)
{
    const auto& o = c.Option();
    auto it = std::find(current_padoptions.begin(), current_padoptions.end(), o);
    if(it == current_padoptions.end()) {
        current_padoptions.emplace_back(o);
    }
    return *this;
}

canvas &canvas::operator<<(const padoption::unset &c)
{
    const auto& o = c.Option();
    current_padoptions.remove(o);
    return *this;
}


canvas &canvas::operator>>(const string &filename)
{
    TCanvas* c = find();
    if(c) {
        c->SaveAs(filename.c_str());
    }
    return *this;
}


//************ hstack ************


hstack::hstack(const string &name, const std::string &title):
    stack(new THStack(name.c_str(),title.c_str())),
    current_option("")
{}

hstack::~hstack()
{}

hstack &hstack::operator<<(TH1D *hist)
{
    stack->Add(hist, current_option.c_str());
    return *this;
}

hstack &hstack::operator<<(const drawoption &c)
{
    current_option = c.Option();
    return *this;
}

TObject *hstack::GetObject()
{
    return stack;
}

void hstack::Draw(const string &option) const
{
    stack->Draw(option.c_str());
}

const std::vector<Color_t> ant::ColorPalette::Colors = {kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kPink+9, kSpring+10, kGray};
