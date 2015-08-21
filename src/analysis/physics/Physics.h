#pragma once

#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"

#include "base/std_ext.h"

#include <list>
#include <string>
#include <memory>
#include <functional>

class TFile;
class TDirectory;

namespace ant {

class TAntHeader;

namespace analysis {

class Physics {
private:
    std::string name_;
protected:
    SmartHistFactory HistFac;
public:
    Physics(const std::string& name);
    virtual ~Physics() {}
    virtual void ProcessEvent(const data::Event& event) =0;
    virtual void Finish() =0;
    virtual void ShowResult() =0;
    std::string GetName() { return name_; }
};

template<class T>
std::unique_ptr<Physics> physics_factory()
{
    return std::move(std_ext::make_unique<T>());
}

using physics_creator = std::function<std::unique_ptr<Physics>()>;

class PhysicsRegistry
{
private:
    using physics_creators_t = std::map<std::string, physics_creator>;
    physics_creators_t physics_creators;

public:
    static PhysicsRegistry& get();

    static std::unique_ptr<Physics> Create(const std::string& name);

    std::vector<std::string> GetList() const;

    void RegisterPhysics(physics_creator c, const std::string& name) {
        physics_creators[name] = c;
    }

    static void PrintRegistry();

};

class PhysicsRegistration
{
public:
    PhysicsRegistration(physics_creator c, const std::string& name);
};

#define AUTO_REGISTER_PHYSICS(physics, name) \
    ant::analysis::PhysicsRegistration _physics_registration_ ## physics(ant::analysis::physics_factory<physics>,name);

}} // nammespace ant::analysis
