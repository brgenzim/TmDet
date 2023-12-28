#ifndef __UNITMP_PDBLIB_TYPES_PDB_METHOD_TYPES__
#define __UNITMP_PDBLIB_TYPES_PDB_METHOD_TYPES__

#include <unordered_map>

namespace UniTmp::PdbLib::Types {

    struct PdbMethodType {
        const char *code;
        const char *name;
    };

    namespace PdbMethodTypes {
        const PdbMethodType Model {"Model","Model"};
        const PdbMethodType EM {"EM","Electron microscopy"};
        const PdbMethodType Xray {"Xray","X-ray diffraction"};
        const PdbMethodType Eray {"Eray", "Electron crystallography"};
        const PdbMethodType NMR {"NMR", "Solution NMR"};
        const PdbMethodType Fiber {"Fiber", "Fiber diffraction"};
        const PdbMethodType IM {"IM", "Infrared microscopy"};
        const PdbMethodType IR {"IR", "Infrared spectroscopy"};
        const PdbMethodType Neutron {"Neutron", "Neutron diffraction"};
        const PdbMethodType Solution {"Solution", "Solution scattering"};
        const PdbMethodType Solid {"Solid", "Solid-state NMR"};
        const PdbMethodType Synchroton {"Synchroton", "Synchroton"};
        const PdbMethodType Other {"Other", "Other"};
        const PdbMethodType Powder {"Powder", "Powder diffraction"};

        const std::unordered_map<const char*, PdbMethodType> all {
            {"Model",Model},
            {"EM",EM},
            {"Xray",Xray},
            {"Eray", Eray},
            {"NMR", NMR},
            {"Fiber", Fiber},
            {"IM", IM},
            {"IR", IR},
            {"Neitron", Neutron},
            {"Solution", Solution},
            {"Solid", Solid},
            {"Synchroton", Synchroton},
            {"Other", Other},
            {"Powder", Powder},
        };

        const std::unordered_map<const char*, PdbMethodType> fromName {
            {"Model",Model},
            {"Electron microscopy",EM},
            {"X-ray diffraction",Xray},
            {"Electron crystallography", Eray},
            {"Solution NMR", NMR},
            {"Fiber diffraction", Fiber},
            {"Infrared microscopy", IM},
            {"Infrared spectroscopy", IR},
            {"Neutron diffraction", Neutron},
            {"Solution scattering", Solution},
            {"Solid-state NMR", Solid},
            {"Synchroton", Synchroton},
            {"Other", Other},
            {"Powder diffraction", Powder},
        };
    }

}

#endif
