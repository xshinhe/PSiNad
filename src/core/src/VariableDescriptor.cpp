#include "kids/VariableDescriptor.h"

#include <iostream>
#include <regex>

#include "kids/debug_utils.h"

namespace PROJECT_NS {

VariableDescriptor::VariableDescriptor(const std::string& token_string, const std::string& save, bool is_output)
    : tokenString{token_string}, save{save}, isOutput{is_output} {
    const std::regex pattern("([^{<:]*)(?:\\{([^@<]*)(?:@(.*?))?\\})?(?:<([^>]*)>)?(?::([^:]*))?");
    std::smatch      match;
    if (std::regex_match(tokenString, match, pattern)) {
        std::tie(name, field, time, index, type) = std::make_tuple(  //
            match[1], match[2], match[3], match[4], match[5]);
    } else {
        throw kids_error(utils::concat("cannot match variable pattern!", tokenString));
    }
    // process field and parse keyRaw/keyRec
    if (field == "" || field == "I") field = "integrator";
    if (field == "M") field = "model";
    if (field == "P") field = "parameter";
    if (field == "R") field = "record";
    if (index == "") {};  // do nothing
    dataType = kids_void_type;
    if (type == "R") dataType = kids_real_type;
    if (type == "C") dataType = kids_complex_type;

    isTabular = isOutput;
    if (isOutput && name[0] == '*') {
        isTabular = false;
        name      = name.substr(1, name.size());
    }

    keyRaw = utils::concat(field, ".", name);
    keyRec = keyRaw, keyRes0 = keyRaw, keyRes1 = keyRaw, keyRes2 = keyRaw;
    if (isOutput) {
        keyRec  = utils::concat("record.", name);
        keyRes0 = utils::concat("_.0.", field, ".", name);
        keyRes1 = utils::concat("_.1.", field, ".", name);
        keyRes2 = utils::concat("_.2.", field, ".", name);
    } else {
        keyRec = (time == "") ? keyRaw : utils::concat("@.", time, ".", field, ".", name);
    }
}

void VariableDescriptor::defineIn(std::shared_ptr<DataSet> DS, kids_dtype data_type,
                                  const std::vector<std::size_t>& cxxshape, std::size_t totalFrameNumber) {
    if (isOutput) {
        // keyRec cannot be defined twice
        // if (DS->haskey(keyRec)) throw kids_error(utils::concat("conflict of : ", keyRec));
        dataPointerRaw = nullptr;

        std::vector<std::size_t> cxxstackedshape;
        cxxstackedshape.push_back(totalFrameNumber);
        for (auto& d : cxxshape) cxxstackedshape.push_back(d);

        switch (data_type) {
            case kids_real_type:
                DS->def_real(keyRec, cxxshape, "CUSTOM DEFINED");
                if (isTabular) {
                    DS->def_real(keyRes0, cxxstackedshape, "CUSTOM DEFINED");
                    DS->def_real(keyRes1, cxxstackedshape, "CUSTOM DEFINED");
                    DS->def_real(keyRes2, cxxstackedshape, "CUSTOM DEFINED");
                }
                break;
            case kids_complex_type:
                DS->def_complex(keyRec, cxxshape, "CUSTOM DEFINED");
                if (isTabular) {
                    DS->def_complex(keyRes0, cxxstackedshape, "CUSTOM DEFINED");
                    DS->def_complex(keyRes1, cxxstackedshape, "CUSTOM DEFINED");
                    DS->def_complex(keyRes2, cxxstackedshape, "CUSTOM DEFINED");
                }
                break;
        }
        std::tie(dataType, dataPointerTrace, shape) = DS->obtain(keyRec);
        if (isTabular) {
            std::tie(dataType, dataPointerRes0, stackedshape) = DS->obtain(keyRes0);
            std::tie(dataType, dataPointerRes1, stackedshape) = DS->obtain(keyRes1);
            std::tie(dataType, dataPointerRes2, stackedshape) = DS->obtain(keyRes2);
        }
    } else {
        // keyRaw must be in DataSet
        if (!DS->haskey(keyRaw)) throw kids_error(utils::concat("loss of : ", keyRaw));
        std::tie(dataType, dataPointerRaw, shape) = DS->obtain(keyRaw);
        switch (dataType) {
            case kids_real_type:
                dataPointerTrace = (void*) DS->def_real(keyRec, *shape, "CUSTOM DEFINED");
                break;
            case kids_complex_type:
                dataPointerTrace = (void*) DS->def_complex(keyRec, *shape, "CUSTOM DEFINED");
                break;
        }
        stackedshape = nullptr;
    }
    // #define LOCAL_DEBUG
    // #undef LOCAL_DEBUG
    // #ifdef LOCAL_DEBUG
    //     std::cout << LOC() << "VariableDescriptor Data After:\n"
    //               << ".tokenString = " << tokenString << "\n"       //
    //               << ".name = " << name << "\n"                     //
    //               << ".field = " << field << "\n"                   //
    //               << ".index = " << index << "\n"                   //
    //               << ".type = " << type << "\n"                     //
    //               << ".time = " << time << "\n"                     //
    //               << ".dataType = " << dataType << "\n"             //
    //               << ".dataPointer = " << dataPointerRaw << "\n"    //
    //               << ".dataPointer = " << dataPointerTrace << "\n"  //
    //               << ".dataPointer = " << dataPointerRes0 << "\n"   //
    //               << ".shape = " << ((shape) ? shape->to_string() : "") << "\n";
    // #endif  // LOCAL_DEBUG
}

void VariableDescriptor::checkTrace(int sampleIndex) {
    if (time == "" || isOutput) return;
    if (stoi(time) == sampleIndex) {
        if (dataType == kids_real_type) {
            kids_real* fromdata = (kids_real*) (dataPointerRaw);
            kids_real* todata   = (kids_real*) (dataPointerTrace);
            for (int i = 0; i < shape->size(); ++i) todata[i] = fromdata[i];
        }
        if (dataType == kids_complex_type) {
            kids_complex* fromdata = (kids_complex*) (dataPointerRaw);
            kids_complex* todata   = (kids_complex*) (dataPointerTrace);
            for (int i = 0; i < shape->size(); ++i) todata[i] = fromdata[i];
        }
    }
};

};  // namespace PROJECT_NS