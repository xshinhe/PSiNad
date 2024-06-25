#include "kids/VariableDescriptor.h"

#include <iostream>
#include <regex>

namespace PROJECT_NS {

VariableDescriptor::VariableDescriptor(const std::string& token_string) : tokenString{token_string} {
    const std::regex pattern("([^{<:]*)(?:\\{([^@<]*)(?:@(.*?))?\\})?(?:<([^>]*)>)?(?::([^:]*))?");
    std::smatch      match;
    if (std::regex_match(tokenString, match, pattern)) {
        std::tie(name, field, time, index, type) = std::make_tuple(  //
            match[1], match[2], match[3], match[4], match[5]);
    } else {
        throw kids_error("cannot match variable pattern!");
    }

    if (field == "0") field = "init";
    if (field == "I") field = "integrator";
    if (field == "M") field = "model";
    if (field == "R") field = "result";
    if (field == "") field = "integrator";

    if (time == "" || time == "T") time = "t";

    if (index == "") {};  // do nothing

    if (type == "R") {
        dataType = kids_real_type;
    } else if (type == "C") {
        dataType = kids_complex_type;
    } else {
        dataType = kids_void_type;
    }
}

void VariableDescriptor::referIn(std::shared_ptr<DataSet>& DS) {
    if (field == "" || name == "") throw kids_error("bad key!");
    std::string key                        = utils::concat(field, ".", name);
    std::tie(dataType, dataPointer, shape) = DS->obtain(key);
    // #define LOCAL_DEBUG
    // #ifdef LOCAL_DEBUG
    //     std::cout << LOC() << "VariableDescriptor Data:\n"
    //               << ".tokenString = " << tokenString << "\n"  //
    //               << ".name = " << name << "\n"                //
    //               << ".field = " << field << "\n"              //
    //               << ".index = " << index << "\n"              //
    //               << ".type = " << type << "\n"                //
    //               << ".time = " << time << "\n"                //
    //               << ".dataType = " << dataType << "\n"        //
    //               << ".dataPointer = " << dataPointer << "\n"  //
    //               << ".shape = " << ((shape) ? shape->to_string() : "") << "\n";
    // #endif  // LOCAL_DEBUG
}

void VariableDescriptor::defineIn(std::shared_ptr<DataSet>& DS, kids_dtype data_type,
                                  const std::vector<std::size_t>& cxxshape) {
    if (field == "" || name == "") throw kids_error("bad key!");
    std::string key = utils::concat(field, ".", name);
    switch (data_type) {
        case kids_real_type:
            DS->def_real(key, cxxshape, "CUSTOM DEFINED");
            break;
        case kids_complex_type:
            DS->def_complex(key, cxxshape, "CUSTOM DEFINED");
            break;
    }
    std::tie(dataType, dataPointer, shape) = DS->obtain(key);
}

};  // namespace PROJECT_NS