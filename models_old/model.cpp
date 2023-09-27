#include "model.h"

Model::Model(const Param& iparm) {
    parm = iparm;  // a deepcopy of iparm

    // try to update param
    if (Param_IfHaveKey(iparm, "set_id")) {
        std::string set_id   = Param_GetT(std::string, iparm, "set_id");
        std::string set_file = Param_GetT(std::string, iparm, "set_file");

        Param_ParseFile(parm, set_file);
        parm = parm[set_id];
        tag  = set_id;
    } else {
        tag = "default";
    }
    iou = parse_unit_parm(parm);
};

Model::~Model(){};

IOUnit Model::parse_unit_parm(const Param& parm) {
    if (!Param_IfHaveKey(parm, "unit")) return IOUnit();

    IOUnit parsed_iou;
    auto uparm = parm["unit"];

    /**
     * @warning you shouldn't change system units in simulations.
     * if you want to use original hartree unit, please refer phys.h and:
     * ```
     *      namespace au = phys::hartree; // only set {hb, ke, c, me} = 1
     * ```
     * else if you want to use generalized hartree units, please also see phys.h
     * ```
     *      namespace au = phys::ghartree; // set {hb, ke, c, me, kB, NA} = 1
     * ```
     * but DO NOT change namespace here!
     */
    namespace used_us = phys::au;

    switch (uparm.type()) {
        /**
         * @first-case, you can pass a list-style node like:
         * ```
         *    "unit": ["fs", "1 Angs", "amu", "J", "K"] // use time unit = 1 fs and so on
         * ```
         */
        case configor::config_value_type::array: {  // list initialization
            for (auto& term : uparm) {
                auto u = phys::us::parse(term.get<std::string>());
                if (u.dim == phys::time_d) parsed_iou.time = phys::us::conv(used_us::unit, u);
                if (u.dim == phys::length_d) parsed_iou.leng = phys::us::conv(used_us::unit, u);
                if (u.dim == phys::mass_d) parsed_iou.mass = phys::us::conv(used_us::unit, u);
                if (u.dim == phys::energy_d) parsed_iou.ener = phys::us::conv(used_us::unit, u);
                if (u.dim == phys::thermodynamic_temperature_d) parsed_iou.temp = phys::us::conv(used_us::unit, u);
            }
            break;
        }
        /**
         * @second-case, you can pass a dict-style node like:
         * ```
         *    "unit": {"time":"1 fs", "leng": "1 Bohr", "mass": "au", "ener": "Hart"]
         * ```
         */
        case configor::config_value_type::object: {  // dict initialization
            parsed_iou.leng = phys::us::conv(used_us::unit, Param_GetT(std::string, uparm, "leng", "1"));
            parsed_iou.time = phys::us::conv(used_us::unit, Param_GetT(std::string, uparm, "time", "1"));
            parsed_iou.mass = phys::us::conv(used_us::unit, Param_GetT(std::string, uparm, "mass", "1"));
            parsed_iou.ener = phys::us::conv(used_us::unit, Param_GetT(std::string, uparm, "ener", "1"));
            parsed_iou.temp = phys::us::conv(used_us::unit, Param_GetT(std::string, uparm, "temp", "1"));
            break;
        }
        default:
            LOG(FATAL);
    }
    return parsed_iou;
}
