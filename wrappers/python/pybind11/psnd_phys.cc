py::module phys_m      = m.def_submodule("phys");
py::module phys_math_m = phys_m.def_submodule("math");

phys_math_m.attr("eu")       = phys::math::eu;
phys_math_m.attr("pi")       = phys::math::pi;
phys_math_m.attr("twopi")    = phys::math::twopi;
phys_math_m.attr("halfpi")   = phys::math::halfpi;
phys_math_m.attr("eps8")     = phys::math::eps8;
phys_math_m.attr("eps16")    = phys::math::eps16;
phys_math_m.attr("eps32")    = phys::math::eps32;
phys_math_m.attr("sqrttwo")  = phys::math::sqrttwo;
phys_math_m.attr("sqrthalf") = phys::math::sqrthalf;
phys_math_m.attr("im")       = phys::math::im;
phys_math_m.attr("iu")       = phys::math::iu;
phys_math_m.attr("iz")       = phys::math::iz;


py::class_<phys::dimension7>(phys_m, "dim")
    .def(py::init<const std::array<double, 7>&>())
    .def("__mul__", [](const phys::dimension7 self, const phys::dimension7 rdim) { return self * rdim; })
    .def("__truediv__", [](const phys::dimension7 self, const phys::dimension7 rdim) { return self / rdim; })
    .def("__pow__", [](const phys::dimension7 self, const double b) { return self.power(b); })
    .def("__repr__", [](const phys::dimension7 self) { return self.to_string(); });

phys_m.attr("dimensionless_d")                = phys::dimensionless_d;
phys_m.attr("length_d")                       = phys::length_d;
phys_m.attr("time_d")                         = phys::time_d;
phys_m.attr("mass_d")                         = phys::mass_d;
phys_m.attr("electric_current_d")             = phys::electric_current_d;
phys_m.attr("thermodynamic_temperature_d")    = phys::thermodynamic_temperature_d;
phys_m.attr("amount_of_substance_d")          = phys::amount_of_substance_d;
phys_m.attr("luminous_intensity_d")           = phys::luminous_intensity_d;
phys_m.attr("current_d")                      = phys::current_d;
phys_m.attr("temperature_d")                  = phys::temperature_d;
phys_m.attr("amount_d")                       = phys::amount_d;
phys_m.attr("none_d")                         = phys::none_d;
phys_m.attr("distance_d")                     = phys::distance_d;
phys_m.attr("wavelength_d")                   = phys::wavelength_d;
phys_m.attr("wave_number_d")                  = phys::wave_number_d;
phys_m.attr("area_d")                         = phys::area_d;
phys_m.attr("volume_d")                       = phys::volume_d;
phys_m.attr("frequency_d")                    = phys::frequency_d;
phys_m.attr("angular_velocity_d")             = phys::angular_velocity_d;
phys_m.attr("angular_acceleration_d")         = phys::angular_acceleration_d;
phys_m.attr("activity_of_a_nuclide_d")        = phys::activity_of_a_nuclide_d;
phys_m.attr("speed_d")                        = phys::speed_d;
phys_m.attr("acceleration_d")                 = phys::acceleration_d;
phys_m.attr("jerk_d")                         = phys::jerk_d;
phys_m.attr("jounce_d")                       = phys::jounce_d;
phys_m.attr("crackle_d")                      = phys::crackle_d;
phys_m.attr("pop_d")                          = phys::pop_d;
phys_m.attr("absement_d")                     = phys::absement_d;
phys_m.attr("area_flow_rate_d")               = phys::area_flow_rate_d;
phys_m.attr("volume_flow_rate_d")             = phys::volume_flow_rate_d;
phys_m.attr("kinematic_viscosity_d")          = phys::kinematic_viscosity_d;
phys_m.attr("thermal_diffusivity_d")          = phys::thermal_diffusivity_d;
phys_m.attr("specific_energy_d")              = phys::specific_energy_d;
phys_m.attr("dose_equivalent_d")              = phys::dose_equivalent_d;
phys_m.attr("absorbed_dose_d")                = phys::absorbed_dose_d;
phys_m.attr("absorbed_dose_rate_d")           = phys::absorbed_dose_rate_d;
phys_m.attr("substance_permeability_d")       = phys::substance_permeability_d;
phys_m.attr("inertia_d")                      = phys::inertia_d;
phys_m.attr("mass_line_density_d")            = phys::mass_line_density_d;
phys_m.attr("mass_area_density_d")            = phys::mass_area_density_d;
phys_m.attr("mass_density_d")                 = phys::mass_density_d;
phys_m.attr("specific_volume_d")              = phys::specific_volume_d;
phys_m.attr("mass_flow_rate_d")               = phys::mass_flow_rate_d;
phys_m.attr("mass_flow_acceleration_d")       = phys::mass_flow_acceleration_d;
phys_m.attr("mass_flow_jerk_d")               = phys::mass_flow_jerk_d;
phys_m.attr("force_d")                        = phys::force_d;
phys_m.attr("momentum_d")                     = phys::momentum_d;
phys_m.attr("energy_d")                       = phys::energy_d;
phys_m.attr("moment_of_force_d")              = phys::moment_of_force_d;
phys_m.attr("torque_d")                       = phys::torque_d;
phys_m.attr("angular_momentum_d")             = phys::angular_momentum_d;
phys_m.attr("action_d")                       = phys::action_d;
phys_m.attr("inv_ener_d")                     = phys::inv_ener_d;
phys_m.attr("power_d")                        = phys::power_d;
phys_m.attr("energy_density_d")               = phys::energy_density_d;
phys_m.attr("pressure_d")                     = phys::pressure_d;
phys_m.attr("surface_tension_d")              = phys::surface_tension_d;
phys_m.attr("energy_line_density_d")          = phys::energy_line_density_d;
phys_m.attr("power_density_d")                = phys::power_density_d;
phys_m.attr("power_area_density_d")           = phys::power_area_density_d;
phys_m.attr("dynamic_viscosity_d")            = phys::dynamic_viscosity_d;
phys_m.attr("heat_flow_rate_d")               = phys::heat_flow_rate_d;
phys_m.attr("heat_density_d")                 = phys::heat_density_d;
phys_m.attr("heat_density_flow_rate_d")       = phys::heat_density_flow_rate_d;
phys_m.attr("heat_flux_density_d")            = phys::heat_flux_density_d;
phys_m.attr("radiant_intensity_d")            = phys::radiant_intensity_d;
phys_m.attr("radiance_d")                     = phys::radiance_d;
phys_m.attr("irradiance_d")                   = phys::irradiance_d;
phys_m.attr("current_density_d")              = phys::current_density_d;
phys_m.attr("electric_charge_d")              = phys::electric_charge_d;
phys_m.attr("electric_charge_density_d")      = phys::electric_charge_density_d;
phys_m.attr("electric_area_charge_density_d") = phys::electric_area_charge_density_d;
phys_m.attr("electric_line_charge_density_d") = phys::electric_line_charge_density_d;
phys_m.attr("electric_dipole_moment_d")       = phys::electric_dipole_moment_d;
phys_m.attr("electric_flux_density_d")        = phys::electric_flux_density_d;
phys_m.attr("electric_displacement_field_d")  = phys::electric_displacement_field_d;
phys_m.attr("electric_polarization_field_d")  = phys::electric_polarization_field_d;
phys_m.attr("magnetic_moment_d")              = phys::magnetic_moment_d;
phys_m.attr("magnetic_field_strength_d")      = phys::magnetic_field_strength_d;
phys_m.attr("magnetization_d")                = phys::magnetization_d;
phys_m.attr("electric_potential_d")           = phys::electric_potential_d;
phys_m.attr("electric_field_strenth_d")       = phys::electric_field_strenth_d;
phys_m.attr("electric_resistance_d")          = phys::electric_resistance_d;
phys_m.attr("electric_conductance_d")         = phys::electric_conductance_d;
phys_m.attr("electric_resistivity_d")         = phys::electric_resistivity_d;
phys_m.attr("electric_conductivity_d")        = phys::electric_conductivity_d;
phys_m.attr("electric_capacitance_d")         = phys::electric_capacitance_d;
phys_m.attr("magnetic_flux_d")                = phys::magnetic_flux_d;
phys_m.attr("magnetic_flux_density_d")        = phys::magnetic_flux_density_d;
phys_m.attr("inductance_d")                   = phys::inductance_d;
phys_m.attr("electric_chargme_mass_ratio_d")  = phys::electric_chargme_mass_ratio_d;
phys_m.attr("magnetic_permeability_d")        = phys::magnetic_permeability_d;
phys_m.attr("permittivity_d")                 = phys::permittivity_d;
phys_m.attr("inv_temp_d")                     = phys::inv_temp_d;
phys_m.attr("heat_capacity_d")                = phys::heat_capacity_d;
phys_m.attr("entropy_d")                      = phys::entropy_d;
phys_m.attr("heat_transfer_coefficient_d")    = phys::heat_transfer_coefficient_d;
phys_m.attr("specific_heat_capacity_d")       = phys::specific_heat_capacity_d;
phys_m.attr("thermal_conductivity_d")         = phys::thermal_conductivity_d;
phys_m.attr("thermal_insulance_d")            = phys::thermal_insulance_d;
phys_m.attr("thermal_resistance_d")           = phys::thermal_resistance_d;
phys_m.attr("thermal_resistivity_d")          = phys::thermal_resistivity_d;
phys_m.attr("concentration_d")                = phys::concentration_d;
phys_m.attr("molar_energy_d")                 = phys::molar_energy_d;
phys_m.attr("molar_entropy_d")                = phys::molar_entropy_d;
phys_m.attr("luminous_flux_d")                = phys::luminous_flux_d;
phys_m.attr("illuminance_d")                  = phys::illuminance_d;
phys_m.attr("luminance_d")                    = phys::luminance_d;

py::class_<phys::uval>(phys_m, "uval")
    .def_readonly("dim", &phys::uval::dim)
    .def_readonly("value", &phys::uval::value)
    .def(py::init<const double&>())
    .def(py::init<const phys::dimension7&, const double&>())
    .def(
        "__add__", [](const phys::uval& lhs, const phys::uval& rhs) { return lhs + rhs; }, py::is_operator())
    .def(
        "__sub__", [](const phys::uval& lhs, const phys::uval& rhs) { return lhs - rhs; }, py::is_operator())
    .def(
        "__mul__", [](const phys::uval& lhs, const phys::uval& rhs) { return lhs * rhs; }, py::is_operator())
    .def(
        "__mul__", [](const phys::uval& lhs, const double& rhs) { return lhs * rhs; }, py::is_operator())
    .def(
        "__rmul__", [](const phys::uval& from_rhs, const double& from_lhs) { return from_lhs * from_rhs; },
        py::is_operator())
    .def(
        "__truediv__", [](const phys::uval& lhs, const phys::uval& rhs) { return lhs / rhs; }, py::is_operator())
    .def(
        "__truediv__", [](const phys::uval& lhs, const double& rhs) { return lhs / rhs; }, py::is_operator())
    .def(
        "__truediv__", [](const double& lhs, const phys::uval& rhs) { return lhs / rhs; }, py::is_operator())
    .def("__pow__", [](const phys::uval& self, const double b) { return phys::power(self, b); })
    .def("__repr__", [](const phys::uval& self) { return phys::to_string(self); });

phys_m.attr("_base_1")    = phys::_base_1;
phys_m.attr("_base_1m")   = phys::_base_1m;
phys_m.attr("_base_1s")   = phys::_base_1s;
phys_m.attr("_base_1kg")  = phys::_base_1kg;
phys_m.attr("_base_1A")   = phys::_base_1A;
phys_m.attr("_base_1K")   = phys::_base_1K;
phys_m.attr("_base_1mol") = phys::_base_1mol;
phys_m.attr("_base_1cd")  = phys::_base_1cd;
phys_m.attr("_base_1Hz")  = phys::_base_1Hz;
phys_m.attr("_base_1N")   = phys::_base_1N;
phys_m.attr("_base_1Pa")  = phys::_base_1Pa;
phys_m.attr("_base_1J")   = phys::_base_1J;
phys_m.attr("_base_1W")   = phys::_base_1W;
phys_m.attr("_base_1C")   = phys::_base_1C;
phys_m.attr("_base_1V")   = phys::_base_1V;
phys_m.attr("_base_1F")   = phys::_base_1F;
phys_m.attr("_base_1S")   = phys::_base_1S;
phys_m.attr("_base_1Om")  = phys::_base_1Om;
phys_m.attr("_base_1Wb")  = phys::_base_1Wb;
phys_m.attr("_base_1T")   = phys::_base_1T;
phys_m.attr("_base_1H")   = phys::_base_1H;

phys_m.attr("G_gravitional_constant") = phys::G_gravitional_constant;
phys_m.attr("c_lightspeed")           = phys::c_lightspeed;
phys_m.attr("ep0_permittivity")       = phys::ep0_permittivity;
phys_m.attr("mu0_permeability")       = phys::mu0_permeability;
phys_m.attr("ke_Comloub")             = phys::ke_Comloub;
phys_m.attr("R_gas_constant")         = phys::R_gas_constant;
phys_m.attr("Rydberg_constant")       = phys::Rydberg_constant;
phys_m.attr("Faraday_constant")       = phys::Faraday_constant;
phys_m.attr("Stefan_constant")        = phys::Stefan_constant;
phys_m.attr("muB_magnetic_moment")    = phys::muB_magnetic_moment;
phys_m.attr("muN_magnetic_moment")    = phys::muN_magnetic_moment;
phys_m.attr("Bohr_length")            = phys::Bohr_length;
phys_m.attr("h_Planck")               = phys::h_Planck;
phys_m.attr("hb_Planck")              = phys::hb_Planck;
phys_m.attr("me_mass")                = phys::me_mass;
phys_m.attr("mp_mass")                = phys::mp_mass;
phys_m.attr("mn_mass")                = phys::mn_mass;
phys_m.attr("amu_mass")               = phys::amu_mass;
phys_m.attr("e_charge")               = phys::e_charge;
phys_m.attr("k_Boltzman")             = phys::k_Boltzman;
phys_m.attr("N_Avagadro")             = phys::N_Avagadro;

py::module phys_si_m = phys_m.def_submodule("si");
phys_si_m.attr("c")  = phys::si::c;
phys_si_m.attr("h")  = phys::si::h;
phys_si_m.attr("hb") = phys::si::hb;
phys_si_m.attr("ke") = phys::si::ke;
phys_si_m.attr("me") = phys::si::me;
phys_si_m.attr("e")  = phys::si::e;
phys_si_m.attr("k")  = phys::si::k;
phys_si_m.attr("N")  = phys::si::N;
phys_si_m.attr("G")  = phys::si::G;
phys_si_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::si::as));
phys_si_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::si::as));

py::module phys_plk_m = phys_m.def_submodule("plk");
phys_plk_m.attr("c")  = phys::plk::c;
phys_plk_m.attr("h")  = phys::plk::h;
phys_plk_m.attr("hb") = phys::plk::hb;
phys_plk_m.attr("ke") = phys::plk::ke;
phys_plk_m.attr("me") = phys::plk::me;
phys_plk_m.attr("e")  = phys::plk::e;
phys_plk_m.attr("k")  = phys::plk::k;
phys_plk_m.attr("N")  = phys::plk::N;
phys_plk_m.attr("G")  = phys::plk::G;
phys_plk_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::plk::as));
phys_plk_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::plk::as));

py::module phys_god_m = phys_m.def_submodule("god");
phys_god_m.attr("c")  = phys::god::c;
phys_god_m.attr("h")  = phys::god::h;
phys_god_m.attr("hb") = phys::god::hb;
phys_god_m.attr("ke") = phys::god::ke;
phys_god_m.attr("me") = phys::god::me;
phys_god_m.attr("e")  = phys::god::e;
phys_god_m.attr("k")  = phys::god::k;
phys_god_m.attr("N")  = phys::god::N;
phys_god_m.attr("G")  = phys::god::G;
phys_god_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::god::as));
phys_god_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::god::as));

py::module phys_sty_m = phys_m.def_submodule("sty");
phys_sty_m.attr("c")  = phys::sty::c;
phys_sty_m.attr("h")  = phys::sty::h;
phys_sty_m.attr("hb") = phys::sty::hb;
phys_sty_m.attr("ke") = phys::sty::ke;
phys_sty_m.attr("me") = phys::sty::me;
phys_sty_m.attr("e")  = phys::sty::e;
phys_sty_m.attr("k")  = phys::sty::k;
phys_sty_m.attr("N")  = phys::sty::N;
phys_sty_m.attr("G")  = phys::sty::G;
phys_sty_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::sty::as));
phys_sty_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::sty::as));

py::module phys_nat_m = phys_m.def_submodule("nat");
phys_nat_m.attr("c")  = phys::nat::c;
phys_nat_m.attr("h")  = phys::nat::h;
phys_nat_m.attr("hb") = phys::nat::hb;
phys_nat_m.attr("ke") = phys::nat::ke;
phys_nat_m.attr("me") = phys::nat::me;
phys_nat_m.attr("e")  = phys::nat::e;
phys_nat_m.attr("k")  = phys::nat::k;
phys_nat_m.attr("N")  = phys::nat::N;
phys_nat_m.attr("G")  = phys::nat::G;
phys_nat_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::nat::as));
phys_nat_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::nat::as));

py::module phys_cgs_m = phys_m.def_submodule("cgs");
phys_cgs_m.attr("c")  = phys::cgs::c;
phys_cgs_m.attr("h")  = phys::cgs::h;
phys_cgs_m.attr("hb") = phys::cgs::hb;
phys_cgs_m.attr("ke") = phys::cgs::ke;
phys_cgs_m.attr("me") = phys::cgs::me;
phys_cgs_m.attr("e")  = phys::cgs::e;
phys_cgs_m.attr("k")  = phys::cgs::k;
phys_cgs_m.attr("N")  = phys::cgs::N;
phys_cgs_m.attr("G")  = phys::cgs::G;
phys_cgs_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::cgs::as));
phys_cgs_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::cgs::as));

py::module phys_ryd_m = phys_m.def_submodule("ryd");
phys_ryd_m.attr("c")  = phys::ryd::c;
phys_ryd_m.attr("h")  = phys::ryd::h;
phys_ryd_m.attr("hb") = phys::ryd::hb;
phys_ryd_m.attr("ke") = phys::ryd::ke;
phys_ryd_m.attr("me") = phys::ryd::me;
phys_ryd_m.attr("e")  = phys::ryd::e;
phys_ryd_m.attr("k")  = phys::ryd::k;
phys_ryd_m.attr("N")  = phys::ryd::N;
phys_ryd_m.attr("G")  = phys::ryd::G;
phys_ryd_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::ryd::as));
phys_ryd_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::ryd::as));

py::module phys_hat_m = phys_m.def_submodule("hat");
phys_hat_m.attr("c")  = phys::hat::c;
phys_hat_m.attr("h")  = phys::hat::h;
phys_hat_m.attr("hb") = phys::hat::hb;
phys_hat_m.attr("ke") = phys::hat::ke;
phys_hat_m.attr("me") = phys::hat::me;
phys_hat_m.attr("e")  = phys::hat::e;
phys_hat_m.attr("k")  = phys::hat::k;
phys_hat_m.attr("N")  = phys::hat::N;
phys_hat_m.attr("G")  = phys::hat::G;
phys_hat_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::hat::as));
phys_hat_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::hat::as));

py::module phys_au_m = phys_m.def_submodule("au");
phys_au_m.attr("c")  = phys::au::c;
phys_au_m.attr("h")  = phys::au::h;
phys_au_m.attr("hb") = phys::au::hb;
phys_au_m.attr("ke") = phys::au::ke;
phys_au_m.attr("me") = phys::au::me;
phys_au_m.attr("e")  = phys::au::e;
phys_au_m.attr("k")  = phys::au::k;
phys_au_m.attr("N")  = phys::au::N;
phys_au_m.attr("G")  = phys::au::G;
phys_au_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::au::as));
phys_au_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::au::as));

py::module phys_qcd_m = phys_m.def_submodule("qcd");
phys_qcd_m.attr("c")  = phys::qcd::c;
phys_qcd_m.attr("h")  = phys::qcd::h;
phys_qcd_m.attr("hb") = phys::qcd::hb;
phys_qcd_m.attr("ke") = phys::qcd::ke;
phys_qcd_m.attr("me") = phys::qcd::me;
phys_qcd_m.attr("e")  = phys::qcd::e;
phys_qcd_m.attr("k")  = phys::qcd::k;
phys_qcd_m.attr("N")  = phys::qcd::N;
phys_qcd_m.attr("G")  = phys::qcd::G;
phys_qcd_m.def("_as", static_cast<double (*)(const phys::dimension7, const std::string&)>(&phys::qcd::as));
phys_qcd_m.def("_as", static_cast<double (*)(const phys::dimension7, const phys::uval&)>(&phys::qcd::as));

phys_m.attr("au_2_amu")       = phys::au_2_amu;
phys_m.attr("au_2_ang")       = phys::au_2_ang;
phys_m.attr("au_2_ev")        = phys::au_2_ev;
phys_m.attr("au_2_J_1mea")    = phys::au_2_J_1mea;
phys_m.attr("au_2_kcal_1mea") = phys::au_2_kcal_1mea;
phys_m.attr("au_2_wn")        = phys::au_2_wn;
phys_m.attr("au_2_fs")        = phys::au_2_fs;
phys_m.attr("au_2_ps")        = phys::au_2_ps;
phys_m.attr("au_2_K")         = phys::au_2_K;
phys_m.attr("au_2_angoverps") = phys::au_2_angoverps;
