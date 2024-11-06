class calcType(Enum):
    opt = 0
    grad = 1
    ts = -1
    freq = 2
    irc = 3
    ci = -2
    mdv = 5

class BindKeywords(TypedDict):
command[1]    calc_type: calcType
command[2]    verbosity: int
command[3]    freeze_medium: int #?
command[4]    freeze_border_atoms: bool
command[5]    memory: str
command[6]    read_wavefunc: str
command[7]    nproc_parallel_runs: int
command[8]    parallel_finite_difference: bool
command[9]    nproc_parallel_freq: int
command[10]    numerical_computation_displacement: int #? not clear
command[11]    optimize_ci_coefficients: bool
command[12]    displacement_step: float
command[13]    target_state: int
command[14]    derivative_coupling_type: dcType
command[15]    save_single_points: bool
command[16]    numerical_differentiation_link_atoms: bool
command[18]    perform_qm_calculations_frequency: bool
command[19]    branching_plane_definition_ci: bpType
command[20]    compute_forces_with_amber_velocities: int #?
command[21]    link_distance_ch: float
command[22]    link_distance_oh: float
command[23]    link_distance_nh: float
command[24]    equilibrium_distance_cc: float
command[25]    equilibrium_distance_co: float
command[26]    equilibrium_distance_cn: float
command[40]    lbound_normal_modes_frequency: float
command[41]    cutoff_value_mm_calculations: float
command[51]    qm_calculation_type: qmType
command[53]    memory_definition_gaussian: str
command[60]    opt_maxiter: int
command[61]    correct_gradient_q1: bool
command[65]    ngridpoints_irc: int
command[66]    stepsize_irc: float
command[67]    thres_first_derivative: float
command[68]    max_opt_step: int
command[71]    mm_calculation_type: mmType
command[80]    random_number: int
command[81]    hopping_scheme_state_limit: int
command[83]    first_timestep_md: float
command[84]    second_timestep_md: float
command[85]    activate_hopping_type: hoppingType
command[86]    thres_hopping_energy: float
command[87]    back_surface_hopping: bool
command[88]    charge_surface_hopping: int #?
command[89]    thres_hopping_energy_gap: float
command[90]    es_gs_hop_type: esgshopType
command[91]    thres_es_gs_hop_gap: float
command[92]    ow_switch_off_weighting_smoothly: int
command[93]    ow_switch_off_weightingAbruptly: int
command[94]    ow_intruder_state_detection: int
command[95]    ow_scale_sigmoid_exponent: float
command[96]    ow_scale_half_height_position: float
command[97]    td_type: tdType
command[98]    thres_wf_expansion: float
command[99]    save_log_file_steps: int
command[100]    save_wavefunction_files_steps: int
command[101]    internal_molcas_numerics: bool
command[102]    ci_vector_rotation_protocol: int #?
command[103]    compute_tddft_after_hopping: bool
command[110]    qm_code_sharc: str
command[120]    use_molcas_ghost_atoms: bool
command[122]    criterion_stopping_traj: trajstopType
command[130]    freeze_gradient_high_layer: bool
command[190]    excited_state_gradient_turbomole: int #?
command[191]    use_ri_rigrad_turbomole: int #?
command[192]    disable_nosym_warning: bool
command[193]    scale_initial_vel_isotopes: bool
command[194]    use_cartesian_functions_molcas: bool
command[195]    ricd_threshold_molcas: float
command[196]    use_ri_approximation_molcas: bool
command[197]    basis_set_molcas: str
command[198]    molpro_dont_bomb_wfconvergence: bool
command[199]    molcas_use_SS_in_GS: bool
command[200]    use_single_point_on_top_of_SS: bool
command[201]    nac_corr_scheme: naccorrType
command[203]    switch_to_mp2_in_gs: bool
command[204]    thres_states_in_computation: float
command[205]    include_states_scheme: incstateType
command[206]    correct_velocity_type: corrvelType
command[207]    thres_es_es_hop: float

def makehard() -> CobraCommand:
    command = CobraCommand(
        calc_type='optxg',
        verbosity=0,
        freeze_medium=0,
        freeze_border_atoms=0,
        memory='500MB',
        read_wavefunc='auto',
        nproc_parallel_runs=1,
        parallel_finite_difference=0,
        nproc_parallel_freq=1,
        numerical_computation_displacement=1,
        optimize_ci_coefficients=1,
        displacement_accuracy=0.001,
        state_to_relax=1,
        derivative_coupling_type=0,
        save_single_points=0,
        numerical_differentiation_link_atoms=1,
        perform_qm_calculations_frequency=0,
        branching_plane_definition_ci=0,
        compute_forces_with_amber_velocities=0,
        link_distance_ch=1.090,
        link_distance_oh=0.947,
        link_distance_nh=1.008,
        equilibrium_distance_cc=1.526,
        equilibrium_distance_co=1.410,
        equilibrium_distance_cn=1.475,
        frequency_normal_modes=0,
        cutoff_value_mm_calculations=999,
        qm_calculation_type=1,
        memory_definition_gaussian='700MB',
        number_optimization_cycles=100,
        correct_gradient_q1=1,
        number_points_irc=10,
        step_size_irc=10,
        convergence_threshold_first_derivative=300,
        max_size_optimization_step=30,
        mm_calculation_type=0,
        random_number=0,
        hopping_scheme_state_limit=0,
        first_time_step_md=1.0,
        second_time_step_md=0.25,
        activate_surface_hopping=0,
        energy_threshold_surface_hopping=1000.0,
        back_surface_hopping=1,
        charge_surface_hopping=0,
        energy_gap_threshold_charge_hopping=15.0,
        es_gs_hop_type=0,
        threshold_es_gs_hop=2.0,
        ow_switch_off_weighting_smoothly=0,
        ow_switch_off_weightingAbruptly=0,
        ow_intruder_state_detection=0,
        ow_scale_sigmoid_exponent=0,
        ow_scale_half_height_position=0,
        td_type=0,
        wf_expansion_threshold=1.0,
        save_log_file_steps=1,
        save_wavefunction_files_steps=-1,
        internal_molcas_numerics=0,
        ci_vector_rotation_protocol=0,
        compute_tddft_after_hopping=1,
        name_qm_code_sharc='',
        use_molcas_ghost_atoms=0,
        criterion_stopping_trajectory=0,
        compute_gradient_high_layer=0,
        excited_state_gradient_turbomole=0,
        use_ri_rigrad_turbomole=0,
        nosym_warning=0,
        scale_initial_vel_isotopes=0,
        use_cartesian_functions_molcas=0,
        ricd_threshold_molcas=1.0E-4,
        use_ri_approximation_molcas=0,
        basis_set_molcas='6-31Gp',
        molproWF_convergence_bomb_job=1,
        use_molcasSS_in_GS=0,
        use_single_point_on_top_of_SS=0,
        massey_tully_scheme=0,
        switch_to_mp2_in_gs=0,
        threshold_states_in_computation=0,
        include_states_below_threshold=0,
        correct_velocity_after_hop=0,
        threshold_es_es_hop=0
    )
    return command