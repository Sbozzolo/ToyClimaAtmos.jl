module ToyClimaAtmos

import ClimaAtmos as CA
import ClimaComms as CC
import CLIMAParameters as CP
import ClimaCore: InputOutput, Meshes, Spaces
using Dates: DateTime, @dateformat_str

import OrdinaryDiffEq as ODE

"""
    ToyAtmosIntegrator(param_dict)

Return a ClimaAtmos simulation object with the given parameters.

Does not directly depend on any file or command-line argument.

"""
struct ToyAtmosIntegrator
    # To run a ClimaAtmos simulation, one needs a CA.AtmosConfig. In turn CA.AtmosConfig has three fields:
    # - toml_dict, which is effectively just a dictionary
    # - parsed_args, which is effectively just a dictionary
    # - comms_ctx, which is the device where the simulation is run
    # This is a lot of structure if one wants to run a simple simulation from an
    # external library (such as ToyClimaAtmos).
    #
    # ToyAtmosIntegrator provides a simplified interface that a dictionary of
    # parameters and sets a simulation.

    param_dict::Dict{String,Any}
    ode::Any
end

function ToyAtmosIntegrator(param_dict::Dict{String,Any})

    # A GPU if available, otherwise a single-threaded CPU
    comms_ctx = CC.context(CC.device())
    CC.init(comms_ctx)

    FT = Float64

    # We retrieve the default configuration provided by ClimaAtmos, then we merge
    # it with the user-provided one
    default_config = CA.default_config_dict()
    config = merge(default_config, param_dict)

    dt = FT(CA.time_to_seconds(config["dt"]))

    # From CA.create_parameter_set()
    overrides = (;
        R_d = 287.0,
        MSLP = 1.0e5,
        grav = 9.80616,
        Omega = 7.29212e-5,
        planet_radius = 6.371229e6,
        ρ_cloud_liq = 1e3,
        τ_precip = dt,
        qc_0 = 5e-6,
    )

    # CA.is_column_* look into specific keys in param_dict. We don't know if
    # these keys are available or not, but the functions access them directly,
    # erroring out if the keys are not found. For this reason, we have to wrap
    # these calls in try-catch.
    try
        if CA.is_column_edmf(config)
            overrides = (; MSLP = 100000.0, τ_precip = dt)
        elseif CA.is_column_without_edmf(param_dict)
            overrides = (; τ_precip = dt)
        end
    catch
    end

    # CA.create_climaatmos_parameter_set takes a CP.AbstractTOMLDict as first entry.
    # AbstractTOMLDict are just dictionaries, so let's wrap our parameter dictionary into
    # a CP.ParamDict
    base_param_CP = CP.create_toml_dict(FT)

    params = CA.create_climaatmos_parameter_set(base_param_CP, config, overrides)

    # From CA.get_integrator

    # get_atmos

    moisture_model = CA.get_moisture_model(config)
    precip_model = CA.get_precipitation_model(config)
    radiation_mode = CA.get_radiation_mode(config, FT)
    forcing_type = CA.get_forcing_type(config)

    diffuse_momentum = !(forcing_type isa CA.HeldSuarezForcing)

    edmfx_adv_test = get(config, "edmfx_adv_test", false)
    @assert edmfx_adv_test in (false, true)

    edmfx_entr_detr = get(config, "edmfx_entr_detr", false)
    @assert edmfx_entr_detr in (false, true)

    edmfx_sgs_flux = get(config, "edmfx_sgs_flux", false)
    @assert edmfx_sgs_flux in (false, true)

    edmfx_nh_pressure = get(config, "edmfx_nh_pressure", false)
    @assert edmfx_nh_pressure in (false, true)

    model_config = CA.get_model_config(config)
    vert_diff = CA.get_vertical_diffusion_model(diffuse_momentum, config, FT)
    atmos = CA.AtmosModel(;
        model_config,
        perf_mode = CA.get_perf_mode(config),
        moisture_model,
        energy_form = CA.get_energy_form(config, vert_diff),
        radiation_mode,
        subsidence = CA.get_subsidence_model(config, radiation_mode, FT),
        ls_adv = CA.get_large_scale_advection_model(config, FT),
        edmf_coriolis = CA.get_edmf_coriolis(config, FT),
        edmfx_entr_detr,
        edmfx_sgs_flux,
        edmfx_nh_pressure,
        precip_model,
        forcing_type,
        turbconv_model = CA.get_turbconv_model(
            FT,
            moisture_model,
            precip_model,
            config,
            params.turbconv_params,
        ),
        non_orographic_gravity_wave = CA.get_non_orographic_gravity_wave_model(
            config,
            model_config,
            FT,
        ),
        orographic_gravity_wave = CA.get_orographic_gravity_wave_model(config, FT),
        hyperdiff = CA.get_hyperdiffusion_model(config, FT),
        vert_diff,
        viscous_sponge = CA.get_viscous_sponge_model(config, FT),
        rayleigh_sponge = CA.get_rayleigh_sponge_model(config, FT),
        sfc_temperature = CA.get_sfc_temperature_form(config),
        surface_model = CA.get_surface_model(config),
    )

    @info "AtmosModel: \n$(summary(atmos))"

    numerics = CA.get_numerics(config)

    job_id = if isnothing(config["job_id"])
        CA.job_id_from_config(config)
    else
        config["job_id"]
    end
    out_dir = config["output_dir"]
    output_dir = isnothing(out_dir) ? joinpath("/tmp/output", job_id) : out_dir

    isdir(output_dir) || @info "Creating output directory: $(mkdir(output_dir))"

    simulation = (;
        comms_ctx,
        is_debugging_tc = config["debugging_tc"],
        output_dir,
        restart = haskey(ENV, "RESTART_FILE"),
        nothing,
        dt = FT(CA.time_to_seconds(config["dt"])),
        start_date = DateTime(config["start_date"], dateformat"yyyymmdd"),
        t_end = FT(CA.time_to_seconds(config["t_end"])),
    )
    n_steps = floor(Int, simulation.t_end / simulation.dt)
    @info(
        "Time info:",
        dt = config["dt"],
        t_end = config["t_end"],
        floor_n_steps = n_steps,
    )
    initial_condition = CA.get_initial_condition(config)
    surface_setup = CA.get_surface_setup(config)

    s = CA.@timed_str begin
        spaces = CA.get_spaces(config, params, comms_ctx)
        Y = CA.InitialConditions.atmos_state(
            initial_condition(params),
            atmos,
            spaces.center_space,
            spaces.face_space,
        )
        t_start = CA.Spaces.undertype(axes(Y.c))(0)
    end

    @info "Allocating Y: $s"

    s = CA.@timed_str begin
        p = CA.get_cache(
            Y,
            config,
            params,
            spaces,
            atmos,
            numerics,
            simulation,
            initial_condition,
            surface_setup,
        )
    end
    @info "Allocating cache (p): $s"

    if get(config, "discrete_hydrostatic_balance", false)
        CA.set_discrete_hydrostatic_balanced_state!(Y, p)
    end

    FT = CA.Spaces.undertype(axes(Y.c))
    s = CA.@timed_str begin
        ode_algo = CA.ode_configuration(FT, config)
    end
    @info "ode_configuration: $s"

    s = CA.@timed_str begin
        callback = CA.get_callbacks(config, simulation, atmos, params)
    end
    @info "get_callbacks: $s"
    @info "n_steps_per_cycle_per_cb: $(CA.n_steps_per_cycle_per_cb(callback, simulation.dt))"
    @info "n_steps_per_cycle: $(CA.n_steps_per_cycle(callback, simulation.dt))"
    tspan = (t_start, simulation.t_end)
    s = CA.@timed_str begin
        integrator_args, integrator_kwargs =
            CA.args_integrator(config, Y, p, tspan, ode_algo, callback)
    end

    s = CA.@timed_str begin
        integrator = ODE.init(integrator_args...; integrator_kwargs...)
    end
    @info "init integrator: $s"
    return ToyAtmosIntegrator(param_dict, integrator)
end

# baroclinic = Dict(
#     "dt_save_to_disk" => "2days",
#     "regression_test" => true,
#     "initial_condition" => "DryBaroclinicWave",
#     "dt" => "580secs",
#     "t_end" => "10days",
#     "job_id" => "sphere_baroclinic_wave_rhoe",
# )

# ToyAtmosIntegrator(baroclinic)

end
