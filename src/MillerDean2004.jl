"""
#######################################
# Miller and Dean (2004) model        #
# Hb: wave breaking height            #
# depthb: depth of breaking           #
# sl: sea level                       #
# Y0: base position                   #
# dt: time step                       #
# D50: median grain size              #
# Hberm: berm height                  #
# kero: erosion coefficient           #
# kacr: accretion coefficient         #
# Yi: initial position                #
# flagP: flag for kero and kacr       #
# Omega: Dean's parameter             #
#######################################
"""

function run_MillerDean()

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = NCDataset(dats*"wav.nc")
    parF = NCDataset(dats*"par.nc")
    slF = NCDataset(dats*"sl.nc")

    configF = NCDataset(dats*"config.nc")

    println("Unpacking datasets...")

    dt = configF["dt"][:][1]
    
    yi = configF["yi"][:][1]

    kacr = parF["kacr"][:][1]
    kero = parF["kero"][:][1]
    Y0 = parF["Y0"][:][1]
    
    brk, angBati, depth, Hberm, D50 = configF["brk"][:][1], configF["angBati"][:][1], configF["depth"][:][1], configF["Hberm"][:][1], configF["D50"][:][1]

    if brk == 1
        
        Hs, Tp, θ_w = collect(skipmissing(wavF["Hs"][:])), collect(skipmissing(wavF["Tp"][:])), collect(skipmissing(wavF["Dir"][:]))

        auxAng, auxDepth = similar(Hs), similar(Hs)
        auxAng .= angBati
        auxDepth .= depth

        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    else
        Hb, Tp, Hs, depthb = wavF["Hb"][:], wavF["Tp"][:], wavF["Hs"][:], wavF["hb"][:]
    end


    Hs = convert(Array{Float64},Hs)
    Tp = convert(Array{Float64},Tp)

    kacr = convert(Array{Float64},kacr)
    kero = convert(Array{Float64},kero)
    Y0 = convert(Array{Float64},Y0)

    Yi = convert(Array{Float64},yi)

    ########## START HERE #############
    w = wMOORE(D50)

    Omega = Hb ./ (w .* Tp)

    println("Starting Jaramillo et al. (2021) - Shoreline Rotation Model...")


    Y_t = MileerDean(Hb, depthb, sl, Y0, dt, D50, Hberm, kero, kacr, Yi, flagP, Omega)

    println("\n\n****************Finished****************\n\n")

    return Y_t
    
end


function MileerDean(Hb, depthb, sl, Y0, dt, D50, Hberm, kero, kacr, Yi, flagP = 1, Omega = 0)
    if flagP == 1
        kero = fill(kero, length(Hb))
        kacr = fill(kacr, length(Hb))
    elseif flagP == 2
        kero = kero .* Hb .^ 2
        kacr = kacr .* Hb .^ 2
    elseif flagP == 3
        kero = kero .* Hb .^ 3
        kacr = kacr .* Hb .^ 3
    elseif flagP == 4
        kero = kero .* Omega
        kacr = kacr .* Omega
    end

    yeq = similar(Hb)
    Y = similar(Hb)
    wl = 0.106 .* Hb .+ sl
    Wast = wast(depthb, D50)
    yeq .= Y0 .- Wast .* wl ./ (Hberm .+ depthb)
    acr_count = 0
    ero_count = 0

    for i in eachindex(sl)
        if i == 1
            r = yeq[i] .- Y[i] > 0
            k = kacr[i] * r + kero[i] * !r
            Y[i] = Yi
        else
            r = yeq[i] .- Y[i-1] > 0
            k = kacr[i] * r + kero[i] * !r
            A = k .* dt .* 0.5
            Y[i] = (Y[i-1] + A .* (yeq[i] + yeq[i-1] - Y[i-1])) ./ (1 + A)
        end
        acr_count = acr_count + r
        ero_count = ero_count + !r
    end

    return Y
end


function cal_MillerDean()

    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data/"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Reading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"
    parF = dats*"par.nc"
    slF = dats*"sl.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]
    
    Yi = ncread(configF, "Yi")[1]

    brk, angBati, depth, Hberm, D50 = ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1], ncread(configF, "Hberm")[1], ncread(configF, "D50")[1]

    flagP = ncread(configF, "flagP")[1]

    MetObj = ncread(configF, "MetObj")[1]

    sl = ncread(slF, "sl")

    if brk == 1
        
        Hs, Tp, θ_w = ncread(wavF, "Hs"), ncread(wavF, "Tp"), ncread(wavF, "Dir")

        auxAng, auxDepth = similar(Hs), similar(Hs)
        auxAng .= angBati
        auxDepth .= depth

        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    else
        Hb, Tp, Hs, depthb = ncread(wavF, "Hb"), ncread(wavF, "Tp"), ncread(wavF, "Hs"), ncread(wavF, "hb")
    end

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")
    
    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)
    t_wav = DateTime.(YY, MM, DD, HH)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, hb, tp, sl, hs, depthb = t_wav[ii], Hb[ii], Tp[ii], sl[ii], Hs[ii], depthb[ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = findall((x)-> x == t_obs[i], t_wav)[1]
    end

    idx_obs = convert(Array{Int64},idx_obs)

    ########## START HERE #############
    w = wMOORE(D50)

    Omega = Hb ./ (w .* Tp)

    println("Starting Miller and Dean (2004) - Shoreline Evolution Model...")

    function Calibra_MDr(Χ)

        Ymd = MileerDean(Hb, depthb, sl, Χ[3], dt, D50, Hberm, exp(Χ[1]), exp(Χ[2]), Yi, flagP, Omega)
        # Ymd, _ = HM.MILLER_DEAN_CSonly(hb,hb./.78,sl,exp(Χ[1]),dt,D50,Hberm, exp(Χ[2]), exp(Χ[3]),Χ[4], flagP, Omega)

        YYsl = Ymd[idx_obs]
        if MetObj == "Pearson"
            return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
        elseif MetObj == "RMSE"
            return abs(sqrt(mean((YYsl .- Y_obs).^2)))
        elseif MetObj == "MSS"
            return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
        elseif MetObj == "BSS"
            return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
        elseif MetObj == "Double"
            return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))))
        end
    end
    
    boundsr = [(log(1e-7), log(1e-1)),
               (log(1e-7), log(1e-1)),
               (minimum(Y_obs), maximum(Y_obs))] 
    
    resr = bboptimize(Calibra_MDr; 
                        # Method = :simultaneous_perturbation_stochastic_approximation,
                        SearchRange = boundsr,
                        NumDimensions = 3,
                        PopulationSize = 5000,
                        MaxSteps = 100000,
                        FitnessTolerance = 1e-6,
                        FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                        TraceMode=:compact,
                        ϵ=0.01,
                        τ = 0.25,
                        MaxStepsWithoutEpsProgress = 20000,
                        TargetFitness=(1, 0.1),
                        Method=:borg_moea)
    
    objr = best_fitness(resr)
    popr = best_candidate(resr)
    
    Ymdr = MileerDean(Hb, depthb, sl, popr[3], dt, D50, Hberm, exp(popr[1]), exp(popr[2]), Yi, flagP, Omega)
    
    Ysl = Ymdr[idx_obs]
    aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
    aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
    aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)
    

    println("\n\n****************Finished****************\n\n")

    Y_atts = Dict("units" => "m",
        "long_name" => "Shoreline position",
        "standard_name" => "Y",
        "data_min" => minimum(Ymdr),
        "data_max" => maximum(Ymdr))
    kacr_atts = Dict("units" => "-",
        "long_name" => "Accretion coefficient",
        "standard_name" => "kacr")
    kero_atts = Dict("units" => "-",
        "long_name" => "Erosion coefficient",
        "standard_name" => "kero")
    Y0_atts = Dict("units" => "m",
        "long_name" => "Base position",
        "standard_name" => "Y0")
    RP_atts = Dict("units" => "-",
        "long_name" => "Pearson correlation coefficient",
        "standard_name" => "RP")
    RMSE_atts = Dict("units" => "m",
        "long_name" => "Root mean square error",
        "standard_name" => "RMSE")
    MSS_atts = Dict("units" => "-",
        "long_name" => "Mielke Skill Score",
        "standard_name" => "MSS")

    nccreate("Calib_MD.nc", "Y", "Time", t_wav, atts = Y_atts)
    ncwrite(Ymdr, "Calib_MD.nc", "Y")
    nccreate("Calib_MD.nc", "kero", atts = kero_atts)
    ncwrite(exp(popr[1]), "Calib_MD.nc", "kero")
    nccreate("Calib_MD.nc", "kacr", atts = kacr_atts)
    ncwrite(exp(popr[2]), "Calib_MD.nc", "kacr")
    nccreate("Calib_MD.nc", "Y0", atts = Y0_atts)
    ncwrite(popr[3], "Calib_MD.nc", "Y0")
    nccreate("Calib_MD.nc", "RP", atts = RP_atts)
    ncwrite(aRP, "Calib_MD.nc", "RP")
    nccreate("Calib_MD.nc", "RMSE", atts = RMSE_atts)
    ncwrite(aRMSE, "Calib_MD.nc", "RMSE")
    nccreate("Calib_MD.nc", "MSS", atts = MSS_atts)
    ncwrite(aMSS, "Calib_MD.nc", "MSS")
        
    
end