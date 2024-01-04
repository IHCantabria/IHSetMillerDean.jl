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

    wavF = dats*"wav.nc"
    parF = dats*"par.nc"
    slF = dats*"sl.nc"

    configF = dats*"config.nc"

    dt = ncread(configF, "dt")[1]
    
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

    t_wav = DateTime.(YY, MM, DD, HH)

    kacr = ncread(parF, "kacr")
    kero = ncread(parF, "kero")
    Y0 = ncread(parF, "Y0")

    ########## START HERE #############
    w = wMOORE(D50)

    Omega = Hb ./ (w .* Tp)

    println("Starting Jaramillo et al. (2021) - Shoreline Rotation Model...")


    Y_t = MileerDean(Hb, depthb, sl, Y0, dt, D50, Hberm, kero, kacr, Yi, flagP, Omega)

    println("\n\n****************Finished****************\n\n")

    return Y_t
    
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

    brk, angBati, depth, Hberm, D50 = ncread(configF, "brk")[1], ncread(configF, "angBati")[1], ncread(configF, "depth")[1], ncread(configF, "Hberm")[1], ncread(configF, "D50")[1]

    if length(ncread(configF, "flagP")) == 1
        flagP = ncread(configF, "flagP")[1]
    else
        flagP = ncread(configF, "flagP")
    end

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

    Ymdr = Dict()
    aRP = Dict()
    aRMSE = Dict()
    aMSS = Dict()
    popr = Dict()
    objr = Dict()
    
    for i in eachindex(flagP)
            
        function Calibra_MDr(Χ)

            # Ymd = MileerDean(Hb, depthb, sl, Χ[3], dt, D50, Hberm, exp(Χ[1]), exp(Χ[2]), Y_obs[1], flagP[i], Omega)
            Ymd = MileerDean(Hb, depthb, sl, Χ[3], dt, D50, Hberm, exp(Χ[1]), exp(Χ[2]), Χ[4], flagP[i], Omega)
            YYsl = Ymd[idx_obs]
            if MetObj == "Pearson"
                return 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl)))
            elseif MetObj == "RMSE"
                return abs(sqrt(mean((YYsl .- Y_obs).^2))/5)
            elseif MetObj == "MSS"
                return sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2)
            elseif MetObj == "BSS"
                return (mean((YYsl .- Y_obs).^2) - mean((YYref .- Y_obs).^2))/mean((YYref .- Y_obs).^2)
            elseif MetObj == "Double"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5))
            elseif MetObj == "Triple"
                return (sum((YYsl .- Y_obs).^2)/length(YYsl)/(var(YYsl)+var(Y_obs)+(mean(YYsl)-mean(Y_obs))^2), abs(sqrt(mean((YYsl .- Y_obs).^2))/5), 1 -  abs(sum((YYsl.-mean(YYsl)).*(Y_obs .- mean(Y_obs)))/(std(YYsl)*std(Y_obs)*length(YYsl))))
            end
        end

        boundsr = [(log(1e-7), log(1e-1)),
                   (log(1e-7), log(1e-1)),
                   (0.25*minimum(Y_obs), 2*maximum(Y_obs)),
                   (0.25*minimum(Y_obs), 2*maximum(Y_obs))] 

        if MetObj == "Double"
            resr = bboptimize(Calibra_MDr; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.01,
                            τ = 0.25,
                            MaxStepsWithoutEpsProgress = 10000,
                            Method=:borg_moea)
        elseif MetObj == "Triple"
            resr = bboptimize(Calibra_MDr; 
                            # Method = :simultaneous_perturbation_stochastic_approximation,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                            TraceMode=:compact,
                            ϵ=0.01,
                            τ = 0.25,
                            MaxStepsWithoutEpsProgress = 10000,
                            Method=:borg_moea)
        else
            resr = bboptimize(Calibra_MDr; 
                            Method = :adaptive_de_rand_1_bin,
                            SearchRange = boundsr,
                            NumDimensions = 4,
                            PopulationSize = 500,
                            MaxSteps = 5000,
                            FitnessTolerance = 1e-6,
                            TraceMode=:compact,
                            ϵ=0.01,
                            τ = 0.25,
                            MaxStepsWithoutEpsProgress = 10000)
        end


        

        objr[string(i)] = best_fitness(resr)
        popr[string(i)] = best_candidate(resr)

        Ymdr[string(i)] = MileerDean(Hb, depthb, sl, popr[string(i)][3], dt, D50, Hberm, exp(popr[string(i)][1]), exp(popr[string(i)][2]), popr[string(i)][4], flagP[i], Omega)

        Ysl = Ymdr[string(i)][idx_obs]
        aRP[string(i)] = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
        aRMSE[string(i)] = sqrt(mean((Ysl .- Y_obs).^2))
        aMSS[string(i)] = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

    end
    println("\n\n****************Finished****************\n\n")

    year_atts = Dict("long_name" => "Year")
    month_atts = Dict("long_name" => "Month")
    day_atts = Dict("long_name" => "Day")
    hour_atts = Dict("long_name" => "Hour")
    println("Writing output...")

    output = wrkDir*"/results/Shoreline_MD.nc"
    nccreate(output, "year",
                "dim", length(YY),
                atts = year_atts)
    ncwrite(YY, output, "year")
    nccreate(output, "month",
                "dim", length(MM),
                atts = month_atts)
    ncwrite(MM, output, "month")
    nccreate(output, "day",
                "dim", length(DD),
                atts = day_atts)
    ncwrite(DD, output, "day")
    nccreate(output, "hour",
                "dim", length(HH),
                atts = hour_atts)
    ncwrite(HH, output, "hour")  

    Y_atts = Dict("units" => "m",
            "long_name" => "Shoreline position",
            "standard_name" => "Y")
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

    for i in eachindex(flagP)
          
        nccreate(output, "kacr_flagP="*string(i),
                    "len", 1,
                    atts = kacr_atts)
        ncwrite([exp(popr[string(i)][1])], output, "kacr_flagP="*string(i))
        nccreate(output, "kero_flagP="*string(i),
                    "len", 1,
                    atts = kero_atts)
        ncwrite([exp(popr[string(i)][2])], output, "kero_flagP="*string(i))
        nccreate(output, "Y0_flagP="*string(i),
                    "len", 1,
                    atts = Y0_atts)
        ncwrite([popr[string(i)][3]], output, "Y0_flagP="*string(i))
        nccreate(output, "RP_flagP="*string(i),
                    "len", 1,
                    atts = RP_atts)
        ncwrite([aRP[string(i)]], output, "RP_flagP="*string(i))
        nccreate(output, "RMSE_flagP="*string(i),
                    "len", 1,
                    atts = RMSE_atts)
        ncwrite([aRMSE[string(i)]], output, "RMSE_flagP="*string(i))
        nccreate(output, "MSS_flagP="*string(i),
                    "len", 1,
                    atts = MSS_atts)
        ncwrite([aMSS[string(i)]], output, "MSS_flagP="*string(i))

        nccreate(output, "Y_flagP="*string(i),
                    "dim", length(Ymdr[string(i)]),
                    atts = Y_atts)
        ncwrite(Ymdr[string(i)], output, "Y_flagP="*string(i))
    end

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

    yeq = zeros(length(Hb))
    Y = zeros(length(Hb))
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

    end

    return Y
end


