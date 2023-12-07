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
    mods = wrkDir*"Modules\\"
    dats = wrkDir*"Data\\"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Loading datasets...")

    wavF = NCDataset(dats*"wav.nc")
    ensF = NCDataset(dats*"ens.nc")
    parF = NCDataset(dats*"par.nc")

    configF = NCDataset(dats*"config.nc")

    println("Unpacking datasets...")

    dt = configF["dt"][:][1]
    
    yi = collect(skipmissing(configF["yi"][:]))

    kacr = collect(skipmissing(parF["kacr"][:]))
    kero = collect(skipmissing(parF["kero"][:]))
    Y0 = collect(skipmissing(parF["Y0"][:]))
    
    brk, angBati, depth, Hberm, D50 = configF["brk"][:][1], configF["angBati"][:][1], configF["depth"][:][1], configF["Hberm"][:][1], configF["D50"][:][1]

    if brk == 1
        
        Hs, Tp, θ_w = collect(skipmissing(wavF["Hs"][:])), collect(skipmissing(wavF["Tp"][:])), collect(skipmissing(wavF["Dir"][:]))

        auxAng, auxDepth = similar(Hs), similar(Hs)
        auxAng .= angBati
        auxDepth .= depth

        println("Breaking waves by linear theory...")
        Hb, θ_b, depthb = WAV.BreakingPropagation(Hs, Tp, θ_w, auxAng, auxDepth, "spectral")
    else
        Hb, tp, hs, depthb = wavF["Hb"][:], wavF["Tp"][:], wavF["Hs"][:], wavF["hb"][:]
    end
    
    close(wavF)
    close(ensF)
    close(configF)
    close(parF)

    println("Datasets closed...")

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


    Y_t = MileerDean(Hb, depthb, sl, Y0, dt, D50, Hberm, kero, kacr, Yi, flagP = 1, Omega = 0)

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

    return Y, yeq
end