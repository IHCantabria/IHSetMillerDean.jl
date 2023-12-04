"""
#######################################
# Miller and Dean (2004) model        #
# hb: wave breaking height            #
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
using IHSetUtils

function MILLER_DEAN_extr(hb, depthb, sl, Y0, dt, D50, Hberm, kero, kacr, Yi, flagP = 1, Omega = 0)
    if flagP == 1
        kero = fill(kero, length(hb))
        kacr = fill(kacr, length(hb))
    elseif flagP == 2
        kero = kero .* hb .^ 2
        kacr = kacr .* hb .^ 2
    elseif flagP == 3
        kero = kero .* hb .^ 3
        kacr = kacr .* hb .^ 3
    elseif flagP == 4
        kero = kero .* Omega
        kacr = kacr .* Omega
    end

    yeq = similar(hb)
    Y = similar(hb)
    wl = 0.106 .* hb .+ sl
    Wast = wast(depthb, D50)
    @views yeq .= Y0 .- Wast .* wl ./ (Hberm .+ depthb)
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

    return Y, yeq, acr_count, ero_count
end

function MILLER_DEAN(hb, yeq, Y0, dt, kero, kacr, Yi, flagP = 1, Omega = 0)
    
    if flagP == 1
        kero = fill(kero, length(yeq))
        kacr = fill(kacr, length(yeq))
    elseif flagP == 2
        kero = fill(kero, length(hb)) .* hb .^ 2
        kacr = fill(kacr, length(hb)) .* hb .^ 2
    elseif flagP == 3
        kero = fill(kero, length(hb)) .* hb .^ 3
        kacr = fill(kacr, length(hb)) .* hb .^ 3
    elseif flagP == 4
        kero = fill(kero, length(yeq)) .* Omega
        kacr = fill(kacr, length(yeq)) .* Omega
    end

    # yeq = similar(hb)
    Y = similar(yeq)
    # wl = 0.106 .* hb .+ sl
    @views yeq .= Y0 .- yeq

    for i in eachindex(yeq)
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
