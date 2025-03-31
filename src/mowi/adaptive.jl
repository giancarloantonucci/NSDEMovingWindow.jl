struct AdaptiveMoWiParameters{δτ⁻_T, δτ⁺_T, δΔτ⁻_T, δΔτ⁺_T, δw⁻_T, δw⁺_T, R_T, fineupdate_T} <: AbstractAdaptiveMoWiParameters
    δτ⁻::δτ⁻_T
    δτ⁺::δτ⁺_T
    δΔτ⁻::δΔτ⁻_T
    δΔτ⁺::δΔτ⁺_T
    δw⁻::δw⁻_T
    δw⁺::δw⁺_T
    R::R_T
    fineupdate::fineupdate_T
end

AdaptiveMoWiParameters(; δτ⁻=0.75, δτ⁺=1.75, δΔτ⁻=0.25, δΔτ⁺=2.0, δw⁻=2.0, δw⁺=0.25, R=10, fineupdate=true) = AdaptiveMoWiParameters(δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, δw⁻, δw⁺, R, fineupdate)

increaseweights() = ()
increasewindowlength() = ()
increaseshifting() = ()

decreaseweights() = ()
decreasewindowlength() = ()
decreaseshifting() = ()


# Stretch

struct StretchParameters{δτ⁻_T, δτ⁺_T, R_T, fineupdate_T} <: AbstractAdaptiveMoWiParameters
    δτ⁻::δτ⁻_T
    δτ⁺::δτ⁺_T
    R::R_T
    fineupdate::fineupdate_T
end

StretchParameters(; δτ⁻=0.5, δτ⁺=1.25, R=10, fineupdate=false) = StretchParameters(δτ⁻, δτ⁺, R, fineupdate)

# Leap

struct LeapParameters{δΔτ⁻_T, δΔτ⁺_T, R_T} <: AbstractAdaptiveMoWiParameters
    δΔτ⁻::δΔτ⁻_T
    δΔτ⁺::δΔτ⁺_T
    R::R_T
end

LeapParameters(; δΔτ⁻=0.25, δΔτ⁺=2.0, R=10) = LeapParameters(δΔτ⁻, δΔτ⁺, R)

struct ZoomParameters{δw⁻_T, δw⁺_T, R_T} <: AbstractAdaptiveMoWiParameters
    δw⁻::δw⁻_T
    δw⁺::δw⁺_T
    R::R_T
end

ZoomParameters(; δw⁻=2.0, δw⁺=0.25, R=10) = ZoomParameters(δw⁻, δw⁺, R)
