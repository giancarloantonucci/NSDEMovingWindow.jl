struct AdaptiveMoWiParameters{δτ⁻_T, δτ⁺_T, δΔτ⁻_T, δΔτ⁺_T, R_T, fineupdate_T} <: AbstractAdaptiveMoWiParameters
    δτ⁻::δτ⁻_T
    δτ⁺::δτ⁺_T
    δΔτ⁻::δΔτ⁻_T
    δΔτ⁺::δΔτ⁺_T
    R::R_T
    fineupdate::fineupdate_T
end

AdaptiveMoWiParameters(; δτ⁻=0.5, δτ⁺=1.25, δΔτ⁻=0.25, δΔτ⁺=2.0, R=10, fineupdate=true) = AdaptiveMoWiParameters(δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, R, fineupdate)

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
