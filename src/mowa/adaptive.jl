struct AdaptiveParameters{δτ⁻_T, δτ⁺_T, δΔτ⁻_T, δΔτ⁺_T, Nᵣ_T, fineupdate_T} <: AbstractAdaptiveParameters
    δτ⁻::δτ⁻_T
    δτ⁺::δτ⁺_T
    δΔτ⁻::δΔτ⁻_T
    δΔτ⁺::δΔτ⁺_T
    Nᵣ::Nᵣ_T
    fineupdate::fineupdate_T
end
AdaptiveParameters(; δτ⁻=0.5, δτ⁺=1.25, δΔτ⁻=0.25, δΔτ⁺=2.0, Nᵣ=10, fineupdate=true) = AdaptiveParameters(δτ⁻, δτ⁺, δΔτ⁻, δΔτ⁺, Nᵣ, fineupdate)
