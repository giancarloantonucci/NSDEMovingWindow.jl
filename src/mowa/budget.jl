mutable struct Budget{Kmax_T} <: AbstractBudget
    Kmax::Kmax_T
end
Budget(; Kmax=10) = Budget(Kmax)
