using SummationByPartsOperators, OrdinaryDiffEq
using LaTeXStrings; using Plots: Plots, plot, plot!, savefig

# general parameters
xmin = -1.
xmax = +1.
tspan = (0., 8.0)
u0_func(x) = exp(-20x^2)
v0_func(x) = zero(x)
# HomogeneousNeumann, HomogeneousDirichlet, and NonReflecting BCs are available
left_bc  = Val(:HomogeneousNeumann)
right_bc = Val(:HomogeneousDirichlet)

# setup spatial semidiscretization
D2 = derivative_operator(MattssonSvärdShoeybi2008(), derivative_order=2,
                         accuracy_order=2, xmin=xmin, xmax=xmax, N=101)
semi = WaveEquationNonperiodicSemidiscretization(D2, left_bc, right_bc)
ode = semidiscretize(v0_func, u0_func, semi, tspan)


# solve second-order ODE using a Runge-Kutta-Nyström method
sol = solve(ode, DPRKN6(), saveat=range(first(tspan), stop=last(tspan), length=200))