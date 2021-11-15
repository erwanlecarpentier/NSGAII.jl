using NSGAII
using CartesianGeneticProgramming
using Cambrian
using UnicodePlots
using LinearAlgebra

function out(plt)
    println(IOContext(stdout, :color=>true), plt)
end

cfg = get_config("cfg/cgp.yaml")
n_population = cfg.n_population
n_gen = cfg.n_gen
p_mut = 1.0 # Probability of mutation at each generation for each individual
plotevery = 100
obj1, obj2 = Vector{Float64}(), Vector{Float64}()
seeding_ind = Vector{Float64}() # Individuals for init

x = [i for i in 0:0.1:1]
f1(x) = 0.25 * (2.0 + cos(2.0*π*x))
f2(x) = 0.25 * (2.0 + cos(2.5*π*x))
y1 = f1.(x)
y2 = f2.(x)
out(lineplot([f1, f2], 0, 1, border=:dotted))

# Fitness function (normalized)
function my_fitness(chromosome::Vector{Float64})
    ind = CGPInd(cfg, chromosome)
    y_hat = zeros(length(y1))
    @inbounds for i in eachindex(y1)
        y_hat[i] = process(ind, [x[i]])[1]
    end
    o1 = norm(y1 - y_hat)
    o2 = norm(y2 - y_hat)
    o1, o2
end

# Mutation function
function my_mutate!(chromosome::Vector{Float64})
    ind = CGPInd(cfg, chromosome)
    child = goldman_mutate(cfg, ind)
    chromosome .= child.chromosome
end

# Crossover function
function my_crossover!(
    indiv_a_chromosome::Vector{Float64},
    indiv_b_chromosome::Vector{Float64},
    child_a_chromosome::Vector{Float64},
    child_b_chromosome::Vector{Float64}
)
    nothing # No crossover
end

# Chromosome initialization function
function my_init()
    ind = CGPInd(cfg)
    ind.chromosome
end

# Display/log some results at each `plotevery` generation
function my_fplot(pop::AbstractArray)
    o1 = [ind.y[1] for ind in pop]
    o2 = [ind.y[2] for ind in pop]
    println(IOContext(stdout, :color=>true),
            scatterplot(o1, o2, title = "Paretto front",
            xlim=[0,1], ylim=[0,1]))
    push!(obj1, minimum(o1))
    push!(obj2, minimum(o2))
end

# Run
paretto = nsga(#_max(
    n_population, n_gen, my_fitness, my_init;
    pmut = p_mut,
    fmut = my_mutate!,#ind -> my_mutate!(ind.x),
    fcross = my_crossover!,#(indivs...) -> my_crossover!(),#getproperty.(indivs, :x)...),
    seed = seeding_ind,
    fplot = my_fplot,
    plotevery = plotevery,
    showprogress = true
)

# Plot
println(IOContext(stdout, :color=>true), lineplot(obj1, title = "Objective 1", ylim=[0,1]))
println(IOContext(stdout, :color=>true), lineplot(obj2, title = "Objective 2", ylim=[0,1]))

# Plot all computed functions
for i in eachindex(paretto)
    ind = CGPInd(cfg, paretto[i].x)
    f(x) = process(ind, [x])[1]
    out(lineplot([f1, f2, f], 0, 1, border=:dotted))
end
