using NSGAII
using CartesianGeneticProgramming
using Cambrian
using UnicodePlots

cfg = get_config("cfg/cgp.yaml")

n_population = cfg.n_population
n_gen = cfg.n_gen
p_mut = 1.0 # Probability of mutation at each generation for each individual
plotevery = 100
obj1, obj2 = Vector{Float64}(), Vector{Float64}()
seeding_ind = Vector{Float64}() # Individuals for init

# Fitness function (normalized)
function my_fitness(chromosome::Vector{Float64})
    o1 = 1.0 - sum(chromosome) / length(chromosome)
    n = Int64(ceil(0.1*length(chromosome)))
    o2 = sum(sort(chromosome)[1:n]) / n
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
nsga(#_max(
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
