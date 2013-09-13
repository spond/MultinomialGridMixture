#!/usr/bin/env julia --color=yes

using ArgParse
using FastaIO

function countmsa(msa)
    counts = nothing
    ncols = 0
    for (name, seq) in FastaReader(msa)
        if counts == nothing
            ncols = length(seq)
            counts = zeros(Int64, (ncols, 4))
        elseif length(seq) != ncols
            error("provided file is not an MSA! length($name) = $(length(seq)), not $ncols!")
        end
        for (i, l) in enumerate(seq)
            L = uppercase(l)
            j = 0
            if L == 'A'
                j = 1
            elseif L == 'C'
                j = 2
            elseif L == 'G'
                j = 3
            elseif L == 'T'
                j = 4
            end
            if j > 0
                counts[i, j] += 1
            end
        end
    end
    counts
end

function rategrid(n::Integer)
    rates = 1 / (2 .^ linspace(2, 10, n))
    rg = Array(Float64, (4 * length(rates) ^ 3, 4))
    i = 1
    for a in rates
        for b in rates
            for c in rates
                M = 1 - (a + b + c)
                rg[i    , :] = [M a b c]
                rg[i + 1, :] = [a M b c]
                rg[i + 2, :] = [a b M c]
                rg[i + 3, :] = [a b c M]
                i += 4
            end
        end
    end
    rg
end

lfact_cache = Float64[log(1)]

function lfact(N::Integer)
    if N == 0
        return 0.0
    end
    n = length(lfact_cache)
    if n < N
        resize!(lfact_cache, N)
        for i in (n + 1):N
            lfact_cache[i] = lfact_cache[i - 1] + log(i)
        end
    end
    lfact_cache[N]
end

function lmc(data::Array{Int64})
    c = lfact(sum(data))
    for d in data
        c -= lfact(d)
    end
    c
end

const smin = log(realmin(Float64))

function gridscores(counts::Array{Int64,2}, rates::Array{Float64,2})
    npoints, _ = size(rates)
    nsites, _ = size(counts)
    conditionals = Array(Float64, (npoints, nsites))
    scalers = Array(Float64, (nsites,))
    for i in 1:nsites
        m = realmin(Float64)
        for j in 1:npoints
            s = lmc(counts[i, :])
            for k in 1:4
                s += counts[i, k] * log(rates[j, k])
            end
            if s > m
                m = s
            end
            conditionals[j, i] = s
        end
        scalers[i] = m
        for j in 1:npoints
            s = conditionals[j, i] - m
            if s < smin
                conditionals[j, i] = 0
            else
                conditionals[j, i] = exp(s)
            end
        end
    end
    (conditionals, scalers)
end

function lddens(weights::Array{Float64,2}, alpha::Float64)
    if min(weights) <= 1e-10
        return 1e-10
    end
    if alpha == 1
        return 0
    end
    n = length(weights)
    r = lgamma(n * alpha) - n * lgamma(alpha)
    for i in 1:n
        r += log(weights[i]) * (alpha - 1)
    end
    r
end

function jll(conditionals::Array{Float64,2}, scalers::Array{Float64,1}, weights::Array{Float64,2}, alpha::Float64)
    ll = sum(log(weights * conditionals) + scalers')
    dir = lddens(weights, alpha)
    (ll, dir)
end

# alpha was once called concentration_parameter
function mcmc(conditionals, scalers, ntotal=10_000_000, nburnin=5_000_000, expected_nsamples=100, alpha=0.5)
    npoints, nsites = size(conditionals)

    normalized_by_site = ones((1, npoints)) * conditionals
    normalized_weights = conditionals * ((1 / normalized_by_site) .* eye(nsites))
    sum_by_site = normalized_weights * ones((nsites, 1))

    weights = transpose((rand((npoints, 1)) .* (1.2 * sum_by_site - 0.8 * sum_by_site)) + 0.8 * sum_by_site)
    weights /= sum(weights)

    stepsize = max(min(1 / npoints, 1 / nsites), median(weights))

    nsample = div(ntotal - nburnin, expected_nsamples)
    sampled_weights = zeros(Float64, (expected_nsamples, npoints))
    sampled_likelihoods = zeros(Float64, (expected_nsamples,))

    current_site_likelihoods = weights * conditionals
    current_site_logsum = sum(log(current_site_likelihoods))
    current_ll, current_dir = jll(conditionals, scalers, weights, alpha)

    step = 1
    accepted_steps = 0
    mean_sampled_ll = 0
    sample_idx = 0

    diffvec = Array(Float64, (1, nsites))
    rng = 1:npoints
    llstr = "(burning in)"

    tic()

    for step in 1:ntotal

        i = rand(rng)
        j = rand(rng)
        while i == j
            j = rand(rng)
        end

        change = rand() * stepsize

        if weights[i] > change
            for k in 1:nsites
                diffvec[k] = (conditionals[j, k] - conditionals[i, k]) * change
            end

            lldiff = 0
            for k in 1:nsites
                lldiff += log(current_site_likelihoods[k] + diffvec[k])
            end
            lldiff -= current_site_logsum
            dirdiff = (alpha - 1) * (log((weights[i] - change) / weights[i]) + log((weights[j] + change) / weights[j]))
            movecost = exp(lldiff + dirdiff)

            if rand() <= movecost
                current_ll += lldiff
                current_dir += dirdiff

                current_site_likelihoods += diffvec
                current_site_logsum += lldiff

                weights[i] -= change
                weights[j] += change

                accepted_steps += 1
            end
        end

        if step > nburnin
            if (step - nburnin) % nsample == 0
                sample_idx += 1
                sampled_weights[sample_idx, :] = weights
                sampled_likelihoods[sample_idx] = current_ll
                mean_sampled_ll += current_ll
                llstr = "$(trunc(mean_sampled_ll / sample_idx, 3))"
            end
        end

        if step % nsample == 0
            progress("running MCMC chain: step = $step/$ntotal, mean logL = $llstr, acceptance rate = $(trunc(accepted_steps / step, 3)) .. ")
        end
    end

    done()

    (sampled_likelihoods, sampled_weights)
end

function loadrates()
    rates3 = readdlm("data/rates.txt", ' ', Float64)
    nrates, _ = size(rates3)
    rates = Array(Float64, (nrates, 4))
    rates[:, 1:3] = rates3
    for i in 1:nrates
        rates[i, 4] = 1.0 - sum(rates[i, 1:3])
    end
    rates
end

function postproc(conditionals::Array{Float64,2}, ws::Array{Float64,2})
    nsites, npoints = size(conditionals)
    nsamples, _ = size(ws)
    priors = zeros(Float64, (1, npoints))
    for i in 1:nsamples
        priors += ws[i, :]
    end
    priors = priors' / nsamples
    priors /= sum(priors)
    normalization = conditionals * priors
    posteriors = conditionals * (priors .* eye(npoints))
    ((1 / normalization) .* eye(nsites)) * posteriors
end

function callvariants(
        counts::Array{Int64,2},
        rates::Array{Float64,2},
        posteriors::Array{Float64,2},
        threshold::Float64=0.1,
        posterior_threshold::Float64=0.999)
    nsites, npoints = size(posteriors)
    posterior_prob = posteriors * map(x -> x > threshold, rates)
    for i in 1:nsites
        s = Set{Char}()
        for (j, c) in enumerate("ACGT")
            if posterior_prob[i, j] > posterior_threshold
                push!(s, c)
            end
        end
        if length(s) > 1
            println("$i: $(join(sort!(collect(s)))),\t[$(join(counts[i, :], '\t'))]")
        end
    end
end

function done()
    print_with_color(:blue, STDERR, "done, took $(toq())s")
    println("")
end
function progress(args...)
    print_with_color(:blue, STDERR, "\rINFO: ", args...)
end
function update(args...)
    print_with_color(:blue, STDERR, "INFO: ", args...)
    tic()
end

function main(args)
    s = ArgParseSettings()
    s.description = "call variants using a multinomial model sampled by MCMC"
    @add_arg_table s begin
        "--grid-density", "-g"
            arg_type = Int64
            default = 14
        "--chain-length", "-c"
            arg_type = Int64
            default = 2_000_000
        "--burnin-length", "-b"
            arg_type = Int64
            default = 1_000_000
        "--target-rate", "-t"
            arg_type = Float64
            default = 0.01
        "--posterior-threshold", "-p"
            arg_type = Float64
            default = 0.95
        "msa"
            required = true
    end
    parsed_args = parse_args(args, s)

    update("loading MSA .. ")
    counts = countmsa(parsed_args["msa"])
    done()

    update("generating rate grid .. ")
    # rates = loadrates()
    rates = rategrid(parsed_args["grid-density"])
    done()

    update("computing grid scores .. ")
    conditionals, scalers = gridscores(counts, rates)
    done()

    lls, ws = mcmc(conditionals, scalers, parsed_args["chain-length"], parsed_args["burnin-length"])
    posteriors = postproc(conditionals', ws)

    callvariants(counts, rates, posteriors, parsed_args["target-rate"], parsed_args["posterior-threshold"])
end

main(ARGS)
