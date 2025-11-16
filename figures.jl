## Script utilisé pour faire les figures illustrant les différents phénomènes

# Pour le fait que la couverture change de façon périodique:
let
    h_km = 800
    eps_deg = 10
    i_deg = 30
    P = 2
    S = 178
    F = 0
    a = Re + h_km*1e3

	sats=walker_delta(P,S,F,i_deg,a)
	
	Ts = 0:100:86800
	covs = [coverage_fraction(sats, T, -90, 90, eps_deg) for T in Ts]
	
	plot = plot(Ts, round.(covs; digits=2),
	    xlabel = "Temps (s)",
	    ylabel = "Couverture (%)",
	    legend = false,
	    title = "Évolution de la couverture en fonction du temps\nWalker-Delta P=$(P), S=$(S), i=$(i_deg)°",
	    markersize = 3,
	    color = :blue)
	savefig(plot, "./coverage_fraction_période")
end

# Pour montrer que la couverture moyenne sur une période converge avec un certain nombre d'échantillons:
let
    h_km = 800
    eps_deg = 10
    i_deg = 30
    P = 2
    S = 5
    F = 0
    a = Re + h_km*1e3

    sats = walker_delta(P, S, F, i_deg, a)

    N = 200
    ns = 2 .* (1:N)
    vals = zeros(N)

    for k in 1:N
        vals[k] = mean_coverage_fraction(sats, -i_deg, i_deg, eps_deg; n=ns[k])
    end

    plot = plot(
        ns, vals;
        xlabel = "Nombre d'échantillons temporels",
        ylabel = "Couverture moyenne (%)",
        title = "Convergence de la couverture moyenne\nWalker-Delta P=$(P), S=$(S), i=$(i_deg)°",
        lw = 2,
        markershape = :circle,
        markerstrokewidth = 0,
        legend = false,
        grid = true,
    )
	savefig(plot,"./convergence_mean_coverage")
end