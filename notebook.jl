### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 9deb23ea-b9c1-11f0-8f30-17ad0f50be39
begin
	using Pkg
	Pkg.activate("/home/kali/Documents/Projet ES313")
	using Plots
end

# ╔═╡ 12648801-3cd3-4060-89ad-ec012fdfd5a7
begin
	using Orbits
	using Unitful, UnitfulAstro
	using UnitfulRecipes
end

# ╔═╡ 973f1398-fa7c-479a-a895-efae3c395458
begin
	using LinearAlgebra

	# Constantes
	μ = 398600.4418       # km^3/s^2
	r_earth = 6371.0      # km
	altitude = 500.0      # km
	r = r_earth + altitude
	v = sqrt(μ / r)
	
	# Nombre de satellites
	N = 20
	
	# Positions des satellites
	positions = [r * [cos(θ), sin(θ), 0.0] for θ in LinRange(0, 2π, N+1)[1:end-1]]
	x_sat = [p[1] for p in positions]
	y_sat = [p[2] for p in positions]
	z_sat = [p[3] for p in positions]
	
	# 🌍 Génération de la Terre comme sphère
	θ = range(0, 2π, length=50)
	ϕ = range(0, π, length=50)
	x_earth = [r_earth * sin(ϕi) * cos(θj) for ϕi in ϕ, θj in θ]
	y_earth = [r_earth * sin(ϕi) * sin(θj) for ϕi in ϕ, θj in θ]
	z_earth = [r_earth * cos(ϕi) for ϕi in ϕ, _ in θ]
	
	# 📊 Visualisation 3D
	#plotlyjs()  # ou gr(), selon ton backend préféré
	surface(x_earth, y_earth, z_earth, color=:lightblue, alpha=0.5, label="Terre")
	scatter3d!(x_sat, y_sat, z_sat, markersize=4, color=:red, label="Satellites LEO")
	plot!(xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", title="Réseau de satellites LEO autour de la Terre", legend=:topright)
end

# ╔═╡ fe06a036-485c-4a2d-8480-8ab3a8f0dfba
begin
	using Geodesics;
	
	# Coordonnées de départ (Trafalgar Square, Londres)
	lon = deg2rad(0.1281);   # longitude en radians
	lat = deg2rad(51.5080);  # latitude en radians
	az = deg2rad(45.0);      # azimut (direction) en radians
	dist = 30_000.0;         # distance en mètres
	
	# Paramètres ellipsoïdaux WGS84
	a = Geodesics.EARTH_R_MAJOR_WGS84;  # rayon équatorial
	f = Geodesics.F_WGS84;              # aplatissement
	
	# Calcul du point d’arrivée
	lon2, lat2 = Geodesics.forward(lon, lat, az, dist, a, f);
	
	# Résultat en degrés
	println("Destination : longitude = $(rad2deg(lon2)), latitude = $(rad2deg(lat2))");
	
end

# ╔═╡ 4f22e607-6967-4028-bc3a-47737d1e4c15
let
	using GeneralAstrodynamics
	
	orbit = rand(R2BPOrbit)
	trajectory = propagate(orbit, orbital_period(orbit))
	
	furnsh(
	    de440s(),                   # position and velocity data for nearby planets
	    latest_leapseconds_tls(),   # timekeeping, parsing epochs
	    gm_de440(),                 # mass parameters for major solar system bodies
	    pck00011(),                 # physical properties of major solar system bodies
	)
	
	μ = reduced_mass(
	  gm("earth"),
	  gm("moon"),
	)
	
	orbit, T = let
	  u, T = halo(μ, 2; amplitude=1e-2)
	
	  CR3BPOrbit(CartesianState(u), CR3BParameters(μ)), T
	end
	
	trajectory = propagate(orbit, T)
end

# ╔═╡ 17749c20-280d-4535-a540-fea79c03690d
html"""
 <! -- this adapts the width of the cells to display its being used on -->
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ 629be12c-8cff-44b1-98cb-6a0d48b4df44
begin

# orbital params for SAO 136799
distance = inv(6.92e-3)u"pc"

orbit = KeplerianOrbit(;
    period = 40.57u"yr",
    ecc = 0.42,
    Omega = 318.6u"°",
    tp = 1972.12u"yr",
    incl = 54.7u"°",
    a = 0.154u"arcsecond" * distance |> u"AU",
    omega = 72.6u"°",
)

# get position at specific time
t = 2022.134u"yr"
pos = relative_position(orbit, t)
ra_off, dec_off = @. pos[1:2] / distance |> u"arcsecond"
end

# ╔═╡ 8c5668cb-94cd-4533-9bbf-aa6e98eac08c
begin
	# plot using Unitful recipes
	plot(orbit; distance, lab="", leg=:topleft)
	scatter!([0u"arcsecond" ra_off], [0u"arcsecond" dec_off],
	          c=[:black 1], m=[:+ :o], lab=["SAO 136799A" "B ($t)"])
end

# ╔═╡ a5cfdd08-9e89-4410-bc2f-a74d26bdd149
let
	# Constantes;
	μ = 398600.4418;       # km^3/s^2;
	r_earth = 6371.0;      # km;
	altitude = 2000.0;      # km;
	r = r_earth + altitude;
	v = sqrt(μ / r);       # vitesse orbitale [km/s];
	
	# Réseau de satellites;
	n_planes = 3;
	sats_per_plane = 4;
	ΔΩ = 2π / n_planes;
	Δν = 2π / sats_per_plane;
	inclination_deg = 53.0;
	i = deg2rad(inclination_deg);
	
	# Orbites;
	orbite_points = 100;
	positions = [];
	orbite_lines = [];
	
	for p in 0:n_planes-1;
	    Ω = p * ΔΩ;
	    for s in 0:sats_per_plane-1;
	        ν_offset = s * Δν;
	        orbite = [];
	        for ν in LinRange(0, 2π, orbite_points);
	            ν_total = ν + ν_offset;
	            x = r * (cos(Ω) * cos(ν_total) - sin(Ω) * cos(i) * sin(ν_total));
	            y = r * (sin(Ω) * cos(ν_total) + cos(Ω) * cos(i) * sin(ν_total));
	            z = r * (sin(i) * sin(ν_total));
	            push!(orbite, [x, y, z]);
	        end;
	        push!(orbite_lines, orbite);
	        push!(positions, orbite[1]);  # position initiale du satellite;
	    end;
	end;
	
	# Coordonnées satellites;
	x_sat = [r[1] for r in positions];
	y_sat = [r[2] for r in positions];
	z_sat = [r[3] for r in positions];
	
	# Coordonnées orbites;
	x_orb = [[r[1] for r in orb] for orb in orbite_lines];
	y_orb = [[r[2] for r in orb] for orb in orbite_lines];
	z_orb = [[r[3] for r in orb] for orb in orbite_lines];
	
	# Sphère Terre;
	θ = range(0, 2π, length=50);
	ϕ = range(0, π, length=50);
	x_earth = [r_earth * sin(ϕi) * cos(θj) for ϕi in ϕ, θj in θ];
	y_earth = [r_earth * sin(ϕi) * sin(θj) for ϕi in ϕ, θj in θ];
	z_earth = [r_earth * cos(ϕi) for ϕi in ϕ, _ in θ];
	
	# Visualisation 3D;
	#plotlyjs();  # ou gr() selon ton backend;
	surface(x_earth, y_earth, z_earth, color=:lightblue, alpha=0.4, label="Terre");
	scatter3d!(x_sat, y_sat, z_sat, markersize=4, color=:red, label="Satellites");
	
	for i in 1:length(x_orb);
	    plot3d!(x_orb[i], y_orb[i], z_orb[i], color=:gray, label=false);
	end;
	
	plot!(axis=false, grid=false, ticks=nothing, legend=false, colorbar=false, background_color=:black, title="", aspect_ratio=:equal, size=(600,600));

end

# ╔═╡ fe273246-b4db-4909-8029-b141b8137ba3
let
	# Coordonnées en radians
	lon1 = deg2rad(0.0);     # Greenwich
	lat1 = deg2rad(51.5);    # Londres
	
	lon2 = deg2rad(2.35);    # Paris
	lat2 = deg2rad(48.85);
	
	# Paramètres WGS84
	a = Geodesics.EARTH_R_MAJOR_WGS84;
	f = Geodesics.F_WGS84;
	
	# Distance et azimut
	s, az1, az2 = Geodesics.inverse(lon1, lat1, lon2, lat2, a, f);
	println("Distance Londres–Paris : $(s/1000) km");
	
end

# ╔═╡ Cell order:
# ╟─17749c20-280d-4535-a540-fea79c03690d
# ╠═9deb23ea-b9c1-11f0-8f30-17ad0f50be39
# ╠═12648801-3cd3-4060-89ad-ec012fdfd5a7
# ╠═629be12c-8cff-44b1-98cb-6a0d48b4df44
# ╟─8c5668cb-94cd-4533-9bbf-aa6e98eac08c
# ╟─973f1398-fa7c-479a-a895-efae3c395458
# ╟─a5cfdd08-9e89-4410-bc2f-a74d26bdd149
# ╠═fe06a036-485c-4a2d-8480-8ab3a8f0dfba
# ╠═fe273246-b4db-4909-8029-b141b8137ba3
# ╠═4f22e607-6967-4028-bc3a-47737d1e4c15
