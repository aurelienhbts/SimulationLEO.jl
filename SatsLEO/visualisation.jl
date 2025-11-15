using LinearAlgebra, Plots

export show_coverage_heatmap, plot_constellation!, plot_earth

"""
show_coverage_heatmap(sats,t,eps_deg)

Affiche les zones couvertes par les satellites. Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol.
"""
function show_coverage_heatmap(sats,t,eps_deg)
    lats = collect(-90:1:90)
    lons = collect(-180:1:180)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]

    M = falses(length(lats), length(lons))

    for (i, lat) in enumerate(lats), (j, lon) in enumerate(lons)
        M[i, j] = any(r -> visible(r, lat, lon, eps_deg), r_ecef)
    end
	p = heatmap(lons, lats, Int.(M))
	p = plot!(p, xlabel = "Longitude [°]", ylabel = "Latitude [°]", title = "Zones couvertes par $(length(sats)) satellites")
	
    return plot!(p, aspect_ratio=1, colorbar=false, framestyle=:none)
end

"""
plot_earth(;Rearth=Re)

Plot une sphère (la Terre) avec des lignes pour les lattitudes et l'équateur en bleu.
"""
function plot_earth(;Rearth=Re)
	# Hémisphère arrière de la Terre (longitudes ~ [-π, 0])
	u = range(-pi/2, pi/2, 60); v1 = range(-pi, 0, 60)
	xs = [Re*cos(ui)*cos(vi) for ui in u, vi in v1]
	ys = [Re*cos(ui)*sin(vi) for ui in u, vi in v1]
	zs = [Re*sin(ui)         for ui in u, vi in v1]
	p = surface(xs, ys, zs, color=:lightblue, opacity=0.03, linecolor=:transparent)
	
	# Hémisphère avant de la Terre (longitudes ~ [0, π])
	v2 = range(0, pi, 60)
	xs = [Re*cos(ui)*cos(vi) for ui in u, vi in v2]
	ys = [Re*cos(ui)*sin(vi) for ui in u, vi in v2]
	zs = [Re*sin(ui)         for ui in u, vi in v2]
	p = surface!(p, xs, ys, zs, color=:lightblue, opacity=0.05, linecolor=:transparent)

	# Ajout de lignes pour les latitudes
	for lat_deg in -80:10:-10
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:black, alpha=0.3)
    end
	for lat_deg = 0
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:blue) # Equateur en bleu
    end
	for lat_deg in 10:10:80
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:black, alpha=0.3)
    end
	
	return plot!(p, aspect_ratio=:equal, legend=false, colorbar=false, size=(600,600), framestyle=:none, ticks=:none, camera=(0,5))
end

"""
plot_constellation!(sats,t;Rearth=Re)

Plot les positions instantanées des satellites en 3D autour de la Terre.
"""
function plot_constellation!(sats,t;Rearth=Re)
	p = plot_earth()
	
	θs = range(0, 2π, 100)
    for s in sats
        R = R3(s.Ω) * R1(s.i)  # orientation du plan orbital
        orb = [R * [s.a * cos(θ), s.a * sin(θ), 0.0] for θ in θs]
        Xorb = [r[1] for r in orb]; Yorb = [r[2] for r in orb]; Zorb = [r[3] for r in orb]
    p = plot3d!(p, Xorb, Yorb, Zorb, lw=0.5, c=:gray, alpha=0.5)  # Tracé des plans orbitaux
    end
	
    r_ecef = [eci_pos(s, t) for s in sats] # Position des satellites
    X = [r[1] for r in r_ecef]
    Y = [r[2] for r in r_ecef]
    Z = [r[3] for r in r_ecef]
    p = scatter3d!(p, X, Y, Z, marker=:circle, ms=3.5, color=:orange, title="Visualisation des satellites en t = $(Int(t))s") # Plot des satellites

	return p
end