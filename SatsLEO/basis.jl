using LinearAlgebra

## Constantes
const μ = 3.986004418e14 	# Paramètre gravitationnel terrestre (m³/s²)
const Re = 6.371e6 		    # Rayon moyen de la Terre (m)
const ωe = 7.2921150e-5 	# Vitesse de rotation de la Terre (rad/s)

## Structure pour stocker les satellites (définis par a, i, Ω et M0).
struct Sat
    a::Float64   # Demi-grand axe (m) → distance moyenne au centre de la Terre
    i::Float64   # Inclinaison orbitale (rad) → angle entre le plan orbital et l’équateur
    Ω::Float64   # Longitude du nœud ascendant (rad) → orientation du plan orbital autour de la Terre
    M0::Float64  # Anomalie moyenne initiale (rad) → position du satellite sur son orbite
end

## deg2rad et rad2deg functions
deg2rad(x)=x*pi/180
rad2deg(x)=180*x/pi

## Matrices de rotation
R1(θ)=[1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] # Autour de l'axe X → Sert pour incliner le plan orbital d’un angle i (l'inclinaison).
R3(θ)=[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] # Autour de l'axe Z → Sert pour tourner le plan orbital d’un angle Ω (ascension du nœud) ou de l’anomalie vraie ν.

"""
ecef_from_eci(r,t)

Convertit un vecteur de position d’un satellite du repère ECI (Earth-Centered Inertial) vers le repère ECEF (Earth-Centered Earth-Fixed).
"""
function ecef_from_eci(r,t)
	return R3(-ωe*t)*r # page 129 - SE216 (v.Aout 2024)
end

"""
latlon_from_ecef(r)

Convertit un vecteur position exprimé dans le repère ECEF (Earth-Centered Earth-Fixed) en latitude et longitude (degrés).
"""
function latlon_from_ecef(r)
    x, y, z = r # Vecteur position
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre
	
    ϕ = asin(z / ρ) 	# Latitude -> angle entre le vecteur position et le plan équatorial.  
    λ = atan(y, x) 		# Longitude -> angle dans le plan équatorial entre l’axe X (Greenwich) et la projection du vecteur sur ce plan.

    return rad2deg(ϕ), rad2deg(λ)
end

"""
eci_pos(sat,t)

Calcule la position d’un satellite dans le repère inertiel (ECI) à un instant t, en supposant une orbite circulaire.
"""
function eci_pos(sat,t;mu=μ)
    n = sqrt(mu / sat.a^3) 			    # Vitesse angulaire du satellite d'après la 3e loi de Kepler
    u = sat.M0 + n * t 				    # Position angulaire du satellite sur son orbite (Position initiale (M₀) + n * t)
    R = R3(sat.Ω) * R1(sat.i) * R3(u) 	# Matrice de transformation (Formule 7.35 - SE216 (v.Aout 2024))
    return R * [sat.a, 0, 0] 			# Position finale du satellite
end

"""
walker_delta(P,S,F,i_deg,a)

Initialise une constellation de type Walker-Delta. Utilisée pour répartir uniformément des satellites sur plusieurs plans orbitaux inclinés.
"""
function walker_delta(P,S,F,i_deg,a)
    sats = Sat[]
    for p in 0:P-1, s in 0:S-1 # On itère sur chaque plan orbital et sur chaque satellite par plan
        
		Ω = 2π * p / P # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l’axe z.
		
        # Répartition uniforme des satellites sur chaque orbite (s/S),
        M0 = 2π * (s / S + (F * p) / (S * P)) # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.

        push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
    end
    return sats
end

"""
myconstellation(vec,i_deg,a)

Initialise une constellation de type Walker-delta mais avec plus de flexibilité. vec est un vecteur comportant le nombre de satellite par plan orbital.
"""
function myconstellation(vec,F,i_deg,a)
	sats = Sat[]
	P = length(vec)
	for p in 1:P
		S = vec[p]
		for s in 1:S
			Ω = 2π * p / P # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l’axe z.
			
			# Répartition uniforme des satellites sur chaque orbite (s/S),
        	M0 = 2π * (s / S + (F * p) / (S * P)) # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.
			
			push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
		end
	end
	return sats
end