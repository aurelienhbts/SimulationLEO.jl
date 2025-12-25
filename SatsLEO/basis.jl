## Constantes
const μ = 3.986004418e14 	# Paramètre gravitationnel terrestre (m³/s²)
const Re = 6.371e6 		    # Rayon moyen de la Terre (m)s
const ωe = 7.2921150e-5 	# Vitesse de rotation de la Terre (rad/s)

"""
  Sat

Structure représentant un satellite sur une orbite (supposée circulaire) autour de la Terre, définie par ses éléments orbitaux principaux : demi-grand axe, inclinaison, longitude du nœud ascendant et anomalie moyenne initiale.

# Champs
- a  : Demi-grand axe (m), distance moyenne au centre de la Terre.
- i  : Inclinaison orbitale (rad), angle entre le plan orbital et le plan équatorial.
- Ω  : Longitude du nœud ascendant (rad), orientation du plan orbital autour de l'axe de rotation terrestre.
- M0 : Anomalie moyenne initiale (rad), position du satellite sur son orbite à l'instant de référence.
"""
struct Sat
    a::Float64
    i::Float64
    Ω::Float64
    M0::Float64
end

"""
  deg2rad(x)

Convertit un angle en degrés vers des radians.

# Argument
- x : Angle en degrés.

# Valeur retournée
- L'angle correspondant en radians.
"""
deg2rad(x) = x * π / 180


"""
  rad2deg(x)

Convertit un angle en radians vers des degrés.

# Argument
- x : Angle en radians.

# Valeur retournée
- L'angle correspondant en degrés.
"""
rad2deg(x) = 180 * x / pi

"""
  R1(θ)

Matrice de rotation autour de l'axe X d'un angle `θ` (en radians).

Utilisation :
- Rotation d'inclinaison d'un plan orbital (angle `i`).

# Argument
- θ : Angle de rotation en radians.

# Valeur retournée
- Matrice 3x3 représentant la rotation autour de l'axe X.
"""
R1(θ) = [1 0 0;
         0 cos(θ) -sin(θ);
         0 sin(θ)  cos(θ)]


"""
  R3(θ)

Matrice de rotation autour de l'axe Z d'un angle `θ` (en radians).

Utilisation :
- Rotation selon l'ascension du nœud ascendant (Ω).
- Rotation selon l'anomalie vraie (ν).
- Rotation horaire du repère dans les conversions ECI/ECEF.

# Argument
- θ : Angle de rotation en radians.

# Valeur retournée
- Matrice 3x3 représentant la rotation autour de l'axe Z.
"""
R3(θ) = [ cos(θ) -sin(θ) 0;
          sin(θ)  cos(θ) 0;
               0       0 1]


"""
  ecef_from_eci(r, t)

Convertit un vecteur de position d'un satellite du repère inertiel ECI (Earth-Centered Inertial) vers le repère tournant ECEF (Earth-Centered Earth-Fixed).  
La transformation applique une rotation autour de l'axe Z d'un angle égal à la rotation terrestre sur l'intervalle de temps `t`.

# Arguments
- r : Vecteur position du satellite dans le repère ECI.
- t : Temps écoulé (en secondes) depuis l'instant de référence.

# Valeur retournée
- Vecteur position transformé dans le repère ECEF.
"""
function ecef_from_eci(r,t)
	return R3(-ωe*t)*r # page 129 - SE216 (v.Aout 2024)
end

"""
  latlon_from_ecef(r)

Convertit un vecteur position exprimé dans le repère ECEF (Earth-Centered Earth-Fixed) en coordonnées géographiques latitude-longitude, exprimées en degrés.

# Arguments
- r : Vecteur position (x, y, z) en coordonnées ECEF.

# Valeurs retournées
- (latitude_deg, longitude_deg) : Latitude et longitude en degrés.
"""
function latlon_from_ecef(r)
    x, y, z = r # Vecteur position
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre
	
    ϕ = asin(z / ρ) 	# Latitude -> angle entre le vecteur position et le plan équatorial.  
    λ = atan(y, x) 		# Longitude -> angle dans le plan équatorial entre l'axe X (Greenwich) et la projection du vecteur sur ce plan.

    return rad2deg(ϕ), rad2deg(λ)
end

"""
  eci_pos(sat, t; mu=μ)

Calcule la position d'un satellite dans le repère inertiel ECI (Earth-Centered Inertial) à l'instant `t`, en supposant une orbite circulaire.

# Arguments
- sat : Structure décrivant le satellite (doit contenir a, i, Ω, M0).
- t   : Temps écoulé (en secondes) depuis l'instant de référence.
- mu  : Paramètre gravitationnel (par défaut `μ`).

# Valeur retournée
- Vecteur position du satellite dans le repère ECI.
"""
function eci_pos(sat,t;mu=μ)
    n = sqrt(mu / sat.a^3) 			    # Vitesse angulaire du satellite d'après la 3e loi de Kepler
    u = sat.M0 + n * t 				    # Position angulaire du satellite sur son orbite (Position initiale (M₀) + n * t)
    R = R3(sat.Ω) * R1(sat.i) * R3(u) 	# Matrice de transformation (Formule 7.35 - SE216 (v.Aout 2024))
    return R * [sat.a, 0, 0] 			# Position finale du satellite
end

"""
  walker_delta(P, S, F, i_deg, a)

Construit une constellation de type Walker-Delta, répartissant uniformément les satellites sur P plans orbitaux inclinés, avec S satellites par plan et un phasage défini par F.

# Arguments
- P      : Nombre total de plans orbitaux.
- S      : Nombre de satellites par plan orbital.
- F      : Facteur de déphasage entre les plans (Walker phasing).
- i_deg  : Inclinaison orbitale (en degrés).
- a      : Demi-grand axe de l'orbite (en mètres).

# Valeur retournée
- Un vecteur de structures `Sat`, chacune définie par ses paramètres
  orbitaux (a, i, Ω, M0).
"""
function walker_delta(P,S,F,i_deg,a)
    sats = Sat[]
    for p in 0:P-1, s in 0:S-1
      # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l'axe z.
      Ω = 2π * p / P
      # Répartition uniforme des satellites sur chaque orbite (s/S)
      # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.
      M0 = 2π * (s / S + (F * p) / (S * P))
      push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
    end
    return sats
end

"""
  myconstellation(vec, F, i_deg, a)

Construit une constellation de type Walker-Delta **généralisée**, où le nombre de satellites par plan orbital est défini par le vecteur `vec`.  
Cela permet une plus grande flexibilité que la version classique `walker_delta(P, S, F, i_deg, a)`.

# Arguments
- vec    : Vecteur de longueur P, où `vec[p]` représente le nombre de satellites dans le plan orbital p.
- F      : Paramètre de phasage (Walker phasing) appliqué entre les plans.
- i_deg  : Inclinaison orbitale en degrés.
- a      : Demi-grand axe de l'orbite (en mètres).

# Valeur retournée
- Un vecteur d'objets `Sat` correspondant aux satellites définis par (a, i, Ω, M0).
"""
function myconstellation(vec,F,i_deg,a)
	sats = Sat[]
  sizehint!(sats, sum(vec))
	P = length(vec)
	for p in 1:P
		S = vec[p]
		for s in 1:S
      # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l'axe z.
			Ω = 2π * p / P
			# Répartition uniforme des satellites sur chaque orbite (s/S)
      # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.
      M0 = 2π * (s / S + (F * p) / (S * P))
			push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
		end
	end
	return sats
end