using Pkg
Pkg.activate(".")

# For the Notebooks
Pkg.add("Pluto")
Pkg.add("PlutoUI")

# Visualisations
Pkg.add(PackageSpec(name="Plots", version="1.35.0")) # Qt6 doesn't work with the .exe files
Pkg.add("FileIO")

# Others
Pkg.add("LinearAlgebra")

