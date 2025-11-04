# activate the environment
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# change to the root directory of the course
cd(joinpath(@__DIR__, ".."))

# start Pluto
using Pluto
Pluto.run()