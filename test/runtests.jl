using ToyClimaAtmos
using Test
using Aqua

@testset "ToyClimaAtmos.jl" begin

    # Overall quality test provided by Aqua
    #
    # Currently, there are problems with ambiguities in the base modules
    # https://github.com/JuliaTesting/Aqua.jl/issues/77
    Aqua.test_all(ToyClimaAtmos, ambiguities = false)
    Aqua.test_ambiguities(ToyClimaAtmos)


    # Test if everything is formatted
end
