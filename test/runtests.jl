using TTVFaster
using Test
using DelimitedFiles
include("../examples/test_ttv.jl") 


function kepler62_test()
    data=readdlm("../examples/kepler62ef_planets.txt",',',Float64)
    @time ttv1,ttv2=test_ttv(5,40,20,data[1:10],false)
    inner=readdlm("../examples/inner_ttv.txt")
    outer=readdlm("../examples/outer_ttv.txt")
    diff1=ttv1[1].-inner[1,2]
    # diff2=ttv2[1].-outer[2,1]
    # println(diff1)
    return diff1
end

@testset "TTVFaster.jl" begin
    @test kepler62_test() == (0.0)
    # @test kepler62_test()[1,1] == kepler62_test()[2,1]
end
