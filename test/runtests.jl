using TTVFaster
using Test
using DelimitedFiles
include("../examples/test_ttv.jl") 

function run_test()
    data=readdlm("../examples/kepler62ef_planets.txt",',',Float64)
    @time ttv1,ttv2=test_ttv(5,40,20,data[1:10]);
    return ttv1[1]
end

@testset "TTVFaster.jl" begin
    # @time ttv1,ttv2=test_ttv(5,40,20,data[1:10]);
    @test run_test() != 0
end
