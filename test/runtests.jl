using TTVFaster
using Test
using DelimitedFiles
include("../examples/test_ttv.jl") 


function kepler62_test()
    data=readdlm("../examples/kepler62ef_planets.txt",',',Float64)
    @time ttv1,ttv2=test_ttv(5,40,20,data[1:10],false)
    inner=readdlm("../examples/inner_ttv.txt")
    outer=readdlm("../examples/outer_ttv.txt")
    inner_ref=readdlm("inner_ttv.txt.ref")
    outer_ref=readdlm("outer_ttv.txt.ref")
    diffs=[maximum(abs.(inner_ref .- inner)) maximum(abs.(outer_ref .- outer))]
    max_diff=maximum(diffs)
    return round(max_diff)
end

@testset "TTVFaster.jl" begin
    @test kepler62_test() == (0.0)
end
