using GridapMHDCalculations.serialTest: Run_test_channel, Run_test_multigrid
using Test

#Test with default options. Future work: implement a test set
@time @test isnothing(Run_test_channel())     #TBD The test could check the pressure drop gradient according to Shercliff correlation
@time @test isnothing(Run_test_multigrid()) 

