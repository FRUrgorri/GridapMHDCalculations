using GridapMHDCalculations: kp_shercliff_cartesian
using GridapMHDCalculations.serialTest: Run_test_channel, Run_test_multigrid
using Test

#Test with default options. 
@time @testset begin 
    @test Run_test_channel(;b=1, Ha=10.0) ≈ kp_shercliff_cartesian(1.0,10.0) atol=0.01 

    kp,u_wall = Run_test_multigrid(;b=1, Ha=10.0)
    @test kp ≈ kp_shercliff_cartesian(1.0,10.0) atol=0.01
    for i in 1:3
        @test u_wall[i] ≈ 0.0 atol=0.001
    end
    
end  

