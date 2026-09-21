using GridapMHDCalculations: kp_shercliff_cartesian
using GridapMHDCalculations.mpiTest: Run_test_channel, Run_test_multigrid
using Test

#Mpi (or sequential) test with default options. 
@testset begin 
    @time @test Run_test_channel((1,1,2);b=1.0, Ha=10.0, backend = :mpi) ≈ kp_shercliff_cartesian(1.0,10.0) atol=0.01 

    @time kp, u_wall = Run_test_multigrid((1,1,2);b=1.0, Ha=10.0b, ackend = :mpi)
    @test kp ≈ kp_shercliff_cartesian(1.0,10.0) atol=0.01
    for i in 1:3
        @test u_wall[i] ≈ 0.0 atol=0.01
    end
end  
