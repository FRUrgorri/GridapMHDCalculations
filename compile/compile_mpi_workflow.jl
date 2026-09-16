using GridapMHDCalculations.mpiTest: Run_test_channel, Run_test_multigrid

Run_test_channel((2,2,2);b=1, Ha=10.0)
Run_test_multigrid((2,2,2);b=1, Ha=10.0)