# GridapMHDCalculations

 Julia Project for launching simulations using GridapMHD app.
 
 Since GridapMHD is currently unregister the dependency needs to be added when instantiating the project. 
 Current is based on:

    add https://github.com/FRUrgorri/GridapMHD.jl.git#mirror-upstream

Note: The current "__precompile(false)__" of GridapMHD prevents the precompilation of GridapMHDCalculations (which returns a ?). 