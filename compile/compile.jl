using PackageCompiler
using Pkg

pkgs = Symbol[]
append!(pkgs, [Symbol(v.name) for v in values(Pkg.dependencies()) if v.is_direct_dep]) #This does not include GridapMHD since it has a __precompile(false)___ which is incompatible with creating the sysimage

create_sysimage(pkgs,sysimage_path="GridapMHDCalculations.so",precompile_execution_file="test/serial/runtests_serial.jl")


