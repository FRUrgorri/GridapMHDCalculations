using PackageCompiler
using Pkg

#This does not include GridapMHD since it has a __precompile(false)___ which is incompatible with creating the sysimage
pkgs = Symbol[]
append!(pkgs, [Symbol(v.name) for v in values(Pkg.dependencies()) if v.is_direct_dep]) 
deleteat!(pkgs,findall(isequal(:GridapMHD),pkgs))

create_sysimage(pkgs,sysimage_path="compile/GridapMHDCalculations.so",precompile_execution_file="compile/compile_workflow.jl", sysimage_build_args=`-O3 --check-bounds=no`)


