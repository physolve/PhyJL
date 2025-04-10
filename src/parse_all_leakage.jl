# parse all adata
# R1 case

using DelimitedFiles
using CSV
using DataFrames
using Plots
using Glob

# original
runT = 5
file="2025-04-10_R1_$(runT).txt"
dir = "resultLeakage"
path = joinpath(@__DIR__, "data", dir, file);

# pack
files = glob("2025-04-*.txt","src/data/resultLeakage/")
files_model = glob("2025-04-*model.txt","src/data/resultLeakage/")
files_real = setdiff(files, files_model)
