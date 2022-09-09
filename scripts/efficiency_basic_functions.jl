using BenchmarkTools,SpecialFunctions

x = 225.3
@btime rand() # 5 ns
@btime randn() # 6 ns
@btime exp($x) # 6 ns
@btime log($x) # 7 ns
@btime cos($x) # 5.5 ns
@btime sin($x) # 5.2 ns
