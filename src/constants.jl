"""Constants"""

const 𝕦 = UnitsConfig()  # 𝕦 is written as \bbu
const G = ustrip( uconvert(𝕦.l*𝕦.v^2/𝕦.m, u"G") )
const 𝕔 = CosmosConfig()
const 𝕤 = SolverConfig()
const 𝕗 = FrictionConfig()
const 𝕛 = JeansConfig()
const sis = SA[1,2,3]
const siss = sis .+ 3
const six = SA[1,2,3,4,5,6]

