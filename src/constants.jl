"""Constants"""

const 𝕦 = UnitsConfig()  # 𝕦 is written as \bbu
const G = ustrip( uconvert(𝕦.l*𝕦.v^2/𝕦.m, u"G") )
const 𝕔 = CosmosConfig()
const 𝕤 = SolverConfig()
