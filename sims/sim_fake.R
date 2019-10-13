
source('./src.r')

load(file="alpha_true.Rda")
load(file='sources.Rda')
load(file='sink.Rda')

result <- em(sink, sources, iters=1000)
