
source('./src.r')

load(file="alpha_true.Rda")
load(file='sources.Rda')
load(file='sink.Rda')

result <- em(sink, sources, unk=0, iters=100, alpha_true=alpha_true)
