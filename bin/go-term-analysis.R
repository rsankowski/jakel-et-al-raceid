source('~/Documents/Single cell analysis/Advanced-plots/20190214-go_term-analysis.R')

date = Sys.Date()
load('data/sc.Robj')
source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")

df <- read_csv('data/unique-up-genes.csv')

enrich_up <- go_term_analysis()

dir.create('GO-terms')
dir.create('GO-terms/bp')
write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms.csv')

#mf terms
enrich_up_mf <- go_term_analysis(ontogeny = 'MF')

dir.create('GO-terms/mf')
write.csv(enrich_up_mf, 'GO-terms/mf/mf_GO_terms.csv')
