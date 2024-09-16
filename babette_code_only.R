knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 12, fig.height = 8)
library(tidyverse)
library(babette)
library(ape)
library(beautier)
library(mauricer)
library(tracerer)

## install.packages("ape")
## install.packages("beastier")
## install.packages("tracerer")
## install.packages("mauricer")
## install.packages("beautier")
## 
## # this installs beast2:
## remotes::install_github("richelbilderbeek/beastierinstall")
## beastierinstall::install_beast2()
## 
## install.packages("babette")
## 
## # and for plotting:
## install.packages("tidyverse")

## install.packages("rJava")

## library(tidyverse)
## library(babette)
## library(ape)
## library(beautier)
## library(mauricer)
## library(tracerer)

nexus_filename <- beastier::get_beast2_example_filename("dna.nex")

nexus_data <- ape::read.nexus.data(nexus_filename)
output_data <- ape::as.DNAbin(nexus_data)

ape::write.FASTA(output_data, file = "dna.fas")

str(nexus_data)

bd_prior <- beautier::create_tree_prior_bd()

rln_clock <- beautier::create_clock_model_rln()

sub_model <- beautier::create_site_model_jc69()

mcmc_settings <- beautier::create_mcmc(chain_length = 3e6, store_every = 1000)

inf_model <- beautier::create_inference_model(tree_prior = bd_prior,
                                              site_model = sub_model,
                                              clock_model = rln_clock,
                                              mcmc = mcmc_settings)

beast_result <- babette::bbt_run_from_model(fasta_filename = "dna.fas",
                                            inference_model = inf_model,
                                            beast2_options = 
                                 beastier::create_beast2_options(rng_seed = 42))

beast_result$estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

beast_result$estimates <- tracerer::remove_burn_ins(beast_result$estimates, 0.1)
beast_result$estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

tracerer::calc_esses(beast_result$estimates, sample_interval = 1000)

all_trees <- beast_result$dna_trees
start <- floor(0.1 * length(all_trees))
end <- length(all_trees)
all_trees <- all_trees[start:end]

babette::plot_densitree(all_trees, alpha = 0.1)

cons_tree <- phytools::consensus.edges(trees = all_trees,
                                       method = "least.squares")
plot(cons_tree)

sub_model <- beautier::create_site_model_gtr()

inf_model <- beautier::create_inference_model(tree_prior = bd_prior,
                                              site_model = sub_model,
                                              clock_model = rln_clock,
                                              mcmc = mcmc_settings)

beast_result_gtr <- babette::bbt_run_from_model(fasta_filename = "dna.fas",
                                            inference_model = inf_model,
                                            beast2_options = 
                              beastier::create_beast2_options(rng_seed = 4321))

beast_result_gtr$estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

gtr_estimates <- tracerer::remove_burn_ins(beast_result_gtr$estimates, 0.1)
gtr_estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

gtr_estimates %>%
  gather(key = "statistic", val = "value", 
         c(rateAC, rateAG, rateAT, rateCG, rateGT)) %>%
  group_by(statistic) %>%
  summarise("median" = median(value))

gtr_trees <- beast_result_gtr$dna_trees
start <- floor(0.1 * length(gtr_trees))
end <- length(gtr_trees)
gtr_trees <- gtr_trees[start:end]
babette::plot_densitree(gtr_trees, alpha = 0.1)

cons_tree_gtr <- phytools::consensus.edges(trees = gtr_trees,
                                       method = "least.squares")

par(mfrow = c(1, 2))
plot(cons_tree, main = "JC69", type = "cladogram")
add.scale.bar()
plot(cons_tree_gtr, main = "GTR", type = "cladogram")
add.scale.bar()

ape::ltt.plot(cons_tree,
              col = "blue", lwd = 3)
ape::ltt.lines(cons_tree_gtr,
               col = "red", lwd = 1)
legend("topleft", legend = c("JC", "GTR"),
       col = c("blue", "red"), lwd = c(3, 1), lty = c(1, 1))

mrca_prior <- beautier::create_mrca_prior(
     taxa_names = beautier::get_taxa_names("dna.fas"),
     is_monophyletic = TRUE,
     mrca_distr = beautier::create_normal_distr(mean = 319, sigma = 1.64))

mcmc_settings <- beautier::create_mcmc(chain_length = 3e6, store_every = 1000)


beauti_options <- beautier::create_beauti_options_v2_6()

inf_model <- beautier::create_inference_model(tree_prior = bd_prior,
                                              clock_model = rln_clock,
                                              site_model = sub_model,
                                              mcmc = mcmc_settings,
                                              mrca_prior = mrca_prior,
                                              beauti_options = beauti_options)

beast_result <- babette::bbt_run_from_model(fasta_filename = "dna.fas",
                                            inference_model = inf_model,
                                            beast2_options = 
                                beastier::create_beast2_options(rng_seed = 42))

beast_result$estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

estimates <- tracerer::remove_burn_ins(beast_result$estimates, 0.2)
estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free")
tracerer::calc_esses(estimates, sample_interval = 1000)

all_trees <- beast_result$dna_trees
start <- floor(0.2 * length(all_trees))
end <- length(all_trees)
all_trees <- all_trees[start:end]

babette::plot_densitree(all_trees, alpha = 0.1)

nexus_filename <- beastier::get_beast2_example_filename("Primates.nex")
fasta_filename <- "primates.fas"

nexus_data <- ape::read.nexus.data(nexus_filename)
output_data <- ape::as.DNAbin(nexus_data)

ape::write.FASTA(output_data, file = fasta_filename)
str(nexus_data)
