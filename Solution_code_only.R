knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 12, fig.height = 8)
library(tidyverse)
library(babette)
library(ape)
library(beautier)
library(mauricer)
library(tracerer)

nexus_filename <- beastier::get_beast2_example_filename("Primates.nex")
fasta_filename <- "primates.fas"

nexus_data <- ape::read.nexus.data(nexus_filename)
output_data <- ape::as.DNAbin(nexus_data)

ape::write.FASTA(output_data, file = fasta_filename)
str(nexus_data)

bd_prior <- beautier::create_tree_prior_bd()

rln_clock <- beautier::create_clock_model_rln()

sub_model <- beautier::create_site_model_jc69()

mcmc_settings <- beautier::create_mcmc(chain_length = 3e6, store_every = 5000)

mrca_prior <- beautier::create_mrca_prior(
     taxa_names = beautier::get_taxa_names("primates.fas"),
     is_monophyletic = TRUE,
     mrca_distr = beautier::create_normal_distr(mean = 67.5, sigma = 0.75))

beauti_options <- beautier::create_beauti_options_v2_6()

inf_model <- beautier::create_inference_model(tree_prior = bd_prior,
                                              clock_model = rln_clock,
                                              site_model = sub_model,
                                              mcmc = mcmc_settings,
                                              mrca_prior = mrca_prior,
                                              beauti_options = beauti_options)

beast_result <- babette::bbt_run_from_model(fasta_filename = "primates.fas",
                                            inference_model = inf_model,
                                            beast2_options = 
                                beastier::create_beast2_options(rng_seed = 421))

beast_result$estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

estimates <- tracerer::remove_burn_ins(beast_result$estimates, 0.3)
estimates %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

tracerer::calc_esses(estimates, sample_interval = 5000)

all_trees <- beast_result$primates_trees
start <- floor(0.3 * length(all_trees))
end <- length(all_trees)
all_trees <- all_trees[start:end]

babette::plot_densitree(all_trees, alpha = 0.1)

cons_tree <- phytools::consensus.edges(trees = all_trees,
                                       method = "least.squares")
plot(cons_tree)

cons_tree$tip.label
human_index <- which(cons_tree$tip.label == "Homo_sapiens")

cons_tree$edge
tip_index <- which(cons_tree$edge[, 2] == human_index)
tip_index

cons_tree$edge.length[tip_index]

plot(cons_tree)
ape::edgelabels(round(cons_tree$edge.length, 2), col = "blue")

get_homo_age <- function(tree) {
  human_index <- which(tree$tip.label == "Homo_sapiens")
  tip_index <- which(tree$edge[, 2] == human_index)
  return(cons_tree$edge.length[tip_index])
}

all_ages <- lapply(all_trees, get_homo_age)
all_ages <- unlist(all_ages)
hist(all_ages)
mean(all_ages)
median(all_ages)
quantile(all_ages, c(0.025, 0.975))

sub_model <- beautier::create_site_model_gtr()

mcmc_settings <- beautier::create_mcmc(chain_length = 5e7, store_every = 5000)

inf_model <- beautier::create_inference_model(tree_prior = bd_prior,
                                              clock_model = rln_clock,
                                              site_model = sub_model,
                                              mcmc = mcmc_settings,
                                              mrca_prior = mrca_prior,
                                              beauti_options = beauti_options)

beast_result_gtr <- babette::bbt_run_from_model(fasta_filename = "primates.fas",
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

estimates_gtr <- tracerer::remove_burn_ins(beast_result_gtr$estimates, 0.1)
estimates_gtr %>%
  mutate("Sample" = Sample / 1000000) %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free") +
  xlab("Sampled trees (Millions)")

tracerer::calc_esses(estimates_gtr, sample_interval = 5000)

all_trees_gtr <- beast_result_gtr$primates_trees
start <- floor(0.1 * length(all_trees_gtr))
end <- length(all_trees_gtr)
all_trees_gtr <- all_trees_gtr[start:end]

babette::plot_densitree(all_trees_gtr, alpha = 0.1)

cons_tree_gtr <- phytools::consensus.edges(trees = all_trees_gtr,
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
