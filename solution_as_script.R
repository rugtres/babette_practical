nexus_filename <- beastier::get_beast2_example_filename("Primates.nex")
fasta_filename <- "primates.fas"

nexus_data <- ape::read.nexus.data(nexus_filename)
output_data <- ape::as.DNAbin(nexus_data)

ape::write.FASTA(output_data, file = fasta_filename)
str(nexus_data)


bd_prior <- beautier::create_tree_prior_bd()

rln_clock <- beautier::create_clock_model_rln()

# this chain is quite long, so will take some time to run! Reduce length if you
# need quick and dirty results
mcmc_settings <- beautier::create_mcmc(chain_length = 1e7, store_every = 5000)

mrca_prior <- beautier::create_mrca_prior(
  taxa_names = beautier::get_taxa_names("primates.fas"),
  is_monophyletic = TRUE,
  mrca_distr = beautier::create_normal_distr(mean = 67.5, sigma = 0.75))

beauti_options <- beautier::create_beauti_options_v2_6()

inf_model <- beautier::create_inference_model(tree_prior = bd_prior,
                                              clock_model = rln_clock,
                                              mcmc = mcmc_settings,
                                              mrca_prior = mrca_prior,
                                              beauti_options = beauti_options)

beast_result <- babette::bbt_run_from_model(fasta_filename = "primates.fas",
                                            inference_model = inf_model,
                                            beast2_options = beastier::create_beast2_options(rng_seed = 42))

beast_result$estimates %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free")



estimates <- tracerer::remove_burn_ins(beast_result$estimates, 0.1)
estimates %>%
  gather(key = "statistic", val = "value", -c(Sample)) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_line() +
  facet_wrap(~statistic, scales = "free")

tracerer::calc_esses(estimates, sample_interval = 5000)



all_trees <- beast_result$primates_trees
start <- floor(0.1 * length(all_trees))
end <- length(all_trees)
all_trees <- all_trees[start:end]

babette::plot_densitree(all_trees, alpha = 0.1)

cons_tree <- phytools::consensus.edges(trees = all_trees)
plot(cons_tree)



## Age of Homo sapiens:
cons_tree$tip.label
human_index <- which(cons_tree$tip.label == "Homo_sapiens")
tip_index <- which(cons_tree$edge[, 2] == human_index)
cons_tree$edge.length[tip_index]


