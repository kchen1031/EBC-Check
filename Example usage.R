# ====================================================================
# Supported combinations (paper, dataset)  -- 23 keys / 9 sources
# --------------------------------------------------------------------
# Willwerscheid et al. 2024 - wOBA (ebnm)
#   willwerscheid_woba_normal       get_data("Willwerscheid", "wOBA_normal")
#   willwerscheid_woba_unimodal     get_data("Willwerscheid", "wOBA_unimodal")
#   willwerscheid_woba_npmle        get_data("Willwerscheid", "wOBA_npmle")
#
# Koenker et al. 2016 (REBayes)
#   koenker_norberg                 get_data("Koenker", "Norberg")
#   koenker_flies                   get_data("Koenker", "Flies")
#
# Narasimhan et al. 2020 (deconvolveR)
#   narasimhan_surgery              get_data("Narasimhan", "Surgery")
#   narasimhan_shakespeare          get_data("Narasimhan", "Shakespeare")
#
# Gu & Koenker - Bayesball (REBayes)
#   gu_baseball                     get_data("Gu", "baseball")
#
# Baek & Park 2022 - species richness (zero-truncated Poisson)
#   baek_brainvessel                get_data("Baek", "brainvessel")
#   baek_hardcandy                  get_data("Baek", "hardcandy")
#   baek_deathnotice                get_data("Baek", "deathnotice")
#   baek_traffic                    get_data("Baek", "traffic")
#   baek_microbial                  get_data("Baek", "microbial")
#   baek_shakespeare                get_data("Baek", "shakespeare")
#
# Jana et al. - hockey (Poisson minimum-distance)
#   jana_hockey_npmle               get_data("Jana", "hockey_npmle")
#   jana_hockey_hellinger           get_data("Jana", "hockey_hellinger")
#   jana_hockey_chisquare           get_data("Jana", "hockey_chisquare")
#
# You 2026 - Project STAR (Gaussian EB)
#   you_star_gaussian               get_data("You", "STAR_gaussian")
#   you_star_npmle                  get_data("You", "STAR_npmle")
#
# Wang et al. 2024 - neural-g
#   wang_arizona                    get_data("Wang", "arizona")
#   wang_framingham                 get_data("Wang", "framingham")
#
# Dicker & Zhao 2016 - Golub leukemia (REBayes::GLmix)
#   dickerzhao_golub_all            get_data("DickerZhao", "golub_all")
#   dickerzhao_golub_aml            get_data("DickerZhao", "golub_aml")
# ====================================================================

remotes::install_github("kchen1031/EBC-Check")

library(ebcalibration)

bundle <- get_data("wang", "arizona")


result <- eb_calibration_check(
  X_obs              = bundle$X_obs,
  prior_sampler      = bundle$prior_sampler,
  likelihood_sampler = bundle$likelihood_sampler,
  feature_transform  = bundle$feature_transform,
  B = 200, seed = 123
)
print(result)
plot(result)


