setwd('/regeneron_effectiveness/code/R')

# Distance method for matching
# 0 = mahalanobis distance,
# 1 = Propensity score
psm = 0


rmarkdown::render(input = "matching_for_parametric_modeling.Rmd",
                  output_dir = "./Results",
                  output_format = "html_document",
                  params = list("psm" = psm),
                  output_file = "matching_for_parametric_modeling.html")


rmarkdown::render(input = "survival_parametric.Rmd",
                  output_dir = "./Results",
                  output_format = "html_document",
                  params = list("sub_grp_anlsys" = "full", "psm" = psm),
                  output_file = "full.html")

rmarkdown::render(input = "survival_parametric.Rmd",
                  output_dir = "./Results",
                  output_format = "html_document",
                  params = list("sub_grp_anlsys" = "below_60", "psm" = psm),
                  output_file = "subgroup_below_60.html")

rmarkdown::render(input = "survival_parametric.Rmd",
                  output_dir = "./Results",
                  output_format = "html_document",
                  params = list("sub_grp_anlsys" = "60_and_above", "psm" = psm),
                  output_file = "subgroup_60_and_above.html")








