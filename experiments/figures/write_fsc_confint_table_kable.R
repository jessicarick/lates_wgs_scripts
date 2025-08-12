library(kableExtra)

confint.table <- bs_CI %>%
  filter(! param %in% c("MaxObsLhood","MaxEstLhood")) %>%
  mutate(param_type = c(rep("Ne change",3),
                        #rep("Lhood",2),
                        rep("Migration rate",12),
                        rep("Starting Ne",6),
                        rep("Ne change",2),
                        rep("Admixture event timing",6),
                        rep("Divergence event timing",4))) %>%
  filter(param_type %in% c("Migration rate","Admixture event timing","Starting Ne","Divergence event timing")) %>%
  relocate(param_type) %>%
  select(param:median,q0.025:q0.975) %>%

  kbl(booktabs = T,caption="Caption here",
      format="latex",
      col.names=c("Parameter","Mean","Median","2.5% Quantile","97.5% Quantile")
      ) %>%
  #add_header_above(c(" ","All Trees" = 4, "Low ILS" = 4, "Medium ILS" = 4,  "High ILS" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Migration Rates", 1, 12) %>%
  pack_rows("Effective Population Sizes", 13, 18) %>%
  pack_rows("Admixture Events (generations ago)", 19,24) %>%
  pack_rows("Divergence Events (generations ago)", 25, 28) 
writeLines(confint.table, '~/Dropbox/Apps/Overleaf/LatesAdmixture_Preprint_PNAS_v1/figures/fbranch_model_confint_table.tex')
