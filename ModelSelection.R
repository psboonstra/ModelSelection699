## ----setup, include=FALSE, cache = FALSE--------------------------
library(MASS);
library(tidyverse); library(broom); library(knitr); library(glue);
library(modelr);
library(ggcorrplot)
knitr::opts_chunk$set(echo = T, warning = F, message = F);
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size);
})


def.chunk.hook  <- knitr::knit_hooks$get("chunk")
knitr::knit_hooks$set(chunk = function(x, options) {
  x <- def.chunk.hook(x, options)
  ifelse(options$size != "normalsize", paste0("\n \\", options$size,"\n\n", x, "\n\n \\normalsize"), x)
})


#knitr::opts_chunk$set(width = 10);
#knitr::opts_chunk$set(tidy.opts=list(width.cutoff=40));
recache = F;
options(digits = 3);
figure_scaler = 1/2;#1/2 for ioslides; ~1/3 for word, pdf
text_scaler = 3/3;#1 for ioslides; 2/3 for word, pdf
fig.x = 16 * figure_scaler;
fig.y = 9 * figure_scaler;


## ---- include = T, echo = T, cache = T, fig.width = fig.x, fig.height = fig.y, message = F----
library(tidyverse); library(broom);
set.seed(1);
n_sim <- 2e3;
#subjects / simulations
n_subj <- 100;
#compound symmetric correlation
rho <- 0.10;
#true generating betas
p_null <- 19;
#one non-zero beta;
beta <- c(1, numeric(p_null));
#true_rsq implies value of sigma
true_rsq <- 0.2;
#remove all estimates with univariable p-value exceeding 
p.value_threshold <- 0.10;
source("stepAICc.R"); # discussed later
source("varselect_sim.R");


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_model_coefs %>% 
  ggplot(mapping = aes(x = term, y = estimate)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), 
             size = 1/3, 
             alpha = 0.05,
             color = "red") + 
  geom_hline(yintercept = 0) + 
  scale_x_discrete(name = "X", limits = x_labels, labels = NULL) +
  scale_y_continuous(name = "") +
  coord_cartesian(ylim = c(-0.75, 1.75)) + 
  theme(text = element_text(size = text_scaler * 26),
        axis.text.y = element_text(size = text_scaler * 20));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
selected_model_coefs %>%
  filter(!is.na(statistic)) %>%
  ggplot(mapping = aes(x = term, y = estimate)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), 
             size = 1/3, 
             alpha = 0.05,
             color = "red") + 
  geom_hline(yintercept = 0) + 
  geom_label(data = selected_model_coefs %>% filter(!is.na(statistic)) %>% count(term), 
             aes(x = term, label = n), 
             y = 1.75) + 
  scale_x_discrete(name = "X", limits = x_labels, labels = NULL) +
  scale_y_continuous(name = "") +
  coord_cartesian(ylim = c(-0.75, 1.75)) + 
  theme(text = element_text(size = text_scaler * 26),
        axis.text.y = element_text(size = text_scaler * 20));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_model_summaries %>% 
  ggplot() + 
  geom_histogram(aes(x = r.squared), bins = 30) +
  geom_vline(xintercept = true_rsq, size = 1.5) + 
  scale_x_continuous(name = expression(R^2),
                     limits = c(0, 1), 
                     expand = c(0,0)) + 
  scale_y_continuous(name = "", labels = NULL) +
  theme(text = element_text(size = text_scaler * 26));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
selected_model_summaries %>%
  ggplot() + 
  geom_histogram(aes(x = r.squared), bins = 30) +
  geom_vline(xintercept = true_rsq, size = 1.5) + 
  scale_x_continuous(name = expression(R^2),
                     limits = c(0, 1), 
                     expand = c(0,0)) + 
  scale_y_continuous(name = "", labels = NULL) +
  theme(text = element_text(size = text_scaler * 26));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
selected_model_coefs %>% 
  filter(!is.na(estimate)) %>%
  group_by(sim_id) %>%
  count() %>%
  ggplot() + 
  geom_bar(aes(x = I(n+2))) + 
  scale_x_continuous(name = "", 
                     breaks = seq(1, 20, by = 2)) +
  scale_y_continuous(name = "Percent of sims", 
                     breaks = seq(0, n_sim, length = 6), 
                     labels = paste0(100 * seq(0, n_sim, length = 6) / n_sim,"%")) +
  theme(text = element_text(size = text_scaler * 26));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_model_summaries %>%
  mutate(aic_compare = logLik_new - logLik + df + 2) %>%
  ggplot() + 
  geom_histogram(aes(x = aic_compare), 
                 bins = 30, 
                 fill = "#d8b365") +
  geom_point(aes(x = mean(aic_compare), 
                 y = -5),
             size = 1.5,
             color = "#d8b365") + 
  scale_x_continuous(name = "" ) + 
  scale_y_continuous(name = "", labels = NULL) +
  theme(text = element_text(size = text_scaler * 26));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----

selected_model_summaries %>%
  mutate(aic_compare = logLik_new - logLik + df + 2) %>%
  ggplot() + 
  geom_histogram(aes(x = aic_compare), 
                 bins = 30, 
                 fill = "#d8b365") +
  geom_point(aes(x = mean(aic_compare),
                 y = -5),
             size = 1.5,
             color = "#d8b365") + 
  scale_x_continuous(name = "" ) + 
  scale_y_continuous(name = "", labels = NULL) +
  theme(text = element_text(size = text_scaler * 26));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_join(full_model_summaries %>% 
            dplyr::select(sim_id, logLik, AIC) %>% 
            mutate(AIC = AIC/(-2)) %>%
            pivot_longer(cols = logLik:AIC, 
                         names_to = "name", 
                         values_to = "full_measure"),
          selected_model_summaries %>%
            dplyr::select(sim_id, logLik, AIC) %>% 
            mutate(AIC = AIC/(-2)) %>%
            pivot_longer(cols = logLik:AIC, 
                         names_to = "name", 
                         values_to = "selected_measure")) %>%
  ggplot() +
  geom_point(aes(x = full_measure, 
                 y = selected_measure, 
                 col = name, 
                 shape = name),
             size = 1.75,
             alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(name = "No variable selection") + 
  scale_y_continuous(name = "Variable selection") +
  scale_color_brewer(name = "Measure", 
                     palette = "Dark2",
                     labels = c("LogLik - Model Size","LogLik")) +
  scale_shape_discrete(name = "Measure", 
                       labels = c("LogLik - Model Size","LogLik")) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) + 
  theme(text = element_text(size = text_scaler * 26),
        legend.position = "top");


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_model_summaries %>%
  mutate(aic_compare =
           logLik_new - logLik + df + 2, 
         aicc_compare =
           logLik_new - logLik + (df + 2) * n_subj / (n_subj - df - 3)) %>%
  ggplot() + 
  geom_histogram(aes(x = aic_compare, 
                     fill = "aic_compare"), 
                 bins = 30) +
  geom_point(aes(x = mean(aic_compare), 
                 y = -5),
             color = "#d8b365", 
             size = 1.5) +
  geom_histogram(aes(x = aicc_compare,
                     fill = "aicc_compare"), 
                 bins = 30) +
  geom_point(aes(x = mean(aicc_compare), 
                 y = -5),
             color = "#01665e", 
             size = 1.5) + 
  scale_x_continuous(name = "") + 
  scale_y_continuous(name = "", labels = NULL) +
  scale_fill_manual(name = "", 
                    breaks = c("aic_compare", "aicc_compare"),
                    values = c("#d8b365", "#01665e"), 
                    labels = c("LogLikNew - LogLik + Model Size", 
                               "LogLikNew - LogLik + Adj. Model Size")) + 
  guides(fill = guide_legend(nrow = 2)) +
  theme(text = element_text(size = text_scaler * 26),
        legend.position = "top");


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
selected_model_summaries %>%
  mutate(aic_compare =
           logLik_new - logLik + df + 2, 
         aicc_compare =
           logLik_new - logLik + (df + 2) * n_subj / (n_subj - df - 3)) %>%
  ggplot() + 
  geom_histogram(aes(x = aic_compare, 
                     fill = "aic_compare"), 
                 bins = 30) +
  geom_point(aes(x = mean(aic_compare), 
                 y = -5),
             color = "#d8b365", 
             size = 1.5) +
  geom_histogram(aes(x = aicc_compare,
                     fill = "aicc_compare"), 
                 bins = 30) +
  geom_point(aes(x = mean(aicc_compare), 
                 y = -5),
             color = "#01665e", 
             size = 1.5) + 
  scale_x_continuous(name = "" ) + 
  scale_y_continuous(name = "", labels = NULL) +
  scale_fill_manual(name = "", 
                    breaks = c("aic_compare", "aicc_compare"),
                    values = c("#d8b365", "#01665e"), 
                    labels = c("LogLikNew - LogLik + Model Size", 
                               "LogLikNew - LogLik + Adj. Model Size")) + 
  guides(fill = guide_legend(nrow = 2)) + 
  theme(text = element_text(size = text_scaler * 26), 
        legend.position = "top");


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_join(full_model_summaries %>% 
            dplyr::select(sim_id, logLik, AIC, df) %>% 
            mutate(AIC = -0.5 * AIC,
                   AICc = logLik - df * (n_subj/(n_subj - df - 3))) %>%
            dplyr::select(sim_id, AIC, AICc) %>%
            pivot_longer(cols = AIC:AICc, names_to = "name", values_to = "full_measure"),
          selected_model_summaries %>%
            dplyr::select(sim_id, logLik, AIC, df) %>% 
            mutate(AIC = -0.5 * AIC,
                   AICc = logLik - df * (n_subj/(n_subj - df - 3))) %>%
            dplyr::select(sim_id, AIC, AICc) %>%
            pivot_longer(cols = AIC:AICc, names_to = "name", values_to = "selected_measure")) %>%
  ggplot() +
  geom_point(aes(x = full_measure, 
                 y = selected_measure, 
                 col = name, 
                 shape = name),
             size = 1.75,
             alpha = 0.6) + 
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_continuous(name = "No variable selection") + 
  scale_y_continuous(name = "Variable selection") +
  scale_color_manual(name = "", 
                     values = c("#1B9E77","#7570B3"),
                     labels = c("LogLik - Model Size","LogLik - Adj. Model Size")) +
  scale_shape_discrete(name = "", 
                       labels = c("LogLik - Model Size","LogLik - Adj. Model Size")) +
  guides(color = guide_legend(nrow = 2,override.aes = list(size = 3.5, alpha = 1))) + 
  theme(text = element_text(size = text_scaler * 26),
        legend.position = "top");


## ----echo=FALSE, out.width='30%'----------------------------------
knitr::include_graphics('anna.png')


## ---- include = T, echo = T, cache = !recache, fig.width = fig.x, fig.height = fig.y----
breast_dx <-
  read_csv("bdiag.csv", show_col_types = FALSE) %>%
  # Translate M/D into 1/0
  mutate(malignant = 1 * (diagnosis == "M")) %>% 
  # Drop errant space in 'concave points_mean' variable name 
  rename_with(~str_replace(string = ., pattern = " ", replacement = "")) %>%
  # Focus only on worst measurements
  select(malignant, 
         #contains("_mean"), 
         #contains("_se"), 
         contains("_worst")) 



## ---- include = T, echo = T, cache = !recache, fig.width = 1.3*fig.x, fig.height = 1.3*fig.y----

ggcorrplot(cor(select(breast_dx, -malignant)), hc.order = T)



## ----bootstrap1, include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
full_fmla <- 
  glue("~",
       glue_collapse(glue("{setdiff(colnames(breast_dx),'malignant')}"),sep = "+")) %>%
  as.formula();
n_boot = 2e3

breast_dx_bootstrap <- 
  breast_dx %>% 
  bootstrap(n = n_boot) 

forward_aic <- 
  breast_dx_bootstrap %>%
  map(.x = .$strap, 
      .f = ~ stepAIC(glm(malignant ~ 1, 
                         data = .,
                         family = "binomial"),
                     scope = list(upper = full_fmla),
                     direction = "forward", 
                     trace = F)) %>%
  map_dfr(tidy, .id = "boot_id") %>%
  filter(term != "(Intercept)") %>% 
  select(boot_id, term) %>%
  arrange(term, boot_id)

forward_aic_counts <- 
  forward_aic %>% 
  pivot_wider(id_cols = boot_id, 
              names_from = term, 
              values_from = term, 
              values_fill = "") %>% 
  select(-boot_id) %>% 
  group_by_all() %>% 
  count(name = "n_aic")

forward_aicc <- 
  breast_dx_bootstrap %>%
  map(.x = .$strap, 
      .f = ~ stepAICc(glm(malignant ~ 1, 
                          data = .,
                          family = "binomial"),
                      scope = list(upper = full_fmla),
                      direction = "forward", 
                      trace = F)) %>%
  map_dfr(tidy, .id = "boot_id") %>%
  filter(term != "(Intercept)") %>% 
  select(boot_id, term) %>%
  arrange(term, boot_id)


forward_aicc_counts <- 
  forward_aicc %>% 
  pivot_wider(id_cols = boot_id, 
              names_from = term, 
              values_from = term, 
              values_fill = "") %>% 
  select(-boot_id) %>% 
  group_by_all() %>% 
  count(name = "n_aicc")



## ----bootstrap2, include = T, echo = F, cache = !recache, fig.width = 0.8 * fig.x, fig.height = 0.8 * fig.y----

bind_rows(
  #forward_aic %>% mutate(criterion = "AIC"),
  forward_aicc %>% mutate(criterion = "AICc")) %>%
  mutate(term = factor(term) %>% fct_infreq()) %>% 
  ggplot() + 
  geom_bar(aes(x = term,
               #fill = criterion
               ),
           position = "dodge") +
  scale_x_discrete(name = NULL) +
  #scale_fill_manual(breaks = c("AIC", "AICc"),
  #                  values = c("#d8b365", "#01665e")) + 
  scale_y_reverse(name = "Percent selected", 
                  breaks = seq(0, n_boot, length = 6), 
                  labels = paste0(100 * seq(0, n_boot, length = 6) / n_boot,"%"), 
                  expand = expansion(0.02)) +
  coord_flip() + 
  theme(text = element_text(size = text_scaler * 19),
        legend.position = "top",
        axis.text.y = element_text(size = text_scaler * 13))



## ----bootstrap3, include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
# Actual selected models
aic_selected <- 
  stepAIC(glm(malignant ~ 1, 
              data = breast_dx,
              family = "binomial"),
          scope = list(upper = full_fmla),
          direction = "forward", 
          trace = F)$coef[-1] %>% names() %>% sort() %>% paste0(collapse = ", ")
aicc_selected <- 
  stepAICc(glm(malignant ~ 1, 
               data = breast_dx,
               family = "binomial"),
           scope = list(upper = full_fmla),
           direction = "forward", 
           trace = F)$coef[-1] %>% names() %>% sort() %>% paste0(collapse = ", ")



most_common_models <- 
  full_join(forward_aic_counts, forward_aicc_counts) %>%
  unite(col = "selected", 
        -c(n_aic, n_aicc),
        sep = ",") %>%
  mutate(selected =
           str_replace_all(selected, ",+", ", ") %>%
           str_replace("^, ","") %>%
           str_replace(",$","")) %>%
  full_join(
    tibble(selected = aic_selected, 
           aic_selected_obs = T)) %>%
  full_join(
    tibble(selected = aicc_selected, 
           aicc_selected_obs = T)) %>%
  mutate(aic_selected_obs = replace_na(aic_selected_obs, FALSE), 
         aicc_selected_obs = replace_na(aicc_selected_obs, FALSE), 
         model = 
           glue("malignant ~ {str_replace_all(selected, ',','+')}") %>%
           as.character()) %>%
  arrange(-(n_aicc*n_aic));

observed_aics <- 
  most_common_models %>%
  nest_by(model) %>%
  mutate(fit = list(glm(formula = model, data = breast_dx, family = "binomial"))) %>%
  mutate(obs_aic = extractAIC(fit)[[2]],
         obs_aicc = extractAICc.glm(fit)[[2]]) %>%
  ungroup() %>%
  select(-fit, -data);


all_results <- 
  full_join(most_common_models, observed_aics) %>%
  mutate(n_aic = replace_na(n_aic, 0), 
         n_aicc = replace_na(n_aicc, 0), 
         obs_aic = obs_aic - obs_aic[which(aic_selected_obs)], 
         obs_aicc = obs_aicc - obs_aicc[which(aicc_selected_obs)])

competing_models <-
  all_results %>% 
  filter(n_aic == max(n_aic) | aic_selected_obs | obs_aic == min(obs_aic))

all_results  %>%
  filter(obs_aic < 0.5 | obs_aicc < 0.5) %>%
  mutate(selected = str_replace_all(selected, "_worst","")) %>%
  select(selected, #n_aic, obs_aic, 
         n_aicc, obs_aicc) %>%
  arrange(obs_aicc) %>%
  mutate(#obs_aic = formatC(obs_aic, format = "f", digits = 2),
         obs_aicc = formatC(obs_aicc, format = "f", digits = 2)) %>%
  add_row(selected = "Subtotal", 
          #n_aic = sum(.$n_aic), 
          #obs_aic = "",
          n_aicc = sum(.$n_aicc), 
          obs_aicc = "") %>%
  mutate(#n_aic = paste0(formatC(100 * n_aic / n_boot, format = "f", digits = 1),"%"), 
         n_aicc = paste0(formatC(100 * n_aicc / n_boot, format = "f", digits = 1),"%")) %>%
  rename(`Variable set` = selected, 
         #`Pct. Sel. AIC` = n_aic, 
         #`Delta(Obs. AIC)` = obs_aic,
         `Pct. Sel. AICc` = n_aicc,
         `Delta(Obs. AICc)` = obs_aicc) %>%
  kable();


## ----bootstrap4, include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
competing_models %>%
  mutate(selected = str_replace_all(selected, "_worst","")) %>%
  select(selected, #n_aic, obs_aic, 
         n_aicc, obs_aicc) %>%
  arrange(obs_aicc) %>%
  mutate(#obs_aic = formatC(obs_aic, format = "f", digits = 2),
         obs_aicc = formatC(obs_aicc, format = "f", digits = 2))  %>%
  mutate(#n_aic = paste0(formatC(100 * n_aic / n_boot, format = "f", digits = 1),"%"), 
         n_aicc = paste0(formatC(100 * n_aicc / n_boot, format = "f", digits = 1),"%")) %>%
  rename(`Variable set` = selected, 
         #`Pct. Sel. AIC` = n_aic, 
         #`Delta(Obs. AIC)` = obs_aic,
         `Pct. Sel. AICc` = n_aicc,
         `Delta(Obs. AICc)` = obs_aicc) %>%
  kable();


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----

bind_rows(
  # Most commonly selected by bootstrap
  competing_models %>% 
    filter(n_aicc == max(n_aicc)) %>%
    pull(model) %>%
    glm(family = "binomial",
        data = breast_dx) %>%
    tidy() %>%
    mutate(Mechanism = "Most common (bootstrap)"),
  # Selected by forward selection on original data
  competing_models %>% 
    filter(aicc_selected_obs) %>%
    pull(model) %>%
    glm(family = "binomial",
        data = breast_dx) %>%
    tidy() %>%
    mutate(Mechanism = "Forward selected"),
  # Selected by forward selection on original data
  competing_models %>% 
    filter(obs_aicc == min(obs_aicc)) %>%
    pull(model) %>%
    glm(family = "binomial",
        data = breast_dx) %>%
    tidy() %>%
    mutate(Mechanism = "Smallest AICc (bootstrap)")) %>% 
  arrange(term, Mechanism) %>%
  filter(term != "(Intercept)") %>%
  left_join(
    apply(breast_dx[-1],2,sd) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "term") %>%
      as_tibble() %>%
      rename(stdev = ".")) %>%
  mutate(std_estimate = paste0("(",formatC(estimate*stdev, format = "f", digits = 1),")"),
         estimate = formatC(estimate, format = "f", digits = 2)) %>%
  unite("pretty_value", c(estimate,std_estimate), sep ="") %>%
  pivot_wider(id_cols = term,
              names_from = Mechanism,
              values_from = pretty_value, 
              values_fill = "") %>%
  rename(` ` = term) %>%
  kable(align = c("lrrr"))



## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
aic_model_coefs %>%
  ggplot(mapping = aes(x = term, y = estimate)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.2), 
             size = 1/3, 
             alpha = 0.05,
             color = "red") + 
  geom_hline(yintercept = 0) + 
  geom_label(data = aic_model_coefs %>% filter(!is.na(statistic)) %>% count(term), 
             aes(x = term, label = n), 
             y = 1.75) + 
  scale_x_discrete(name = "X", limits = x_labels, labels = NULL) +
  scale_y_continuous(name = "") +
  coord_cartesian(ylim = c(-0.75, 1.75)) + 
  theme(text = element_text(size = text_scaler * 26),
        axis.text.y = element_text(size = text_scaler * 20));


## ---- include = T, echo = F, cache = !recache, fig.width = fig.x, fig.height = fig.y----
aic_model_summaries %>%
  ggplot() + 
  geom_histogram(aes(x = r.squared), bins = 30) +
  geom_vline(xintercept = true_rsq, size = 1.5) + 
  scale_x_continuous(name = expression(R^2),
                     limits = c(0, 1), 
                     expand = c(0,0)) + 
  scale_y_continuous(name = "", labels = NULL) +
  theme(text = element_text(size = text_scaler * 26));

