library(tidyverse); library(broom); 

# Generate data ----
p <- length(beta);
x <- #matrix(replicate(p, rep(rnorm(n_subj),times = n_sim)), ncol = p) %*% 
  matrix(rnorm(n_subj * p * n_sim), ncol = p) %*% 
  chol(diag(1 - rho, p) + rho);

#Add enough noise to get desired true_rsq
true_sigma <- 
  sqrt(drop(t(beta)%*%(diag(1 - rho, p) + rho)%*%beta) *
         (1 / true_rsq - 1));
y <- drop(x %*% beta) + 
  true_sigma * rnorm(nrow(x));

#one row per simulation, per subject
wide_dat <- 
  data.frame(sim_id = rep(1:n_sim, each = n_subj), 
             subj_id = rep(1:n_subj, time = n_sim), 
             y = y,
             x = x);

#one row per simulation, per subject, per variable
tall_dat <- 
  wide_dat %>%
  pivot_longer(cols = starts_with("x"), values_to = "x") %>%
  arrange(sim_id, name, subj_id)


# No variable selection ----

#Fit full model with no variable selection
full_model_fmla <- as.formula(paste0("y ~ ", paste0("x.", 1:p, collapse = "+")));
full_model_fits <- 
  wide_dat %>%
  nest_by(sim_id) %>%
  mutate(model = list(lm(formula = full_model_fmla, data = data))) 

#Store coefficient estimates
full_model_coefs <- 
  full_model_fits %>%
  summarize(broom::tidy(model), 
            .groups = "drop") %>%
  #Remove rows for intercept term
  filter(term != "(Intercept)") %>%
  mutate(term = factor(term, levels = paste0("x.", 1:p), ordered = T))

#Store model summaries
full_model_summaries <- 
  full_model_fits %>%
  summarize(broom::glance(model), 
            .groups = "drop")

# p-value-based variable selection ----

univariable_model_fits <- 
  tall_dat %>%
  nest_by(sim_id, name) %>%
  mutate(model = list(lm(formula = y ~ x , data = data))) %>%
  summarize(broom::tidy(model), 
            .groups = "drop") %>%
  #Remove rows for intercept term
  filter(term != "(Intercept)") 

#Do model selection by setting. This will be identical to 'wide_dat'
#but with all 'non-significant' variables to zero
wide_dat_selected <- 
  univariable_model_fits %>%
  mutate(selected = p.value < p.value_threshold) %>%
  dplyr::select(sim_id, name, selected) %>%
  full_join(x = tall_dat) %>%
  mutate(x = x * selected) %>%
  dplyr::select(-selected) %>%
  pivot_wider(id_cols = sim_id:y, names_from = name, values_from = x) %>%
  dplyr::select(sim_id, subj_id, y, paste0("x.", 1:p))

#Now refit with revised data
selected_model_fits <- 
  wide_dat_selected %>%
  nest_by(sim_id) %>%
  mutate(model = list(lm(formula = full_model_fmla, data = data))) 

selected_model_coefs <- 
  selected_model_fits %>%
  summarize(broom::tidy(model), 
            .groups = "drop") %>%
  #Remove rows for intercept term, unselected terms
  filter(term != "(Intercept)") %>%
  mutate(p.value = replace_na(p.value, 1),
         term = factor(term, levels = paste0("x.", 1:p), ordered = T)) %>%
  arrange(sim_id, term)

selected_model_summaries <- 
  selected_model_fits %>%
  summarize(broom::glance(model), 
            .groups = "drop") %>%
  mutate(df = replace_na(df, 0))

#Calculate loglikelihoods on independent data given same covariates
y_new <- drop(x %*% beta) + true_sigma * rnorm(nrow(x))

yhat_full <- 
  pull(full_model_fits, model) %>%
  map(predict) %>%
  unlist();

full_model_summaries <-
  full_model_summaries %>%
  mutate(logLik_new = 
           dnorm(x = y_new - yhat_full, sd = rep(pull(full_model_summaries, sigma) * sqrt((n_subj - pull(full_model_summaries,df) - 1) / n_subj), each = n_subj), log = T) %>%
           #dnorm(x = y_new - yhat_full, sd = rep(pull(full_model_summaries, sigma), each = n_subj), log = T) %>%
           matrix(nrow = n_subj) %>%
           colSums());

yhat_selected <- 
  pull(selected_model_fits , model) %>% 
  map(predict) %>% 
  unlist();

selected_model_summaries <-
  selected_model_summaries %>%
  mutate(logLik_new = 
           dnorm(x = y_new - yhat_selected, sd = rep(pull(selected_model_summaries,sigma) * sqrt((n_subj - pull(selected_model_summaries,df) - 1) / n_subj), each = n_subj), log = T) %>%
           #dnorm(x = y_new - yhat_selected, sd = rep(pull(selected_model_summaries,sigma), each = n_subj), log = T) %>%
           matrix(nrow = n_subj) %>%
           colSums()) 

x_labels <- 
  #setdiff(levels(full_model_coefs$term), "x.1")
  levels(full_model_coefs$term)

# forward selection using AIC ----

# selected models
aic_model_fits <- 
  wide_dat %>%
  nest_by(sim_id) %>%
  mutate(model = list(stepAICc(object = lm(formula = y ~ 1, data = data), 
                               scope = list(upper = full_model_fmla),
                               direction = "forward", 
                               trace = FALSE))) 

#Store coefficient estimates
aic_model_coefs <- 
  aic_model_fits %>%
  summarize(broom::tidy(model), 
            .groups = "drop") %>%
  #Remove rows for intercept term
  filter(term != "(Intercept)") %>%
  mutate(term = factor(term, levels = paste0("x.", 1:p), ordered = T))

#Store model summaries
aic_model_summaries <- 
  aic_model_fits %>%
  summarize(broom::glance(model), 
            .groups = "drop")
