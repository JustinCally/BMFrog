#### Posterior predictive checks ####
# Check number of zeros (first)
# Check number of obs pred (site by group)
# prop zero
prop_zero<- function(x) mean(x == 0, na.rm = T)
prop_low <- function(x, bound = 20) mean(x < bound, na.rm = T)
mean_narm <- function(x) mean(x, na.rm = T)
sd_narm <- function(x) sd(x, na.rm = T)
max_narm <- function(x) max(x, na.rm = T)

#' Run posterior checks for CTDS distance sampling models
#'
#' @param model cmdstanr model
#' @param model_data model data generated from prepare_model_data()
#' @param stat stat passed to ppc_stat
#' @param title title for plot
#' @param ... additional arguments passed to ppc_stat
#'
#' @return plot grid
#' @export

posterior_checks_multispecies <- function(model, model_data, species_index,
                                          stat = "mean_narm", title, hits = NULL,
                                          integrated = F, only_det = F, ...) {

  # get n_obs from data
  n_obs_totals_0 <- model_data[[species_index]] %>%
    filter(Detected == 0) %>%
    pull(Prob_f)

  n_obs_totals_1 <- model_data[[species_index]] %>%
    filter(Detected == 1) %>%
    pull(Prob_f)

  which_inc_0 <- seq(1:length(n_obs_totals_0))
  which_inc_1 <- seq(1:length(n_obs_totals_1))

  model_draws <- model$draws("score_pred", format = "matrix")

   model_z <- 1-model$draws("det_z", format = "matrix")
   model_z[model_z == 0] <- NA
   model_draws_0 <- model_draws * model_z

   model_z <- model$draws("det_z", format = "matrix")
   model_z[model_z == 0] <- NA
   model_draws_1 <- model_draws * model_z

  which_sp <- which(stringr::str_detect(string = colnames(model_draws_0),
                                        pattern = paste0(species_index, "\\]")))

  model_draws_0_df <- data.frame(variable = "No detection",
                              value = model_draws_0[,which_sp] %>%
                                apply(MARGIN = 2, FUN = stat) %>%
                                unname())

  model_draws_1_df <- data.frame(variable = "Detection",
                                 value = model_draws_1[,which_sp] %>%
                                   apply(MARGIN = 2, FUN = stat) %>%
                                   unname())

  model_draws_df <- bind_rows(model_draws_0_df, model_draws_1_df)
  y_df <- data.frame(variable = c("No detection",
                                  "Detection"),
                     value = c(do.call(stat, args = list(x = n_obs_totals_0)),
                               do.call(stat, args = list(x = n_obs_totals_1))))

    ppc_plots <- ggplot(data = model_draws_df) +
      geom_histogram(aes(fill = variable, x = value, colour = variable),
                     linewidth = 0.15, na.rm = TRUE, binwidth = 0.005) +
      geom_vline(data = y_df, aes(xintercept = value, colour = variable), linewidth = 1) +
      bayesplot:::scale_color_ppc(values = c("#005954", "#894400"),
                                  labels = c(expression(italic(T(italic(`y=1`)))),
                                             expression(italic(T(italic(`y=0`)))))) +
      bayesplot:::scale_fill_ppc(values = c("#00B2A9", "#E57200"),
                      labels = c(expression(italic(T)(italic(`y=1`)[rep])),
                                 expression(italic(T)(italic(`y=0`)[rep])))) +
      guides(color = guide_legend(title = NULL),
             fill = guide_legend(order = 1, title = bayesplot:::stat_legend_title(stat,
                                                                      deparse(substitute(stat))))) +
      bayesplot:::dont_expand_y_axis() +
      bayesplot:::bayesplot_theme_get() +
      bayesplot:::no_legend_spacing() +
      bayesplot:::xaxis_title(FALSE) +
      bayesplot:::yaxis_text(FALSE) +
      bayesplot:::yaxis_ticks(FALSE) +
      bayesplot:::yaxis_title(FALSE) +
      ggplot2::ggtitle(label = title) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top")

  return(ppc_plots)
}
