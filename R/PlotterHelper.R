#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import zoo
#' @import xts
#' @import lubridate
#' @import ggplot2
#' @import stats
#' @import utils
NULL

# SPDX-License-Identifier: GPL-3.0-or-later
#
#  Copyright (C) 2024-2026 Domenic Franjic
#
#  This file is part of TwoStepSDFM.
#
#  TwoStepSDFM is free software: you can redistribute
#  it and/or modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  TwoStepSDFM is distributed in the hope that it
#  will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.

#' Helper for plotting factor time series
#' @keywords internal
plotFactorEstimates <- function(factors, smoothed_state_variance, no_of_factors,
                                axis_text_size) {
  pot_list_factors <- list()
  seq_along_dates <- 1:length(as.Date(time(factors)))
  for(factor in 1:no_of_factors){
    correction_factor <- 1.96 * sqrt(pmax(smoothed_state_variance[factor, seq(1, length(as.Date(time(factors))) * no_of_factors, by = no_of_factors) + (factor - 1)], 1e-15))
    current_factor <- data.frame(
      Date = as.Date(time(factors)),
      Value = factors[seq_along_dates, factor],
      `Upper 95%-CI` = factors[seq_along_dates, factor] + correction_factor,
      `Lower 95%-CI` = factors[seq_along_dates, factor] - correction_factor,
      check.names = FALSE
    )
    
    current_line_plot <- ggplot(current_factor, aes(x = Date, y = Value)) +
      geom_ribbon(aes(ymin = `Lower 95%-CI`, ymax = `Upper 95%-CI`), fill = "#88ccee", alpha = 0.4) +
      geom_line(colour = "black") +
      labs(title = paste0("Factor ", factor), y = "") +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
            text = element_text(size = axis_text_size))
    
    pot_list_factors[[factor]] <- current_line_plot
  }
  
  return( patchwork::wrap_plots(pot_list_factors, ncol = 1))
}

#' Helper for plotting loading matrix heat maps
#' @keywords internal
plotLoadingHeatMap <- function(loading_matrix_estim, series_names, no_of_factors, axis_text_size,
                                legend_title_text_size) {
  lambda_df <- as.data.frame(loading_matrix_estim)
  colnames(lambda_df) <- paste0("Factor ", 1:no_of_factors)
  if(dim(loading_matrix_estim)[2] == 1){
    stacked_loadings <- lambda_df
    stacked_loadings$Factor <- "Factor 1"
  }else{
    stacked_loadings <- stack(lambda_df[, ]) 
  }
  colnames(stacked_loadings) <- c("Loading", "Factor")
  stacked_loadings$Variable <- factor(rep(series_names, no_of_factors), levels = rev(series_names))
  heat_map_plot <- ggplot(stacked_loadings, aes(x = Factor, y = Variable)) +
    geom_tile(data = stacked_loadings, aes(fill = Loading), width = 0.9, height = 0.8) +
    geom_tile(data = subset(stacked_loadings, Loading == 0), fill = "black", width = 0.9, height = 0.8) +
    scale_fill_gradient2(low = "#88ccee", high = "#117733", na.value = "#882255", mid = "#FFFFFF") +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          text = element_text(size = axis_text_size),
          legend.title = element_text(size = legend_title_text_size),
          strip.text.y = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(heat_map_plot)
}

#' Helper for plotting measurement error var.cov. heatmap
#' @keywords internal
plotMeasVarCovHeatmap <- function(measurement_error_var_cov_df, series_names, axis_text_size,
                               legend_title_text_size) {
  colnames(measurement_error_var_cov_df) <- series_names
  stacked_measurement_error_var_cov <- stack(measurement_error_var_cov_df[, series_names])
  colnames(stacked_measurement_error_var_cov) <- c("(Co-)Variance", "Variable")
  stacked_measurement_error_var_cov$`(Co-)Variable` <- factor(rep(series_names, length(series_names)), levels = rev(series_names))
  heat_map_plot <- 
    ggplot(stacked_measurement_error_var_cov, aes(x = Variable, y = `(Co-)Variable`)) +
    geom_tile(data = stacked_measurement_error_var_cov, aes(fill = `(Co-)Variance`), width = 0.8, height = 0.8) +
    scale_fill_gradient2(low = "#88ccee", high = "#117733", na.value = "#882255", mid = "#FFFFFF") +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          text = element_text(size = axis_text_size),
          legend.title = element_text(size = legend_title_text_size),
          strip.text.y = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(heat_map_plot)
}

#' Helper for plotting measurement error var.cov. eigenvalues
#' @keywords internal
plotMeasVarCovEigenvalues <- function(eigen_values, no_of_factors, axis_text_size,
                                  legend_title_text_size) {
  eig_val_df <- data.frame("Value" = eigen_values,
                           "cut_off" = "normal",
                           "Eigen Value" = paste0("E.V. ", 1:length(eigen_values)),
                           check.names = FALSE)
  eig_val_df$`Eigen Value` <- factor(eig_val_df$`Eigen Value`, levels = eig_val_df$`Eigen Value`) 
  eig_val_df$cut_off[no_of_factors] <- "highlight"
  eig_val_plot <- ggplot(eig_val_df, aes(x = `Eigen Value`, y = Value, , fill = cut_off)) +
    geom_col() +
    geom_text(data = eig_val_df[no_of_factors, , drop = FALSE], 
              aes(x = `Eigen Value`, y = 0, label = "No. of factors chosen"),
              angle = 90, vjust = 0.5, hjust = 0, color = "white", size = axis_text_size * 0.33,
              inherit.aes = FALSE) + 
    scale_fill_manual(values = c(normal   = "#88ccee", highlight = "#117733")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          text = element_text(size = axis_text_size),
          legend.position = "none")
  
  return(eig_val_plot)
}




