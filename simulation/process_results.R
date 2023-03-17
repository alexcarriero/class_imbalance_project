# ------------------------ view results ----------------------------------------

calibration_plot <- function(dataframe){
  
  dataframe %>%
    filter(problematic != 1) %>% 
    drop_na() %>% 
    ggplot(aes(x = pred, y = class, group = as.factor(iter))) +
    geom_line(stat = "smooth",
              method = stats::loess, 
              formula = y ~ x,
              se = F, 
              linewidth = 0.5, 
              alpha = 0.1,
              color = "#0d0887") +
    geom_abline(slope = 1, intercept = 0, size = 1) +
    theme_minimal() +
    xlab(" ") + 
    ylab(" ") +
    xlim(0,1) + 
    ylim(0,1)  
}