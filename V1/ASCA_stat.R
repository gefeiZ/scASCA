library(glmmTMB)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(dplyr)
library(tidyr)
library(dplyr)
library(tidyr)
library(ggplot2)


df <- read.table("asca_counts_input_cd4.csv", header=TRUE, sep="\t")


df <- df %>%
  filter(!grepl("chrY|chrM", SNP)) %>%
  mutate(
    feature = paste(SNP, Peak, sep=":"),
    Total_Count = Count_Pos + Count_Neg,
    AF = Count_Pos / (Total_Count + 0.00001), 
    AF_dev = AF - 0.5 
  ) %>%
  filter(Total_Count >= 3) # 简单的深度过滤


feature_stats <- df %>%
  group_by(feature) %>%
  summarise(
    n_valid_donors = n_distinct(Donor),
    .groups = "drop"
  ) %>%
  filter(n_valid_donors >= 3) 

df_testing <- df %>% semi_join(feature_stats, by="feature")

run_test_robust_tracking <- function(sub_df) {

  if(nrow(sub_df) < 3) return(list(p_val = NA_real_, method = "TooFewDonors"))
  
  pos <- sub_df$Count_Pos
  neg <- sub_df$Count_Neg
  sum_pos <- sum(pos)
  sum_total <- sum(pos + neg)
  
  # Perfect Separation

  if (sum_pos == 0 || sum_pos == sum_total) {
    res <- binom.test(sum_pos, sum_total, p = 0.5)
    return(list(p_val = res$p.value, method = "ExactBinomial_Perfect"))
  }
  
  # MAIN：Beta-Binomial ( Donor Aware)

  tryCatch({
    fit <- glmmTMB(cbind(pos, neg) ~ 1, family = betabinomial(link = "logit"))
    return(list(p_val = summary(fit)$coefficients$cond[1, "Pr(>|z|)"], 
                method = "BetaBinomial"))
    
  }, error = function(e) {
    
    # Quasi-Binomial (outline backup)
    tryCatch({
      fit_q <- glm(cbind(pos, neg) ~ 1, family = quasibinomial(link = "logit"))
      

      if (!fit_q$converged) {
        res <- binom.test(sum_pos, sum_total, p = 0.5)
        return(list(p_val = res$p.value, method = "ExactBinomial_NonConverged"))
      }
      coef_val <- coef(fit_q)[1]
      if(abs(coef_val) > 10) { 
        res <- binom.test(sum_pos, sum_total, p = 0.5)
        return(list(p_val = res$p.value, method = "ExactBinomial_HighCoef"))
      }
      return(list(p_val = summary(fit_q)$coefficients[1, "Pr(>|t|)"], 
                  method = "QuasiBinomial"))
      
    }, error = function(e2) {
      return(list(p_val = NA_real_, method = "Failed"))
    })
  })
}

# ==========================================
# RUN PIPELINE
# ==========================================
results_final <- df_testing %>%
  group_by(feature) %>%
  summarise(
    SNP = first(SNP),
    Peak = first(Peak),
    n_donors = n(),
    mean_AF = sum(Count_Pos) / sum(Total_Count),
    total_depth = sum(Total_Count),
    test_result = list(run_test_robust_tracking(cur_data())),
    .groups = "drop"
  ) %>%
  unnest_wider(test_result) %>%
  filter(!is.na(p_val)) %>%
  mutate(fdr = p.adjust(p_val, method="BH"))

n_sig <- sum(results_final$fdr < 0.05)
message("High-Confidence ASCA Hits: ", n_sig)

hits <- results_final %>% filter(fdr < 0.05)
#save hits

########################

library(ggplot2)
library(dplyr)


plot_data_main_only <- results_final %>%
  filter(method == "BetaBinomial") %>%
  mutate(
    delta_AF = mean_AF - 0.5,
    sig_category = ifelse(fdr < 0.05, "Significant", "Not Significant")
  )


p_beta <- ggplot(plot_data_main_only, aes(x = delta_AF, y = -log10(fdr))) +
  

  geom_point(data = subset(plot_data_main_only, sig_category == "Not Significant"),
             color = "grey80", size = 1, alpha = 0.5) +
  
  geom_point(data = subset(plot_data_main_only, sig_category == "Significant"),
             color = "#E41A1C", size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey90") +
  
  theme_classic() +
  labs(
    title = "Volcano Plot",
    subtitle = "Excluding perfectly mono-allelic sites to focus on donor-variable imbalance.",
    x = "Allelic Imbalance (Mean AF - 0.5)",
    y = "-log10(FDR)",
    color = "Significance"
  )

print(p_beta)
# ggsave("ASCA_Volcano.pdf", p_quasi, width = 7, height = 6)


