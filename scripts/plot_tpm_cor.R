# install.packages("tidyverse")

library(tidyverse)
library(ggplot2)


data_path <- "/Users/haydenji/Desktop/Classes/ela/results"
it <- "50it"
gt_path <- file.path(data_path, "na12878_dRNA_transcriptome_quantification.tsv")
#pred_path <- file.path(data_path, it, "na12878_sm_both.abundance.modified.tsv")
pred_path <- file.path(data_path, "50it_confCo", "na12878_sm_both.abundance.modified.tsv")
gt_quant <- read_tsv(gt_path)
gt_quant <- gt_quant %>% rename(tpm = TPM) %>% rename(target_id = ID) %>% 
  select(target_id, tpm) %>% filter(tpm > 0) %>% mutate(tpm = tpm + 1e-5) %>% 
  mutate(tpm = log2(tpm))
head(gt_quant) 
pred_quant <- read_tsv(pred_path)
head(pred_quant)

pred_quant <- pred_quant %>% rename(target_id = transcript_id) %>% 
  filter(target_id != "None") %>% filter(tpm > 0) %>% mutate(tpm = log2(tpm))
# mutate(tpm = tpm + 1e-5) 
head(pred_quant)
quant_df <- inner_join(pred_quant, gt_quant, by = "target_id")

quant_df <- quant_df %>% rename(predicted= tpm.x) %>% rename(actual = tpm.y)

head(quant_df)

g <- ggplot(quant_df, aes(x=actual, y=predicted)) + geom_point(alpha=0.05) + 
  ggtitle("TranSigner vs. Ground Truth")
#  xlim(0, 15) + ylim(0, 15)
plot(g)

