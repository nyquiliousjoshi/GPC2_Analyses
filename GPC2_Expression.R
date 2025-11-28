# GPC2 Expression Analyses
# Author: Nikhil Joshi

library(tidyverse)
library(survminer)
library(ggsurvfit)

# Solid Tissue ----

# This data was abstracted from OpenPedCan V15 accessioned via pedcbioportal and the OpenPedCan Project
ClinicalData_GPC2 <- read.delim("~/Downloads/GPC2 Analyses/ClinicalData.BOXPLOT.tsv")
GPC2_Expression <- read.delim("~/Downloads/GPC2 Analyses/Gpc2_Expression.Boxplot.txt")

SummaryData = ClinicalData_GPC2 %>%
  dplyr::filter(SAMPLE_TYPE == "Solid Tissue") %>%
  dplyr::filter(TUMOR_TYPE == "primary") %>%
  inner_join(GPC2_Expression, by = join_by(Sample.ID == SAMPLE_ID)) %>%
  group_by(CANCER_TYPE_DETAILED) %>% summarise(count = n())

# Pulling primary solid tissue samples
ClinicalData_GPC2 <- ClinicalData_GPC2 %>%
  dplyr::filter(SAMPLE_TYPE == "Solid Tissue") %>%
  dplyr::filter(TUMOR_TYPE == "primary" | TUMOR_TYPE == "Primary Tissue" | is.na(TUMOR_TYPE)) %>%
mutate(HISTOLOGY_short = case_when( CANCER_TYPE_DETAILED == "Diffuse midline glioma, H3 K28-altered" ~ "DMG (n = 128)", 
                                    CANCER_TYPE_DETAILED == "Osteosarcoma" ~ "OS (n = 72)", 
                                    CANCER_TYPE_DETAILED == "Embryonal tumor with multilayer rosettes, NOS" ~ "ETMR (n = 4)",
                                    CANCER_TYPE_DETAILED == "Diffuse hemispheric glioma, H3 G35-mutant" ~ "DHG (n = 26)",
                                    CANCER_TYPE_DETAILED == "Medulloblastoma, group 4" ~ "MB, GRP4 (n = 86)",
                                    CANCER_TYPE_DETAILED == "Medulloblastoma, group 3" ~ "MB, GRP3 (n = 47)",
                                    CANCER_TYPE_DETAILED == "Medulloblastoma, SHH-activated" ~ "MB, SHH (n = 48)",
                                    CANCER_TYPE_DETAILED == "Medulloblastoma, WNT-activated" ~ "MB, WNT (n = 22)",
                                    HISTOLOGY == "Low-grade glioma" ~ "LGG (n = 325)",
                                     HISTOLOGY == "Atypical Teratoid Rhabdoid Tumor" ~ "AT/RT (n = 56)",
                                     HISTOLOGY == "High-grade glioma" ~ "HGG (n = 139)",
                                     HISTOLOGY == "Ependymoma" ~ "Ependymoma (n = 112)", 
                                     HISTOLOGY == "Acute Lymphoblastic Leukemia" ~ "ALL (n = 646)",
                                     HISTOLOGY == "Wilms tumor" ~ "Nephroblastoma (n = 124)",
                                     HISTOLOGY == "Rhabdoid tumor of the kidney" ~ "Extracranial MRT (n = 64)",
                                     HISTOLOGY == "Neuroblastoma" ~ "Neuroblastoma (n = 300)",
                                     HISTOLOGY == "Acute Myeloid Leukemia" ~ "AML (n = 215)",
                                     HISTOLOGY == "Chordoma" ~ "Chordoma (n = 4)",
                                     HISTOLOGY == "Brain (GTEx)" ~ "GTEx Brain (n = 2642)"
  )) %>%
  drop_na(HISTOLOGY_short) %>%
  dplyr::select(Sample.ID, CANCER_TYPE_DETAILED, HISTOLOGY_short, OS_STATUS, OS_MONTHS)

GPC2_Plot <- ClinicalData_GPC2 %>% 
  inner_join(GPC2_Expression, by = join_by(Sample.ID == SAMPLE_ID)) %>%
  drop_na(GPC2)
  

GPC2.plot <- GPC2_Plot %>%
  ggplot() +
  geom_boxplot(aes(x = fct_rev(reorder(HISTOLOGY_short, log2(GPC2+1), mean)), y = log2(GPC2+1)), outliers = F, fatten = NULL) +
  geom_jitter(aes(x = reorder(HISTOLOGY_short, -log2(GPC2 + 1), mean), y = log2(GPC2 + 1)), alpha = 0.2, width = 0.1) +
  stat_summary(aes(x = reorder(HISTOLOGY_short, log2(GPC2+1), mean), y = log2(GPC2+1)), fun=median, geom="crossbar", linetype = "solid", width = 0.75, size = 0.3, color = "red") +
  #scale_x_discrete(labels= labels_STUFF)+
  theme_classic() +
  labs(x = "", y = "GPC2 Expression (LogTPM)", title = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(hjust = 1, angle = 90, size = 13))

ggsave(GPC2.plot, filename = "GPC2.boxplot.solidtissue.pdf", width = 10, height = 8, dpi = 1000)

