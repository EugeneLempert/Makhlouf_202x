# Import the raw data
Input_list <-
  list.files(path = "./Input/",
             pattern = "*.txt", 
             full.names = TRUE) %>% 
  lapply(read_plus)

Input_list_names <-   
  list.files(path = "./Input/",
             pattern = "*.txt") %>%
  str_replace(pattern = ".txt", 
              replacement = "")

names(Input_list) <- Input_list_names

# Import the primer efficiency file
AA10 <- read_csv("AA10.csv")

# Plate 1-1----
OR1A1_V1_vs_V2.2_1 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220202`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Acsm4", "Acss2", "Slc25a35"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("One", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_1 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_1, 
                                            MeltFile = Input_list$`Exp-Melt-220202`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_1 <- OR1A1_V1_vs_V2.2_QE_1$QE_table %>% Plot_QE(MutationColumn = "Mutation")

# Plate 1-2----

OR1A1_V1_vs_V2.2_2 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220203_1`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Acsm4", "Acss2", "Slc25a35"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Two", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_2 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_2, 
                                            MeltFile = Input_list$`Exp-Melt-220203_1`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_2 <- OR1A1_V1_vs_V2.2_QE_2$QE_table %>% Plot_QE(MutationColumn = "Mutation")

# Plate 1-3----

OR1A1_V1_vs_V2.2_3 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220203_2`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr358", "Olfr609", "Olfr690", "Acsm4", "Acss2", "Slc25a35"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Three", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_3 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_3, 
                                            MeltFile = Input_list$`Exp-Melt-220203_2`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_3 <- OR1A1_V1_vs_V2.2_QE_3$QE_table %>% Plot_QE(MutationColumn = "Mutation")

# Plate 2-1----
OR1A1_V1_vs_V2.2_4 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220218_1`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Four", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_4 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_4, 
                                            MeltFile = Input_list$`Exp-Melt-220218_1`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_4 <- OR1A1_V1_vs_V2.2_QE_4$QE_table %>% Plot_QE(MutationColumn = "Mutation")

# Plate 2-2----
OR1A1_V1_vs_V2.2_5 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220218_2`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Five", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_5 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_5, 
                                            MeltFile = Input_list$`Exp-Melt-220218_2`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_5 <- OR1A1_V1_vs_V2.2_QE_5$QE_table %>% Plot_QE(MutationColumn = "Mutation")

# Plate 2-3----
OR1A1_V1_vs_V2.2_6 <- Label_qPCR_data(CqFile = Input_list$`Exp-Cq-220218_3`, 
                                      Style = "Simple", 
                                      SampleNames = c("MSS-3184_V2.1", "MSS-3254_V1.1.2", "MSS-3255_V1.1.2", "MSS-3256_V1.1.2", "MSS-3257_V1.1.2", "MSS-3258_V1.1.2",
                                                      "MSS-3337_V1.1.2", "MSS-3338_V1.1.2", "MSS-3339_V1.1.2", "MSS-3340_V1.1.2", 
                                                      "MSS-3372_V2.2", "MSS-3373_V2.2", "MSS-3374_V2.2", "MSS-3388_V2.2", "MSS-3389_V2.2", "MSS-3390_V2.2", 
                                                      "WT_gDNA", "UPwater"), 
                                      GeneNames = c("Olfr1154", "Acsm4", "Acss2", "Slc25a35", "Olfr596_603", "Olfr390", "Olfr510"),
                                      Mutation = c("WT_2.1", "WT_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", "OR1A1_1.1.2", 
                                                   "OR1A1_1.1.2", "WT_1.1.2", "OR1A1_1.1.2", "WT_1.1.2", 
                                                   "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "WT_2.2", "OR1A1_2.2", "OR1A1_2.2", "WT_gDNA", "Blank"),
                                      Plate = rep("Six", times = 18),
                                      Replicates = 3)

OR1A1_V1_vs_V2.2_QE_6 <- Quality_Evaluation(OR1A1_V1_vs_V2.2_6, 
                                            MeltFile = Input_list$`Exp-Melt-220218_3`,
                                            RepDiffAllow = 0.6)

OR1A1_V1_vs_V2.2_QE_Plot_6 <- OR1A1_V1_vs_V2.2_QE_6$QE_table %>% Plot_QE(MutationColumn = "Mutation")

# Generating the relative expression and related graphs----
Total_OR1A1_V1_vs_V2.2_QE <- OR1A1_V1_vs_V2.2_QE_1$QE_table %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_2$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_3$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_4$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_5$QE_table) %>%
  bind_rows(OR1A1_V1_vs_V2.2_QE_6$QE_table)

# Filtering out poor performing wells
Total1A1_Input <- Total_OR1A1_V1_vs_V2.2_QE %>%
  filter(Gene != "Olfr609") %>%
  filter(MaxCq == FALSE) %>%
  filter(FalseCq == FALSE) %>%
  filter(Sample %notin% c("WT_gDNA", "UPwater")) %>% 
  filter(Cq <= 30) %>% 
  filter(Cq != 29.8) %>% 
  filter(Peak == "1") %>%
  filter(!(Pos == "J16" & Plate == "Six")) %>%
  filter(!(Pos == "A1" & Plate == "Three")) %>%
  filter(!(Pos == "D23" & Plate == "Six")) %>%
  filter(!(Pos == "A19" & Plate == "Four")) %>%
  filter(!(Pos == "A24" & Plate == "Five")) %>%
  filter(!(Pos == "B18" & Plate == "One")) %>%
  filter(!(Pos == "D19" & Plate == "Four")) %>%
  filter(!(Pos == "M12" & Plate == "Five")) %>%
  filter(!(Pos == "A23" & Plate == "Six"))

Raw_values <- Total1A1_Input %>%
  select(Sample, Mutation, Gene, Plate, Replicate, Cq) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", 
                                        "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154")),
         Plate = factor(Plate, levels = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")),
         Genotype = case_when(Mutation == "OR1A1_1.1.2" ~ "OR1A1-Cherry",
                              Mutation == "OR1A1_2.2" ~ "OR1A1-GCaMP",
                              TRUE ~ "WT")) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "OR1A1-Cherry", "OR1A1-GCaMP"))) %>%
  mutate(Mouse_Line = case_when(Mutation %in% c("WT_1.1.2", "OR1A1_1.1.2") ~ "5x21 OR1A1 Cherry V1.1.2",
                                Mutation %in% c("WT_2.2", "OR1A1_2.2") ~ "9x21 OR1A1 GCaMP V2.2",
                                TRUE ~ "5x21 OR1A1 GCaMP V2.1")) %>%
  rename(Ct = Cq) %>%
  select(Sample, Mouse_Line, Genotype, Gene, Plate, Replicate, Ct)

# The first data compaction step, taking the mean of triplicates on each plate
Intermediate_Data_1A1 <- Raw_values %>%
  group_by(Sample, Mouse_Line, Genotype, Gene, Plate) %>%
  summarize(MeanSampleCt = mean(Ct),
            SampleSD = sd(Ct))

# Calculating the Efficiency-corrected DCt values for each sample/gene from each plate
Interruption1 <- Intermediate_Data_1A1 %>%
  select(-SampleSD) %>%
  left_join(AA10, by = "Gene") %>%
  ungroup() %>%
  mutate(MeanSCtxlog2Eff = MeanSampleCt * log2(Primer_Efficiency))

Interruption2 <- Interruption1 %>%
  filter(Gene %in% c("Acsm4", "Acss2", "Slc25a35")) %>%
  select(Sample, Gene, Plate, MeanSCtxlog2Eff) %>%
  pivot_wider(names_from = Gene, values_from = MeanSCtxlog2Eff) %>%
  rename(MeanSCtxlog2Eff_Acsm4 = Acsm4, 
         MeanSCtxlog2Eff_Acss2 = Acss2,
         MeanSCtxlog2Eff_Slc25a35 = Slc25a35)

Plate_Sample_DCt_1A1 <- Interruption1 %>% 
  left_join(Interruption2, by = c("Sample", "Plate")) %>%
  ungroup() %>%
  mutate(DCtSample = (MeanSCtxlog2Eff_Acsm4 + MeanSCtxlog2Eff_Acss2 + MeanSCtxlog2Eff_Slc25a35)/3 - MeanSCtxlog2Eff)

# Taking the mean of the DCtSample values from each plate to unify across plates
Sample_MeanDCt_SEM <- Plate_Sample_DCt_1A1 %>%
  ungroup() %>% 
  group_by(Sample, Mouse_Line, Genotype, Gene) %>% 
  summarize(MeanDCtSample = mean(DCtSample), 
            DCtSampleSEM = sd(DCtSample)/sqrt(length(DCtSample)))

# Calculating a DCt for each genotype/gene pair
Genotype_MeanDCt_SEM <- Sample_MeanDCt_SEM %>%
  ungroup() %>% 
  group_by(Genotype, Gene) %>% 
  summarize(MeanDCtGenotype = mean(MeanDCtSample),
            DCtGenotypeSize = length(MeanDCtSample),
            DCtGenotypeSEM = sd(MeanDCtSample)/sqrt(DCtGenotypeSize))

# Calculating the DDCt for each genotype/gene pair
WTGenotype_MeanDCt_SEM <- Genotype_MeanDCt_SEM %>%
  filter(Genotype == "WT") %>%
  rename(DCtWTSize = DCtGenotypeSize, 
         MeanDCtWT = MeanDCtGenotype, 
         DCtWTSEM = DCtGenotypeSEM) %>%
  ungroup() %>%
  select(-Genotype)

DDCt_Genotype_SEM_CI <- Genotype_MeanDCt_SEM %>%
  left_join(WTGenotype_MeanDCt_SEM, by = "Gene") %>%
  mutate(MeanDDCt = MeanDCtGenotype - MeanDCtWT, 
         DDCtSEM = sqrt(DCtGenotypeSEM^2 + DCtWTSEM^2),
         DDCtCI = qt(0.975, df = DCtGenotypeSize + DCtWTSize - 2) * DDCtSEM)

# DDCt for each sample
DDCt_Sample <- Sample_MeanDCt_SEM %>%
  ungroup() %>%
  left_join(WTGenotype_MeanDCt_SEM, by = "Gene") %>%
  select(Sample, Gene, Genotype, MeanDCtSample, MeanDCtWT) %>%
  mutate(SampleDDCt = MeanDCtSample - MeanDCtWT)

# Plots
DDCtSampleGenotype <- DDCt_Sample %>%
  left_join(DDCt_Genotype_SEM_CI, by = c("Gene", "Genotype")) %>%
  mutate(Gene = factor(Gene, levels = c("Acsm4", "Acss2", "Slc25a35", "Olfr358", "Olfr390", "Olfr510", "Olfr596_603", "Olfr609", "Olfr690", "Olfr1154"))) %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "OR1A1-Cherry", "OR1A1-GCaMP"))) %>%
  mutate(Class = case_when(Gene %in% c("Olfr596_603", "Olfr690") ~ "Class I", 
                           TRUE ~ "Class II")) %>%
  left_join(DVI2, by = c("Gene" = "Symbol")) %>%
  mutate(DVI = ifelse(is.na(DVI), 1.05, DVI)) %>%
  mutate(Gene = factor(Gene, levels = c("Olfr596_603", "Olfr690", "Olfr510", "Olfr1154", "Olfr358", "Olfr390")))

ggplot(data = DDCtSampleGenotype) +
  geom_point(aes(Gene, SampleDDCt, color = Genotype, group = Genotype), position = position_dodge(width = 0.6), size = 3) +
  scale_color_manual(values=c("WT" = "yellow2", "OR1A1-Cherry" = "red3", "OR1A1-GCaMP" = "skyblue1")) +
  geom_errorbar(aes(x = Gene, group = Genotype, ymin = MeanDDCt - DDCtCI, ymax = MeanDDCt + DDCtCI), 
                width = 0.4, position = position_dodge(0.6), color = "black") +
  scale_y_continuous(breaks=seq(-3,1,1)) +
  labs(title = "Log2 Fold Change in Ct values for Olfr mRNA in OR1A1-Cherry and OR1A1-GCaMP", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(size = 17),
        legend.position = c(0.1, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank()) +
  geom_vline(xintercept = 2.5, linetype="dotted", color = "grey40", size=1.5) +
  annotate("text", label = "Class I", x = 1.5, y = 1, size = 8, colour = "blue", fontface = 2) +
  annotate("text", label = "Class II", x = 4.5, y = 1, size = 8, colour = "red", fontface = 2) +
  annotate("text", label = "DVI: 1.05", x = 1, y = -2.7, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 2, y = -2.7, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 3, y = -2.7, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.05", x = 4, y = -2.7, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.5", x = 5, y = -2.7, size = 4, colour = "black") +
  annotate("text", label = "DVI: 3.6", x = 6, y = -2.7, size = 4, colour = "black")

DDCtSampleGenotype358_690 <- DDCtSampleGenotype %>% filter(Gene %in% c("Olfr358", "Olfr690"))

ggplot(data = DDCtSampleGenotype358_690) +
  geom_point(aes(Gene, SampleDDCt, color = Genotype, group = Genotype), position = position_dodge(width = 0.6), size = 3) +
  scale_color_manual(values=c("WT" = "black", "OR1A1-Cherry" = "red2", "OR1A1-GCaMP" = "greenyellow")) +
  geom_errorbar(aes(x = Gene, group = Genotype, ymin = MeanDDCt - DDCtCI, ymax = MeanDDCt + DDCtCI), 
                width = 0.4, position = position_dodge(0.6), color = "black") +
  scale_y_continuous(breaks=seq(-3,1,1)) +
  labs(title = "Olfr Log2 Fold Change in OR1A1-Cherry and OR1A1-GCaMP", 
       y = expression("Log2 Fold Change (-" ~ Delta*Delta ~ "Ct)") ) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(face="bold", size = 15),
        axis.title.y = element_text(size = 17),
        legend.position = c(0.1, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(face="bold"),
        legend.title = element_blank()) +
  geom_vline(xintercept = 1.5, linetype="dotted", color = "grey40", size=1.5) +
  annotate("text", label = "Class I", x = 1, y = 0.8, size = 8, colour = "blue", fontface = 2) +
  annotate("text", label = "Class II", x = 2, y = 0.8, size = 8, colour = "red", fontface = 2) +
  annotate("text", label = "DVI: 1.05", x = 1, y = -2.7, size = 4, colour = "black") +
  annotate("text", label = "DVI: 1.5", x = 2, y = -2.7, size = 4, colour = "black")

# Statistics----
ForStats <- DDCt_Sample %>%
  filter(Gene %notin% c("Acsm4", "Acss2", "Slc25a35")) #Removing the Reference Gene values from the table

#We will only run the Gene-specific analyses
ggplot(ForStats, aes(SampleDDCt)) + geom_histogram(binwidth = 0.25) +  facet_grid(rows = vars(Gene), cols = vars(Genotype))
ggplot(ForStats, aes(SampleDDCt, Genotype)) + geom_boxplot() + facet_wrap(~Gene)

#A quick ANOVA and krusal test for each gene group
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "anova", group.by = c("Gene"))
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "kruskal.test", group.by = c("Gene"))

#Making subgroups and running the toughest homoscedasticity test 
ForStats_358 <- ForStats %>% filter(Gene == "Olfr358")
ForStats_390 <- ForStats %>% filter(Gene == "Olfr390")
ForStats_510 <- ForStats %>% filter(Gene == "Olfr510")
ForStats_596_603 <- ForStats %>% filter(Gene == "Olfr596_603")
ForStats_690 <- ForStats %>% filter(Gene == "Olfr690")
ForStats_1154 <- ForStats %>% filter(Gene == "Olfr1154")

fligner.test(SampleDDCt ~ Genotype, data = ForStats_358)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_390)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_510)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_596_603)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_690)
fligner.test(SampleDDCt ~ Genotype, data = ForStats_1154)

bartlett.test(SampleDDCt ~ Genotype, data = ForStats_358)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_390)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_510)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_596_603)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_690)
bartlett.test(SampleDDCt ~ Genotype, data = ForStats_1154)

leveneTest(SampleDDCt ~ Genotype, data = ForStats_358)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_390)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_510)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_596_603)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_690)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_1154)

leveneTest(SampleDDCt ~ Genotype, data = ForStats_358, center = mean)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_390, center = mean)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_510, center = mean)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_596_603, center = mean)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_690, center = mean)
leveneTest(SampleDDCt ~ Genotype, data = ForStats_1154, center = mean)

#ANOVA for residuals for normality
AOV358 <- aov(SampleDDCt ~ Genotype, data = ForStats_358)
shapiro.test(AOV358$residuals)
AOV390 <- aov(SampleDDCt ~ Genotype, data = ForStats_390)
shapiro.test(AOV390$residuals)
AOV510 <- aov(SampleDDCt ~ Genotype, data = ForStats_510)
shapiro.test(AOV510$residuals)
AOV596_603 <- aov(SampleDDCt ~ Genotype, data = ForStats_596_603)
shapiro.test(AOV596_603$residuals)
AOV690 <- aov(SampleDDCt ~ Genotype, data = ForStats_690)
shapiro.test(AOV690$residuals)
AOV1154 <- aov(SampleDDCt ~ Genotype, data = ForStats_1154)
shapiro.test(AOV1154$residuals)

#Performing parametric and non-parametric pairwise comparisons
compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "wilcox.test", group.by = c("Gene"))

compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "t.test", group.by = c("Gene")) %>%
  select(Gene, group1, group2, p) %>%
  arrange(p) %>%
  mutate(row = row_number(),
         HolmAdjustment = 0.05/(18 - row + 1), 
         PassHA = p < HolmAdjustment, 
         HolmCorrection = p * (18 - row + 1),
         PassHC = HolmCorrection < 0.05, 
         Pass_DunnSidakCorrection = p < 1 - (1 - 0.05)^(1/18)) %>%
  mutate(NEWp = p.adjust(p, method = "holm"))

Table_1A1_2 <- compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "t.test", group.by = c("Gene")) %>%
  mutate(Comparison = paste(group1, "vs", group2)) %>%
  select(Gene, Comparison, p.adj) %>% 
  mutate(p.adj = case_when(p.adj < 0.001 ~ "***",
                           p.adj < 0.01 ~ "**",
                           p.adj < 0.05 ~ "*", 
                           TRUE ~ "NS")) %>%
  pivot_wider(names_from = Gene, values_from = p.adj) %>%
  select(Comparison, Olfr596_603, Olfr690, Olfr510, Olfr1154, Olfr358, Olfr390)

# Making the combined Table/Text objects
TextGrob1A1 <- text_grob(paste("Eighteen independent t.tests were performed. With alpha = 0.05, Holm-Bonferroni and Dunn-Šidák corrections were calculated. ", 
                               "Dunn-Šidák corrected alpha evaluated to 0.002845571. Boxes in Red are NS comparing raw p-values to the DS alpha value.",
                               "Corrected p-value < 0.001 = ***, < 0.01 = **, < 0.05 = *, Not Significant = NS", sep = "\n"), 
                         face = "italic")

TableGrob_1A1_2 <- tableGrob(Table_1A1_2, rows=NULL) 

find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

ind1 <- find_cell(TableGrob_1A1_2, 4, 3, "core-bg")
TableGrob_1A1_2$grobs[ind1][[1]][["gp"]] <- gpar(fill="red2", col = "red4", lwd=5)

TableGrob_1A1_2 %>%
  grid.arrange(TextGrob1A1, heights = c(0.8, 0.5), as.table = TRUE)

Table_1A1_3 <- compare_means(SampleDDCt ~ Genotype, data = ForStats, method = "t.test", group.by = c("Gene")) %>%
  mutate(Comparison = paste(group1, "vs", group2)) %>%
  select(Gene, Comparison, p.adj) %>% 
  pivot_wider(names_from = Gene, values_from = p.adj) %>%
  select(Comparison, Olfr596_603, Olfr690, Olfr510, Olfr1154, Olfr358, Olfr390)

TableGrob_1A1_3 <- tableGrob(Table_1A1_3, rows=NULL) 

# Code used to color in a cell in a table 
ind1 <- find_cell(TableGrob_1A1_3, 4, 3, "core-bg")
TableGrob_1A1_3$grobs[ind1][[1]][["gp"]] <- gpar(fill="red2", col = "red4", lwd=5)

TableGrob_1A1_3 %>%
  grid.arrange(as.table = TRUE)

# Making other tables with RNASeq data
Comp_Table <- tibble(Gene = c("Olfr690", "Olfr358"), 
                     "RNASeq OR1A1-Cherry " = c("-2.577 (***)", "-2.568 (***)"), 
                     "RTqPCR OR1A1-Cherry" = c("-1.861 (***)", "-1.872 (***)"),
                     "RTqPCR OR1A1-GCaMP" = c("-1.202 (***)", "-2.009 (***)"))

Comp_Table <- Comp_Table %>% 
  pivot_longer(cols = -Gene, names_to = "Analysis", values_to = "Log2 Fold Change (sig.)") %>% 
  pivot_wider(names_from = Gene, values_from = `Log2 Fold Change (sig.)`)

CompGrob_1A1 <- tableGrob(Comp_Table, rows=NULL) 

TitleGrob_1A1 <- text_grob("Log2 Fold Change Comparisons", face = "bold", size = 20)

SubPlotGrob_1A1 <- text_grob(paste("All comparisons relative to WT Samples.",
                                   "Adjusted p-value < 0.001 = (***)", sep = "\n"), face = "italic")

grid.arrange(TitleGrob_1A1, CompGrob_1A1, SubPlotGrob_1A1, heights = c(0.3, 0.4, 0.3), as.table = TRUE)
# CV calculations----
##Raw value Intra-Assay CV: 1.71
Total_OR1A1_V1_vs_V2.2_QE %>%
  ungroup() %>%
  select(Plate, Gene, Sample, RepMean, RepSd) %>%
  distinct() %>%
  mutate(CV = RepSd/RepMean * 100) %>%
  group_by(Plate) %>%
  summarize(MeanCV = mean(CV, na.rm = TRUE)) %>%
  summarize(mean = mean(MeanCV))

##Filtered Raw value Intra-Assay CV: 0.726
Total1A1_Input %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, Cq) %>%
  distinct() %>%
  group_by(Plate, Sample, Gene) %>%
  summarize(Plate = Plate, 
            Sample = Sample, 
            Gene = Gene, 
            RepMean = mean(Cq, na.rm = TRUE),
            RepSd = sd(Cq, na.rm = TRUE)) %>% 
  ungroup() %>% 
  distinct() %>%
  select(Plate, Sample, Gene, RepMean, RepSd) %>% 
  distinct() %>%
  mutate(CV = RepSd/RepMean * 100) %>%
  group_by(Plate) %>%
  summarize(MeanCV = mean(CV, na.rm = TRUE)) %>%
  summarize(mean = mean(MeanCV))

Final2 <- Final_OR1A1 %>%
  filter(Sample != "MSS-3184_V2.1") %>%
  filter(Gene != "Olfr609")

##Raw value Inter-Assay CV: 2.36
Total_OR1A1_V1_vs_V2.2_QE %>%
  ungroup() %>%
  select(Plate, Gene, Sample, RepMean, RepSd) %>%
  distinct() %>%
  group_by(Sample, Gene) %>%
  summarize(PlateMeans = mean(RepMean, na.rm = TRUE), 
            PlateSd = sd(RepMean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  summarise(CV = PlateSd/PlateMeans * 100) %>%
  summarize(meanCV = mean(CV, na.rm = TRUE))

##Filtered Raw value Inter-Assay CV: 0.738
Total1A1_Input %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, Cq) %>%
  distinct() %>%
  group_by(Plate, Sample, Gene) %>%
  summarize(Plate = Plate, 
            Sample = Sample, 
            Gene = Gene, 
            RepMean = mean(Cq, na.rm = TRUE),
            RepSd = sd(Cq, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(Plate, Sample, Gene, RepMean, RepSd) %>% 
  distinct() %>%
  group_by(Sample, Gene) %>%
  summarize(PlateMeans = mean(RepMean, na.rm = TRUE), 
            PlateSd = sd(RepMean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  summarise(CV = PlateSd/PlateMeans * 100) %>%
  summarize(meanCV = mean(CV))
