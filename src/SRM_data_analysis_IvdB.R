library(here)
library(tidyr)
library(dplyr)
library(DescTools)
library(rlang)
library(magrittr)

# define outliers based on expected value and Threshold: 
Is.outlier <- function(i, expected.value, Threshold){
  return(i >= expected.value + expected.value * 0.01 * abs(Threshold) | 
           i < expected.value - expected.value * 0.01 * abs(Threshold) | 
           is.na(i) == TRUE)
}
# define extreme outliers: 
Is.Outlier.or.Extreme <- function(i, low.or.high, outlier.or.extreme){
  if(low.or.high == "low") {
    low <- 0.25
    high <- 1
  } else { if(low.or.high == "high") {
    low <- 0
    high <- 0.75
    } else {
    low <- 0.25
    high <- 0.75
    }
  }
  if(outlier.or.extreme == "outlier") {j <- 1.5
  } else {j <- 3
  }
  return(i < quantile(i, low, na.rm = TRUE) - j * IQR(i, na.rm = TRUE) | 
           i > quantile(i, high, na.rm = TRUE) + j * IQR(i, na.rm = TRUE))
}
# transform results table: 
transform_data <- function(df){
  transformed.df <- df %>%
    mutate_at(vars(one_of(c("Actual.Concentration",
                            "Calculated.Concentration", 
                            "Accuracy",
                            "Area", "IS.Area", "Area.Ratio", "Height", 
                            "Quality", "IS.Quality", 
                            "Retention.Time",
                            "Total.Width", "IS.Total.Width",
                            "Width.at.50.", "IS.Width.at.50.",
                            "Signal...Noise","IS.Signal...Noise",
                            "Asymmetry.Factor",
                            "Ion.Ratio"))), 
              as.numeric, na.rm = FALSE) %>%
    separate(Component.Name, 
             c("sp","Accession", "Protein", "Species","Peptide", "Transition", "Type", "Label"),
             sep ="([\\|\\_\\.])",
             remove = FALSE) %>%
    separate(Type,c("Type", "No"), sep=4) %>% 
    separate(Sample.Name, c("Patient.ID", "Timepoint"),
             sep="([\\M])",
             remove=FALSE,
             fill="right") %>%
    mutate(Peptide_Day = paste(Peptide, Batch.Day, sep="_"),
           Protein_Day = paste(Protein, Batch.Day, sep="_"),
           Peptide_Injection = paste(Peptide, Acquisition.Date...Time, sep="_"),
           Peptide_Sample = paste(Peptide, Sample.Name, sep="_"),
           Protein_Injection = paste(Protein, Acquisition.Date...Time,sep="_"),
           Protein_Sample = paste(Protein, Sample.Name, sep="_"),
           Transition_Sample = paste(Component.Name, Sample.Name, sep="_")) 
  return(transformed.df)    
}
# calculate ion ratio for all combinations of transitions (only works for 3 transitions):
calculate.ion.ratio <- function(df){
  df.with.ion.ratio <- df %>%
    arrange(Batch.Day, Sample.Name, Protein, desc(Type), No) %>%
    group_by(Peptide_Injection) %>%
    mutate(Ion.ratio.1_2 = nth(Calculated.Concentration, 2) / first(Calculated.Concentration),
           Ion.ratio.1_3 = nth(Calculated.Concentration, 3) / first(Calculated.Concentration),
           Ion.ratio.2_3 = nth(Calculated.Concentration, 3) / nth(Calculated.Concentration, 2), na.rm = FALSE) %>%
    ungroup()
  return(df.with.ion.ratio)
}
# calculate the % difference between duplicates for unknowns (only works when 2 replicates):
calculate.duplicate.per.transition <- function(df) {
  df.with.duplicate.per.transition <- df %>%
    arrange(Batch.Day, Sample.Name, Protein, desc(Type), No) %>%
    group_by (Transition_Sample) %>%
    mutate(duplicate.per.no = ifelse(Sample.Type == "Unknown",
                                     abs(first(Calculated.Concentration) - nth(Calculated.Concentration, 2))
                                     / mean(Calculated.Concentration) * 100, NA)) %>%
    ungroup()
  return(df.with.duplicate.per.transition)
}
# keep only duplicates, e,g having 2 variables in common for a shared column (e.g., Sample Name) while differentiated by another column (e.g., Sample ID):
only.duplicate <- function(df, shared.column, diff.column){
  if(nrow(df) == 0) {df.with.duplicates = df
  }
  else{
  shared.column <- enquo(shared.column)
  diff.column <- enquo(diff.column)
    Duplicated <- df %>%
      select(!!shared.column, Sample.Type, !!diff.column) %>%
      unique() %>%
      select(!!shared.column, Sample.Type) %>% 
      table() %>%
      data.frame() %>%
      filter(Freq == 2)
  
  df.with.duplicates <- df %>%
    filter(!!shared.column %in% Duplicated[[as_name(shared.column)]])
  }
  return(df.with.duplicates)
}
# Select data (remove bad injections / bad duplicates and bad sample quality)
select.data <- function(df, sample.matrix, remove.SampleQuality) {
  #1: remove/identify samples with bad quality: for tips incorrectly filled, for other matrices when sample comment present
  df.No1 <- df %>%
    filter(No == 1,
           Sample.Type == "Unknown")
  if(remove.SampleQuality == TRUE){
    if(sample.matrix == "tip"){
      df.GoodQuality <- df.No1 %>%
          filter(grepl("both", Sample.Comment, ignore.case = TRUE) == TRUE |
                   (Sample.Comment == "B" & Sample.ID == "B") | 
                   (Sample.Comment == "A" & Sample.ID == "A")) 
    }
    else{
      df.GoodQuality <- df.No1 %>%
          filter(Sample.Comment == "")
    }
  }
  else{
    df.GoodQuality <- df.No1 
  }
# 2 identify extreme low outliers for area or quality in more than 6 quantifiers to indicate bad injection - needs confirmation
  low.extremes <- df.GoodQuality %>%
    group_by(Peptide_Day) %>%
    filter(Is.Outlier.or.Extreme(Area, "low", "extreme") == TRUE |
             is.na(Area) == TRUE |
             Is.Outlier.or.Extreme(IS.Area, "low", "extreme") == TRUE |
             is.na(IS.Area) == TRUE |
             (Is.Outlier.or.Extreme(Quality, "low", "extreme") == TRUE & Quality < 0.7) |
             is.na(Quality) == TRUE |
             (Is.Outlier.or.Extreme(IS.Quality, "low", "extreme") == TRUE & IS.Quality < 0.7) |
             is.na(IS.Quality) == TRUE) %>%
    ungroup()
  Bad.Injection <- low.extremes %>%
    filter(Type == "QUAN") %>%
    select(Acquisition.Date...Time, Sample.Name) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 8) 
  
 # 3: remove / identify both or single replicate with bad duplication
  #a identify samples with 6 or more quantifiers with % difference > 40%
  Bad.duplicate <- df.GoodQuality %>%
    filter(!Acquisition.Date...Time %in% Bad.Injection$Acquisition.Date...Time,
           Type == "QUAN") %>%
    only.duplicate(Sample.Name, Sample.ID) %>%
    filter(duplicate.per.no > 40) %>%
    select(Sample.Name, Batch.Day, Protein_Sample, duplicate.per.no) %>%
    unique() %>%
    select(Sample.Name, Batch.Day) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 5) %>%
    mutate_at(1, as.character) %>%
    mutate_at(2, as.integer) %>%
    left_join(df.No1[, -1], by = "Sample.Name") %>%
    select(Batch.Day, Acquisition.Date...Time, Sample.Name, Sample.ID, Freq) %>%
    unique()
  
  #b identify samples with 6 or more quantifiers with outliers (1.5 * IQR) for Area or Quality
  Outliers <- df.GoodQuality %>%
    filter(!Acquisition.Date...Time %in% Bad.Injection$Acquisition.Date...Time,
           Type == "QUAN") %>%
    group_by(Peptide_Day) %>%
    filter(Is.Outlier.or.Extreme(Area, "both", "outlier") == TRUE |
             is.na(Area) == TRUE |
             Is.Outlier.or.Extreme(IS.Area, "both", "outlier") == TRUE |
             is.na(IS.Area) == TRUE |
             (Is.Outlier.or.Extreme(Quality, "both", "outlier") == TRUE & Quality < 0.7) |
             is.na(Quality) == TRUE |
             (Is.Outlier.or.Extreme(IS.Quality, "both", "outlier") == TRUE & IS.Quality < 0.7) |
             is.na(IS.Quality) == TRUE) %>%
    ungroup() %>%
    select(Sample.ID, Sample.Name) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 5) %>%
    mutate_at(c(1, 2), as.character) %>%
    left_join(df.No1, by = c("Sample.Name", "Sample.ID")) %>%
    select(Batch.Day, Acquisition.Date...Time, Sample.Name, Sample.ID, Freq) %>%
    unique()
    
  #c identify bad duplicate with outlier in 1 replicate
  remove.single.replicate <- intersect(Bad.duplicate[-5], Outliers[-5]) %>%
    setdiff(intersect(Bad.duplicate[-5], Outliers[-5]) %>%
              mutate(Sample.Type = "Unknown") %>%
              only.duplicate(Sample.Name, Sample.ID) %>%
              select(-5))
  
  #d identify bad duplicate without outliers or where both have outliers
  remove.both.replicates <- setdiff(Bad.duplicate[-5], Outliers[-5]) %>%
    filter(!Sample.Name %in% remove.single.replicate$Sample.Name) %>%
    bind_rows(intersect(Bad.duplicate[-5], Outliers[-5]) %>%
                mutate(Sample.Type = "Unknown") %>%
                only.duplicate(Sample.Name, Sample.ID)) %>%
    select(-5)
  
  # 4: ONLY IDENTIFY (do not remove) peptides with high duplicate and extreme outliers for Area: 
  Bad.Peptide <- low.extremes %>%
    filter(!Acquisition.Date...Time %in% Bad.Injection$Acquisition.Date...Time,
           !Sample.Name %in% remove.both.replicates$Sample.Name,
           !Acquisition.Date...Time %in% remove.single.replicate$Acquisition.Date...Time,
           duplicate.per.no > 50)

  # 5: Output is a list with each subgroup identified for evaluation/visual inspection 
  Select.data <- list()
  Select.data$Bad.Quality <- df %>%
    filter(!Acquisition.Date...Time %in% df.GoodQuality$Acquisition.Date...Time,
           Sample.Type == "Unknown")
  Select.data$Bad.Injection <- df %>%
    filter(Acquisition.Date...Time %in% Bad.Injection$Acquisition.Date...Time)
  Select.data$Bad.Duplicate.with.outlier <- df %>%
    filter(Acquisition.Date...Time %in% remove.single.replicate$Acquisition.Date...Time)
  Select.data$Bad.Duplicate <- df %>%
    filter(Acquisition.Date...Time %in% remove.both.replicates$Acquisition.Date...Time)
  Select.data$Bad.Peptide <- df %>%
    filter(Peptide_Injection %in% Bad.Peptide$Peptide_Injection)
  Select.data$standards.and.QC <- df %>%
    filter(Sample.Type != "Unknown")
  Select.data$check.integrations <- Select.data$Bad.Injection %>% 
    #simplified overview for easy check in MultiQuant
    select(Batch.Day, Acquisition.Date...Time, Sample.Name, Sample.ID) %>%
    unique() %>%
    mutate(Peptide = "multiple") %>%
    bind_rows(Select.data$Bad.Peptide %>%
                select(Batch.Day, Acquisition.Date...Time, Sample.Name, Sample.ID, Peptide) %>%
                unique())
  Select.data$Select.df <- setdiff(df, bind_rows(
    Select.data$Bad.Quality, Select.data$Bad.Injection, Select.data$Bad.Duplicate.with.outlier,
    Select.data$Bad.Duplicate)) %>%
    filter(Sample.Type == "Unknown")
  #Select.df keeps BAD PEPTIDE, these are only IDENTIFIED to ALLOW early RE-CHECK of CHROMATOGRAMS

  return(Select.data)
}
# Create table with IR outlier thresholds based on specified dataframe including at least n = x samples from n = y different individuals
define.IR.thresholds <- function(df){
  IR.thresholds <- df %>%
    filter (No == 1,
            !Peptide_Day %in% subset(
              df, grepl("fail", Component.Comment, ignore.case = TRUE) == TRUE)$Peptide_Day) %>%
    group_by(Protein, Peptide, Type) %>%
    summarise(n.sample = n(),
              mean1_2 = mean(Ion.ratio.1_2, na.rm = TRUE),
              CV1_2 = sd(Ion.ratio.1_2, na.rm = TRUE) / mean1_2 * 100,
              mean1_3 = mean(Ion.ratio.1_3, na.rm = TRUE),
              CV1_3 = sd(Ion.ratio.1_3, na.rm = TRUE) / mean1_3 * 100,
              mean2_3 = mean(Ion.ratio.2_3, na.rm = TRUE),
              CV2_3 = sd(Ion.ratio.2_3, na.rm = TRUE) / mean2_3 * 100) %>%
    ungroup() %>%
    arrange(Protein, desc(Type)) %>%
    transmute(Protein = Protein,
              Peptide = Peptide,
              Type = Type,
              n = n.sample,
              TH1_2 = ifelse(Type == "QUAN", ifelse(2 * CV1_2 <= 20, 20, ifelse(2 * CV1_2 >= 50, 50,
                                                                                RoundTo(2 * CV1_2, 5, ceiling))),
                             ifelse(3 * CV1_2 <= 30, 30, ifelse (3* CV1_2 >= 50, 50,
                                                                 RoundTo(3 * CV1_2, 5, ceiling)))),
              TH1_3 = ifelse(Type == "QUAN", ifelse(2 * CV1_3 <= 20, 20, ifelse(2 * CV1_3 >= 50, 50,
                                                                                RoundTo(2 * CV1_3, 5, ceiling))),
                             ifelse(3 * CV1_3 <= 30, 30, ifelse(3 * CV1_3 >= 50, 50,
                                                                RoundTo(3 * CV1_3, 5, ceiling)))),
              TH2_3 = ifelse(Type == "QUAN", ifelse(2 * CV2_3 <= 20, 20, ifelse(2 * CV2_3 >= 50, 50,
                                                                                RoundTo(2 * CV2_3, 5, ceiling))),
                             ifelse(3 * CV2_3 <= 30, 30, ifelse(3 * CV2_3 >= 50, 50,
                                                                RoundTo(3 * CV2_3, 5, ceiling)))))
    return(IR.thresholds)
}
# Create column to define type of outliers for QC and IR:
identify.outliers.QC.IR <- function(df, Thresholds.IR){ 
  # subset where one or more transitions fail QC
  FailQC.1ormore <- df %>% 
    filter(grepl("fail", Component.Comment, ignore.case = TRUE) == TRUE) %>%
    distinct(Peptide_Day, No, .keep_all = TRUE) 
  # subset where two or more transitions fail QC
  FailQC <- FailQC.1ormore %>%
    group_by(Peptide_Day) %>%
    mutate(n.failQC = n()) %>%
    ungroup() %>%
    filter(n.failQC > 1) %>%
    distinct(Peptide_Day, .keep_all = TRUE) %>%
    select(-n.failQC)

  # subsets with IR outliers:
  Outlier.IR <- lapply(unique(df$Peptide), function(Pep) {
      bind_rows(df %>%
                  filter(Peptide == Pep,
                  Peptide_Day %in% FailQC$Peptide_Day |
                    !Peptide_Day %in% FailQC.1ormore$Peptide_Day,
                  No == 1,
                  Is.outlier(Ion.ratio.1_2, 1, subset(
                    Thresholds.IR, Peptide == Pep)$TH1_2) == TRUE,
                  Is.outlier(Ion.ratio.1_3, 1, subset(
                    Thresholds.IR, Peptide == Pep)$TH1_3) == TRUE),
                df %>%
                  filter(Peptide == Pep,
                  !Peptide_Day %in% FailQC$Peptide_Day,
                  Peptide_Day %in% FailQC.1ormore$Peptide_Day,
                  grepl("fail", Component.Comment, ignore.case = TRUE) == TRUE, 
                  (No == 3 & Is.outlier(Ion.ratio.1_2, 1, subset(
                      Thresholds.IR, Peptide == Pep)$TH1_2) == TRUE)|
                  (No == 2 & Is.outlier(Ion.ratio.1_3, 1, subset(
                      Thresholds.IR, Peptide == Pep)$TH1_3) == TRUE)|
                  (No == 1 & Is.outlier(Ion.ratio.2_3, 1, subset(
                      Thresholds.IR, Peptide == Pep)$TH2_3) == TRUE)))
    }) %>%
    bind_rows()
  
  IR.Q2 <- lapply(unique(df$Peptide), function(Pep){
    Outlier.IR %>%
      filter(Peptide == Pep,
             No == 1,
             Is.outlier(Ion.ratio.2_3, 1,
                         subset(Thresholds.IR, Peptide == Pep)$TH2_3) == FALSE)
  }) %>%
    bind_rows()
  
  IR.outlier.x2 <- setdiff(Outlier.IR, IR.Q2) %>%
    only.duplicate(Peptide_Sample, Sample.ID)

  IR.outlier.x1 <- setdiff(setdiff(Outlier.IR, IR.Q2), IR.outlier.x2)
# create 2 additional columns to differentiate type of QC and IR outlier in input df
# some comments, e.g, "Q2.replace" are used for creation of Figures but are treated as "pass" or "no.outliers" in data-analysis process
  df.QC.IR <- df %>%
    mutate(QC = ifelse(!Peptide_Day %in% FailQC.1ormore$Peptide_Day, "Pass",
                     ifelse(Peptide_Day %in% FailQC$Peptide_Day, "Fail", 
                            ifelse(Peptide_Day %in% subset(FailQC.1ormore, No == "1")$Peptide_Day, 
                                   "Q2.replace", "Fail1"))),
            IR = ifelse(Peptide_Injection %in% IR.Q2$Peptide_Injection, "Q2.replace",
                     ifelse(Peptide_Injection %in% IR.outlier.x2$Peptide_Injection, "IRx2",
                            ifelse(Peptide_Injection %in% IR.outlier.x1$Peptide_Injection, "IRx1",
                                   "no.outlier"))))
  return(df.QC.IR)
}
# select transitions (Q1 or Q2):
select.transition <- function(df) {
  df.1transition <- df %>%
    filter(ifelse(QC == "Q2.replace", No == "2",
                  ifelse(IR == "Q2.replace", No == "2",
                         No == "1")))
}
# calculate peptide ratio, mean, and % difference based on selected transitions:
calculate.peptide.ratio.and.mean <- function(df) {
  if(nrow(df) > nrow(distinct(df, Peptide_Injection))) {
    print("Need to select 1 transition per peptide to compute Peptide Ratio")
    } 
  else {
    df.PR.mean <- df %>%
      group_by(Protein_Injection) %>%
      mutate(Peptide.Ratio = ifelse(Type=="QUAL", NA, 
                                Calculated.Concentration[Type=="QUAN"]/Calculated.Concentration[Type=="QUAL"])) %>%
      ungroup() %>%
      group_by(Peptide_Sample) %>%
      mutate(Mean = mean(Calculated.Concentration),
             Diff.Mean = as.numeric(abs(first(Calculated.Concentration) -
                                          nth(Calculated.Concentration, 2)) / 
                                                                       Mean * 100)) %>%
      ungroup()
    
    return(df.PR.mean)
  }
}
# Create table with PR outlier thresholds based on specified dataset including at least n = x samples from n = y different individuals
define.PR.thresholds <- function(df){
  PR.thresholds <- df %>%
    filter(QC != "Fail", 
           QC != "Q2.replace",
           IR == "no.outlier") %>%
    only.duplicate(Protein_Injection, Peptide_Injection) %>%
    filter(Type == "QUAN") %>%
    group_by(Protein, Patient.ID) %>%
    summarise(n.ID = n(),
              mean = mean(Peptide.Ratio),
              sd.square = sd(Peptide.Ratio)^2) %>%
    ungroup() %>%
    drop_na() %>%
    group_by(Protein) %>%
    summarise(n.total = sum(n.ID),
              n.patient = n(),
              median = median(mean),
              CVinter = sd(mean)/median(mean) * 100,
              CVintra = mean(sd.square)^0.5 * 100) %>%
    mutate(median.PR = RoundTo(median, 0.05, round),
           Percent.threshold = ifelse(2 * CVintra <= 20, 20, 
                                      ifelse(2 * CVintra > 40, 40, 
                                             RoundTo(2 * CVintra, 5, ceiling))))
  return(PR.thresholds)
}
# Create column to define outliers for Peptide Ratio
identify.outliers.PR <- function(df, Thresholds.PR) {
  PR.outlier <- bind_rows(lapply(unique(df$Protein), function(Prot){
    df %>%
      filter(Protein == Prot,
             Type == "QUAN",
             Is.outlier(Peptide.Ratio, subset(Thresholds.PR, Protein == Prot)$median.PR,
                        subset(Thresholds.PR, Protein == Prot)$Percent.threshold) == TRUE)
  })) 
  PRx2 <- only.duplicate(PR.outlier, Protein_Sample, Sample.ID)
  PRx1 <- setdiff(PR.outlier, PRx2)
  
  df.PR <- df %>%
    mutate(PR = ifelse(Protein_Injection %in% PRx1$Protein_Injection, "PRx1",
                       ifelse(Protein_Injection %in% PRx2$Protein_Injection, "PRx2",
                              "no.outlier")))
  return(df.PR)
}
# Create table with duplicate Tier1 and Tier2 threshold
define.duplicate.thresholds <- function(df){
  Summary.dup <- df %>%
    filter(QC != "Fail",
           QC != "Q2.replace",
           IR == "no.outlier", 
           PR == "no.outlier") %>%
    only.duplicate(Peptide_Sample, Sample.ID) %>% 
    filter(Sample.ID == "A") %>%
    group_by(Protein, Peptide, Type) %>%
    summarise(n = n(),
              Tier1.threshold = ifelse(median(Diff.Mean) * 2 >= 25, 25, 
                                       RoundTo(median(Diff.Mean) * 2, 5, ceiling)),
              Tier2.threshold = ifelse(median(Diff.Mean) * 3 >= 40, 40, 
                                       RoundTo(median(Diff.Mean) * 3, 5, ceiling))) %>%
    arrange(Protein, desc(Type))
    return(Summary.dup)
}
# Identify outliers for duplicates:
identify.outliers.duplicate <- function(df, Thresholds.Duplicate) {
  Duplicate.outlier.Tier1 <- bind_rows(lapply(unique(df$Peptide), function(Pep){
    df %>%
      filter(Peptide == Pep,
             Diff.Mean >= subset(Thresholds.Duplicate, Peptide == Pep)$Tier1.threshold)
      }))
  Duplicate.outlier.Tier2 <- bind_rows(lapply(unique(df$Peptide), function(Pep){
    df %>%
      filter(Peptide == Pep,
             Diff.Mean >= subset(Thresholds.Duplicate, Peptide == Pep)$Tier2.threshold)
  }))
  df.dup <- df %>%
    mutate(Duplicate = ifelse(!Peptide_Sample %in% only.duplicate(
      df, Peptide_Sample, Sample.ID)$Peptide_Sample, "Single",
      ifelse(Peptide_Sample %in% Duplicate.outlier.Tier2$Peptide_Sample, "Fail.Tier2",
             ifelse(Peptide_Sample %in% Duplicate.outlier.Tier1$Peptide_Sample, "Fail.Tier1",
                                     "Pass"))))
  return(df.dup)
}
# Create final results:
create.report <- function(df, AllowSingles){

remove.Tier1 <- df %>%
  filter(QC == "Fail" | IR == "IRx2" | PR == "PRx2") %>%
  filter(Type == "QUAN")

pass.Tier1 <- setdiff(df, remove.Tier1) %>%
  filter(Duplicate == "Pass",
         Type == "QUAN")
  
# remove QC fail QUAN/QUAL, remove IRx2 QUAN/QUAL remove IRx1 QUAN or QUAL(in singles only)
remove.Tier2 <- setdiff(df, union(remove.Tier1, pass.Tier1)) %>%
  filter(IR == "IRx1" | PR == "PRx1" |
           Protein_Sample %in% subset(df, QC == "Fail" | IR == "IRx2")$Protein_Sample |
           (Duplicate == "Single" & Protein_Sample %in% subset(df, IR == "IRx1")$Protein_Sample)) %>%
  filter(Type == "QUAN")

widows.Tier2 <- df %>% 
  filter(Peptide_Sample %in% remove.Tier2$Peptide_Sample,
         !Peptide_Injection %in% remove.Tier2$Peptide_Injection,
         Duplicate != "Single")

pass.Tier2 <- setdiff(df, union(union(remove.Tier1, remove.Tier2), widows.Tier2)) %>%
  filter(Duplicate == "Fail.Tier1",
         Type == "QUAN")

Extremes <- df %>%
  group_by(Peptide) %>%
  filter(Is.Outlier.or.Extreme(Calculated.Concentration, "both", "extreme") == TRUE) %>%
  ungroup()

bad.duplicates <- setdiff(df, union(union(remove.Tier1, remove.Tier2), widows.Tier2)) %>%
  filter(Duplicate == "Fail.Tier2")

Both.replicates.not.extreme <- bad.duplicates %>%
  setdiff(Extremes) %>%
  only.duplicate(Peptide_Sample, Sample.ID) #both replicates are not extremes

Any.replicate.extreme <- bad.duplicates %>%
  intersect(Extremes)

Widow.with.outlier.IR.QUAL <- setdiff (bad.duplicates, Both.replicates.not.extreme) %>%
  filter(Protein_Injection %in% subset(df, IR == "IRx1")$Protein_Injection)
# IRx1 outliers for QUAL are allowed when duplicates are OK but are removed for single results 

remove.replicate <- union(union(Both.replicates.not.extreme, Any.replicate.extreme),
                          Widow.with.outlier.IR.QUAL) %>%
  filter(Type == "QUAN")

widow.replicate <- setdiff(bad.duplicates, remove.replicate) %>%
  filter(Type == "QUAN")
         
singles <- setdiff(df, union(remove.Tier1, remove.Tier2)) %>%
  filter(Duplicate == "Single", 
         Type == "QUAN")

if(AllowSingles == FALSE) {
  Results <- list(
    remove.Tier1 = remove.Tier1,
    remove.Tier2 = union(remove.Tier2, widows.Tier2),
    remove.replicates = union(remove.replicate, widow.replicate),
    remove.singles = singles,
    results.Duplicate.Tier1 = pass.Tier1,
    results.Duplicate.Tier2 = pass.Tier2)
}
if(AllowSingles == TRUE) {
  Results <- list(
    remove.Tier1 = remove.Tier1,
    remove.Tier2 = remove.Tier2,
    remove.replicates = remove.replicate,
    Extremes.singles =  intersect(singles, Extremes),
    Extremes.widows = intersect(widows.Tier2, Extremes),
    results.singles = setdiff(singles, Extremes),
    results.widows.Tier2 = setdiff(widows.Tier2, Extremes),
    results.widows = widow.replicate,
    results.Duplicate.Tier1 = pass.Tier1,
    results.Duplicate.Tier2 = pass.Tier2)
}
return(Results)
}

# combine above functions into 2 functions for data-analysis process:
# 1) function to define Thresholds for Ion Ratio, Peptide Ratio and Duplicates:
define.QC.thresholds <- function(df.for.thresholds, sample.matrix, remove.SampleQuality) {
  df.select <- df.for.thresholds %>%
    transform_data %>%
    calculate.ion.ratio() %>%
    calculate.duplicate.per.transition() %>%
    select.data(., sample.matrix, remove.SampleQuality)
  
  Thresholds.IR <- define.IR.thresholds(df.select$Select.df)

  df.1transition <- df.select$Select.df %>%
    identify.outliers.QC.IR(., Thresholds.IR) %>%
    select.transition() %>%
    calculate.peptide.ratio.and.mean()

  Thresholds.PR <- define.PR.thresholds(df.1transition)
  
  df.1peptide <- df.1transition %>%
    identify.outliers.PR(., Thresholds.PR) 

  Thresholds.Duplicate <- define.duplicate.thresholds(df.1peptide)
  
  Thresholds <- list(
    IR = Thresholds.IR,
    PR = Thresholds.PR,
    Duplicate = Thresholds.Duplicate)
  
  return(Thresholds)
}
# 2) function to create final results:
create.results <- function(df.for.QC, sample.matrix, remove.SampleQuality, AllowSingles){
  df <- df.for.QC %>%
    transform_data() %>%
    calculate.ion.ratio() %>%
    calculate.duplicate.per.transition() %>%
    select.data(., sample.matrix, remove.SampleQuality) %$% # note the $ to extract element from list, magrittr
    Select.df %>%
    identify.outliers.QC.IR(., Thresholds$IR)%>%
    select.transition() %>%
    calculate.peptide.ratio.and.mean() %>%
    identify.outliers.PR(., Thresholds$PR) %>%
    identify.outliers.duplicate(.,Thresholds$Duplicate) %>%
    create.report(., AllowSingles)
  return(df)
}

#### LOAD DATA - DEFINE THRESHOLDS - CREATE REPORT ----
df.for.thresholds <- read.csv(here("data/Mitra_all.csv"), sep=",", na.strings = "NA", 
                              strip.white = TRUE, header = TRUE, stringsAsFactors = FALSE)
Thresholds <- define.QC.thresholds(df.for.thresholds, "tip", TRUE)
# "tip for Mitra tips with sample comment regarding fill quality)
# Remove.SampleQuality = TRUE will remove identified bad quality samples from data

df.for.QC <- read.csv(here("data/Mitra_all.csv"), sep=",", na.strings = "NA", 
                      strip.white = TRUE, header = TRUE, stringsAsFactors = FALSE)

# Select Data: remove bad samples/injections/duplicates and visually inspect:
Select.Data <- df.for.QC %>%
  transform_data() %>%
  calculate.ion.ratio() %>%
  calculate.duplicate.per.transition() %>%
  select.data(., "tip", TRUE) #a)

# Create Result (AllowSingle=True allows inclusion of single results):
Results.withsingle <- Select.Data$Select.df %>%
  identify.outliers.QC.IR(., Thresholds$IR) %>% #a)
  select.transition() %>%
  calculate.peptide.ratio.and.mean() %>%
  identify.outliers.PR(., Thresholds$PR) %>% #a)
  identify.outliers.duplicate(., Thresholds$Duplicate) %>% #a)
  create.report(., AllowSingles = TRUE) #b)

# Write included data to csv with identifier for groups
write.csv(bind_rows(Results.withsingle[["results.Duplicate.Tier1"]],
                    Results.withsingle[["results.Duplicate.Tier2"]],
                    Results.withsingle[["results.singles"]],
                    Results.withsingle[["results.widows"]],
                    Results.withsingle[["results.widows.Tier2"]],
                    .id = "Tier"),
          here("data/ResultswithSingles.csv"), row.names = FALSE)