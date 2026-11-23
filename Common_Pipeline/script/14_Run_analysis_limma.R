
# LOAD LIBRARY
suppressMessages({
  library(readxl)
})


### DECIDE WHICH SCRIPT TO RUN DEPENDING ON THE NR OF MALES AND FEMALES AFTER QC
samplesheet <- as.data.frame(read_xlsx("output/Tables/samplesheet_after_QC.xlsx"))

if (multiple_samples == "yes"){
  samplesheet <- samplesheet[which(samplesheet$Baseline_sample == "Y"),]
}

nr_males <- nrow(samplesheet[which(samplesheet$Sex == "M"), ])
nr_females <- nrow(samplesheet[which(samplesheet$Sex == "F"), ])



  if (tissue_type == "blood"){
  
    if (baseline_severity_analysis == "yes"){
      print("Running sex-adjusted baseline severity analysis.")
      
      source("script/14_BaselineSeverity_Analysis_Females_Males_limma_blood.R")
      
      # CREATE A REPORT
      print("Create report")
      source("2_Paths.R")
      rmarkdown::render("script/Report_Part4_BaselineSeverity_sex-adjusted.Rmd", 
                        output_file = paste0(dir_gen,"Report_Part5_BaselineSeverity_Analysis_sex_adjusted.html"),
                        params = list(dir_gen = dir_gen)
      )
    }
  }
  
  if (tissue_type == "saliva"){
    if (baseline_severity_analysis == "yes"){
      print("Running sex-adjusted baseline severity analysis.")
      
      source("script/14_BaselineSeverity_Analysis_Females_Males_limma_saliva.R")
      
      # CREATE A REPORT
      print("Create report")
      source("2_Paths.R")
      rmarkdown::render("script/Report_Part4_BaselineSeverity_sex-adjusted.Rmd", 
                        output_file = paste0(dir_gen,"Report_Part5_BaselineSeverity_Analysis_sex_adjusted.html"),
                        params = list(dir_gen = dir_gen)
      )
    }
  }




