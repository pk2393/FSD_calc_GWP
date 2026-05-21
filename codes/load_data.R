library(dplyr)
library(tidyr)
library(magrittr)

data_file <- "data_mers/table_from_cauchemez_lancetID.csv"
if (!file.exists(data_file)) {
  data_file <- file.path("..", "data_mers", "table_from_cauchemez_lancetID.csv")
}

df <- readr::read_csv(file = data_file)

df %<>%
  separate(`Cluster ID`, into = c("Cluster_num", "Star"), sep = "(?=\\*)", remove = FALSE)%>%
  mutate(
    age_num = suppressWarnings(as.numeric(age)), 
    category = case_when(
      is.na(age_num) ~ "NR",              
      age_num < 20   ~ "child",           
      TRUE           ~ "adult"
    )
  )

df_summ <- df  %>%
  group_by(Cluster_num) %>%
  summarise(
    n_reports = n(), 
    n_child = sum(category == "child", na.rm = TRUE),
    n_adult = sum(category == "adult", na.rm = TRUE),
    has_NR  = any(category == "NR"), 
    first_report_date = min(`Date reported`, na.rm = TRUE)
  )

df_summ$index_category <- NA

for(n in 1:nrow(df_summ)){
  clusterID = df_summ$Cluster_num[n]
  date = df_summ$first_report_date[n]
  
  
  (df %>% filter(
    Cluster_num == clusterID, 
    `Date reported` == date))$category ->   df_summ$index_category[n]
  
}

