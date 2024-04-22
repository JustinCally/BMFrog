library(tidyverse)

# species we are interested in: add to this list if there are others
species_list <- c(`Barking Marsh Frog` = 13059,
                  `Eastern Sign-bearing Froglet` = 13131,
                  `Common Froglet` = 13134,
                  # `Sloane's Froglet` = 13135,
                  `Common Spadefoot Toad` = 13086,
                  `Peron's Tree Frog` = 13204,
                  `Pobblebonk Frog` = 63913,
                  `Spotted Marsh Frog (race unknown)` = 13063)

nspecies <- length(species_list)

sp_l_col <- paste0("sp_", species_list)

# site_records <- read_csv("data/TLM_2022-23_Barmah-Millewa_FrogsEnsemble1_2023.csv") %>%
#   rename(Site_prefix = Site) %>%
#   tidyr::separate_wider_delim(cols = "Filename", delim = "/",
#                               names = c(NA,NA,NA,"Site", NA))

site_records <- read_csv("data/model_results/TLM_Barmah_Frogs_March2024.csv") %>%
  # rename(Site_prefix = Site) %>%
  tidyr::separate_wider_delim(cols = "Filename", delim = "/",
                              names = c(NA,NA,NA,NA, "Site", "HexDate"), too_many = "debug")

site_records_cleaned <- site_records %>%
  filter(HexDate != "NOISE") %>%
  mutate(HexDate = stringr::str_remove_all(HexDate, ".WAV"),
         Site = case_when(stringr::str_detect(Site, "_AM") ~ Site,
                          TRUE ~ paste0(Site, "_AM1"))) #%>%
  # sample_n(10000)

underscores <- stringr::str_detect(site_records_cleaned$HexDate, "_")

site_records_cleaned$DateTime <- NA_POSIXct_

site_records_cleaned$DateTime[underscores] <- as.POSIXct(site_records_cleaned$HexDate[underscores],
                                                            format = c("%Y%m%d_%H%M%OS"), tz = Sys.timezone())

site_records_cleaned$DateTime[!underscores] <- as.POSIXct(as.numeric(as.hexmode(site_records_cleaned$HexDate[!underscores])),
                                                             origin = "1970-01-01", tz = "UTC")

site_records_cleaned$DateTime <- as.POSIXct(site_records_cleaned$DateTime, tz = "Australia/Melbourne")

times <- lubridate::hour(site_records_cleaned$DateTime)
hist(times[!underscores])
#list all of species codes and names (Is this all the classes in the model currently?)
class_list<-read_csv("data/ClassList.csv")

#data with column containing string of category probabilities from the AI
site_records_cleaned$probvec <- site_records_cleaned[[13]]
probcolname <- names(site_records_cleaned)[13]

#clean and split the probabilities
probs <- site_records_cleaned %>% select(probvec) %>%
  mutate(probvec=str_remove(probvec, "\n")) %>%
  mutate(probvec=str_remove(probvec, "\\[")) %>%
  mutate(probvec=str_remove(probvec, "\\]")) %>%
  mutate(probvec=str_squish(probvec)) %>%
  separate(probvec, into=paste0("sp_",class_list$TaxonId), sep="\\s") %>%
  mutate(across(`sp_-2`:`sp_525739`, as.numeric))

# filter into the species we are interested in
frog_probs <- probs %>%
  dplyr::select(all_of(sp_l_col)) %>%
  `colnames<-`(names(species_list))

site_records_date <- site_records_cleaned %>%
  mutate(Date = as.Date(DateTime - hours(12))) %>%
  bind_cols(frog_probs)

daily_max_record <- list()

for(i in 1:nspecies) {
  froggo_name <- names(species_list)[i]
  daily_max_record[[froggo_name]] <- site_records_date %>%
    group_by(Site, Date) %>%
    filter(!!sym(froggo_name) == max(!!sym(froggo_name))) %>%
    sample_n(1) %>% # randomly sample if multiple cases of there being a 'top' score
    select(FileId,
           Filename,
           Site,
           Date,
           DateTime,
           Start_Time, End_Time,
           Orig_Start = Start_Time, Orig_End = End_Time,
           Observed_TaxonId = Predict_TaxonId,
           Observed_Label = Predict_Label,
           Predict_TaxonId,
           Predict_Label,
           Probability,
           all_of(names(species_list)),
           !!sym(probcolname),
           )

}

records_to_validate <- list()

for(i in 1:nspecies) {
  froggo_name <- names(species_list)[i]
  records_to_validate[[froggo_name]] <- daily_max_record[[froggo_name]] %>%
    mutate(ProbGrp = cut(!!sym(froggo_name),
                         breaks = c(0,0.2,0.4,0.6,0.8,1), include.lowest = T)) %>%
    # group_by(Site) %>%
    # mutate(MinDate = min(Date),
    #        MaxDate = max(Date),
    #        DaysDep = MaxDate-MinDate,
    #        TopPredictionCounts = sum(Predict_TaxonId == species_list[i]),
    #        CountsPerDay = TopPredictionCounts/as.numeric(DaysDep)) %>%
    group_by(ProbGrp) %>%
    slice(sample(n(), min(50, n()))) %>%
    ungroup() %>%
    select(-ProbGrp) %>%
    distinct()
}

dir_write_path <- "data/for_validation/"
dir_daily_data <- "data/daily_data/"
year <- "OtherYears/"

for(i in 1:nspecies) {
 write.csv(records_to_validate[[i]],
           file = paste0(dir_write_path, year, "to_validate_", names(species_list)[i], ".csv"))
}

for(i in 1:nspecies) {
  write.csv(daily_max_record[[i]],
            file = paste0(dir_daily_data, year, "daily_max_", names(species_list)[i], ".csv"))
}


#### Check validated records ####
files_to_read <- list.files("data/validated_calls")
sp_read_order <- stringr::str_remove(files_to_read, "_validation.csv")

validated_calls <- list()
for(i in 1:nspecies) {
  validated_calls[[sp_read_order[i]]] <- readr::read_csv(paste0("data/validated_calls/", files_to_read[i])) %>%
    mutate(Validated_Grp = as.character(species_list[[sp_read_order[i]]]),
           Validation_Species = sp_read_order[i],
           Validated_Grp_score = !!sym(sp_read_order[i]))
}

validated_calls_combined <- bind_rows(validated_calls) %>%
  rowwise() %>%
  # mutate(Validation_Labels = str_replace_all(Validation_Labels, "Perons", "Peron's"),
  #        Validation_Labels = str_remove_all(Validation_Labels, " [(]race unknown[)]"),
  #        Validated_Grp_score = str_remove_all(Validated_Grp_score, " [(]race unknown[)]")) %>%
  mutate(Detected = stringr::str_detect(Validation_TaxonIDs, Validated_Grp))

validated_calls_combined %>%
  ggplot() +
  geom_density(aes(x = Validated_Grp_score, fill = Detected, colour = Detected), position = "jitter", alpha = 0.25) +
  facet_wrap(~Validation_Species)
