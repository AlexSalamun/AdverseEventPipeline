# Side Effect Pipeline Part 5
# This section will compare the resultant drug safety profiles to that of the current label.

library(tidyr)

CurrentLabel<- label.entities %>%
  mutate(low_AE = str_to_lower(AdverseEvent)) %>%
  inner_join(AllDrugVocab, by = c("low_AE" = "Name")) %>%
  filter(grepl(drug_name, `Drug Name`, ignore.case = T),
         Association == "Side Effect") %>%
  select(AdverseEvent) %>%
  unique() %>%
  mutate(Source = "DrugPackage")

RedditLabel <- reddit.dpa %>%
  filter(`Signal Detected` == "*") %>%
  select(AdverseEvent) %>%
  mutate(Source = "Reddit")

openFDA.Label <- openFDA.dpa %>%
  filter(`Signal Detected` == "*") %>%
  select(AdverseEvent) %>%
  mutate(Source = "openFDA")

Combined.Label <- Total.dpa %>%
  filter(`Signal Detected` == "*") %>%
  select(AdverseEvent) %>%
  mutate(Source = "Combined")

# this generates the tables that are seen in the excel files  for full drug safety profiles
FullLabel <- CurrentLabel %>%
  full_join(y = RedditLabel) %>%
  full_join(y = openFDA.Label) %>%
  full_join(y = Combined.Label) %>%
  mutate(value = 1) %>%
  spread(key = Source, value) %>%
  mutate(Drug = drug_name) %>%
  select(c(Drug, AdverseEvent, DrugPackage, Reddit, openFDA, Combined)) %>%
  replace_na(list(Reddit = 0, openFDA = 0, Combined = 0, DrugPackage = 0)) %>%
  mutate(Commonality = Reddit + openFDA + Combined + DrugPackage) %>%
  arrange(desc(Commonality))
