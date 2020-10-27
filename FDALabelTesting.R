# Side Effect Pipeline Part 4
# This section will evaluate the accuracy of the results obtained from the web based sources
# These results are compared to the official drug labels
# It requires a manual input screening words that are typically used as indications across the designated drugs.
# The screens for Tigecycline, Lisinopril, and Aspirin are already created as examples

#install.packages("tm")
library(pdftools)
library(caret)

# Insert the chosen labels
pdf.labels <- c("Tigecycline.pdf","Lisinopril.pdf","Aspirin.pdf")

# Read in the pdf as a character vector
clean.label <- pdf_text(pdf = pdf.labels[grepl(drug_name,pdf.labels, ignore.case = T)])

# screen out the indications
lisinopril.screen <- "lisinopril|hypertension|myocardial|heart"
aspirin.screen <- "aspirin|pain|headache|tooth|creamps|menstrual|premenstrual|arthritis"
tigecycline.screen <- "tigecycline|pneumonia|infection|bacter|diabetic|resist|kleb|staph"

TotalDrugVocab <- DrugIndications %>%
  union(DrugSideEffects)
TotalDrugVocab <- strsplit(x = TotalDrugVocab$Name, split = " ")
TotalDrugVocab <- data.frame("low_AE" = str_to_lower(unlist(TotalDrugVocab))) %>% unique()
  
screens <- c(lisinopril.screen, aspirin.screen, tigecycline.screen)
current.screen <- screens[grepl(drug_name, x = screens, ignore.case = T)]

reddit.dpa1 <- reddit.dpa %>%
  filter(!grepl(current.screen, AdverseEvent, ignore.case = T))

openFDA.dpa1 <- openFDA.dpa %>%
  filter(!grepl(current.screen, AdverseEvent, ignore.case = T)) %>%
  mutate(AdverseEvent = gsub("oe", replacement = "E", x = AdverseEvent, ignore.case = T)) %>%
  mutate(AdverseEvent = gsub("ae", replacement = "E", x = AdverseEvent,ignore.case = T))

Total.dpa1 <- Total.dpa %>%
  filter(!grepl(current.screen, AdverseEvent, ignore.case = T)) %>%
  mutate(AdverseEvent = gsub("oe", replacement = "E", x = AdverseEvent,ignore.case = T)) %>%
  mutate(AdverseEvent = gsub("ae", replacement = "E", x = AdverseEvent,ignore.case = T))

# clean and parse the label for comparison to results
clean.label <- gsub(pattern = "(f|ht)tp(s?)://(.*)[.][a-z]+",replacement = "",x = clean.label)
clean.label <- gsub(pattern = "[[:punct:]]", replacement = " ", x = clean.label)
clean.label <- gsub(pattern = "NA", replacement = "", x = clean.label)
clean.label <- gsub("(?=[0-9])[a-z0-9.-]+","",clean.label, ignore.case = T, perl = T)
clean.label <- gsub(pattern = "\\s+", replacement = " ", x = clean.label)
clean.label <- trimws(clean.label, "both")

clean.df <- clean.label %>% data.frame()

label.entities <- data.frame()
for (i in 1:6) {
  windows.ngram <- unnest_tokens(tbl = clean.df,
                                 output = word, 
                                 input = ".", 
                                 token ="ngrams", 
                                 n = i, 
                                 drop = T,
                                 collapse = F) %>%
    mutate(word = tolower(word))
  label.entities <- rbind(label.entities, windows.ngram)
}
names(label.entities) <- c("AdverseEvent")
label.entities$AdverseEvent <- str_to_upper(label.entities$AdverseEvent)

CurrentLabel<- label.entities %>%
  mutate(low_AE = str_to_lower(AdverseEvent)) %>%
  inner_join(AllDrugVocab, by = c("low_AE" = "Name")) %>%
  filter(grepl(drug_name, `Drug Name`, ignore.case = T),
         Association == "Side Effect") %>%
  select(AdverseEvent) %>%
  unique() %>%
  mutate(Source = "DrugPackage")

# Create accuracy tables
# openFDA
openFDA.accuracy <- openFDA.dpa1 %>%
  mutate(AdverseEvent = str_to_upper(AdverseEvent)) %>%
  full_join(CurrentLabel) %>%
  mutate(signal = as.factor(ifelse(`Signal Detected` == "*","TRUE","FALSE")),
         labelmatch = as.factor(case_when(Source == "DrugPackage" ~ "TRUE",
                                          TRUE ~ "FALSE"))) %>%
  filter(N >= B_A[1]*0.001)

# Reddit
reddit.accuracy <- reddit.dpa1 %>%
  mutate(AdverseEvent = str_to_upper(AdverseEvent)) %>%
  full_join(CurrentLabel) %>%
  mutate(signal = as.factor(ifelse(`Signal Detected` == "*","TRUE","FALSE")),
         labelmatch = as.factor(case_when(Source == "DrugPackage" ~ "TRUE",
                                          TRUE ~ "FALSE"))) %>%
  filter(N >= length(post_ids)*0.001)

# Combined
Total.accuracy <- Total.dpa1 %>%
  mutate(AdverseEvent = str_to_upper(AdverseEvent)) %>%
  full_join(CurrentLabel) %>%
  mutate(signal = as.factor(ifelse(`Signal Detected` == "*","TRUE","FALSE")),
         labelmatch = as.factor(case_when(Source == "DrugPackage" ~ "TRUE",
                                          TRUE ~ "FALSE"))) %>%
  filter(N >= (B_A[1]+length(post_ids))*0.001)


# use a confusion matrix to get the sensitivity and specificity of each
openFDA.cm <- confusionMatrix(openFDA.accuracy$signal, openFDA.accuracy$labelmatch)
reddit.cm <- confusionMatrix(reddit.accuracy$signal, reddit.accuracy$labelmatch)
Total.cm <- confusionMatrix(Total.accuracy$signal, Total.accuracy$labelmatch)

openFDA.cm
reddit.cm
Total.cm

## The final testing step is to use the dpa measures to create a logistic regression model
# A successful item is one that was matched in the label
# Based on this it can determine the appropriate cutoff for determining signal detection
# the goal is to raise specificity as high as possible, so the cutoff will be lowered until it is maximized

# Logit Reddit
logit.reg <- glm(data = reddit.accuracy, formula = labelmatch ~ N + EB05, family = "binomial")
summary(logit.reg)

glm.probs <- predict(logit.reg, type = "response")
logit.df <- data.frame(actual = reddit.accuracy$labelmatch, predicted = glm.probs)
logit.cm <- confusionMatrix(factor(ifelse(glm.probs > .9,"TRUE","FALSE")), reference = factor(reddit.accuracy$labelmatch))
Spec <- logit.cm$byClass[2]
cutoff <- 0.9
while (Spec < 0.9) {
  cutoff = cutoff - 0.01
  cm.iter <- confusionMatrix(factor(ifelse(glm.probs > cutoff,"TRUE","FALSE")), reference = factor(reddit.accuracy$labelmatch))
  Spec <- cm.iter$byClass[2]
  print(cutoff)
}
if (Spec == 1) {
  cutoff = cutoff + 0.01
}
reddit.cutoff <- cutoff
reddit.LR.cm <- confusionMatrix(factor(ifelse(glm.probs > cutoff,"TRUE","FALSE")), reference = factor(reddit.accuracy$labelmatch))
reddit.LR.cm 

reddit.accuracy <- reddit.accuracy %>%
  mutate(P_Logit = 1/(1+exp(-(logit.reg$coefficients[1]+logit.reg$coefficients[2]*reddit.accuracy$EBGM+logit.reg$coefficients[3]*reddit.accuracy$EB05))))

# #compute with logistic regression cutoff values
reddit.LR <- reddit.accuracy %>%
  filter(P_Logit > reddit.cutoff)

# Logit openFDA
logit.reg <- glm(data = openFDA.accuracy, formula = labelmatch ~ EBGM + EB05, family = "binomial")
summary(logit.reg)

glm.probs <- predict(logit.reg, type = "response")
logit.df <- data.frame(actual = openFDA.accuracy$labelmatch, predicted = glm.probs)
logit.cm <- confusionMatrix(factor(ifelse(glm.probs > .9,"TRUE","FALSE")), reference = factor(openFDA.accuracy$labelmatch))
Spec <- logit.cm$byClass[2]
cutoff <- 0.9
while (Spec < 0.9) {
  cutoff = cutoff - 0.01
  cm.iter <- confusionMatrix(factor(ifelse(glm.probs > cutoff,"TRUE","FALSE")), reference = factor(openFDA.accuracy$labelmatch))
  Spec <- cm.iter$byClass[2]
  print(cutoff)
}
if (Spec == 1) {
  cutoff = cutoff + 0.01
}
openFDA.cutoff <- cutoff
openFDA.LR.cm <- confusionMatrix(factor(ifelse(glm.probs > cutoff,"TRUE","FALSE")), reference = factor(openFDA.accuracy$labelmatch))
openFDA.LR.cm 

openFDA.accuracy <- openFDA.accuracy %>%
  mutate(P_Logit = 1/(1+exp(-(logit.reg$coefficients[1]+logit.reg$coefficients[2]*openFDA.accuracy$EBGM+logit.reg$coefficients[3]*openFDA.accuracy$EB05))))

#compute with logistic regression cutoff values
openFDA.LR <- openFDA.accuracy %>%
  filter(P_Logit >= openFDA.cutoff)
# openFDA.LR.cm <- confusionMatrix(openFDA.LR$signal, openFDA.LR$labelmatch)
# min(openFDA.LR$EBGM)

# Logit Total
#Total.accuracy$labelmatch <- as.numeric(as.logical(Total.accuracy$labelmatch))
logit.reg <- glm(data = Total.accuracy, formula = labelmatch ~ EBGM + EB05, family = "binomial")
summary(logit.reg)

glm.probs <- predict(logit.reg, type = "response")
logit.df <- data.frame(actual = Total.accuracy$labelmatch, predicted = glm.probs)
logit.cm <- confusionMatrix(factor(ifelse(glm.probs > .9,"TRUE","FALSE")), reference = factor(Total.accuracy$labelmatch))
Spec <- logit.cm$byClass[2]
cutoff <- 0.9
while (Spec < 0.9) {
  cutoff = cutoff - 0.01
  cm.iter <- confusionMatrix(factor(ifelse(glm.probs > cutoff,"TRUE","FALSE")), reference = factor(Total.accuracy$labelmatch))
  Spec <- cm.iter$byClass[2]
}
if (Spec == 1) {
  cutoff = cutoff + 0.01
}
Total.cutoff <- cutoff
Total.LR.cm <- confusionMatrix(factor(ifelse(glm.probs > cutoff,"TRUE","FALSE")), reference = factor(Total.accuracy$labelmatch))
Total.LR.cm 

Total.accuracy <- Total.accuracy %>%
  mutate(P_Logit = 1/(1+exp(-(logit.reg$coefficients[1]+logit.reg$coefficients[2]*Total.accuracy$EBGM+logit.reg$coefficients[3]*Total.accuracy$EB05)))) %>%
  arrange(desc(P_Logit))

# #compute with logistic regression cutoff values
Total.LR <- Total.accuracy %>%
  filter(P_Logit >= Total.cutoff)
# Total.LR.cm <- confusionMatrix(Total.LR$signal, Total.LR$labelmatch)
# min(Total.LR$EBGM)

