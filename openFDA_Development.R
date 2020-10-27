# Side Effect Pipeline Part 3
# openFDA Development
# This section includes all of the api calls that need to be made to the openFDA
# it will also calculate the DPA values and aggregate with the pre-existing RedditData


############################### Connect to openFDA ##########################################

# creating the api call for the overall adverse event profile of Aspirin
# Create [A] Fields in the contingency table
url <- paste0("https://api.fda.gov/drug/event.json?",
              "api_key=",api_key,
              "&search=patient.drug.openfda.generic_name.exact:\"",drug_name,"\"",
              "&count=patient.reaction.reactionmeddrapt.exact&limit=1000")
query <- GET(url = url)
queryContent <- checkQuery(query)
query.df <- fromJSON(queryContent)
results <- query.df$results


############################### Connect to openFDA ##########################################

## specific drug name [B] in the table]
url <- paste0("https://api.fda.gov/drug/event.json?",
              "api_key=",api_key,
              "&search=patient.drug.openfda.generic_name.exact:\"",drug_name,"\"",
              "&limit=1")
query <- GET(url = url)
queryContent <- checkQuery(query)
query.df <- fromJSON(queryContent)
B_A <- rep(query.df$meta$results$total,length(results$term))
B <- B_A - results$count
results <-  cbind(results, B)


# For the [C] value I'll need to create a loop to run through all the different adverse events
# this is the rate limiting step because the openFDA only allows 240 queries per minute so it'll take up to five minutes

C_A <- as.numeric()
for (i in 1:nrow(results)) {
  if (i %% 200 == 0) {
    Sys.sleep(65)
  }
  reaction <- gsub(" ", "+", results$term[239])
  reaction <- gsub("\\^", "%27", reaction)
  url <- paste0("https://api.fda.gov/drug/event.json?",
                "api_key=",api_key,
                "&search=patient.reaction.reactionmeddrapt.exact:\"",reaction,"\"",
                "&limit=1")
  query <- GET(url = url)
  queryContent <- checkQuery(query)
  query.df <-  fromJSON(queryContent)
  C_A[i] <- as.numeric(query.df$meta$results$total)
}
C <- C_A - results$count
results <-  cbind(results, C)


## all events in the table [D in the table]
url <- paste0("https://api.fda.gov/drug/event.json?",
              "api_key=",api_key,
              "&search=_exists_:(patient.reaction.reactionmeddrapt.exact)&limit=1")
query <- GET(url = url)
queryContent <- checkQuery(query)
query.df <- fromJSON(queryContent)
D_A <- rep(query.df$meta$results$total,length(results$term))
D <- D_A - B_A - C_A
results <-  cbind(results, D)


openFDA.draft <- results %>%
  mutate(Drug = rep(drug_name, length(results$term))) %>%
  select(c(Drug, "AE" = term, "A" = count, B,C, D))

##########################################################3

o.df <- data.frame(N = as.integer(),
                   E = as.numeric(),
                   RR = as.numeric(),
                   PRR = as.numeric(),
                   EBGM = as.numeric(),
                   EB05 = as.numeric(),
                   EB95 = as.numeric())



for (i in 1:length(openFDA.draft$AE)) {
  o.df <- rbind(o.df, EBGM(openFDA.draft$A[i],openFDA.draft$B[i],openFDA.draft$C[i],openFDA.draft$D[i]))
}

openFDA.dpa <- data.frame(Drug = openFDA.draft$Drug,
                          AdverseEvent = openFDA.draft$AE,
                          o.df) %>% arrange(desc(N)) %>%
  mutate("Signal Detected" = ifelse(EB05 >= 2.000,"*","")) %>%
  #filter(EB05 >= 2.000) %>%
  arrange(desc(EB05))

openFDA.dpa <- openFDA.dpa[,c(1:4, 7:8, 10)]


####################################################3
# Parts 2 and 3 have now created the Reddit and openFDA only datasets
# The next step is to merge these together and calculate the combined signal

Total.draft <- openFDA.draft %>%
  merge(y = reddit.draft, by = c("Drug","AE"), all = T)

Total.draft[is.na(Total.draft)] <- 0
Total.draft <- Total.draft %>%
  mutate(Total.A = A.x + A.y,
         Total.B = B.x + B.y,
         Total.C = C.x + C.y,
         Total.D = D.x + D.y)

t.df <- data.frame(N = as.integer(),
                   E = as.numeric(),
                   RR = as.numeric(),
                   PRR = as.numeric(),
                   EBGM = as.numeric(),
                   EB05 = as.numeric(),
                   EB95 = as.numeric())

for (i in 1:length(Total.draft$Drug)) {
  t.df <- rbind(t.df, EBGM(Total.draft$Total.A[i],
                           Total.draft$Total.B[i],
                           Total.draft$Total.C[i],
                           Total.draft$Total.D[i]))
}

Total.dpa <- data.frame(Drug = Total.draft$Drug,
                        AdverseEvent = Total.draft$AE,
                        t.df) %>% arrange(desc(N)) %>%
  mutate("Signal Detected" = ifelse(EB05 >= 2.000,"*","")) %>%
  #filter(EB05 >= 2.000) %>%
  arrange(desc(EB05))

Total.dpa <- Total.dpa[,c(1:4, 7:8, 10)]






