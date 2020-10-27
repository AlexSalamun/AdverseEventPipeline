# Side Effect Pipeline Part 2
#Reddit Development
# This section includes all of the web based resources needed for performing reddit analysis
# It also has all of the text cleaning and other processes run
# The end results will be a DPA for Reddit only data

library(tidyverse)
library(httr)
library(jsonlite)
library(rvest)
library(tidytext)
library(qdapDictionaries)

############################### Download the MedDRA data #############################################

# side effects
se.url <- "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
tmp <- tempfile()
download.file(se.url, tmp)
data <- read.csv(
  gzfile(tmp),
  sep="\t",
  header=F,
  stringsAsFactors=FALSE)

SideEffects <- data
names(SideEffects) <- c("Flat_ID","Stereo_ID", "UMLS ID for Label",
                        "MedDRA Concept Type", "UMLS ID for MedDRA", "Side Effect Name")


# indications
ind.url <- "http://sideeffects.embl.de/media/download/meddra_all_indications.tsv.gz"
tmp <- tempfile()
download.file(ind.url, tmp)
data <- read.csv(
  gzfile(tmp),
  sep="\t",
  header=F,
  stringsAsFactors=FALSE)

Indications <- data
names(Indications) <- c("Flat_ID","UMLS ID for Label", "Method of Detection",
                        "Concept Name", "MedDRA Concept Type", "UMLS ID for MedDRA", 
                        "MedDRA Concept Name")


# Drug names
drugs.url <- "http://sideeffects.embl.de/media/download/drug_names.tsv"
tmp <- tempfile()
download.file(drugs.url, tmp)
data <- read.csv(
  gzfile(tmp),
  sep="\t",
  header=F,
  stringsAsFactors=FALSE)

DrugNames <- data
names(DrugNames) <- c("Flat_ID", "Drug Name")

DrugIndications <- DrugNames %>%
  inner_join(y = Indications, by = "Flat_ID") %>%
  select(c(`Drug Name`, "Name" = `MedDRA Concept Name`)) %>%
  mutate(Association = "Indication") %>%
  unique()

DrugSideEffects <- DrugNames %>%
  inner_join(y = SideEffects, by = "Flat_ID") %>%
  select(c(`Drug Name`, "Name" = `Side Effect Name`)) %>%
  mutate(Association = "Side Effect") %>%
  unique() %>%
  arrange(`Drug Name`)

DrugVocab <- DrugIndications %>%
  union(DrugSideEffects) %>%
  filter(grepl(drug_name, x = `Drug Name`, ignore.case = T)) %>%
  arrange(Name)

AllDrugVocab <- DrugIndications %>%
  union(DrugSideEffects) %>%
  mutate(Name = str_to_lower(Name)) %>%
  arrange(Name)
oe_words <- AllDrugVocab %>%
  filter(grepl(pattern = "oe", x = Name, ignore.case = T))
oe_words$Name <- gsub(pattern = "oe", replacement = "e", x = oe_words$Name)
AllDrugVocab <- rbind(AllDrugVocab, oe_words)

ae_words <- AllDrugVocab %>%
  filter(grepl(pattern = "ae", x = Name, ignore.case = T)) 
ae_words$Name <- gsub(pattern = "ae", replacement = "e", x = ae_words$Name)
AllDrugVocab <- rbind(AllDrugVocab, ae_words)

n_occurs <- data.frame(table(DrugVocab$Name))
DrugVocab$Association[DrugVocab$Name %in% n_occurs$Var1[n_occurs$Freq > 1]] <- "Mix"
DrugVocab <- DrugVocab %>% unique()

DrugVocab$Name <- str_to_lower(DrugVocab$Name)
oe_words <- DrugVocab %>%
  filter(grepl(pattern = "oe", x = Name, ignore.case = T))
oe_words$Name <- gsub(pattern = "oe", replacement = "e", x = oe_words$Name)
DrugVocab <- rbind(DrugVocab, oe_words)

ae_words <- DrugVocab %>%
  filter(grepl(pattern = "ae", x = Name, ignore.case = T)) 
ae_words$Name <- gsub(pattern = "ae", replacement = "e", x = ae_words$Name)
DrugVocab <- rbind(DrugVocab, ae_words)

DrugVocab$Name <- gsub(pattern = "diarrhoea",replacement = "diarrhea", x = DrugVocab$Name)
DrugVocab$Name <- gsub(pattern = "foetal",replacement = "fetal", x = DrugVocab$Name)

################################# Connect to RxNorm from the UMLS to identify all brand names of the drug ###############################
url <- paste0("https://rxnav.nlm.nih.gov/REST/drugs?name=",drug_name)
query <- GET(url = url,
             user_agent(user_agent))
queryContent <- checkQuery(query)
query.df <- fromJSON(queryContent)
if (drug_name == "ASPIRIN") {
  druglist <- query.df$drugGroup$conceptGroup$conceptProperties[[3]]$name
  druglist <- druglist[grepl("325", x = druglist, ignore.case = T)]
  drugnames <- c(unique(str_extract(druglist, pattern = "(?<=\\[).+?(?=\\])")),drug_name)
  drugnames <- gsub(" ", replacement = "%20", x = drugnames)
  drugnames <- str_to_upper(drugnames)
} else {
  druglist <- query.df$drugGroup$conceptGroup$conceptProperties[[3]]$name
  drugnames <- c(unique(str_extract(druglist, pattern = "(?<=\\[).+?(?=\\])")),drug_name)
  drugnames <- gsub(" ", replacement = "%20", x = drugnames)
  drugnames <- str_to_upper(drugnames)
}


################################# Next step is to create the reddit data and add it to this list ########################################
# scrape the list of medical related subreddits
# this list will be looped through as a targeted way to reduce query time
url <- 'https://www.reddit.com/r/ListOfSubreddits/wiki/health'
subreddits <- read_html(url) %>%
  html_nodes(xpath = '//ul') %>%
  html_text() %>%
  tail(1) %>%
  strsplit(split = '\n') %>%
  unlist() %>%
  str_remove("/r/|'") %>%
  c('cancer',
    'MultipleSclerosis',
    'rheumatoid',
    'CrohnsDisease',
    'Asthma',
    'testicularcancer',
    'Parkinsons',
    'Hashimotos',
    'Alzheimers',
    'breastcancer',
    'braincancer',
    'pancreaticcancer',
    'lymphoma',
    'leukemia',
    'kidney',
    'multiplemyeloma',
    'thyroidcancer',
    'lungcancer',
    'skincancer',
    'adverseeffects',
    'braincancer'
  ) %>%
  toupper() %>%
  unique() %>%
  data.frame() %>%
  filter(!grepl("fit|supp|run|lbs|veg|sport|iron|beaut|lift|paleo|nat|juic|form|protein|border|fat|cong",.,ignore.case = TRUE)) %>%
  select(Subreddits = ".")
subreddits <- gsub(pattern = " ", "", x = subreddits$Subreddits)

#############################
# Now loop through every medical subreddit and grab all the posts that mention any of the drug names
post_ids <- character()
products <- character()
thread_poster <- character()
reddit_searches <- character()

for (j in 1:length(drugnames)) {
  drug_name <- drugnames[j]
  for (i in 1:length(subreddits)){
    url <- paste0("https://www.reddit.com/r/",subreddits[i],"/search.json?q=",drug_name,"&restrict_sr=1&limit=100&sort=new")
    query <- GET(url = url, user_agent(user_agent))
    queryContent <- checkQuery(query)
    query.df <- fromJSON(queryContent)
    if(length(query.df$data$children) == 0) {
      next
    }
    if(query.df$data$dist == 100) {
      querySize <- as.numeric(query.df$data$dist)
      iter <- 0
      while (querySize == 100) {
        print(querySize)
        iter = iter + 1
        print(iter)
        url <- paste0("https://www.reddit.com/r/",subreddits[i],"/search.json?q=",drug_name,
                      "&restrict_sr=1&limit=100&sort=new&after=",query.df$data$after)
        query <- GET(url = url, user_agent(user_agent))
        queryContent <- checkQuery(query)
        query.df <- fromJSON(queryContent)
        querySize = as.numeric(query.df$data$dist)
        post_ids <- c(post_ids, query.df$data$children$data$id)
        thread_poster <- c(thread_poster, query.df$data$children$data$author)
        products <- c(products, rep(drug_name, length(query.df$data$children$data$name)))
        reddit_searches <- c(reddit_searches, rep(subreddits[i], length(query.df$data$children$data$name)))
      }
      next
    }
    post_ids <- c(post_ids, query.df$data$children$data$id)
    thread_poster <- c(thread_poster, query.df$data$children$data$author)
    products <- c(products, rep(drug_name, length(query.df$data$children$data$name)))
    reddit_searches <- c(reddit_searches, rep(subreddits[i], length(query.df$data$children$data$name)))
  }
}

threads <- character()
posttitles <- character()

# loop through every post and get the content out of it
# length(post_ids)
for (i in 1:length(post_ids)) {
  comments.url <- paste0("https://www.reddit.com/comments/",post_ids[i],".json")
  comments.query <- GET(url = comments.url, user_agent(user_agent))
  comments.queryContent <- checkQuery(comments.query)
  comments.query.df <- fromJSON(comments.queryContent)
  postTitle <- comments.query.df$data$children[[1]]$data$title
  postbody <- comments.query.df$data$children[[1]]$data$selftext
  comments <- comments.query.df$data$children[[2]]$data$body
  comments.true <- comments[!grepl("I am a bot", comments)]
  comments.responses <- comments.query.df$data$children[[2]]$data$replies
  responses <- comments.responses[lengths(comments.responses)>1]
  if (is.null(responses)) {
    Thread <- paste(postbody)
    posttitles <- c(posttitles, postTitle)
    threads <- c(threads, Thread)
    next
  }
  responseData <- unlist(responses) %>% data.frame() %>% select(responses = ".")
  UserResponseBody <- responseData["data.children.data.body",]
  OtherResponses <- responseData["data.children.data.replies.data.children.data.body",]
  UserThread <- paste(postbody, UserResponseBody, collapse = "")
  OtherThread <- paste(comments.true, OtherResponses, collapse = "")
  Thread <- paste(UserThread, OtherThread, collapse = "")
  threads <- c(threads, Thread)
  posttitles <- c(posttitles, postTitle)
}


# # separate out threads into their own data frame
Reddit.df <- data.frame(Post_ID = post_ids,
                        Product = products,
                        Subreddit = reddit_searches, 
                        Username = thread_poster, 
                        PostTitle = posttitles, 
                        Thread = threads)

# This portion will perform text cleaning on the reddit threads to prepare them for topic modeling
Reddit.df$Thread <- gsub(pattern = "(f|ht)tp(s?)://(.*)[.][a-z]+",replacement = "",x = Reddit.df$Thread)
Reddit.df$Thread <- gsub(pattern = "[[:punct:]]", replacement = "", x = Reddit.df$Thread)
Reddit.df$Thread <- gsub(pattern = "NA", replacement = "", x = Reddit.df$Thread)
Reddit.df$Thread <- gsub(pattern = "\\s+", replacement = " ", x = Reddit.df$Thread)
Reddit.df$Thread <- trimws(Reddit.df$Thread, "both")

# make a dataframe and vector version of the cleaned threads
thread.df <- Reddit.df$Thread %>% data.frame()

# Building a loop to check through different levels of ngrams 1 thru 6
entities <- data.frame()
for (i in 1:6) {
  windows.ngram <- unnest_tokens(tbl = thread.df,
                                 output = word, 
                                 input = ".", 
                                 token ="ngrams", 
                                 n = i, 
                                 drop = F,
                                 collapse = F) %>%
    mutate(word = tolower(word))
  entities <- rbind(entities, windows.ngram)
}
names(entities) <- c("thread", "word")

# create the stopwords
extra_stopwords <- data.frame(word = c("thank", "you", "character0",
                                       "app","get","think","will","name","can","dont",
                                       "list()","hi","i","i'm","pain"))
stopwords <- stop_words %>%
  select(-lexicon) %>%
  rbind(extra_stopwords) %>%
  distinct()

# create a word list of the GradyAugmented english dictionary
`%notin%` <- Negate(`%in%`)

grady <- data.frame(GradyAugmented) %>%
  mutate(Name = as.character(GradyAugmented)) %>%
  mutate(Association = "English") %>%
  subset(select = -GradyAugmented) %>%
  filter(Name %notin% DrugVocab$Name)

# combine the English dictionary with the MedDRA dictionary for the specific drug
english <- DrugVocab %>% 
  select(-`Drug Name`) %>%
  union(grady) %>%
  distinct()

# Remove stopwords and filter to only the domain vocabulary
entity.match <- entities %>%
  anti_join(stopwords) %>%
  inner_join(english, by = c("word" = "Name")) %>%
  #filter(Association %in% c("Side Effect","Mix","Indication","English")) %>%
  filter(nchar(word) >= 3)


# create tables of the different associations
Indications <- entity.match %>%
  filter(Association == "Indication") %>%
  group_by(word) %>%
  tally() %>%
  arrange()

# Only record one unique reaction per thread
# one thread can have multiple associated reactions but each one can only be associated once
Effects <- Reddit.df %>%
  left_join(entity.match, by = c("Thread" = "thread")) %>%
  filter(Association %in% c("Side Effect","Mix")) %>%
  group_by(Post_ID, word) %>%
  tally() %>%
  mutate(n = 1) %>%
  group_by(word) %>%
  tally()


# # build the contingency table for each Side Effect from Reddit
reaction_names <- gsub(" ", replacement = "%20", x = Effects$word)
reaction_names <- str_to_upper(reaction_names)
reaction_counts <- numeric()
endpoint <- "https://api.pushshift.io/reddit/search/submission/?"

subreddit.string <- paste(subreddits, collapse = ",")

for (j in 1:length(reaction_names)) {
  if (j %% 25 == 0) {
    Sys.sleep(60)
  }
  reaction <- reaction_names[j]
  url <- paste0(endpoint,"q=",reaction,"&subreddit=",subreddit.string,"&aggs=subreddit&size=0")
  query <- GET(url = url)
  queryContent <- checkQuery(query)
  query.df <- fromJSON(queryContent)
  reaction_counts <- c(reaction_counts, sum(query.df$aggs$subreddit$doc_count))
}

####################################################
# build D, contingency table value for no drug exposure and no adverse event

url <- paste0(endpoint,"subreddit=",subreddit.string,"&aggs=subreddit&size=0")
query <- GET(url = url, user_agent("alexander.salamun@marquette.edu"))
queryContent <- checkQuery(query)
query.df <- fromJSON(queryContent)
totals <- sum(query.df$aggs$subreddit$doc_count)

reddit.draft <- data.frame(Drug = rep(drug_name, length(reaction_names)),
                           AE = str_to_upper(Effects$word),
                           A = Effects$n,
                           B = reaction_counts - Effects$n,
                           C = rep(length(post_ids), length(reaction_names)) - Effects$n,
                           D = rep(totals, length(reaction_names))- reaction_counts - rep(length(post_ids), length(reaction_names))
)

# now calculate the EBGM values for the Reddit data
r.df <- data.frame(N = as.integer(),
                   E = as.numeric(),
                   RR = as.numeric(),
                   PRR = as.numeric(),
                   EBGM = as.numeric(),
                   EB05 = as.numeric(),
                   EB95 = as.numeric())

for (i in 1:length(reddit.draft$Drug)) {
  r.df <- rbind(r.df, EBGM(reddit.draft$A[i],reddit.draft$B[i],reddit.draft$C[i],reddit.draft$D[i]))
}

options(digits = 3, scipen = 999)
reddit.dpa <- data.frame(Drug = reddit.draft$Drug, 
                         AdverseEvent = reddit.draft$AE,
                         r.df) %>% arrange(desc(N)) %>%
  mutate("Signal Detected" = ifelse(EB05 >= 2.000,"*","")) %>%
  #filter(EB05 >= 2.000) %>%
  arrange(desc(EB05))
reddit.dpa <- reddit.dpa[,c(1:4,7:8,10)]

