# data source: http://www.dt.fee.unicamp.br/~tiago/smsspamcollection/

sms <- read.csv("sms_spam.csv", stringsAsFactors = FALSE)
sms$type <- factor(sms$type)

install.packages("tm")
library(tm)

#data cleaning

corpus <- Corpus(VectorSource(sms$text))
corpus_clean <- tm_map(corpus, content_transformer(tolower))
corpus_clean <- tm_map(corpus_clean, removeNumbers)
corpus_clean <- tm_map(corpus_clean, removeWords, stopwords())
corpus_clean <- tm_map(corpus_clean, removePunctuation)
corpus_clean <- tm_map(corpus_clean, stripWhitespace)

dtm <- DocumentTermMatrix(corpus_clean)
dtm

#train and test datasets

sms_train <- sms[1:4169, ]
sms_test  <- sms[4170:5559, ]

dtm_train <- dtm[1:4169, ]
dtm_test  <- dtm[4170:5559, ]

sms_corpus_train <- corpus_clean[1:4169]
sms_corpus_test  <- corpus_clean[4170:5559]

#dataViz wordcloud

install.packages("wordcloud")
library(wordcloud)

wordcloud(sms_corpus_train, min.freq = 20, random.order = FALSE)
wordcloud(spam$text, max.words = 40, scale = c(3, 0.5))
wordcloud(ham$text, max.words = 40, scale = c(3, 0.5))

