# Analysing some data

library(readr)
library(HistDAWass)
library(tidyverse)
# genres_v2 <- read_csv("genres_v2.csv")
# genres_v2$genre<-factor(genres_v2$genre)

tibble2MATH<-function(data, #a dataframe
                      factorv=1, # the factor variable used for the new row names
                      num_vars=c(1:3), # numerical variables to translate into distributions
                      qua=100 #a number of quantiles in wich divide the distribution
                      ){
  
  np=c(0:qua)/qua
  
  data[[factorv]] <- factor(data[[factorv]])
  names_of_units <- as.character(levels(data[[factorv]]))
  names_of_vars <- colnames(data)[num_vars]

  data2<-data %>% select_at(c(factorv,num_vars))
  MyMath<-MatH(nrows=length(names_of_units),ncols=length(num_vars),
               rownames = names_of_units, varnames=names_of_vars)
  for (i in 1:length(names_of_units)){
    
    for (j in 1:length(names_of_vars)){
    #  browser()
      tmp<-data2 %>% filter(.[[1]]==names_of_units[i]) %>% select_at(j+1)
      
      tmpx<-unname(quantile(unlist(tmp),p=np))
      tmpd <-distributionH(x=tmpx,p=np)
      MyMath@M[i,j][[1]]<-tmpd
    }
  }
  return(MyMath)
}

my_math<-tibble2MATH(genres_v2,factorv = 19,qua=20)

######## spotify variables
# Danceability: Describes how suitable a track is for dancing based on a combination of musical elements including tempo, rhythm stability, beat strength, and overall regularity.
# Valence: Describes the musical positiveness conveyed by a track. Tracks with high valence sound more positive (e.g. happy, cheerful, euphoric), while tracks with low valence sound more negative (e.g. sad, depressed, angry).
# Energy: Represents a perceptual measure of intensity and activity. Typically, energetic tracks feel fast, loud, and noisy. For example, death metal has high energy, while a Bach prelude scores low on the scale.
# Tempo: The overall estimated tempo of a track in beats per minute (BPM). In musical terminology, tempo is the speed or pace of a given piece, and derives directly from the average beat duration.
# Loudness: The overall loudness of a track in decibels (dB). Loudness values are averaged across the entire track and are useful for comparing relative loudness of tracks.
# Speechiness: This detects the presence of spoken words in a track. The more exclusively speech-like the recording (e.g. talk show, audio book, poetry), the closer to 1.0 the attribute value.
# Instrumentalness: Predicts whether a track contains no vocals. “Ooh” and “aah” sounds are treated as instrumental in this context. Rap or spoken word tracks are clearly “vocal”.
# Liveness: Detects the presence of an audience in the recording. Higher liveness values represent an increased probability that the track was performed live.
# Acousticness: A confidence measure from 0.0 to 1.0 of whether the track is acoustic.
# Key: The estimated overall key of the track. Integers map to pitches using standard Pitch Class notation . E.g. 0 = C, 1 = C♯/D♭, 2 = D, and so on.
# Mode: Indicates the modality (major or minor) of a track, the type of scale from which its melodic content is derived. Major is represented by 1 and minor is 0.
# Duration: The duration of the track in milliseconds.
# Time Signature: An estimated overall time signature of a track. The time signature (meter) is a notational convention to specify how many beats are in each bar (or measure).
