# Load the IGT data pool into R
# Data: Steingroever et al. (2015), journal of Open Psychology Data, https://osf.io/8t7rm

load("data-raw/IGTdata.rdata")
# Select Horstmann participants
is_horst <- index_100[, "Study"] == "Horstmann"
sum(is_horst)

# Keep Horstmann rows only

choice_h<-choice_100[is_horst, ]
wi_h<-wi_100[is_horst, ]
lo_h<-lo_100[is_horst, ]
dim(choice_h)
dim(wi_h)
dim(lo_h)

# net outcome per trial = win + loss (losses stored as -#'s)
# model learns from net outcome (reward)

net_h <- wi_h + lo_h
dim(net_h)
range(net_h)











