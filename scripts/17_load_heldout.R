# Load the held-out sample (SteingroverInPrep) for the confirmatory test of the ORL variant.
# Same steps as 01_load_data.R, different study. Check structure only before fitting.
load("data-raw/IGTdata.rdata")
table(index_100[, "Study"])             # study labels; confirms the exact spelling
is_ho <- index_100[, "Study"] == "SteingroverInPrep"
sum(is_ho)  # should be 70

choice_ho <- choice_100[is_ho, ]
wi_ho <- wi_100[is_ho, ]
lo_ho <- lo_100[is_ho, ]
net_ho <- wi_ho + lo_ho
dim(choice_ho)  # should be 70 100

# Payoff check: each deck's distinct net outcomes (should match Horstmann's, scheme 2)
for (d in 1:4) print(sort(unique(net_ho[choice_ho == d])))




