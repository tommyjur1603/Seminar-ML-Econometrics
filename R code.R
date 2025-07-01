
# R code
load("nsw_dw.rda")
data <- nsw_dw
data


url <- "https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/master/datasets/IHDP/csv/ihdp_npci_1.csv"
ihdp <- read.csv(url, header = FALSE)
colnames(ihdp)[1:5] <- c("treatment", "y_factual", "y_cfactual", "mu0", "mu1")
colnames(ihdp)

colnames(ihdp)[6:length(colnames(ihdp))] <- paste0("x", seq_len(ncol(ihdp) - 5))
colnames(ihdp)
