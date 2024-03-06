cnvLogs <- read.table("BRCA_snp_segmented_scna_minus_germline_cnv_hg19.txt", header = T, fill = T)
dim(cnvLogs)
head(cnvLogs)
summary(cnvLogs)

#for one patient
nrow(subset(cnvLogs, Sample == "TCGA-3C-AAAU-10A-01D-A41E-01"))
nrow(subset(cnvLogs, Sample == "TCGA-3C-AAAU-10A-01D-A41E-01" & Segment_Mean > 0))
nrow(subset(cnvLogs, Sample == "TCGA-3C-AAAU-10A-01D-A41E-01" & Segment_Mean < 0))
mean(subset(cnvLogs, Sample == "TCGA-3C-AAAU-10A-01D-A41E-01")[["Segment_Mean"]])