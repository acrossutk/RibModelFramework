t <- read.table("S.cerevisiae.PNAS.elongation.rates.csv", sep=",", header=T)
ref <- t
for (aa in unique(t[,1])) {
  ref[(t[,1] == aa),3] <- t[(t[,1] == aa),3] - t[max(which(t[,1] == aa)),3]
}

ref <- ref[-which(ref[,3] == 0), ]
write.table(ref, file = "S.cerevisiae.PNAS.elongation.rates.ref.csv", sep=",", row.names = F, col.names = F, quote = F)
