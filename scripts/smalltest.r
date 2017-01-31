x <- matrix(1:12, 4)
x_unique <- matrix(1:4, 4)

dfx <- data.frame(x, row.names = LETTERS[1:nrow(x)])
colnames(dfx) <- LETTERS[(nrow(x) +1):(nrow(x) + ncol(x))]

dfx_unique <- data.frame(x_unique, row.names = LETTERS[1:nrow(x_unique)])
colnames(dfx_unique) <- "MEANS"

print(dfx_unique)
print(dfx)

rowMeans(dfx)

dfx_unique[,"MEANS"] <- rowMeans(dfx)
print(dfx_unique)