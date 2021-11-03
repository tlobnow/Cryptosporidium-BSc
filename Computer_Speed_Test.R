# How fast is your Computer

install.packages("benchmarkme")
library(benchmarkme)

res <- benchmark_std(runs = 3)
plot(res)

upload_results(res)