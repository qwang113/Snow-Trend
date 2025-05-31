dat_yisu <- readRDS(here::here("snow_cleaned.Rda"))
colnames(dat_yisu)

dates <- as.Date(sub("^X", "", colnames(dat_yisu)[-c(1,2)]), format = "%m.%d.%Y")
colnames(dat_yisu)[-c(1,2)] <- as.character(dates)

full_dat <- readRDS(here::here("snow_cleaned_full.Rda"))

full_selected <- full_dat[,colnames(dat_yisu)[263:ncol(dat_yisu)]]
yisu_selected <- dat_yisu[,263:ncol(dat_yisu)]
plot(as.matrix(full_selected)[1,])
lines(as.matrix(yisu_selected)[1,])

     