setwd("/Users/ryotaro/R勉強会/rates")

#install.packages("stats")
#install.packages("termstrc")
#library(termstrc)
library(stats)

# termstrc のアーカイブ URL
#url <- "https://cran.r-project.org/src/contrib/Archive/termstrc/termstrc_1.3.2.tar.gz"

# インストール
#install.packages(url, repos = NULL, type = "source")

#install.packages("YieldCurve")
library(YieldCurve)
library(xts)
library(zoo)

Data <- read.table("kinri.csv", header=T, sep= ",")
Data.Rate <- as.matrix(Data[, -1])
#Data.Date <- as.Date(Data[, 1], format = "%Y.%m.%d")
Data.Date <- as.Date(Data[, 1], format = "%Y/%m/%d")

Term <- c(1:10, 15, 20, 25, 30, 40)

Sp.term <- seq(0.5, 40, length=80)
WK1 <- matrix(0, nrow=nrow(Data.Rate),ncol=length(Sp.term))
colnames(WK1) <- c(paste(Sp.term,"Y"))

for (i in 1:nrow(Data.Rate)) {
  sp <- smooth.spline(Term, Data.Rate[i,])
  WK1[i,] <- round(predict(sp, Sp.term)[[2]], 4)
}

WK2 <- matrix(0, nrow=nrow(WK1), ncol=ncol(WK1))
for (k in 1:nrow(WK2)){
  Cr <- WK1[k,]/100
  vt <- matrix(0, nrow=1, ncol=ncol(WK1))
  for (i in 1:length(Cr)){
    vt[1,i]<-(1-sum(Cr[i]*vt)/2)/(1+Cr[i]/2)
    WK2[k,i]<-round((1/(vt[1,i]^(1/Sp.term[i]))-1)*100, 6)
  }
}

#ネルソンシーゲル
yield <- WK2[,Term*2]
yield_xts <- xts(yield, order.by=Data.Date) #金利データをxts型データに変換
ns_res <- Nelson.Siegel(rate=yield_xts, maturity=Term)
ns_sp <- NSrates(ns_res, 1:30)

write.zoo(ns_res, file = "ns_parm.csv", sep = ",")
write.zoo(ns_sp, file = "ns_res.csv", sep = ",")


library(ggplot2)
#install.packages("reshape2")
library(reshape2)


# xts → data.frame（Date列付きに変換）
df_param <- data.frame(Date = index(ns_res), coredata(ns_res))

# 縦持ちに変換（ggplotで使いやすく）
df_param_long <- melt(df_param, id.vars = "Date")

# プロット
ggplot(df_param_long, aes(x = Date, y = value, color = variable)) +
  geom_line() +
  labs(title = "Nelson-Siegel Parameters Over Time", y = "Parameter Value") +
  theme_minimal()

#イールドカーブプロット
dates <- as.Date(c("2007-11-06", "2007-12-03", "2007-12-14"))

# xtsデータのインデックス（Date型）
all_dates <- index(ns_sp)
# 年を文字列で抽出
years <- as.character(2007:2015)
# 各年の最初の日付を1つずつ取得
#library(dplyr)
#library(lubridate)
dates_selected <- sapply(years, function(y) {
  min(all_dates[format(all_dates, "%Y") == y])
})
# 結果はPOSIXct型で返るのでDate型に変換（必要なら）
dates_selected <- as.Date(dates_selected)

curve_df <- data.frame(
  Maturity = as.numeric(gsub("X", "", colnames(ns_sp)))
)

for (d in dates_selected) {
  d_date <- as.Date(d)
  curve_df[[as.character(d_date)]] <- as.numeric(ns_sp[d_date, ])
}

matplot(curve_df$Maturity, curve_df[, -1], type = "l", lty = 1, col = 1:length(dates),
        xlab = "Maturity (Years)", ylab = "Yield (%)", main = "Yield Curves Comparison")
legend("bottomright", legend = as.character(dates_selected), col = 1:length(dates), lty = 1,bty = "n")


Data <- read.table("ns_res.csv", header = T, sep = ",")
input <- Data[, -1]

#対数変化率のデータに変換
Delt <- log((1/(1 + input[-1, ]/100))/(1/(1 + input[-nrow(input), ]/100)))
colnames(Delt) <- colnames(Data[2:ncol(Data)])

#各年限の平均・偏差を出力
mData <- sapply(Delt, mean)
sData <- sapply(Delt, sd)
