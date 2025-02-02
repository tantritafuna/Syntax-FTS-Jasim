library(dplyr)
library(lubridate)

#Input Data
data <- read.csv(file.choose(), header = TRUE)
ihk.ts <- ts(data$IHK, start = c(2017, 1), frequency = 12)

#Membagi Data Training dan Testing
n_total <- length(ihk.ts);n_total
n_train <- floor(0.95 * n_total);n_train
n_test <- n_total - n_train;n_test
train_data <- window(ihk.ts, end = c(2017, n_train))
test_data <- window(ihk.ts, start = c(2017, n_train + 1))

#Mengubah Data Menjadi Numerik
data <- as.numeric(train_data)

#Menentukan Himpunan Semesta Pembicaraan
D_min <- min(data);D_min
D_max <- max(data);D_max
D1 <- abs(floor(D_min) - D_min);D1
D2 <- ceiling(D_max) - D_max;D2
U <- c(D_min - D1, D_max + D2);U

#Menghitung Panjang Interval
n <- length(data);n
differences <- abs(diff(data));differences
av <- sum(differences)/(n-1);av
B <- av/2;B
get_base <- function(B) {
  if (B > 0.1 && B <= 1) return(0.1)
  if (B > 1 && B <= 10) return(1)
  if (B > 10 && B <= 100) return(10)
  if (B > 100 && B <= 1000) return(100)
  if (B > 1000) return(1000)
  return(NULL)
}
base <- get_base(B);base
I <- round(B, 1);I

#Menentukan Banyak Interval
m <- ceiling((D_max + D1 - D_min + D2) / I);m 
##Apabila Menggunakan Metode Sturges
m <- ceiling(1 + 3.322 * log10(n));m
I <- 3.88
##Mencoba Beberapa Jumlah Interval
m <- 20
I <- 1.55
m <- 40
I <- 0.775

#Membentuk Himpunan Fuzzy
fuzzy_sets <- lapply(1:m, function(i) {
  c(U[1] + (i - 1) * I, U[1] + i * I)
})
names(fuzzy_sets) <- paste0("A", 1:m)
intervals <- do.call(rbind, fuzzy_sets)
intervals

#Fuzzifikasi Data 
fuzzify <- function(value, fuzzy_sets) {
  for (i in seq_along(fuzzy_sets)) {
    if (value >= fuzzy_sets[[i]][1] && value <= fuzzy_sets[[i]][2]) {
      return(i)
    }
  }
  return(NA)
}
data_fuzzy <- sapply(data, fuzzify, fuzzy_sets)
data_fuzzy

#Membentuk FLR
FLR <- cbind(
  1:length(data_fuzzy),
  data_fuzzy,
  c(data_fuzzy[-1], NA)
)
FLR

#Membentuk FLRG
# 6. Create FLRG as a list
FLRG <- lapply(1:m, function(i) {
  next_states <- FLR[FLR[,2] == i, 3]
  next_states[!is.na(next_states)]
})
FLRG

#Peramalan Menggunakan FTS Jasim
##Menghitung Y
calculate_Y <- function(Dt, Dt_1, Dt_2) {
  return((Dt - Dt_1) - (Dt_1 - Dt_2))
}
##Aturan 1
rule_1 <- function(Dt, Dt_1, interval_Aj, I) {
  x <- abs(Dt - Dt_1) / 2
  Aj_mid <- (interval_Aj[2] - interval_Aj[1]) / 2
  
  if(x > Aj_mid) {
    # Cenderung menaik
    Ft <- 0.75 * I + interval_Aj[1]
  } else if(x == Aj_mid) {
    # Nilai tengah
    Ft <- (interval_Aj[1] + interval_Aj[2]) / 2
  } else {
    # Cenderung menurun
    Ft <- 0.25 * I + interval_Aj[1]
  }
  return(Ft)
}
##Aturan 2
rule_2 <- function(Y, Dt_1, interval_Aj, I) {
  # Kondisi 1: x = |Y| × 2 + Dt-1 ??? Aj atau x = Dt-1 - |Y| × 2 ??? Aj
  x1 <- abs(Y) * 2 + Dt_1
  x2 <- Dt_1 - abs(Y) * 2
  
  if(x1 >= interval_Aj[1] && x1 <= interval_Aj[2] || 
     x2 >= interval_Aj[1] && x2 <= interval_Aj[2]) {
    # Cenderung menaik
    Ft <- 0.75 * I + interval_Aj[1]
  } else {
    # Kondisi 2: x = |Y|/2 + Dt-1 ??? Aj atau x = Dt-1 - |Y|/2 ??? Aj
    x3 <- abs(Y)/2 + Dt_1
    x4 <- Dt_1 - abs(Y)/2
    
    if(x3 >= interval_Aj[1] && x3 <= interval_Aj[2] || 
       x4 >= interval_Aj[1] && x4 <= interval_Aj[2]) {
      # Cenderung menurun
      Ft <- 0.25 * I + interval_Aj[1]
    } else {
      # Jika bukan keduanya
      Ft <- (interval_Aj[1] + interval_Aj[2]) / 2
    }
  }
  return(Ft)
}
##Aturan 3
rule_3 <- function(Y, Dt_1, interval_Aj, I) {
  # Kondisi 1: x = |Y|/2 + Dt-1 ??? Aj atau x = Dt-1 - |Y|/2 ??? Aj
  x1 <- abs(Y)/2 + Dt_1
  x2 <- Dt_1 - abs(Y)/2
  
  if(x1 >= interval_Aj[1] && x1 <= interval_Aj[2] || 
     x2 >= interval_Aj[1] && x2 <= interval_Aj[2]) {
    Ft <- 0.25 * I + interval_Aj[1]
  } else {
    # Kondisi 2: x = |Y| × 2 + Dt-1 ??? Aj atau x = Dt-1 - |Y| × 2 ??? Aj
    x3 <- abs(Y) * 2 + Dt_1
    x4 <- Dt_1 - abs(Y) * 2
    
    if(x3 >= interval_Aj[1] && x3 <= interval_Aj[2] || 
       x4 >= interval_Aj[1] && x4 <= interval_Aj[2]) {
      Ft <- 0.75 * I + interval_Aj[1]
    } else {
      Ft <- (interval_Aj[1] + interval_Aj[2]) / 2
    }
  }
  return(Ft)
}
##Peramalan berdasarkan FLRG
forecast_FLRG <- function(FLRG, i, j, data, intervals, I) {
  # Cek nilai NA pada input
  if (any(is.na(c(i, j, data[length(data)]))) || length(data) < 3) {
    return((intervals[i, 1] + intervals[i, 2]) / 2)
  }
  
  if (length(FLRG) == 0) {
    #FLRG kosong
    return((intervals[i, 1] + intervals[i, 2]) / 2)
  } else if (length(unique(FLRG)) == 1) {
    #One to one
    Y <- calculate_Y(data[length(data)], data[length(data) - 1], data[length(data) - 2])
    # Tangani nilai NA pada Y
    if (is.na(Y)) {
      return((intervals[i, 1] + intervals[i, 2]) / 2)
    }
    if ((j > i && Y > 0) || (j < i && Y > 0) || (j == i && Y > 0)) {
      return(rule_2(Y, data[length(data) - 1], intervals[j,], I))
    } else if ((j > i && Y < 0) || (j < i && Y < 0) || (j == i && Y < 0)) {
      return(rule_3(Y, data[length(data) - 1], intervals[j,], I))
    }
  } else {
    # One to many
    unique_states <- unique(FLRG)
    if (max(diff(sort(unique_states))) <= 2) {
      # Jika perbedaan <= 2
      interval_means <- sapply(unique_states, function(idx) {
        (intervals[idx, 1] + intervals[idx, 2]) / 2
      })
      return(mean(interval_means))
    } else {
      # Jika perbedaan > 2
      return(mean(sapply(unique_states, function(idx) intervals[idx, 2])))
    }
  }
}

#Peramalan Data Training
train_forecast <- sapply(1:length(data), function(t) {
  if (t == 1) return(NA)
  current_idx <- data_fuzzy[t]
  if (is.na(current_idx)) return(NA)
  next_idx <- ifelse(!is.na(FLR[t, 3]), FLR[t, 3], current_idx)
  if (is.na(next_idx)) next_idx <- current_idx
  forecast_FLRG(FLRG[[current_idx]], current_idx, next_idx, data[1:t], intervals, I)
})

#Peramalan Data Testing
test_data_numeric <- as.numeric(test_data)
test_forecast <- numeric(length(test_data_numeric))
for(i in 1:length(test_data_numeric)) {
  # Fuzzify the test value
  test_value <- test_data_numeric[i]
  current_idx <- fuzzify(test_value, fuzzy_sets)
  if(is.na(current_idx)) {
    test_forecast[i] <- mean(c(intervals[1,1], intervals[m,2]))
    next
  }
  if(length(FLRG[[current_idx]]) > 0) {
    next_idx <- FLRG[[current_idx]][1]
  } else {
    next_idx <- current_idx
  }
  test_forecast[i] <- forecast_FLRG(FLRG[[current_idx]], current_idx, next_idx, 
                                    c(data, test_data_numeric[1:i]), intervals, I)
}

#Menghitung MAPE
calculate_mape <- function(actual, forecast) {
  return(mean(abs((actual - forecast) / actual)) * 100)
}
##MAPE Data Training
train_mape <- calculate_mape(data[-1], train_forecast[-1])
cat("Mean Absolute Percentage Error (MAPE):", round(train_mape, 4), "%\n")
##MAPE Data Testing
test_mape <- calculate_mape(test_data_numeric, test_forecast)
cat("Mean Absolute Percentage Error (MAPE):", round(test_mape, 4), "%\n")

#Hasil Peramalan
train_results <- data.frame(
  Set = "Training",
  TimePoint = format(seq.Date(as.Date("2017-01-01"), by = "month", length.out = length(data)), "%b %Y"),
  Actual = data,
  Forecast = train_forecast,
  Error = data - train_forecast,
  APE = abs((data - train_forecast) / data) * 100
)
test_results <- data.frame(
  Set = "Testing",
  TimePoint = format(seq.Date(as.Date("2017-01-01") + months(length(data)), by = "month", length.out = length(test_data_numeric)), "%b %Y"),
  Actual = test_data_numeric,
  Forecast = test_forecast,
  Error = test_data_numeric - test_forecast,
  APE = abs((test_data_numeric - test_forecast) / test_data_numeric) * 100
)
results <- rbind(train_results, test_results);results

#Plot Data Peramalan
plot(1:nrow(results), results$Actual, type = "l", col = "blue",
     xlab = "Bulan", ylab = "IHK", xaxt = "n")  # Matikan sumbu X default
lines(1:nrow(results), results$Forecast, col = "red", lty = 2)
abline(v = length(data), col = "gray", lty = 3)
legend("bottomleft", 
       legend = c("Aktual", "Peramalan", "Train/Test Split"),
       col = c("blue", "red", "gray"), 
       lty = c(1, 2, 3))
axis(1, at = 1:nrow(results), labels = results$TimePoint, las = 2, cex.axis = 0.53)

#Menghitung Peramalan Periode Selanjutnya
last_value <- tail(ihk.ts, 1)
last_idx <- fuzzify(last_value, fuzzy_sets)
if (!is.na(last_idx)) {
  if (length(FLRG[[last_idx]]) > 0) {
    next_idx <- FLRG[[last_idx]][1]
  } else {
    next_idx <- last_idx
  }
  next_period_forecast <- forecast_FLRG(FLRG[[last_idx]], last_idx, next_idx, 
                                        as.numeric(ihk.ts), intervals, I)
} else {
  next_period_forecast <- mean(c(intervals[1,1], intervals[m,2]))
}
cat("\nForecast for the next period:", round(next_period_forecast, 2), "\n")

#Plot Peramalan Periode Selanjutnya
final_results <- data.frame(
  Period = c(time(ihk.ts), max(time(ihk.ts)) + 1/12),
  Value = c(as.numeric(ihk.ts), next_period_forecast),
  Type = c(rep("Actual", length(ihk.ts)), "Forecast")
)
plot(final_results$Period[final_results$Type == "Actual"], 
     final_results$Value[final_results$Type == "Actual"],
     type = "l", col = "blue",
     xlab = "Tahun", ylab = "IHK")
points(tail(final_results$Period[final_results$Type == "Actual"], 1),
       tail(final_results$Value[final_results$Type == "Actual"], 1),
       col = "blue", pch = 19)
points(tail(final_results$Period, 1),
       tail(final_results$Value, 1),
       col = "red", pch = 19)
lines(tail(final_results$Period, 2),
      tail(final_results$Value, 2),
      col = "red", lty = 2)
legend("topright", 
       legend = c("Data Historis", "Nilai Aktual Terakhir", "Peramalan Periode Berikutnya"),
       col = c("blue", "blue", "red"),
       pch = c(NA, 19, 19),
       lty = c(1, NA, 2))
