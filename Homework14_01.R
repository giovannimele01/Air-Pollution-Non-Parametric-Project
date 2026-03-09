library(data.table)
library(progress)
library(parallel)
setwd("~/Primo Semestre AA 2024-2025/Nonparametric Statistics/Homework")
rm(list = ls())
graphics.off()

#----PRE-PROCESSING----
data <- readRDS("C:/Users/giovi/OneDrive/Documenti/Primo Semestre AA 2024-2025/Nonparametric Statistics/Homework/data.homework.rds")
#write.csv(data, "data_homework.csv", row.names = FALSE)

data <- as.data.table(data)
# Estrazione dell'anno
data[, year := year(Date)]

# Estrazione del numero della settimana
data[, week := week(Date)]

dataset <- as.data.frame(data)
data <- as.data.frame(data)#NO, FAI TUTTO IN DATA TABLE 
#E POI IN FONDO PASSA A DAT FRAME

dataset <- data[,c(1,6,7,13,14)]
colSums(is.na(dataset))
sum(apply(dataset, 1, function(row) any(is.na(row))))

dataset_clean <- na.omit(dataset)

any(is.na(dataset_clean))

setDT(dataset_clean)

dataset_grouped <- dataset_clean[, .(
  max_PM10 = max(PM10, na.rm = TRUE),
  max_Precip = max(Precip, na.rm = TRUE)
), by = .(year, week)]

final_dataset <- dataset_grouped[, .(max_PM10, max_Precip)]
#fwrite(final_dataset, "dataset.csv")
#final_dataset <- as.data.frame(final_dataset)
#----PIECE-WISE POLYNOMIAL REGRESSION----
cutoff <- quantile(final_dataset$max_Precip, 0.5)
#cutoff <- as.numeric(cutoff)
#or
#node <- median(dataset_grouped$max_Precip)

final_dataset$max_Precip_cut <- final_dataset$max_Precip>cutoff

final_dataset$max_Precip_cut_model <- (final_dataset$max_Precip - cutoff)*final_dataset$max_Precip_cut

model_cut_quad <- lm(max_PM10 ~ poly(max_Precip,degree = 2) + max_Precip_cut_model + I(max_Precip>cutoff),  data=final_dataset)

new_data <-
  with(final_dataset, data.frame(
    max_Precip = seq(range(max_Precip)[1], range(max_Precip)[2], by = 0.5)
  ))
new_data$max_Precip_cut_model = (new_data$max_Precip - cutoff) * (new_data$max_Precip > cutoff)
new_data$max_Precip_cut = new_data$max_Precip > cutoff

preds=predict(model_cut_quad,new_data,se=T)
se.bands=cbind(preds$fit +2* preds$se.fit ,preds$fit -2* preds$se.fit)

with(final_dataset, plot(max_Precip ,max_PM10 ,xlim=range(new_data$max_Precip) ,cex =.5, col =" darkgrey " ))
lines(new_data$max_Precip,preds$fit ,lwd =2, col =" blue")
matlines(new_data$max_Precip, se.bands ,lwd =1, col =" blue",lty =3)
summary(model_cut_quad)
abline(v=cutoff)


#seconda prova
final_dataset$max_Precip_cut_model_quadratic <- (final_dataset$max_Precip - cutoff)^2*final_dataset$max_Precip_cut

model_cut_quad_2 <- lm(max_PM10 ~ poly(max_Precip,degree = 2) + max_Precip_cut_model + max_Precip_cut_model_quadratic + I(max_Precip>cutoff),  data=final_dataset)
summary(model_cut_quad_2)
new_data2 <-
  with(final_dataset, data.frame(
    max_Precip = seq(range(max_Precip)[1], range(max_Precip)[2], by = 0.5)
  ))
new_data2$max_Precip_cut_model = (new_data2$max_Precip - cutoff) * (new_data2$max_Precip > cutoff)
new_data2$max_Precip_cut_model_quadratic = (new_data2$max_Precip - cutoff)^2 * (new_data2$max_Precip > cutoff)
new_data2$max_Precip_cut = new_data2$max_Precip > cutoff

preds_2=predict(model_cut_quad_2,new_data2,se=T)
se.bands_2=cbind(preds_2$fit +2* preds_2$se.fit ,preds_2$fit -2* preds_2$se.fit)

with(final_dataset, plot(max_Precip ,max_PM10 ,xlim=range(new_data2$max_Precip) ,cex =.5, col =" darkgrey " ))
lines(new_data2$max_Precip,preds_2$fit ,lwd =2, col =" blue")
matlines(new_data2$max_Precip, se.bands_2 ,lwd =1, col =" blue",lty =3)

all.equal(preds, preds_2) #false --> corretto

#terzo modello
model_cut_quad_2_alt <- lm(max_PM10 ~ poly(max_Precip,degree = 2) * max_Precip_cut,  data=final_dataset)
preds_2_alt=predict(model_cut_quad_2_alt,new_data2,se=T)
se.bands_2_alt=cbind(preds_2_alt$fit +2* preds_2_alt$se.fit ,preds_2_alt$fit -2* preds_2_alt$se.fit)
all.equal(preds_2, preds_2_alt)

summary(model_cut_quad_2_alt)

with(final_dataset, plot(max_Precip ,max_PM10 ,xlim=range(new_data2$max_Precip) ,cex =.5, col =" darkgrey " ))
lines(new_data2$max_Precip,preds_2_alt$fit ,lwd =2, col =" blue")
matlines(new_data2$max_Precip, se.bands_2_alt ,lwd =1, col =" blue",lty =3)

#----BOOTSTRAPPING----

#original_model <- lm(max_PM10 ~ poly(max_Precip,degree = 2) * max_Precip_cut,  data=final_dataset)
original_model <- lm(max_PM10 ~ poly(max_Precip,degree = 2) + max_Precip_cut_model + I(max_Precip>cutoff),  data=final_dataset)
original_coeff <- coef(original_model)

set.seed(1401) 
B <- 1e04

coeff_matrix <- matrix(NA, nrow = B, ncol = length(coef(original_model)))

for (i in 1:B) {
  boot_resid <- sample(residuals(original_model), replace = TRUE)
  
  boot_y <- fitted(original_model) + boot_resid
  
  boot_data <- final_dataset
  boot_data$max_PM10 <- boot_y
  boot_model <- lm(max_PM10 ~ poly(max_Precip,degree = 2) + max_Precip_cut_model + I(max_Precip>cutoff),  data=boot_data)
  
  coeff_matrix[i, ] <- coef(boot_model)
}

coeff_df <- as.data.frame(coeff_matrix)
#colnames(coeff_df) <- names(coef(original_model))

#colnames(coeff_df) <- c("Intercept", "linear term", "quadratic term", "cutoff", "interaction(linear)", "interaction(quadratic)")
colnames(coeff_df) <- c("Intercept", "linear term", "quadratic term", "cutoff", "interaction(linear)")
coeff_df

#reverse percentile intervals

n_cores <- detectCores() / 2  # Use half of the available cores
cl <- makeCluster(n_cores)

clusterExport(cl, varlist = c("coeff_df", "original_coeff"))  

compute_reverse_percentile <- function(j) {
  boot_dist <- coeff_df[[j]]
  
  original_value <- original_coeff[j]
  
  lower_bound <- original_value - (quantile(boot_dist, 0.975) - original_value)
  upper_bound <- original_value - (quantile(boot_dist, 0.025) - original_value)
  
  c(Lower = lower_bound, Upper = upper_bound)
}

results <- parLapply(cl, 1:ncol(coeff_df), compute_reverse_percentile)

reverse_percentile_intervals <- data.frame(
  Coefficient = colnames(coeff_df),
  do.call(rbind, results)
)

stopCluster(cl)

reverse_percentile_intervals




