library(data.table)
library(dplyr)
library(progress)
library(pbapply)
library(parallel)
setwd("~/Primo Semestre AA 2024-2025/Nonparametric Statistics/Homework")
rm(list = ls())
graphics.off()

data <- readRDS("C:/Users/giovi/OneDrive/Documenti/Primo Semestre AA 2024-2025/Nonparametric Statistics/Homework/data.homework.rds")

#Preprocessing
#For each observation, "highlighting" the week of the year and the year
data <- as.data.table(data)
data[, year := year(Date)]
data[, week := week(Date)]

#Selecting the variables for my exercise
dataset <- data[,.(Date,PM10,Precip,year,week)]

#I remove all the rows with at least one NA
colSums(is.na(dataset))
sum(apply(dataset, 1, function(row) any(is.na(row))))
dataset_clean <- na.omit(dataset)
any(is.na(dataset_clean)) 

#For the variables of interest, I select the maximum 
#over each week of each year and added the associated new variables
dataset_grouped <- dataset_clean[, .(
  max_PM10 = max(PM10),
  max_Precip = max(Precip)
), by = .(year, week)]

#I pick the two variables (one covariate and one response) that I will use
final_dataset <- dataset_grouped[, .(max_PM10, max_Precip)]

#Piece-wise Polynomial Model

#Creating the variables for the model: cutoff variable and interaction between the precipitation and the cutoff
cutoff <- quantile(final_dataset$max_Precip, 0.5)
final_dataset$max_Precip_cut <- final_dataset$max_Precip>cutoff
final_dataset$max_Precip_cut_model <- (final_dataset$max_Precip - cutoff)*final_dataset$max_Precip_cut

#Fitting the model
model_cut_quad <- lm(max_PM10 ~ poly(max_Precip,degree = 2) + max_Precip_cut_model + I(max_Precip>cutoff),  data=final_dataset)

#I'm looking at the residuals because after 
#I will have to implement bootstrap
plot(model_cut_quad)

#Summary of the model
summary(model_cut_quad)

#Re-naming the coefficient for a better visualization:
#-Creation of an appropriate map for name changing, 
#-creation of an object easier to handle (tibble, kind of a data.frame) but still containing the summary information
#-changing the names
newnames_map <- c(
  "(Intercept)" = "Intercept",
  "poly(max_Precip, degree = 2)1" = "linear term",
  "poly(max_Precip, degree = 2)2" = "quadratic term",
  "max_Precip_cut_model" = "linear interaction",
  "I(max_Precip > cutoff)TRUE" = "x > cutoff"
)
renamed_coefficients <- broom::tidy(summary(model_cut_quad)) %>%
  dplyr::mutate(term = dplyr::recode(term, !!!newnames_map))
renamed_coefficients

#Plot of the predicted values of the model: 
#-At first the data.frame contains new observations equally spaced in the range of our observations
#-Then I add the other variables (based on the first one)
new_data <-
  with(final_dataset, data.frame(
    max_Precip = seq(range(max_Precip)[1], range(max_Precip)[2], by = 0.01)
  ))
new_data$max_Precip_cut_model = (new_data$max_Precip - cutoff) * (new_data$max_Precip > cutoff)
new_data$max_Precip_cut = new_data$max_Precip > cutoff

#Fitted values of the new observations and their standard errors
preds=predict(model_cut_quad,new_data,se=T)
se.bands=cbind(preds$fit +2* preds$se.fit ,preds$fit -2* preds$se.fit)

#Plot of the new predicted observations, with the standard error bands
with(final_dataset, plot(max_Precip ,max_PM10 ,xlim=range(new_data$max_Precip) ,cex =.5, col =" darkgrey " ))
lines(new_data$max_Precip,preds$fit ,lwd =2, col =" blue")
matlines(new_data$max_Precip, se.bands ,lwd =1, col =" blue",lty =3)

#Bootstrapping 
n <- dim(final_dataset)[1]

original_model <- lm(max_PM10 ~ poly(max_Precip,degree = 2) + max_Precip_cut_model + I(max_Precip>cutoff),  data=final_dataset)
original_coeff <- as.data.frame(coef(original_model))
rownames(original_coeff) <- c("Intercept", "linear term", "quadratic term", "cutoff", "interaction(linear)")
original_coeff

#Setting for parallel computing
n_cores <- detectCores()
n_cores
cl = parallel::makeCluster(min(c(n_cores, 4)))

#Bootstrap iteration on the observations and not on the residuals
bootstrap_iteration <- function(final_dataset) {
  boot <- sample(1:n, replace = TRUE)
  boot_data <- final_dataset[boot,]
  boot_model <- lm(max_PM10 ~ poly(max_Precip, degree = 2) + max_Precip_cut_model + I(max_Precip > cutoff), data = boot_data)
  coef(boot_model)

}

clusterExport(cl, varlist = c("final_dataset","bootstrap_iteration","cutoff","n"))
set.seed(1401)
B <- 1e04 
coeff_matrix <- pbreplicate(B, bootstrap_iteration(final_dataset), cl = cl)
stopCluster(cl)

#I manipulate a little the result of the parallel computation 
dim(coeff_matrix)
coeff_matrix <- t(coeff_matrix)
coeff_df <- as.data.frame(coeff_matrix)
colnames(coeff_df) <- c("Intercept", "linear term", "quadratic term", "cutoff", "interaction(linear)")
head(coeff_df)

#Inverse quantiles confidence intervals

#I create a function that will help me later to collect all the quantiles together
extract_quantiles <- function(column,alpha){
  lower <- quantile(column, alpha/2)
  upper <- quantile(column, 1 - alpha/2)
  return(list(lower = lower, upper = upper))
}
alpha = 0.05

#Using sapply(), I apply extract_quantiles() to each column of the object
#containing the bootstrapped coefficients of the model
quantiles <- sapply(coeff_df, extract_quantiles, alpha = alpha)
quantiles <- t(quantiles)
quantiles <- as.data.frame(quantiles)

original_coeff_n <- as.numeric(original_coeff[, 1])
lower_q <- sapply(quantiles$lower, as.numeric)
upper_q <- sapply(quantiles$upper, as.numeric)

intervals <- data.frame(
  lower.bound = original_coeff_n - (upper_q - original_coeff_n) ,
  estimated = original_coeff_n,
  upper.bound = original_coeff_n - (lower_q - original_coeff_n),
  row.names = rownames(original_coeff)
)
intervals

# bootstrap_iteration <- function(original_model, final_dataset) {
#   
#   boot_resid <- sample(residuals(original_model), replace = TRUE)
#   boot_y <- fitted(original_model) + boot_resid
#   boot_data <- final_dataset
#   boot_data$max_PM10 <- boot_y
#   boot_model <- lm(max_PM10 ~ poly(max_Precip, degree = 2) + max_Precip_cut_model + I(max_Precip > cutoff), data = boot_data)
#   coef(boot_model)
#   
# }

