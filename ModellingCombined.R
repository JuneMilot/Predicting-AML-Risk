library(tidyverse)
library(glmnet)
library(caret)
library(pROC)
set.seed(55)

cat("Dataset: Combined Dataset\n")

target <- readRDS("CleanedDatasets/combined_all_genes.rds")

# set features
x <- as.matrix(target %>% select(-patient, -risk))
y <- target$risk

# split data
train_idx <- createDataPartition(y, p = 0.75, list = FALSE)
x_train <- x[train_idx, ]
y_train <- y[train_idx]
x_test  <- x[-train_idx, ]
y_test  <- y[-train_idx]

# balance classes
y_train <- as.factor(y_train)
train_data <- as.data.frame(x_train)
train_data$risk <- y_train
train_data <- train_data %>% select(where(~ !is.list(.)))
train_balanced <- upSample(x = train_data %>% select(-risk), y = train_data$risk)
x_train <- as.matrix(train_balanced %>% select(-Class))
y_train <- train_balanced$Class

# normalize using training data stats
preProc <- preProcess(x_train, method = c("center", "scale"))
x_train <- predict(preProc, x_train)
x_test  <- predict(preProc, x_test)

# fit logistic with lasso
fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, standardize = FALSE)

# 1 standard error rule
lambda <- fit$lambda.1se

# predict probs
probs <- predict(fit, newx = x_test, s = lambda, type = "response")
preds <- ifelse(probs > 0.5, 1, 0)

# print accuracy
accuracy <- mean(preds == y_test)
print(round(accuracy, 4))

# create ROC curve and print AUC
roc <- roc(y_test, as.numeric(probs))
print(auc(roc))
plot(roc, main = "ROC Curve - Test Set")

# predict probs for training set
train_probs <- predict(fit, newx = x_train, s = lambda, type = "response")
train_preds <- ifelse(train_probs > 0.5, 1, 0)

# create ROC curve and print AUC for training set
train_roc <- roc(y_train, as.numeric(train_probs))
print(auc(train_roc))
plot(train_roc, main = "ROC Curve - Training Set")

# get coefs with 1se rule that are not zero
coefs <- coef(fit, s = lambda)
coef_df <- subset(data.frame(gene=rownames(coefs),
                             coefficient=as.numeric(coefs)),
                  coefficient != 0)
coef_df <- coef_df[order(-abs(coef_df$coefficient)),]

# print coefs
print(nrow(coef_df))
print(coef_df)

# create confusion matrix
cm <- confusionMatrix(factor(preds, levels = c(0, 1)),
                      factor(y_test, levels = c(0, 1)),
                      positive = "1")
print(cm)

