  library(tidyverse)
  library(glmnet)
  library(caret)
  library(pROC)
  set.seed(42)
  files <- c(
    "target_cleaned",
    "gse6891_1_cleaned",
    "gse6891_2_cleaned",
    "gse37642_1_cleaned",
    "gse37642_2_cleaned",
    "tcga_cleaned",
    "BEAT_BMA_cleaned",
    "BEAT_PB_cleaned"
  )
  
  #dataset l  oop
  for (name in files) {
    
    cat("Dataset:", name, "\n")
    
    target <- readRDS(paste0("CleanedDatasets/", name, ".rds"))
    
    #set features
    x <- as.matrix(target %>% select(-patient, -risk))
    y <- target$risk
    
    #split data
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
    
    
    #fit logistic with lasso
    cvfit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
    
    #1 standard error rule
    lambda <- cvfit$lambda.1se
    
    #predict probs
    probs <- predict(cvfit, newx = x_test, s = lambda, type = "response")
    preds <- ifelse(probs > 0.5, 1, 0)
    #print accuracy
    accuracy <- mean(preds == y_test)
    print(round(accuracy, 4))
    
    
    #create roc curve and print auc (area under curve)
    roc <- roc(y_test, as.numeric(probs))
    print(auc(roc))
    plot(roc, main = paste("ROC Curve -", name))
    
    # Predict probabilities on training set
    train_probs <- predict(cvfit, newx = x_train, s = lambda, type = "response")
    train_preds <- ifelse(train_probs > 0.5, 1, 0)
    
    # Create ROC curve and print AUC for training set
    train_roc <- roc(y_train, as.numeric(train_probs))
    print(paste("Training AUC for", name, ":", auc(train_roc)))
    plot(train_roc, main = paste("ROC Curve - Training Set -", name))
    
    
    #gets coef with 1se rule that are not zero
    coefs <- coef(cvfit, s = lambda)
    coef_df <- subset(data.frame(gene=rownames(coefs), coefficient=as.numeric(coefs)), coefficient!=0)
    #order coef
    coef_df <- coef_df[order(-abs(coef_df$coefficient)),]
  
    # print coefs
    print(nrow(coef_df))
    print(coef_df)
    
    #create confusion matrix with caret and print
    cm <- confusionMatrix(factor(preds, levels = c(0, 1)),
                          factor(y_test, levels = c(0, 1)),
                          positive = "1")
    print(cm)
    
  }
