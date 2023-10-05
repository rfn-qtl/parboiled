#####################################
# PC and VT analytic pipeline - Function
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: May 9 2023
#####################################

VT_PC <- function(data){
  
data <- data[,1:(ncol(data)-2)]
# study overall Name
study <- unique(data$studyDescription)

data$locationName <- gsub(paste0(unique(data$studyYear) %% 2000, "_", study, "_"), "", data$studyName, fixed = T)
trials <- unique(data$locationName)
traits <- colnames(data)[31:ncol(data)]

# coerce col to factors or numeric
data[,1:30] <- lapply(data[,1:30], factor)
data[,31:ncol(data)] <- lapply(data[,31:ncol(data)], function(x){as.numeric(x)})
data$row <- as.numeric(data$rowNumber)
data$col <- as.numeric(data$colNumber)

# reorganize the data
data <- reshape2::melt(data[,c(1:30,    which(colnames(data) %in% traits), 
                                        which(colnames(data) == "row"), 
                                        which(colnames(data) == "col"))], 
                                        measure.vars = traits)

#################################
# fitting a model for each trial
#################################
require(foreach)
require(doParallel)
require(doMC)
library(SpATS)
library(car)
library(ggplot2)

# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1) # type the number of cores you want to use
getDoParWorkers()

#running all the single-trials in parallel 
# preparing a grid to run the analysis in parallel
grid <- expand.grid(trials, traits)

system.time(
  results.st <- foreach(i = 1:nrow(grid), 
                        .packages = c("SpATS", "car"), 
                        .combine = "rbind",
                        .export = c("SpATS", "predict.SpATS", "getHeritability", "outlierTest"),
                        .multicombine = TRUE, 
                        .errorhandling = "remove",
                        .verbose = TRUE    
 ) %dopar% {
      # subset the data  
    sample <- droplevels.data.frame(data[data$locationName == grid[i,1] & data$variable == grid[i,2] ,])

    # outlier detection and elimination
    fit <- lm(value ~ replicate + rowNumber + colNumber + germplasmName, data = sample)
    outlier <- names(outlierTest(fit)$p)
    sample[outlier, "value"] <- NA
      
      nrow <- max(sample$row)
      ncol <- max(sample$col)
      nseg.row <- nrow
      nseg.col <- ncol
      
      fitF <- SpATS(response = "value", 
                                fixed = ~ 1, 
                                random = ~ replicate + rowNumber + colNumber, 
                                spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                                genotype = "germplasmName", 
                                genotype.as.random = FALSE, 
                                data = sample)
      
      # obtaining the spatial trends - from raw data to BLUES and BLUPS
      #plot.SpATS(fitF) 
    
      # Estimate BLUEs
      blues <- predict.SpATS(fitF, which = "germplasmName")
      blues <- blues[,c(1,7,8)]
      colnames(blues)[2:3] <- c("BLUE", "sep_BLUE")
      
      # Now, run as random
      fitR <- SpATS(response = "value", 
                    fixed = ~ 1, 
                    random = ~ replicate + rowNumber + colNumber, 
                    spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                    genotype = "germplasmName", 
                    genotype.as.random = TRUE, 
                    data = sample)
      
      # obtaining the spatial trends - from raw data to BLUES and BLUPS
      #plot.SpATS(fitR)       
      
      # to obtain the heritability via the package function we can use 
      h2g <- getHeritability(fitR)
      # Broad-sense heritability based on Cullis method
      Vg <- fitR$var.comp["germplasmName"]
      ng <- length(unique(sample$germplasmName))
      C11_g <- fitR$vcov$C11_inv
      trC11_g <-sum(diag(C11_g))
      av2 <- 2/ng * (trC11_g - (sum(C11_g) - trC11_g) / ng-1) # mean var of a difference between genotypic BLUPS
      H2.Cullis <- 1 - av2 / (2 * Vg)
      # Estimate BLUPs for Grain Yield
      blups <- predict.SpATS(fitR, which = "germplasmName")
      blups <- blups[,c(7,8)]
      colnames(blups)[1:2] <- c("BLUP", "sep_BLUP")
      # Reliability
      rel <- mean(1 - blups$sep_BLUP^2 / fitR$var.comp["germplasmName"])
      # weights for ID's - adjust residual for further analysis
      vcov.mme <- fitR$vcov$C11_inv
      w <- diag(vcov.mme)
      
      output <- data.frame(blues,
                           w = w,
                           blups,
                           Location = as.character(unique(sample$locationName)), 
                           h2g = h2g,
                           r = rel,
                           H.cullis = H2.Cullis,
                           trait = unique(sample$variable)
                          )
        }
)

output <- results.st
# saving the output file for single-trials analysis
write.csv(output, paste0("output/", study, "_", "output_single-trials.csv"))

# and a summary of the single trials
sum.st <- unique.data.frame(output[,7:11])
write.csv(sum.st, paste0("output/", study, "_", "summary_single-trials.csv"))

##################################################
# second step - joint analysis MET
##################################################
library(sommer)
library(statgenGxE)
library(metan)

out.MET <- data.frame()

for(i in 1:length(traits)){

  sample <- output[output$trait == traits[i],] 
  
# Fitting genotype by environment model - joint analysis

  # Fitting genotype by environment model - joint analysis
  fitMET <- mmer(BLUE ~ 1,
                 random= ~ Location + germplasmName + Location:germplasmName,
                 #weights = w,
                 rcov= ~ units,
                 data = sample, 
                 tolParInv = 1e-01,
                 verbose = FALSE)
  
  # Broad-sense heritability
  h2g.MET <- vpredict(fitMET, h2 ~ V2 / ( V2 + V3)) # MET level

# predicting the BLUP - main effect
BLUPs.main <- data.frame(germplasmName = gsub("germplasmName", "", names(fitMET$U$germplasmName$BLUE), fixed = T, ),
                         predicted.value = fitMET$U$germplasmName$BLUE + fitMET$Beta$Estimate)

# predicting the BLUP per environment
BLUPs.MET <- sample[,c(1,2,7)]
colnames(BLUPs.MET)[2] <- "predicted.value"

rice.handbook.table <- round(reshape2::acast(BLUPs.MET, germplasmName ~  Location, value.var = "predicted.value"))
write.csv(rice.handbook.table, paste0("output/", study, "_", traits[i], "_rice.handbook.table.csv"))

# barplot graph with confidence interval using main
data.plot <- BLUPs.main

limits <- aes(ymax = data.plot$predicted.value + sqrt(as.numeric(fitMET$sigma$units))*1.96,
              ymin = data.plot$predicted.value - sqrt(as.numeric(fitMET$sigma$units))*1.96)

p <- ggplot(data = data.plot, aes(x = reorder(germplasmName, -predicted.value), y = predicted.value, 
)) + 
  geom_bar(position = position_dodge(), stat="identity", fill = "purple") +
  scale_fill_brewer(palette="Set1") +
  geom_errorbar(limits, position = position_dodge(),
                width = 0.5) +
  labs(x = "Variety", y = traits[i]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  annotate("text", x=length(unique(data.plot$germplasmName))/2, 
           y=max(data.plot$predicted.value), 
           label = paste("Reliability = ", round(h2g.MET*100), "%")) +
  coord_cartesian(ylim=c(min(data.plot$predicted.value - sqrt(as.numeric(fitMET$sigma$units))*1.96),
                         max(data.plot$predicted.value + sqrt(as.numeric(fitMET$sigma$units))*1.96)))

ggsave(filename = paste0("output/", study, "_", traits[i], '_overall_performances.tiff'),
       plot = p,
       device = 'tiff',
       width = 280,
       height = 140,
       units = 'mm',
       dpi = 300)

data.plot$std.error <- sqrt(as.numeric(fitMET$sigma$units))*1.96
colnames(data.plot)[1:3] <- c("Variety", traits[i], "Confidence_interval")

# saving the output file for single-trials analysis
write.csv(data.plot, paste0("output/", study, "_", traits[i], "_", "overall_performances.csv"))

out.MET <- rbind(out.MET, data.frame(
  Study = study,
  Date = date(),
  Number_of_lines = length(unique(sample$germplasmName)),
  Trait = unique(sample$trait),
  Number_of_locations = length(unique(sample$Location)),
  H2g = as.numeric(h2g.MET[1])
  ))

#######################################
# third step - GGE MET analysis
#######################################

data.MET <- BLUPs.MET
colnames(data.MET)[c(1,2)] <- c("Variety", traits[i]) 

## Create a TD object
dropsTD <- statgenSTA::createTD(data = data.MET, genotype = "Variety", trial = "Location")
## Fit a model where trials are nested within scenarios.
dropsVarComp <- gxeVarComp(TD = dropsTD, trait = traits[i])
# Finlay-Wilkinson Analysis
dropsFW <- gxeFw(TD = dropsTD, trait = traits[i])

# reorganizing the data
test <- matrix(dropsFW$TD[[1]]$beta) %*% matrix(sort(dropsFW$envEffs[,2]), nrow = 1) + dropsFW$estimates[match(dropsFW$TD[[1]]$genotype, dropsFW$estimates$Genotype) ,2]
locations <- c(trials)
locations <- as.character(locations[match(sort(dropsFW$envEffs[,2]), dropsFW$envEffs[,2])])
colnames(test) <- sort(dropsFW$envEffs[,2])
rownames(test) <- droplevels(dropsFW$TD[[1]]$genotype)
test <- data.frame(Variety = rownames(test), test)
colnames(test)[2:ncol(test)] <- sort(dropsFW$envEffs[,2])
test <- reshape2::melt(test)
test$variable <- as.numeric(as.character(test$variable))

## Create line plot for Finlay Wilkinson analysis.
q <- ggplot(data = test, aes(x = variable, 
                               y = value, 
                               group = Variety, 
                               colour = Variety)) + 
  geom_line() + geom_point() + 
  labs(x = "Location", y = traits[i]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(xintercept = meanEnV <- mean(dropsFW$envEffs[,2]), colour = "grey") +
  scale_x_continuous(breaks = round(dropsFW$envEffs[,2]), sec.axis = dup_axis(labels = locations))

ggsave(filename = paste0("output/", study, "_", traits[i], '_stability_adaptability_across_locations.tiff'),
       plot = q,
       device = 'tiff',
       width = 300,
       height = 400,
       units = 'mm',
       dpi = 300)

# GGE Biplot
model <- gge(data.MET, Location, Variety, traits[i], svp = "symmetrical")
a <- plot(model, type = 1)
b <- plot(model, type = 2)
c <- plot(model, type = 3)
d <- arrange_ggplot(a, b, c, tag_levels = "a")

ggsave(filename = paste0("output/", study,"_", traits[i], '_GGE_Biplots.tiff'),
       plot = d,
       device = 'tiff',
       width = 400,
       height = 250,
       units = 'mm',
       dpi = 300)

}

# saving the output file for single-trials analysis
write.csv(out.MET, paste0("output/", study, "_", "MET_KPI_report.csv"))

}
