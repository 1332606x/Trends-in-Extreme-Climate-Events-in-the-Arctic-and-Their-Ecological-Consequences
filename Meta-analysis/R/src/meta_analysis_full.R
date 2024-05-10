## Calculating effect sizes, meta-analysis, and producing orchaRd plots

## load package
library(metafor) # for effect size estimation
library(altmeta) # outliers diagnose
library(dplyr)
library(knitr)
library(kableExtra)
library(orchaRd)
library(ape)
library(rotl)
library(stringr)

## calculate Hedges'd
meta_extracted_data <- read.csv("data/Meta_Analysis_Screening_Data.csv", sep = ";", header = TRUE, row.names = NULL, fileEncoding = "ISO-8859-1")

meta_extracted_data <- subset(meta_extracted_data, !(n_1 == 0 | n_2 == 0))

ehedged <- escalc(measure= "SMDH", n1i = n_2, n2i = n_1, 
                  m1i = x_2, m2i = x_1, 
                  sd1i = sd_2, sd2i = sd_1, data = meta_extracted_data, append = TRUE)

ehedged$final_yi <- ehedged$yi * ehedged$Multiplier

ehedged <- with(ehedged, ehedged[complete.cases(final_yi, vi), ])

## finding the outliers

# in my case, my kept being unable to find the metaoutliers function despite loading the altmeta package
# the functions pasted below are from the package 'altmeta' with github link:
# https://github.com/cran/altmeta/tree/master/R
# I used the 'metaoutliers' function and any other functions it called including:
# metahet.base, tau2.r.solver, and tau2.m.solver

metahet.base <- function(y, s2){
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  n <- length(y)
  
  w <- 1/s2
  mu.bar <- sum(w*y)/sum(w)
  
  out <- NULL
  
  out$weighted.mean <- mu.bar
  
  # the conventional methods
  Q <- sum(w*(y - mu.bar)^2)
  H <- sqrt(Q/(n - 1))
  I2 <- (Q - n + 1)/Q
  tau2.DL <- (Q - n + 1)/(sum(w) - sum(w^2)/sum(w))
  tau2.DL <- max(c(0, tau2.DL))
  out$Q <- Q
  out$H <- H
  out$I2 <- I2
  out$tau2.DL <- tau2.DL
  
  # absolute deviation based on weighted mean
  Qr <- sum(sqrt(w)*abs(y - mu.bar))
  Hr <- sqrt((3.14159*Qr^2)/(2*n*(n - 1)))
  Ir2 <- (Qr^2 - 2*n*(n - 1)/3.14159)/(Qr^2)
  tau2.r <- tau2.r.solver(w, Qr)
  out$Qr <- Qr
  out$Hr <- Hr
  out$Ir2 <- Ir2
  out$tau2.r <- tau2.r
  
  # absolute deviation based on weighted median
  
  expit <- function(x) {ifelse(x >= 0, 1/(1 + exp(-x/0.0001)), exp(x/0.0001)/(1 + exp(x/0.0001)))}
  psi <- function(x) {sum(w*(expit(x - y) - 0.5))}
  mu.med <- uniroot(psi, c(min(y) - 0.001, max(y) + 0.001))$root
  out$weighted.median <- mu.med
  Qm <- sum(sqrt(w)*abs(y - mu.med))
  Hm <- sqrt(3.14159/2)*Qm/n
  Im2 <- (Qm^2 - 2*n^2/3.14159)/Qm^2
  tau2.m <- tau2.m.solver(w, Qm)
  out$Qm <- Qm
  out$Hm <- Hm
  out$Im2 <- Im2
  out$tau2.m <- tau2.m
  
  return(out)
}
tau2.r.solver <- function(w, Qr){
  f <- function(tau2){
    out <- sum(sqrt(1 - w/sum(w) + tau2*(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))) - Qr*sqrt(3.14159/2)
    return(out)
  }
  f <- Vectorize(f)
  tau.upp <- Qr*sqrt(3.14159/2)/sum(sqrt(w - 2*w^2/sum(w) + w*sum(w^2)/(sum(w))^2))
  tau2.upp <- tau.upp^2
  f.low <- f(0)
  f.upp <- f(tau2.upp)
  if(f.low*f.upp > 0){
    tau2.r <- 0
  }else{
    tau2.r <- uniroot(f, interval = c(0, tau2.upp))
    tau2.r <- tau2.r$root
  }
  return(tau2.r)
}
tau2.m.solver <- function(w, Qm){
  f <- function(tau2){
    out <- sum(sqrt(1 + w*tau2)) - Qm*sqrt(3.14159/2)
    return(out)
  }
  f <- Vectorize(f)
  n <- length(w)
  tau2.upp <- sum(1/w)*(Qm^2/n*2/3.14159 - 1)
  tau2.upp <- max(c(tau2.upp, 0.01))
  f.low <- f(0)
  f.upp <- f(tau2.upp)
  if(f.low*f.upp > 0){
    tau2.m <- 0
  }else{
    tau2.m <- uniroot(f, interval = c(0, tau2.upp))
    tau2.m <- tau2.m$root
  }
  return(tau2.m)
}

metaoutliers <- function(y, s2, data, model){
  if(missing(y)) stop("please specify effect size.")
  if(missing(s2)) stop("please specify within-study variance.")
  if(!missing(data)){
    y <- eval(substitute(y), data, parent.frame())
    s2 <- eval(substitute(s2), data, parent.frame())
  }
  if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
  w <- 1/s2
  y.p <- sum(y*w)/sum(w)
  n <- length(y)
  
  if(missing(model)){
    hetmeasure <- metahet.base(y, s2)
    Ir2 <- hetmeasure$Ir2
    if(Ir2 < 0.3){
      model <- "FE"
      cat("This function uses fixed-effect meta-analysis because Ir2 < 30%.\n")
    }else{
      model <- "RE"
      cat("This function uses random-effects meta-analysis because Ir2 >= 30%.\n")
    }
  }
  
  if(!is.element(model, c("FE", "RE"))) stop("wrong input for the argument model.")
  
  y.p.i <- res <- std.res <- numeric(n)
  if(model == "FE"){
    for(i in 1:n){
      w.temp <- w[-i]
      y.temp <- y[-i]
      y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1/sum(w.temp) + s2[i]
      std.res[i] <- res[i]/sqrt(var.res.i)
    }
  }else{
    for(i in 1:n){
      s2.temp <- s2[-i]
      y.temp <- y[-i]
      tau2.temp <- metahet.base(y.temp, s2.temp)$tau2.DL
      w.temp <- 1/(s2.temp + tau2.temp)
      y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
      res[i] <- y[i] - y.p.i[i]
      var.res.i <- 1/sum(w.temp) + s2[i] + tau2.temp
      std.res[i] <- res[i]/sqrt(var.res.i)
    }
  }
  
  outliers <- which(abs(std.res) >= 3)
  if(length(outliers) == 0) outliers <- "All the standardized residuals are smaller than 3"
  
  out <- NULL
  out$model <- model
  out$std.res <- std.res
  out$outliers <- outliers
  
  class(out) <- "metaoutliers"
  return(out)
}

outliers <- metaoutliers(final_yi, vi, ehedged, model = "RE")

ehedged <- ehedged[-outliers$outliers, ]

study_distribution <- ehedged %>%
  group_by(Study_ID) %>%
  summarize(Lat = list(unique(Lat)), Lon = list(unique(Lon))) %>%
  group_by(Lat, Lon) %>%
  summarize(Num_Studies = n()) %>%
  unnest(Lat, Lon) 

file_name <- "output/study_locations.csv"
write.csv(study_distribution, file = file_name, row.names = FALSE)


## Controlling for phylogeny
ehedged$Species <- str_replace(ehedged$Species, "Empetrum hermaphroditum", "Empetrum nigrum") 
ehedged$Species <- str_replace(ehedged$Species, "VacciniumÊmyrtillus", "Vaccinium myrtillus") 
ehedged$Species <- str_replace(ehedged$Species, "MesaphoruraÊmacrochaeta", "Mesaphorura macrochaeta") 
ehedged$Species <- str_replace(ehedged$Species, "ArrhopalitesÊprincipalis", "Arrhopalites principalis") 
ehedged$Species <- str_replace(ehedged$Species, "ProtaphoruraÊgisini", "Protaphorura gisini") 

# Hermaphroditum is a subspecies (Empetrum nigrum sbsp. hermaphroditum) but including both
# prevents the inclusion of Empetrum nigrum so I removed the subspecies
# also re-formatting the names of all species which import with mistakes despite being correct on excel sheet

species_full <- unique(ehedged$Species)
# Match animal species names
name_matches <- tnrs_match_names(species_full, context_name = "All life")

name_matches <- na.omit(name_matches)

rawtree <- tol_induced_subtree(ott_ids = name_matches$ott_id, label_format = "name")

tips <- gsub("_", " ", rawtree$tip.label)
tips <- gsub(" *\\(.*?\\) *", "", tips)
tips

# Which tips are extra? i.e., not in species list
extra <- tips[(tips %in% species_full) == F] 
extra

# Which species missed?
miss <- species_full[(species_full %in% tips) == F] 
miss

tip_to_rename <- "Pygmarrhopalites principalis"
# Find the index of the tip in the tree
tip_index <- which(tips == tip_to_rename)
# Rename the tip label
new_label <- "Arrhopalites principalis"
rawtree$tip.label[tip_index] <- new_label

#plot(rawtree, cex=0.4)
is.binary(rawtree) # tree not binary

# Make tree binary with the bifurcation being "random"
bintree_full <- multi2di(rawtree, random = T) 
is.binary(bintree_full)
#plot(bintree_full,cex=0.2)

# Computing branch lengths
phylo_branch_full <- compute.brlen(bintree_full, method = "Grafen", power = 1)
is.ultrametric(phylo_branch_full)

plot(phylo_branch_full, "f", cex=0.6)

# Making correlation matrix
phylo_branch_full$tip.label <- gsub("\\_", " ", phylo_branch_full$tip.label)
unique(ehedged$Species)[!unique(ehedged$Species) %in% phylo_branch_full$tip.label] 
phylo_branch_full$tip.label[!phylo_branch_full$tip.label %in% unique(ehedged$Species)]

phylo_cor_species <- vcv(phylo_branch_full, cor = T)
phylo_branch_full$node.label <- NULL
phylo_MCMC_species <- MCMCglmm::inverseA(phylo_branch_full, nodes = "ALL", scale = TRUE)$Ainv

## save files to use for analyses 
save(phylo_cor_species, file="./output/phylo_cor_species.Rdata")
save(phylo_branch_full, file="./output/phylotree.Rdata")
save(phylo_MCMC_species, file = "./output/phyloMCMC.Rdata")

jpeg("./output/phylo_branch_full_plot.jpg", quality = 100, width = 3000, height = 3000, units = "px", res = 300)
plot(phylo_branch_full, "f", cex=0.8)
dev.off()

## Models

## We used the `rma.mv` function from the package `metafor` to run all meta-analytic models and meta-regressions. 

ehedged$Species2 <- ehedged$Species

# Changing all categorical variables into factors

ehedged$Kingdom <- as.factor(ehedged$Kingdom) # Kingdom
ehedged$Experimental <- as.factor(ehedged$Experimental)
ehedged$Extreme_event <- as.factor(ehedged$Extreme_event)
ehedged$Experimental <- as.factor(ehedged$Experimental)
ehedged$Measurement_season <- as.factor(ehedged$Measurement_season)
ehedged$Broad_category <- as.factor(ehedged$Broad_category)
ehedged$Species2 <- as.factor(ehedged$Species2)
ehedged$Phylum <- as.factor(ehedged$Phylum)
ehedged$First_author <- as.factor(ehedged$First_author)
ehedged$Kingdom <- as.factor(ehedged$Kingdom) # Kingdom
ehedged$Response_measure <- as.factor(ehedged$Response_measure) # Kingdom


ehedged$Is_Stef_Bokhorst <- as.factor(ifelse(ehedged$First_author == "Stef Bokhorst", TRUE, FALSE))

model.stef <- rma.mv(yi = final_yi, V = vi,
                          mods = ~ Is_Stef_Bokhorst,
                          random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                          R = list(Species = phylo_cor_species),
                          sparse = TRUE,
                          method = "REML", verbose=FALSE, data = ehedged)
model.stef

p.stef <- orchard_plot(model.stef, 
                          mod="Is_Stef_Bokhorst", 
                          xlab = "SMDH", 
                          transfm = "none", 
                          alpha = 0.3, 
                          group = "Study_ID", 
                          angle = 0)
p.stef

## Finding outliers

model_null <- rma.mv(yi = final_yi, V = vi,
                 random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species, ~1|Species2),  
                 R = list(Species = phylo_cor_species),
                 sparse = TRUE,
                 method = "REML", verbose=FALSE, data = ehedged)
model_null

I2.null <- i2_ml(model_null)
round(I2.null,4)

p.null<- orchard_plot(model_null, 
                          xlab = "SMDH", 
                          transfm = "none", 
                          alpha = 0.3, 
                          group = "Study_ID", 
                          angle = 0)
p.null

model.allmods <- rma.mv(yi = final_yi, V = vi,
                 mods = ~ Experimental + Extreme_event + Broad_category + Measurement_season + Kingdom + Phylum, 
                 random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                 R = list(Species = phylo_cor_species),
                 sparse = TRUE,
                 method = "REML", verbose=FALSE, data = ehedged)

model.allmods

p.allmods<- orchard_plot(model.allmods, 
                          xlab = "SMDH", 
                          transfm = "none", 
                          alpha = 0.3, 
                          group = "Study_ID", 
                          angle = 0) 
p.allmods

R2.all <- r2_ml(model.allmods)
round(R2.all,5)*100

# individual mods

model.tax.class <- rma.mv(yi = final_yi, V = vi,
                        mods = ~ Taxonomic_class - 1,
                        random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                        R = list(Species = phylo_cor_species),
                        sparse = TRUE,
                        method = "REML", verbose=FALSE, data = ehedged)
model.tax.class


p.tax.class<- orchard_plot(model.tax.class, 
                              mod="Taxonomic_class", 
                              xlab = "SMDH", 
                              transfm = "none", 
                              alpha = 0.3, 
                              group = "Study_ID", 
                              angle = 0)
p.tax.class

R2.tax.class <- r2_ml(model.tax.class)
round(R2.tax.class,4)*100

model.tax.phylum <- rma.mv(yi = final_yi, V = vi,
                          mods = ~ Phylum - 1,
                          random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                          R = list(Species = phylo_cor_species),
                          sparse = TRUE,
                          method = "REML", verbose=FALSE, data = ehedged)
model.tax.phylum

desired_order <- c("Chordata", "Bryophyta", "Ascomycota", "Arthropoda", "Angiosperms")

desired_order <- c("Angiosperms", "Arthropoda", "Ascomycota", "Bryophyta", "Chordata")

p.tax.phylum<- orchard_plot(model.tax.phylum, 
                           mod="Phylum", 
                           xlab = "SMDH", 
                           transfm = "none", 
                           tree.order = desired_order,
                           alpha = 0.3, 
                           group = "Study_ID", 
                           angle = 0)
p.tax.phylum

R2.tax.phylum <- r2_ml(model.tax.phylum)
round(R2.tax.phylum,5)*100

model.king.grp <- rma.mv(yi = final_yi, V = vi,
                         mods = ~ Kingdom - 1,
                         random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                         R = list(Species = phylo_cor_species),
                         sparse = TRUE,
                         method = "REML", verbose=FALSE, data = ehedged)
model.king.grp

p.king.grp<- orchard_plot(model.king.grp, 
                          mod="Kingdom", 
                          xlab = "SMDH", 
                          transfm = "none", 
                          alpha = 0.3, 
                          group = "Study_ID", 
                          angle = 0)
p.king.grp

R2.king.grp <- r2_ml(model.king.grp)
round(R2.king.grp,4)*100

model.event.grp <- rma.mv(yi = final_yi, V = vi,
                         mods = ~ Extreme_event - 1,
                         random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                         R = list(Species = phylo_cor_species),
                         sparse = TRUE,
                         method = "REML", verbose=FALSE, data = ehedged)
model.event.grp


p.event.grp<- orchard_plot(model.event.grp, 
                          mod="Extreme_event", 
                          xlab = "SMDH", 
                          transfm = "none", 
                          alpha = 0.3, 
                          group = "Study_ID", 
                          angle = 0)
p.event.grp

R2.event.grp <- r2_ml(model.event.grp)
round(R2.event.grp,4)*100


model.exp.grp <- rma.mv(yi = final_yi, V = vi,
                        mods = ~ Experimental - 1,
                        random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                        R = list(Species = phylo_cor_species),
                        sparse = TRUE,
                        method = "REML", verbose=FALSE, data = ehedged)
model.exp.grp

p.exp.grp<- orchard_plot(model.exp.grp, 
                           mod="Experimental", 
                           xlab = "SMDH", 
                           transfm = "none", 
                           alpha = 0.3, 
                           group = "Study_ID", 
                           angle = 0)
p.exp.grp

R2.exp.grp <- r2_ml(model.exp.grp)
round(R2.exp.grp,4)*100

model.cat.grp <- rma.mv(yi = final_yi, V = vi,
                        mods = ~ Broad_category - 1,
                        random = list(~1 | EffectSize_ID, ~1 | Study_ID, ~1|Species2, ~1|Species),  
                        R = list(Species = phylo_cor_species),
                        sparse = TRUE,
                        method = "REML", verbose=FALSE, data = ehedged)
model.cat.grp

p.cat.grp<- orchard_plot(model.cat.grp, 
                         mod="Broad_category", 
                         xlab = "SMDH", 
                         transfm = "none", 
                         alpha = 0.3, 
                         group = "Study_ID", 
                         angle = 0)
p.cat.grp

R2.cat.grp <- r2_ml(model.cat.grp)
round(R2.cat.grp,5)*100


## Summary tables I2 and R2

I2.models <- data.frame("Model" = "Null",
                        "I2_total" = I2.null[1],
                        "I2_EffectSize_ID" = I2.null[2],
                        "I2_study_ID" = I2.null[3],
                        "I2_species" = I2.null[4],
                        "I2_phylo" = I2.null[5])

rownames(I2.models) <- NULL
colnames(I2.models)[2] <- "$I^2_{total}$"
colnames(I2.models)[3] <- "$I^2_{esID}$"
colnames(I2.models)[4] <- "$I^2_{studyID}$"
colnames(I2.models)[5] <- "$I^2_{species}$"
colnames(I2.models)[6] <- "$I^2_{phylo}$"

I2.models %>% 
  kable(format="html", escape=F, caption="Heterogeneity explained in the null model") %>% 
  collapse_rows(columns = 1, valign = "top") %>% 
  kable_styling(c("hover", "condensed"), full_width = F)

R2.models <- data.frame("Model" = c("Full",
                                    "Experimental Y/N",
                                    "Event Type",
                                    "Kingdom",
                                    "Phylum",
                                    "Taxonomic class",
                                    "Broad Category"),
                        
                        "R2" = c(R2.all[1],
                                 R2.exp.grp[1],
                                 R2.event.grp[1],
                                 R2.king.grp[1],
                                 R2.tax.phylum[1],
                                 R2.tax.class[1],
                                 R2.cat.grp[1]))
R2.models2 <- R2.models %>% 
  transmute(Model = Model
            ,R2= round(R2*100,2))

colnames(R2.models2)[2] <- "$R^2(\\%)$"

R2.models2 %>% 
  kable(format="html", escape=F, caption="Heterogeneity explained in the models with moderators") %>% 
  collapse_rows(columns = 1, valign = "top") %>% 
  kable_styling(c("hover", "condensed"), full_width = F)


model.cat.grp

## Funnel plots!

f.null <- funnel(model_null, level=c(90, 95, 99), 
                 shade=c("white", "gray55", "gray75"), 
                 yaxis="seinv", refline=0, 
                 pch=1, 
                 cex= 0.4, bg= "white", 
                 legend=FALSE) 

f.full <- funnel(model.allmods, level=c(90, 95, 99), 
                 shade=c("white", "gray55", "gray75"), 
                 yaxis="seinv", refline=0, 
                 pch=1, 
                 legend=FALSE, ylim= c(0.8, 4))
