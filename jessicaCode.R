#install.packages("vqtl")
#install.packages("qtl")
library(qtl)
library(vqtl)
#we also got rid of "(" in SNP ID rz44bd and rz574bc
#setwd("/work/06156/jg1994/stampede2/Stapleton_Lab")
fr <-read.cross(file = "https://raw.githubusercontent.com/jgray07/Stapleton_Lab/master/familyandrandom.csv")
fr <- drop.nullmarkers(fr)
#scan with variance
fr <- calc.genoprob(fr)
routv <- scanonevar(cross = fr,
                    mean.formula = PlantHeight ~ plotstructure * (mean.QTL.add + mean.QTL.dom),
                    var.formula = ~ plotstructure * (var.QTL.add + var.QTL.dom),
                    return.covar.effects = TRUE)
library("dplyr")
effect.sizes = function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL, 
                         genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL, 
                         title = paste(phenotype.name, "by", paste(focal.groups, 
                                                                   collapse = ", ")), draw_ribbons = TRUE, se_line_size = 1, 
                         point_size = 1) 
{
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))], 
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))], 
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[, 
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups, 
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form), 
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df), 
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim - 
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim), 
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim + 
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>% 
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim), 
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub), 
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb), 
                     group.sd.ub = mean(indiv.sd.ub))
  return(group.prediction.tbl)
}

y = 1:length(routv$result$loc.name)
#effect sizes can not be computed for these 3 SNPs
rsizedf = sapply(y, function(x){
  tryCatch({
    print(x)
    tempm =  effect.sizes(cross = fr,
                          phenotype.name = "PlantHeight",
                          genotype.names = c("AA","BB"),
                          focal.groups = routv$result$loc.name[x])
  }, error = function(e) message(e),
  finally = function(tempm){
    tempv = c(tempm[1,2:7],tempm[2,2:7])
    return(unlist(tempv))
  }
  )
})
rsizedf1 <- as.data.frame(matrix(rep(0,length(y)*12), ncol= 12))
sapply(1:length(rsizedf), function(x){
  print(x)
  if(!is.null(rsizedf[[x]])){
    rsizedf1[x,] <<- c(rsizedf[[x]][1,2:7],rsizedf[[x]][2:7])
  }
})
nall0 <-sapply(1:dim(rsizedf1)[1], function(x){
  !all(rsizedf1[x,] == 0)
})  
ditch <- which(nall0 == F)
rsizedf1 <- rsizedf1[-ditch,]
keep <- 1:3235; keep<- keep[-ditch]
routvdf<- data.frame(routv$result$loc.name[keep],
                     routv$result$pos[keep],
                     routv$result$mQTL.lod[keep],
                     routv$result$mQTL.asymp.p[keep],
                     routv$result$vQTL.lod[keep],
                     routv$result$vQTL.asymp.p[keep],
                     routv$result$mvQTL.lod[keep],
                     routv$result$mvQTL.asymp.p[keep],
                     routv$result$`m__est__(Intercept)`[keep],
                     routv$result$m__est__plotstructure[keep],
                     routv$result$`m__se__(Intercept)`[keep],
                     routv$result$m__se__plotstructure[keep],
                     routv$result$`v__est__(Intercept)`[keep],
                     routv$result$v__est__plotstructure[keep],
                     routv$result$`v__se__(Intercept)`[keep],
                     routv$result$v__se__plotstructure[keep])
routvdf = cbind(routvdf,rsizedf1)
colnames(routvdf) = c("SNP Name",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "Mean Est Intercept",
                      "Mean Est Plot Structure",
                      "Mean SE Intercept",
                      "Mean SE Plot Structure",
                      "Variance Est Intercept",
                      "Variance Est Plot Structure",
                      "Variance SE Intercept",
                      "Variance SE Plot Structure",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")


write.csv(routvdf, file = "FamilyRandom0528.csv")
