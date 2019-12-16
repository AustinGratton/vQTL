install.packages("vqtl")
install.packages("qtl")
library("qtl")
library("vqtl")
library("dplyr")
#we also got rid of "(" in SNP ID rz44bd and rz574bc
#now with the Family dataset
family <-read.cross(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Family.csv"))
family <- drop.nullmarkers(family)
#scan with variance
family <- calc.genoprob(family)
foutv <- scanonevar(cross = family,
                    mean.formula = PlantHeight ~ mean.QTL.add,
                    var.formula = ~ var.QTL.add)

foutvdf<- data.frame(foutv$result$loc.name,
                     foutv$result$mean.lod,
                     foutv$result$mean.asymp.p,
                     foutv$result$var.lod,
                     foutv$result$var.asymp.p,
                     foutv$result$joint.lod,
                     foutv$result$joint.asymp.p)
colnames(foutvdf) = c("SNP Names",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value")
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
#this works
#this is an example of the first one
sizes = effect.sizes(cross = family,
                     phenotype.name = "height.in.",
                     genotype.names = c("AA","BB"),
                     focal.groups = scrub[1])
sizes

fsizedf <- data.frame(NULL)

y = 1:length(foutv$result$loc.name)
y = y[-c(826)]
#effect sizes can not be computed for these 3 SNPs
fsizedf = sapply(y, function(x){
  tempm =  effect.sizes(cross = family,
                        phenotype.name = "PlantHeight",
                        genotype.names = c("AA","BB"),
                        focal.groups = foutv$result$loc.name[x])
  tempv = c(tempm[1,2:7],tempm[2,2:7])
  return(tempv)
})
undebug(effect.sizes)
effect.sizes(cross = random,
             phenotype.name = "height.in.",
             genotype.names = c("AA","BB"),
             focal.groups = foutv$result$loc.name[470])


foutvdf<- data.frame(foutv$result$loc.name,
                     foutv$result$pos,
                     foutv$result$mean.lod,
                     foutv$result$mean.asymp.p,
                     foutv$result$var.lod,
                     foutv$result$var.asymp.p,
                     foutv$result$joint.lod,
                     foutv$result$joint.asymp.p)
#dropping the SNPs whose effect sizes could not be computed
foutvdf = foutvdf[-c(826),]
foutvdf = cbind(foutvdf,fsizedf)
colnames(foutvdf) = c("SNP Names",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
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

write.csv(foutvdf, file = "FamilyvQTL_LOD,Pvals,EffectSizes.csv")
