##----------------
## Fit the time series models for regional trends
## CM: 29/08/2019
## 
##----------------
library(ggplot2); theme_set(theme_bw())
library(TMB)
library(plyr)

## data
source("1_make_data.R")

regions <- as.character(unique(BvB.df$region))
regions <- regions[regions != "Other"] ## no BdivBmsypref
regions <- regions[regions != "Antarctic"] ## one stock
variables <- c("BvB", "UvU", "MCatch")

## remove missing values - will be included again later
BvB.df <- na.omit(subset(BvB.df, year >= 1970 & !region %in% c("Antarctic", "Other")))
UvU.df <- na.omit(subset(UvU.df, year >= 1970 & !region %in% c("Antarctic", "Other")))
MCatch.df <- na.omit(subset(MCatch.df, year >= 1970 & !region %in% c("Antarctic", "Other")))

## NB scaling catch to longterm average
## ddply re-sorts the stockid, so order first
MCatch.df <- MCatch.df[with(MCatch.df, order(stockid, year)), ]
tmp <- ddply(MCatch.df, .(stockid), function(x) {
    cbar <- mean(x$MCatch, na.rm = TRUE)
    z <- data.frame(year = x$year, MCatch = x$MCatch / cbar)
    return(z)
})

all(MCatch.df$year == tmp$year & MCatch.df$stockid == tmp$stockid)

MCatch.df$MCatch <- tmp$MCatch

## get coverage (percentage of stocks present per year)
for(i in variables){
    ## subset for given variable
    var.df <- get(paste0(i, ".df"))
    ## count by region and year
    stock.count <- with(var.df, aggregate(as.formula(paste(i, "~ region + year")), FUN = length))
    stock.count <- stock.count[with(stock.count, order(region, year)), ]
    names(stock.count)[names(stock.count) == i] <- "N"
    ## total number of stocks by region
    stock.total <- aggregate(stockid ~ region, data = unique(var.df[, c("region", "stockid")]), FUN = length)
    names(stock.total)[names(stock.total) == "stockid"] <- "Ntotal"
    stock.count <- merge(stock.count, stock.total)
    stock.count$Coverage <- with(stock.count, N / Ntotal)
    assign(paste0(i, ".coverage"), stock.count)
}

## include others here
ggplot(BvB.coverage, aes(x = year, y = Coverage)) + geom_point() + facet_wrap(~region)+
    geom_point(data= UvU.coverage, colour = "blue") +
    geom_point(data= MCatch.coverage, colour = "red")

## compile TMB code
## with AR(1) on the residuals
compile("dlm_ar1.cpp")

## load the function
dyn.load(dynlib("dlm_ar1"))

get.trend <- function(region.val, variable){
    print(variable)
    dat <- get(paste0(variable, ".df"))
    dat <- na.omit(dat)
    dat <- subset(dat, region == region.val)
    dat <- dat[order(dat$stockid, dat$year), ]
    dat$fyear <- factor(dat$year)
    dat$lnvar <- log(dat[, variable] + 0.001)
    dat <- droplevels(dat)
    ## coverage
    coverage.df <- get(paste0(variable, ".coverage"))
    coverage.df <- subset(coverage.df, region == region.val)
    ##-----
    ## LM
    ##-----
    dat$stockid2 <- dat$stockid
    lm.fit <- lm(lnvar ~ stockid + factor(year), data = dat)
    coef.lm <- coef(lm.fit)
    lm.effects <- c(coef.lm["(Intercept)"], coef.lm["(Intercept)"] + coef.lm[grep("year", names(coef.lm))])
    ##-----
    ## TMB
    ##-----
    ## reshape the data for tmb
    dat.wide <- reshape(dat[, c("stockid", "year", "lnvar")], idvar = "stockid", timevar = "year", direct = "wide")
    rownames(dat.wide) <- dat.wide[, "stockid"]
    dat.wide <- dat.wide[, names(dat.wide) != "stockid"]
    y <- as.matrix(dat.wide)
    ## make sure ordered
    y <- y[, order(as.numeric(gsub("lnvar.", "", colnames(y))))]
    if(!all(order(colnames(y)) == 1:ncol(y))){
        stop("y matrix not ordered correctly")
    }
    n <- ncol(y)
    m <- nrow(y)
    ypresent <- ifelse(is.na(y), 0, 1)
    first.obs <- apply(ypresent, 1, which.max) - 1 ## -1 for start at zero in TMB
    ## create the AD object
    obj <- MakeADFun(
        data = list(
            y = y,
            ypresent = ypresent,
            first_obs = first.obs
        ),
        parameters = list(
            lnsde = log(0.1),
            lnsdx = log(0.1),
            logitrho = -log(2/(1 + 0.5) - 1), ## for AR(1) = 0.5
            x = rep(0, n),
            Apar = rep(0, m - 1)
        ),
        random = c("x"),
        DLL = "dlm_ar1",
        silent = TRUE)    
    ## fit the model
    opt <- nlminb(objective = obj$fn,
                  gradient = obj$gr,
                  start = obj$par,
                  lower = c(lnsde = log(0.05), lnsdx  = log(0.05)),
                  control = list(iter.max = 1e3, eval.max = 1e3))
    if(opt$convergence == 0){
        ## report
        rep <- sdreport(obj)
        srep <- summary(rep)
        xhat <- srep[rownames(srep) == "x", ]
        rownames(xhat) <- NULL
        xhat <- as.data.frame(xhat)
        ## finite population correction
        xhat$region <- region.val
        xhat$year <- as.numeric(gsub("lnvar.", "", colnames(y)))
        xhat <- merge(xhat, coverage.df)
        xhat <- xhat[order(xhat$year), ]
        years <- xhat$year
        coverage <- xhat$Coverage
        ## finite-population corrected
        xhat$fpc.se <- xhat[, "Std. Error"] * with(xhat, sqrt((Ntotal - N)/(Ntotal - 1)))
        dlm.geomean <- exp(xhat[, "Estimate"])
        dlm.upper <- exp(xhat[, "Estimate"] + 1.96 * xhat[, "fpc.se"])
        dlm.lower <- exp(xhat[, "Estimate"] - 1.96 * xhat[, "fpc.se"])
    }else{
        na.vec <- rep(NA, n)
        years <- coverage <- dlm.geomean <- dlm.upper <- dlm.lower <- na.vec
    }
    ## predictions
    pred.df <- data.frame(region = region.val,
                          variable = variable,
                          year = years,
                          Coverage = coverage,
                          dlm.geomean = dlm.geomean,
                          dlm.lower = dlm.lower,
                          dlm.upper = dlm.upper,
                          fixed.effects = exp(as.numeric(lm.effects)),
                          stringsAsFactors = FALSE)    
    pred.df <- pred.df[order(pred.df$year), ]
    ## re-scale to median of full coverage years
    box <- boxplot(as.formula(paste(variable, "~ region + year")), data = dat, plot = FALSE)
    stats.df <- as.data.frame(t(box$stats))
    names(stats.df) <- c("lower.whisker", "q.25", "median", "q.75", "upper.whisker")
    stats.df$year <- as.numeric(gsub(paste0(region.val, "."), "", box$name))
    ## link to coverage to get scaling
    median.df <- stats.df[, c("year", "median")]
    median.df <- merge(median.df, coverage.df)
    ## note 90% because of Med
    median.df <- subset(median.df, region == region.val & Coverage > 0.9)
    tmp <- merge(pred.df, median.df)
    ## scale
    scale.dlm <- with(tmp, sum(median) / sum(dlm.geomean))
    pred.df[, c("dlm.geomean", "dlm.lower", "dlm.upper")] <- scale.dlm * pred.df[, c("dlm.geomean", "dlm.lower", "dlm.upper")]
    ## scale fixed effects
    scale.fixed <- with(tmp, sum(median) / sum(fixed.effects))
    pred.df[, "fixed.effects"] <- scale.fixed * pred.df[, "fixed.effects"]    
    pred.df <- merge(pred.df, stats.df)
    pred.df <- pred.df[order(pred.df$year), ]
    return(pred.df)
}

## container for estimates
est.df <- NULL

for(i in 1:length(variables)){
    for(j in 1:length(regions)){
        tmp <- get.trend(region.val = regions[j], variable = variables[i])
        est.df <- rbind(est.df, tmp)
        }
}

## barplots
BvB.df <- merge(BvB.df, BvB.coverage)
UvU.df <- merge(UvU.df, UvU.coverage)
MCatch.df <- merge(MCatch.df, MCatch.coverage)

## to long format for legend
library(reshape)
est.long.df <- melt(est.df[, c("year", "region", "variable", "dlm.geomean", "fixed.effects", "median")], id.vars = c("year", "region", "variable"), variable_name = "method")
est.long.df$Method <- NA
est.long.df$Method[est.long.df$method == "dlm.geomean"] <- "State space model"
est.long.df$Method[est.long.df$method == "fixed.effects"] <- "Fixed effects"
est.long.df$Method[est.long.df$method == "median"] <- "Median"
names(est.long.df)[names(est.long.df) == "value"] <- "Estimate"

## subset for median and dlm
est.long.df <- subset(est.long.df, method %in% c("dlm.geomean", "median"))
est.long.df <- droplevels(est.long.df)

## BvB
date <- format(Sys.Date(), "%m_%d_%Y")

pdf(paste0("../tex/figures/BvB_trends_boxplots_", date, ".pdf"), height = 12, width = 8)
ggplot(BvB.df, aes(x = year, y = BvB)) +
    geom_boxplot(aes(group = cut_width(year, 1), fill = Coverage, colour = Coverage), outlier.shape = NA, size = 0.05, fatten = NULL, colour = "grey", width = 1) +
    geom_ribbon(data = subset(est.df, variable == "BvB"), aes(ymin = dlm.lower, ymax = dlm.upper, y = dlm.geomean), fill = "darkorange", alpha = 0.4) +
    facet_wrap(~ region, ncol = 3) +
    scale_fill_gradient(name = "Coverage", low = "white", high = "cadetblue", limits=c(0, 1)) +
    coord_cartesian(ylim = c(0, 4)) +
    geom_point(data = subset(est.long.df, variable == "BvB" & !method %in% c("dlm.geomean")), aes(y = Estimate, group = Method, colour = Method), size = 0.5) +
    geom_line(data = droplevels(subset(est.long.df, variable == "BvB" & method %in% c("dlm.geomean"))), aes(y = Estimate, group = Method, colour = Method)) +
    scale_colour_manual(name = "Method", values = c("red", "darkorange"),##values = c("gold", "red", "darkorange"),
                        guide = guide_legend(override.aes = list(
                                                 linetype = c("blank", "solid"),
                                                 shape = c(19, NA),
                                                 size = c(1, 1)))
                        ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")+
    geom_hline(yintercept = 1, linetype = 2, size = 0.4) +
    xlab("Year") +
    ylab(expression("Relative biomass " * "(B/" * B[MSY] * ")"))
    ##scale_fill_manual(values = c(1,2))
dev.off()

## UvU
pdf(paste0("../tex/figures/UvU_trends_boxplots_", date, ".pdf"), height = 12, width = 8)
ggplot(UvU.df, aes(x = year, y = UvU)) +
    geom_boxplot(aes(group = cut_width(year, 1), fill = Coverage, colour = Coverage), outlier.shape = NA, size = 0.05, fatten = NULL, colour = "grey", width = 1) +
    geom_ribbon(data = subset(est.df, variable == "UvU"), aes(ymin = dlm.lower, ymax = dlm.upper, y = dlm.geomean), fill = "darkorange", alpha = 0.4) +
    facet_wrap(~ region, ncol = 3) +
    scale_fill_gradient(name = "Coverage", low = "white", high = "cadetblue", limits=c(0, 1)) +
    coord_cartesian(ylim = c(0, 4)) +
    geom_point(data = subset(est.long.df, variable == "UvU" & !method %in% c("dlm.geomean")), aes(y = Estimate, group = Method, colour = Method), size = 0.5) +
    geom_line(data = droplevels(subset(est.long.df, variable == "UvU" & method %in% c("dlm.geomean"))), aes(y = Estimate, group = Method, colour = Method)) +
    scale_colour_manual(name = "Method", values = c("red", "darkorange"), ## values = c("gold", "red", "darkorange"),
                        guide = guide_legend(override.aes = list(
                                                 linetype = c("blank", "solid"),
                                                 shape = c(19, NA),
                                                 size = c(1, 1)))
                        ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")+
    geom_hline(yintercept = 1, linetype = 2, size = 0.4) +
    xlab("Year") +
    ylab(expression("Relative fishing pressure " * "(U/" * U[MSY] * ")"))
    ##scale_fill_manual(values = c(1,2))
dev.off()

## MCatch
pdf(paste0("../tex/figures/MCatch_trends_boxplots_", date, ".pdf"), height = 12, width = 8)
ggplot(MCatch.df, aes(x = year, y = MCatch)) +
    geom_boxplot(aes(group = cut_width(year, 1), fill = Coverage, colour = Coverage), outlier.shape = NA, size = 0.05, fatten = NULL, colour = "grey", width = 1) +
    geom_ribbon(data = subset(est.df, variable == "MCatch"), aes(ymin = dlm.lower, ymax = dlm.upper, y = dlm.geomean), fill = "darkorange", alpha = 0.4) +
    facet_wrap(~ region, ncol = 3) +
    scale_fill_gradient(name = "Coverage", low = "white", high = "cadetblue", limits=c(0, 1)) +
    coord_cartesian(ylim = c(0, 3)) +
    geom_point(data = subset(est.long.df, variable == "MCatch" & !method %in% c("dlm.geomean")), aes(y = Estimate, group = Method, colour = Method), size = 0.5) +
    geom_line(data = droplevels(subset(est.long.df, variable == "MCatch" & method %in% c("dlm.geomean"))), aes(y = Estimate, group = Method, colour = Method)) +
    scale_colour_manual(name = "Method", values = c("red", "darkorange"), ##values = c("gold", "red", "darkorange"),
                        guide = guide_legend(override.aes = list(
                                                 linetype = c("blank", "solid"),
                                                 shape = c(19, NA),
                                                 size = c(1, 1)))
                        ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")+
    geom_hline(yintercept = 1, linetype = 2, size = 0.4) +
    xlab("Year") +
    ylab("Catch / (Mean Catch)")    
    ##scale_fill_manual(values = c(1,2))
dev.off()

est.df <- subset(est.df, !variable %in% c("BTrend", "UTrend"))

## output a csv
write.csv(est.df, file = paste0("../data/state_space_results_", format(Sys.Date(), "%m_%d_%Y"), ".csv"), row.names = FALSE)

