
# Author: triffe
###############################################################################

source("/data/commons/triffe/git/HMD_Rlifetables_git/RLifeTables/R/Exposures_per.R")
HFDmonthly <- list.files( "/data/commons/triffe/Desktop/HFDmonthly")

# HMD countries minus BEL
HMDcntries <- c("AUS","AUT","BGR","BLR","CAN","CHE","CHL","CZE","DNK","DEUTE","DEUTW","DEUTNP",
  "ESP","EST","FIN","FRATNP","GBR","GBRCENW","GBR_NIR","GBR_NP","GBR_SCO",
  "HUN","IRL","ISL","ITA","JPN","LTU","LVA","NLD","NOR","POL","PRT",
  "NZL_MA","NZL_NM","NZL_NP","SWE","USA","SVK","SVN","UKR","FRACNP","GBRTENW","ISR","LUX","RUS","TWN")
# 30 countries to be tested, not bad
(TestCountries <- HMDcntries[paste0(HMDcntries, "monthly.txt") %in% HFDmonthly])
TestCountries  <- TestCountries[!TestCountries %in% c("TWN")] # permission probs
working.path   <- "/data/commons/hmd/HMDWORK"
TestOutput     <- list()

working.path   <- "/data/commons/hmd/HMDWORK"
TestCountries <- "ITA"

for (.ctry in TestCountries){
  ctry.list <- list()
  for (.sex in c("m","f")){
    ctry.list[[.sex]] <- list()
    for (.v in c(5,6)){
      ctry.list[[.sex]][[paste0("v",.v)]] <-
      try(Exposures_per(country.folder = file.path(working.path, .ctry), 
        sex = .sex, 
        OPENAGE = 110, 
        save.bin = FALSE, 
        MPVERSION = .v, #MPversion = 6
        test = TRUE,
        Monthly.folder = "/data/commons/triffe/Desktop/HFDmonthly"
      ))
    }
  }
  TestOutput[[.ctry]] <- ctry.list 
}
# save, to have a copy
save(TestOutput, file = "/data/commons/triffe/Desktop/TestOutput.Rdata")
TO <- local(get(load("/data/commons/triffe/Desktop/TestOutput.Rdata")))

# if you want Mx estimates by version and sex
Mxmv5 <- lapply(TO, function(x){ x[["m"]][["v5"]][["Deaths"]] / x[["m"]][["v5"]][["Exp"]]})
Mxfv5 <- lapply(TO, function(x){ x[["f"]][["v5"]][["Deaths"]] / x[["f"]][["v5"]][["Exp"]]})
Mxmv6 <- lapply(TO, function(x){ x[["m"]][["v6"]][["Deaths"]] / x[["m"]][["v6"]][["Exp"]]})
Mxfv6 <- lapply(TO, function(x){ x[["f"]][["v6"]][["Deaths"]] / x[["f"]][["v6"]][["Exp"]]})


# look at relative difference:
# males
RExpm <- lapply(TO , function(x){
  v6Exp <- x[["m"]][["v6"]][["Exp"]]
  v5Exp <- x[["m"]][["v5"]][["Exp"]]
  (v6Exp - v5Exp) / v5Exp
})
# females
RExpf <- lapply(TO , function(x){
  v6Exp <- x[["f"]][["v6"]][["Exp"]]
  v5Exp <- x[["f"]][["v5"]][["Exp"]]
  (v6Exp - v5Exp) / v5Exp
})
# a color ramp
colfun          <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9, "RdBu"), space = "Lab")

# where 'path' is the folder to save diagnostics to
plot.Rdiff <- function(R, ctry, sex = "m", path = "/data/commons/triffe/Desktop/ExpDiagnostics"){
  years           <- as.integer(colnames(R))
  ages            <- as.integer(rownames(R))
  Width           <- 7 * (length(years) + 1) / 131 + 3
  Height          <- 7 + 2 # fixed ( 1 inch margins plus 7 inch tall lexis surface)
  brks            <- seq(-.05,.05,by=.01)
  R[R < -.05]      <- -.05
  R[R > .05]       <- .05
  R[zapsmall(R) == 0] <- NA
  # dev.new( height = Height, width = Width)
  pdf(file.path(path,paste0("Rdiff",ctry,sex,".pdf")), height = Height, width = Width)
  fields:::image.plot(x = years + .5, y = ages + .5, t(R), 
                      asp = 1,
                      breaks = brks,
                      col = rev(colfun(length(brks) - 1)),
                      useRaster = TRUE, 
                      zlim = c(-.05, .05),
                      xlim = range(years) + c(0, 1),
                      ylim = range(ages) + c(0, 1),
                      xlab = "",
                      ylab = "",
                      axes = FALSE,
                      main = paste0(ctry, " ", ifelse(sex == "m","Males","Females"),"\nRelative change in exposures, (v6 - v5) / v5"))
  # fill NAs grey- need to superimpose
  image(x = years + .5, y = ages + .5, t(is.na(R)), col = c(NA, gray(.92)), add = TRUE)
  # nifty lexis reference lines
  source("/data/commons/triffe/git/HMD_Rlifetables_git/Diagnostics/utils_LexRef5.R")
  LexRef5(ages, years, col = "#44444420", lwd = .5)
  
  # bounding box
  rect(min(years), 0, max(years) + 1,max(ages) + 1, border = gray(.2), lwd=.5)
  # year axis
  xyrs <- years[years %% 10 == 0]
  segments(xyrs,0,xyrs,-1, xpd =TRUE)
  text(xyrs, -3, xyrs, xpd =TRUE)
  # age axis
  yages <- ages[ages %% 10 == 0]
  segments(min(years),yages,min(years)-1,yages, xpd =TRUE)
  text(min(years),yages, yages, xpd =TRUE,pos=2)
  
  dev.off()
}

# this fills the folder with diagnostic plots. grey are 0 difference values
for (i in 1:length(RExpm)){
  plot.Rdiff(RExpm[[i]],ctry = names(RExpm)[[i]], sex = "m") 
  plot.Rdiff(RExpf[[i]],ctry = names(RExpf)[[i]], sex = "f")
}


# -------------------------------------------------------------------
# Lisa scanned the French monthly births data. Try to work with that now.
#BM              <- read.table("/data/commons/triffe/Desktop/HFDmonthly/FRATNPmonthly.txt",
#                              header = TRUE, 
#                              sep = ",", 
#                              stringsAsFactors = FALSE, 
#                              na.strings = ".")
#head(BM)
#BO <- read.csv("/data/commons/triffe/Desktop/HFDmonthly/FRA.births.by-1861-1945.csv", stringsAsFactors = FALSE, )              
#rownames(BO) <- as.character(1861:1945)
#BO <- as.matrix(BO[, c(2, 5:16)])
#colnames(BO) <- c(13, 1:12)
#BOlong <- reshape2:::melt(BO, 
#                          varnames = c("Year", "Month"), 
#                          value.name = "Births")
#B.add.cols <- c("PopName", "Area", "Year", "YearReg", "Month", "Vital",
#                          "Births", "Access", "Note1", "Note2", "Note3", "RefCode", "LDB")
#B.addmonthly      <- as.data.frame(matrix( ncol = length(B.add.cols), 
#                            nrow = nrow(BOlong),
#                            dimnames = list(NULL, B.add.cols)))
## fill in data
#B.addmonthly[, colnames(BOlong)]    <- BOlong
#                        
## remaining columns:
#B.addmonthly$PopName                    <- "FRATNP"
#B.addmonthly$YearReg                    <- B.addmonthly$Year
#B.addmonthly$Access                     <- "O"
#B.addmonthly$LDB                        <- 1
#                        
## resort (Month within year, TOT (13) at bottom)
#B.addmonthly <- B.addmonthly[with(B.addmonthly,order(Year, Month)), ]
#B.addmonthly$Month[B.addmonthly$Month == 13] <- "TOT"
#
#FRATNPmonthly <- rbind(B.addmonthly, BM)
#FRATNPmonthly[is.na(FRATNPmonthly)]   <- "."
#head(FRATNPmonthly)
#source("/data/commons/triffe/git/HMD_CS/HMDwork/Input_Templates/removeDuplicates.R")
#FRATNPmonthly <- removeDuplicates(FRATNPmonthly)
#
#write.table(FRATNPmonthly, 
#  file = "/data/commons/triffe/Desktop/HFDmonthly/FRATNPmonthly_add.txt", 
#  sep = ",", 
#  col.names = colnames(FRATNPmonthly), 
#  row.names = FALSE,
#  quote = FALSE)

# now just read in /data/commons/triffe/Desktop/HFDmonthly/FRATNPmonthly_add.txt for the France testing.
working.path   <- "/data/commons/hmd/HMDWORK"
fem5 <- try(Exposures_per(WORKING = file.path(working.path, "ITA"), 
    sex = "f", 
    OPENAGE = 110, 
    save.bin = FALSE, 
    MPVERSION = 5, #MPversion = 6
    test = TRUE,
    Monthly.folder = "/data/commons/triffe/Desktop/HFDmonthly"
  ))
fem6 <- try(Exposures_per(country.folder = file.path(working.path, "ITA"), 
    sex = "f", 
    openage = 110, 
    save.bin = FALSE, 
    MPversion = 6, #MPversion = 6
    test = TRUE,
    Monthly.folder = "/data/commons/triffe/Desktop/HFDmonthly"
  ))
mal5 <- try(Exposures_per(country.folder = file.path(working.path, "FRATNP"), 
    sex = "m", 
    openage = 110, 
    save.bin = FALSE, 
    MPversion = 5, #MPversion = 6
    test = TRUE,
    Monthly.folder = "/data/commons/triffe/Desktop/HFDmonthly"
  ))
mal6 <- try(Exposures_per(country.folder = file.path(working.path, "FRATNP"), 
    sex = "m", 
    openage = 110, 
    save.bin = FALSE, 
    MPversion = 6, #MPversion = 6
    test = TRUE,
    Monthly.folder = "/data/commons/triffe/Desktop/HFDmonthly"
  ))
fR <- (fem6$Exp - fem5$Exp) / fem5$Exp
mR <- (mal6$Exp - mal5$Exp) / mal5$Exp

plot.Rdiff(fR,ctry = "FRATNP", sex = "f") 
plot.Rdiff(mR,ctry = "FRATNP", sex = "m")
source("/data/commons/triffe/git/HMD_Rlifetables_git/R/ltper_AxN.R")

# 1915 and 1919 cohorts have biggest variations

Mxm5 <- (mal5$Deaths / mal5$Exp)[1:100,]
Mxf5 <- (fem5$Deaths / fem5$Exp)[1:100,]
Mxm6 <- (mal6$Deaths / mal6$Exp)[1:100,]
Mxf6 <- (fem6$Deaths / fem6$Exp)[1:100,]
#Mx <- Mxm5
# let's look at before and after log Mx surfaces:
colfun          <- grDevices:::colorRampPalette(RColorBrewer:::brewer.pal(9, "Spectral"), space = "Lab")

# where 'path' is the folder to save diagnostics to
plot.MxSurf <- function(Mx, ctry, v = 5, sex = "m", path = "/data/commons/triffe/Desktop/ExpDiagnostics"){
  years           <- as.integer(colnames(Mx))
  ages            <- as.integer(rownames(Mx))
  Width           <- 7 * (length(years) + 1) / 131 + 3
  Height          <- 7 + 2 # fixed ( 1 inch margins plus 7 inch tall lexis surface)
  Mx              <- log(Mx)
  brks            <- seq(-10, .2, by = .2)
  ticks           <- 10 ^ (-5:1) %o% seq(1:9)
  labs            <- ticks
  labs[, c(3,4,6,7,8,9)] <- "" 
  labs            <- c(t(labs))
  ticks           <- c(t(ticks))
  keep            <- log(ticks) >= min(brks) & log(ticks) <= max(brks)
  labs            <- labs[keep]
  ticks           <- ticks[keep]
  
  # dev.new( height = Height, width = Width)
  pdf(file.path(path, paste0(ctry, "Mx", sex, "v", v, ".pdf")), height = Height, width = Width)
  fields:::image.plot(x = years + .5, y = ages + .5, t(Mx), 
    asp = 1,
    breaks = brks,
    col = rev(colfun(length(brks) - 1)),
    useRaster = TRUE, 
    zlim = range(brks),
    xlim = range(years) + c(0, 1),
    ylim = range(ages) + c(0, 1),
    xlab = "",
    ylab = "",
    axes = FALSE,
    main = paste0(ctry, " ", ifelse(sex == "m","Males","Females"),"\nlog Mx surface (version ",v,")"),
    axis.args = list(at = log(ticks), labels = labs))

  # nifty lexis reference lines
  source("/data/commons/triffe/git/HMD_Rlifetables_git/Diagnostics/utils_LexRef5.R")
  LexRef5(ages, years, col = "#85858520", lwd = .5)
  
  # bounding box
  rect(min(years), 0, max(years) + 1,max(ages) + 1, border = gray(.2), lwd=.5)
  # year axis
  xyrs <- years[years %% 10 == 0]
  segments(xyrs,0,xyrs,-1, xpd =TRUE)
  text(xyrs, -3, xyrs, xpd =TRUE)
  # age axis
  yages <- ages[ages %% 10 == 0]
  segments(min(years),yages,min(years)-1,yages, xpd =TRUE)
  text(min(years),yages, yages, xpd =TRUE,pos=2)
  
  dev.off()
}

plot.MxSurf(Mxm5,"FRATNP",5,"m")
plot.MxSurf(Mxf5,"FRATNP",5,"f")
plot.MxSurf(Mxm6,"FRATNP",6,"m")
plot.MxSurf(Mxf6,"FRATNP",6,"f")




