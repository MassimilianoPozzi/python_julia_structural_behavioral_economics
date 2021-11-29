### This script reads in the data from ExpertsCohens.csv.
### and creates Figure 8 and Online Appendix Figure 7 in the paper.

rm(list=ls())
library('Hmisc')

### Set the names for the outputfiles
fig8.outputfile = '~/Desktop/MetaAnalysis/ExpertLitPlot_Fig8.png'
appfig7.outputfile = '~/Desktop/MetaAnalysis/ExpertLitPlot_AppFig7.png'

############################### Figure 8 ###############################
landscape.a4.width = 11.69
landscape.a4.height = 8.27
png.scale = 1

png(filename = fig8.outputfile,
    height = landscape.a4.height*png.scale, width = landscape.a4.width*png.scale,
    units = "in", res = 300)

# Load in data for the plot
df = read.csv('~/Desktop/MetaAnalysis/ExpertsCohens.csv')

### Set some plotting parameters (e.g. colors, symbols, sizes, and limits for x- and y-axes)
my.pch = c(19,17)
my.colors = c("blue","red")
my.xlim = c(1550,1950)
my.ylim = c(1400,2100)

my.cex.main = 1.7
my.cex.lab = 1.4
my.cex.axis = 1.6
my.cex.text = 1.2
my.cex = 1.3

my.arrow.length = 0.1

########### Plot of Predictions of Experts/Literature versus Actual Effort #############
matplot(x=df$mean_t,
        y=cbind(df$forecast_t, df$lit_forecast_t),
        xlim = my.xlim, ylim = my.ylim,
        main = "Expert Prediction/Literature versus Actual Effort",
        ylab = "Effort Predicted by Experts or Implied by Literature",
        xlab = "Actual Treatment Group Effort",
        cex = my.cex, cex.main = my.cex.main, cex.axis = my.cex.axis, cex.lab = my.cex.lab,
        pch = my.pch, col = my.colors)
abline(a=0,b=1,lty=2)

### "Manually" add arrows and labels for each of the different treatments
arrows(x0 = df$mean_t[1]+20, x1=df$mean_t[1],
       y0 = 0.5*df$forecast_t[1] + 0.5*df$lit_forecast_t[1], y1 = df$forecast_t[1],
       length = my.arrow.length)
arrows(x0 = df$mean_t[1]+20, x1=df$mean_t[1],
       y0 = 0.5*df$forecast_t[1] + 0.5*df$lit_forecast_t[1], y1 = df$lit_forecast_t[1],
       length = my.arrow.length)
text(x = df$mean_t[1]+20,
     y = 0.5*df$forecast_t[1] + 0.5*df$lit_forecast_t[1],
     labels = as.character(df$Treatment[1]), pos = 4, cex = my.cex.text)

arrows(x0 = df$mean_t[2]-20, x1=df$mean_t[2],
       y0 = 0.5*df$forecast_t[2] + 0.5*df$lit_forecast_t[2]-150, y1 = df$forecast_t[2],
       length = my.arrow.length)
arrows(x0 = df$mean_t[2]-20, x1=df$mean_t[2],
       y0 = 0.5*df$forecast_t[2] + 0.5*df$lit_forecast_t[2]-150, y1 = df$lit_forecast_t[2],
       length = my.arrow.length)
text(x = df$mean_t[2]-20,
     y = 0.5*df$forecast_t[2] + 0.5*df$lit_forecast_t[2]-150,
     labels = as.character(df$Treatment[2]), pos = 1, cex = my.cex.text)

arrows(x0 = df$mean_t[3]-100, x1=df$mean_t[3],
       y0 = 0.5*df$forecast_t[3] + 0.5*df$lit_forecast_t[3]+20, y1 = df$forecast_t[3],
       length = my.arrow.length)
arrows(x0 = df$mean_t[3]-100, x1=df$mean_t[3],
       y0 = 0.5*df$forecast_t[3] + 0.5*df$lit_forecast_t[3]+20, y1 = df$lit_forecast_t[3],
       length = my.arrow.length)
text(x = df$mean_t[3]-100,
     y = 0.5*df$forecast_t[3] + 0.5*df$lit_forecast_t[3]+20,
     labels = as.character(df$Treatment[3]), pos = 2, cex = my.cex.text)

arrows(x0 = df$mean_t[4]+20, x1=df$mean_t[4],
       y0 = 0.5*df$forecast_t[4] + 0.5*df$lit_forecast_t[4], y1 = df$forecast_t[4],
       length = my.arrow.length)
arrows(x0 = df$mean_t[4]+20, x1=df$mean_t[4],
       y0 = 0.5*df$forecast_t[4] + 0.5*df$lit_forecast_t[4], y1 = df$lit_forecast_t[4],
       length = my.arrow.length)
text(x = df$mean_t[4]+20,
     y = 0.5*df$forecast_t[4] + 0.5*df$lit_forecast_t[4],
     labels = as.character(df$Treatment[4]), pos = 4, cex = my.cex.text)

arrows(x0 = df$mean_t[5]-50, x1=df$mean_t[5],
       y0 = 0.5*df$forecast_t[5] + 0.5*df$lit_forecast_t[5]-50, y1 = 0.5*df$forecast_t[5] + 0.5*df$lit_forecast_t[5],
       length = my.arrow.length)
text(x = df$mean_t[5]-50,
     y = 0.5*df$forecast_t[5] + 0.5*df$lit_forecast_t[5]-50,
     labels = as.character(df$Treatment[5]), pos = 2, cex = my.cex.text)

arrows(x0 = df$mean_t[6]-20, x1=df$mean_t[6],
       y0 = 0.75*df$forecast_t[6] + 0.25*df$lit_forecast_t[6], y1 = df$forecast_t[6],
       length = my.arrow.length)
arrows(x0 = df$mean_t[6]-20, x1=df$mean_t[6],
       y0 = 0.75*df$forecast_t[6] + 0.25*df$lit_forecast_t[6], y1 = df$lit_forecast_t[6],
       length = my.arrow.length)
text(x = df$mean_t[6]-20,
     y = 0.75*df$forecast_t[6] + 0.25*df$lit_forecast_t[6],
     labels = as.character(df$Treatment[6]), pos = 2, cex = my.cex.text)

arrows(x0 = df$mean_t[7]+20, x1=df$mean_t[7],
       y0 = 0.5*df$forecast_t[7] + 0.5*df$lit_forecast_t[7], y1 = df$forecast_t[7],
       length = my.arrow.length)
arrows(x0 = df$mean_t[7]+20, x1=df$mean_t[7],
       y0 = 0.5*df$forecast_t[7] + 0.5*df$lit_forecast_t[7], y1 = df$lit_forecast_t[7],
       length = my.arrow.length)
text(x = df$mean_t[7]+20,
     y = 0.5*df$forecast_t[7] + 0.5*df$lit_forecast_t[7],
     labels = as.character(df$Treatment[7]), pos = 4, cex = my.cex.text)

arrows(x0 = df$mean_t[8]-50, x1=df$mean_t[8],
       y0 = 0.5*df$forecast_t[8] + 0.5*df$lit_forecast_t[8]+100, y1 = df$forecast_t[8],
       length = my.arrow.length)
arrows(x0 = df$mean_t[8]-50, x1=df$mean_t[8],
       y0 = 0.5*df$forecast_t[8] + 0.5*df$lit_forecast_t[8]+100, y1 = df$lit_forecast_t[8],
       length = my.arrow.length)
text(x = df$mean_t[8]-50,
     y = 0.5*df$forecast_t[8] + 0.5*df$lit_forecast_t[8]+100,
     labels = as.character(df$Treatment[8]), pos = 2, cex = my.cex.text)

legend("topleft", pch=my.pch, col=my.colors, cex = my.cex,
       legend=c("Expert Forecasts", "Effort Implied by Literature"))
errbar(x = df$mean_t, y = df$forecast_t, add=TRUE, type="n", errbar.col=my.colors[1],
       yplus = df$forecast_t + 2*df$se_forecast_t,
       yminus = df$forecast_t - 2*df$se_forecast_t)
errbar(x = df$mean_t, y = df$lit_forecast_t, add=TRUE, type="n", errbar.col=my.colors[2],
       yplus = df$lit_forecast_t + 2*df$se_lit_forecast_t,
       yminus = df$lit_forecast_t - 2*df$se_lit_forecast_t)

# "Manually" add text comparing experts' and "literature's" predictions to the truth
text(x = 1550, y = 1900, pos = 4, cex = my.cex.text + 0.1,
     labels = paste0("Mean Absolute Difference of\nExpert Forecast from Effort = ",
                     round(mean(abs(df$mean_t - df$forecast_t))),
                     "\n\nMean Absolute Difference of\nLiterature Forecast from Effort = ",
                     round(mean(abs(df$mean_t - df$lit_forecast_t)))))

dev.off()

############################### Online Appendix Figure 7 ###############################
png(filename = appfig7.outputfile,
    height = landscape.a4.height*png.scale, width = landscape.a4.width*png.scale,
    units = "in", res = 300)

### Set some plotting parameters (e.g. colors, symbols, sizes, and limits for x- and y-axes)
my.pch = c(19,17)
my.colors = c("blue","red")
my.xlim = c(1550,1950)
my.ylim = c(1000,2600)

my.cex.main = 1.7
my.cex.lab = 1.4
my.cex.axis = 1.6
my.cex.text = 1.2
my.cex = 1.3

my.arrow.length = 0.1

########### Expert/Literature versus Actual Effort #############

### Plot with Error Bars
matplot(x=df$mean_t,
        y=cbind(df$forecast_t, df$lit_citewt_forecast_t),
        xlim = my.xlim, ylim = my.ylim,
        main = "Expert Prediction/Literature versus Actual Effort\n(Literature weighted by Citations)",
        ylab = "Effort Predicted by Experts or Implied by Literature",
        xlab = "Actual Treatment Group Effort",
        cex = my.cex, cex.main = my.cex.main, cex.axis = my.cex.axis, cex.lab = my.cex.lab,
        pch = my.pch, col = my.colors)
abline(a=0,b=1,lty=2)

### "Manually" add arrows and labels for each of the different treatments
arrows(x0 = df$mean_t[1]+20, x1=df$mean_t[1],
       y0 = 0.5*df$forecast_t[1] + 0.5*df$lit_citewt_forecast_t[1], y1 = df$forecast_t[1],
       length = my.arrow.length)
arrows(x0 = df$mean_t[1]+20, x1=df$mean_t[1],
       y0 = 0.5*df$forecast_t[1] + 0.5*df$lit_citewt_forecast_t[1], y1 = df$lit_citewt_forecast_t[1],
       length = my.arrow.length)
text(x = df$mean_t[1]+20,
     y = 0.5*df$forecast_t[1] + 0.5*df$lit_citewt_forecast_t[1],
     labels = as.character(df$Treatment[1]), pos = 4, cex = my.cex.text)

arrows(x0 = df$mean_t[2]-20, x1=df$mean_t[2],
       y0 = 0.5*df$forecast_t[2] + 0.5*df$lit_citewt_forecast_t[2]-150, y1 = df$forecast_t[2],
       length = my.arrow.length)
arrows(x0 = df$mean_t[2]-20, x1=df$mean_t[2],
       y0 = 0.5*df$forecast_t[2] + 0.5*df$lit_citewt_forecast_t[2]-150, y1 = df$lit_citewt_forecast_t[2],
       length = my.arrow.length)
text(x = df$mean_t[2]-20,
     y = 0.5*df$forecast_t[2] + 0.5*df$lit_citewt_forecast_t[2]-150,
     labels = as.character(df$Treatment[2]), pos = 1, cex = my.cex.text)

arrows(x0 = df$mean_t[3]-100, x1=df$mean_t[3],
       y0 = 0.5*df$forecast_t[3] + 0.5*df$lit_citewt_forecast_t[3]+200, y1 = df$forecast_t[3],
       length = my.arrow.length)
arrows(x0 = df$mean_t[3]-100, x1=df$mean_t[3],
       y0 = 0.5*df$forecast_t[3] + 0.5*df$lit_citewt_forecast_t[3]+200, y1 = df$lit_citewt_forecast_t[3],
       length = my.arrow.length)
text(x = df$mean_t[3]-100,
     y = 0.5*df$forecast_t[3] + 0.5*df$lit_citewt_forecast_t[3]+200,
     labels = as.character(df$Treatment[3]), pos = 2, cex = my.cex.text)

arrows(x0 = df$mean_t[4]+20, x1=df$mean_t[4],
       y0 = 0.5*df$forecast_t[4] + 0.5*df$lit_citewt_forecast_t[4], y1 = df$forecast_t[4],
       length = my.arrow.length)
arrows(x0 = df$mean_t[4]+20, x1=df$mean_t[4],
       y0 = 0.5*df$forecast_t[4] + 0.5*df$lit_citewt_forecast_t[4], y1 = df$lit_citewt_forecast_t[4],
       length = my.arrow.length)
text(x = df$mean_t[4]+20,
     y = 0.5*df$forecast_t[4] + 0.5*df$lit_citewt_forecast_t[4],
     labels = as.character(df$Treatment[4]), pos = 4, cex = my.cex.text)

arrows(x0 = df$mean_t[5]-50, x1=df$mean_t[5],
       y0 = 0.5*df$forecast_t[5] + 0.5*df$lit_citewt_forecast_t[5] - 50, y1 = df$forecast_t[5],
       length = my.arrow.length)
arrows(x0 = df$mean_t[5]-50, x1=df$mean_t[5],
       y0 = 0.5*df$forecast_t[5] + 0.5*df$lit_citewt_forecast_t[5] - 50, y1 = df$lit_citewt_forecast_t[5],
       length = my.arrow.length)
text(x = df$mean_t[5]-50,
     y = 0.5*df$forecast_t[5] + 0.5*df$lit_citewt_forecast_t[5] - 50,
     labels = as.character(df$Treatment[5]), pos = 2, cex = my.cex.text)

arrows(x0 = df$mean_t[6]-20, x1=df$mean_t[6],
       y0 = 0.75*df$forecast_t[6] + 0.25*df$lit_citewt_forecast_t[6], y1 = df$forecast_t[6],
       length = my.arrow.length)
arrows(x0 = df$mean_t[6]-20, x1=df$mean_t[6],
       y0 = 0.75*df$forecast_t[6] + 0.25*df$lit_citewt_forecast_t[6], y1 = df$lit_citewt_forecast_t[6],
       length = my.arrow.length)
text(x = df$mean_t[6]-20,
     y = 0.75*df$forecast_t[6] + 0.25*df$lit_citewt_forecast_t[6],
     labels = as.character(df$Treatment[6]), pos = 2, cex = my.cex.text)

arrows(x0 = df$mean_t[7]+20, x1=df$mean_t[7],
       y0 = 0.5*df$forecast_t[7] + 0.5*df$lit_citewt_forecast_t[7], y1 = df$forecast_t[7],
       length = my.arrow.length)
arrows(x0 = df$mean_t[7]+20, x1=df$mean_t[7],
       y0 = 0.5*df$forecast_t[7] + 0.5*df$lit_citewt_forecast_t[7], y1 = df$lit_citewt_forecast_t[7],
       length = my.arrow.length)
text(x = df$mean_t[7]+20,
     y = 0.5*df$forecast_t[7] + 0.5*df$lit_citewt_forecast_t[7],
     labels = as.character(df$Treatment[7]), pos = 4, cex = my.cex.text)

arrows(x0 = df$mean_t[8]-50, x1=df$mean_t[8],
       y0 = 0.5*df$forecast_t[8] + 0.5*df$lit_citewt_forecast_t[8]+100, y1 = df$forecast_t[8],
       length = my.arrow.length)
arrows(x0 = df$mean_t[8]-50, x1=df$mean_t[8],
       y0 = 0.5*df$forecast_t[8] + 0.5*df$lit_citewt_forecast_t[8]+100, y1 = df$lit_citewt_forecast_t[8],
       length = my.arrow.length)
text(x = df$mean_t[8]-50,
     y = 0.5*df$forecast_t[8] + 0.5*df$lit_citewt_forecast_t[8]+100,
     labels = as.character(df$Treatment[8]), pos = 2, cex = my.cex.text)

legend("topleft", pch=my.pch, col=my.colors, cex = my.cex,
       legend=c("Expert Forecasts", "Effort Implied by Literature"))
errbar(x = df$mean_t, y = df$forecast_t, add=TRUE, type="n", errbar.col=my.colors[1],
       yplus = df$forecast_t + 2*df$se_forecast_t,
       yminus = df$forecast_t - 2*df$se_forecast_t)
errbar(x = df$mean_t, y = df$lit_citewt_forecast_t, add=TRUE, type="n", errbar.col=my.colors[2],
       yplus = df$lit_citewt_forecast_t + 2*df$se_lit_citewt_forecast_t,
       yminus = df$lit_citewt_forecast_t - 2*df$se_lit_citewt_forecast_t)

# "Manually" add text comparing experts' and "literature's" predictions to the truth
text(x = 1600, y = 1250, pos = 4, cex = my.cex.text + 0.1,
     labels = paste0("Mean Absolute Difference of\nExpert Forecast from Effort = ",
                     round(mean(abs(df$mean_t - df$forecast_t))),
                     "\n\nMean Absolute Difference of\nLiterature Forecast from Effort = ",
                     round(mean(abs(df$mean_t - df$lit_citewt_forecast_t)))))

dev.off()
