remove(list=ls())

library(RODBC)
library(dplyr)
library(forecast)
library(lubridate)
library("xts")
library(ggplot2)
library(hrbrthemes)
library(tseries)
library(fpp2)
library(forecast)
library("gridExtra")


setwd("d:/Studying/4 семестр/Аналіз даних/Курсова робота")

connStr<- "Driver={SQL Server};Server=LT-ALEXEI\\MSSQLSERVER01;Database=RestaurantVisitingDB"
db <- odbcDriverConnect(connStr, case="nochange")


import.to.Dim_DayOfTheWeek<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("dayOfTheWeekID"= rawData$dateID,
                       "dayOfTheWeekName"= rawData$value,stringsAsFactors = F )
  sqlSave(connection, rawData, tablename = "dbo.Dim_DayOfTheWeek", rownames = F, append=T)
}
import.to.Dim_DayOfTheWeek("day_ofTheWeek_info.csv", db)


import.to.Dim_Date<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("dateID"=rawData$dateID,
                       "dayValue"= rawData$calendar_date,
                       "dayOfTheWeekID"=rawData$day_of_the_week_id,
                       "isHoliday"=rawData$holiday_flg,
                       stringsAsFactors = F )
  sqlSave(connection, rawData, tablename = "dbo.Dim_Date", rownames = F, append=T)
}
update.to.Dim_Date<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("dateID"=rawData$dateID,
                       "dayValue"= rawData$calendar_date,
                       "dayOfTheWeekID"=rawData$day_of_week_id,
                       "isHoliday"=rawData$holiday_flg,
                       stringsAsFactors = F )
  sqlUpdate(connection, rawData, tablename = "dbo.Dim_Date")
}
import.to.Dim_Date("date_info_1.csv", db)
update.to.Dim_Date("UPD_date_info.csv", db)


import.to.Dim_Reservation<- function(path="air_reserve_1.csv", connection=db){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("reserveID"=rawData$reserveID,
                       "dateID"= rawData$dateID,
                       "timeValue"=rawData$time_value,
                       "resrvAmntOfVisitors"=rawData$reserve_visitors,
                       "placeID"=rawData$placeID,
                       stringsAsFactors = F )
  sqlSave(connection, rawData, tablename = "dbo.Dim_Reservation", rownames = F, append=T)
}
update.to.Dim_Reservation<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("reserveID"=rawData$reserveID,
                       "dateID"= rawData$dateID,
                       "timeValue"=rawData$time_value,
                       "resrvAmntOfVisitors"=rawData$reserve_visitors,
                       "placeID"=rawData$placeID,
                       stringsAsFactors = F )
  sqlUpdate(connection, rawData, tablename = "dbo.Dim_Reservation")
}
import.to.Dim_Reservation("air_reserve_1.csv", db)
update.to.Dim_Reservation("UPD_air_reserve.csv", db)


import.to.Dim_Cuisine<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("cuisineID"=rawData$cuisine_id,
                       "cuisineName"= rawData$air_genre_name,
                       "loadDate"=rep(Sys.time(), times= nrow(rawData)),
                       "loadEndDate"=rep(as.Date('9999-01-01'), times= nrow(rawData)),
                       stringsAsFactors = F )
  sqlSave(connection, rawData, tablename = "dbo.Dim_Cuisine", rownames = F, append=T)
}
update.to.dim_Cuisine <-function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("cuisineID"= rawData$cuisine_id,
                       "loadEndDate"=rep((Sys.time()), times= nrow(rawData)),
                       stringsAsFactors = F )
  sqlUpdate(connection, rawData, tablename = "dbo.Dim_Cuisine")
  
  x<-paste("SELECT TOP 1 cuisineID  FROM dbo.Dim_Cuisine ORDER BY cuisineID DESC",
           sep = "", collapse=NULL)
  x1<-sqlQuery(connection, x)
  x1<- as.integer(x1$cuisineID)
  
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("cuisineID"= rawData$cuisine_id,
                       "cuisineName"= rawData$air_genre_name,
                       "loadDate"= rep(Sys.time(), times= nrow(rawData)),
                       "loadEndDate"=rep(as.Date('9999-01-01'), times= nrow(rawData)),
                       stringsAsFactors = F )
  for (n in 1:nrow(rawData)) {
    x1<-x1+1
    rawData$cuisineID[n]<-x1
  }
  sqlSave(connection, rawData, tablename = "dbo.Dim_Cuisine", rownames = F, append=T)
  
  
}
import.to.Dim_Cuisine("air_cuisine_info.csv", db)
update.to.dim_Cuisine("UPD_air_cuisine_info.csv", db)


import.to.Dim_Place<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("placeID"=rawData$place_id,
                       "restaurantName"= rawData$air_area_name,
                       "cuisineID"=rawData$cuisine_id,
                       "latitude"=rawData$latitude,
                       "longitude"=rawData$longitude,
                       stringsAsFactors = F )
  sqlSave(connection, rawData, tablename = "dbo.Dim_Place", rownames = F, append=T)
}
update.to.Dim_Place<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("placeID"=rawData$place_id,
                       "restaurantName"= rawData$air_area_name,
                       "cuisineID"=rawData$cuisine_id,
                       "latitude"=rawData$latitude,
                       "longitude"=rawData$longitude,
                       stringsAsFactors = F )
  sqlUpdate(connection, rawData, tablename = "dbo.Dim_Place")
}
import.to.Dim_Place("air_place_info.csv", db)
update.to.Dim_Place("UPD_air_place_info.csv", db)


import.to.Fact_RestaurantVisiting<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("dateID"=rawData$dateID,
                       "placeID"= rawData$placeID,
                       "amountOfVisitors"=rawData$amountOfVisitors,
                       stringsAsFactors = F )
  sqlSave(connection, rawData, tablename = "dbo.Fact_RestaurantVisiting", rownames = F, append=T)
}
update.to.Fact_RestaurantVisiting<- function(path, connection){
  rawData <- read.csv(file= path, stringsAsFactors = F)
  rawData<- data.frame("dateID"=rawData$dateID,
                       "placeID"= rawData$placeID,
                       "amountOfVisitors"=rawData$amountOfVisitors,
                       stringsAsFactors = F )
  sqlUpdate(connection, rawData, tablename = "dbo.Fact_RestaurantVisiting")
}
import.to.Fact_RestaurantVisiting("air_visit_data_1.csv", db)
update.to.Fact_RestaurantVisiting("UPD_air_visit_data.csv", db)


#import from SQL SERVER
dataReserve<-sqlQuery(channel=db, query = "SELECT * FROM dbo.getFactAmountResrvVisitors('Hyogo-ken Kakogawa Kakogawacho Kitazaike')")
dataWholeAmount<-sqlQuery(channel=db, query = "SELECT * FROM dbo.getFactAmountVisitors('Hyogo-ken Kakogawa Kakogawacho Kitazaike')")




#All Functions

#Show By Week
getPlotByAllWeeks<-function(data,dayValue,sumValue,colnamesLabel, plotLabel){
  x<-data.frame(Date=as.Date(dayValue), Reserved=sumValue)
  x<- x %>% group_by(month=floor_date(Date, "week")) %>%
    summarize(Reserved=sum(Reserved))
  dfReserve<-data.frame(month=as.Date(x$month,
                                      origin = "1899-12-30"), Reserved=x$Reserved, stringsAsFactors = T)
  
  dfReserve$month<-as.factor( format(dfReserve$month,"%d/%m/%Y" ))
  dfReserve.xts<-as.xts(dfReserve[,-1], order.by = as.yearmon(as.character(dfReserve[,1]), 
                                                              "%d/%m/%Y"))
  colnames(dfReserve.xts)<-colnamesLabel
  attr(dfReserve.xts, 'frequency') <- 7
  df.ts <- as.ts(dfReserve.xts, start=1900 + .indexyear(dfReserve.xts)[1] +
                   .indexmon(dfReserve.xts)[1] / frequency(dfReserve.xts), names=colnames(dfReserve.xts))
  forecast::ggmonthplot(df.ts, xlab="Weeks", ylab=colnames(as.vector(dfReserve.xts)), main=plotLabel)
}

#Show By All Years
getPlotByAllYears<-function(data, dayValue, sumValue,colnamesLabel, plotLabel){
  x<-data.frame(Date=as.Date(dayValue), Reserved=sumValue)
  x<- x %>% group_by(month=floor_date(Date, "month")) %>%
    summarize(Reserved=sum(Reserved))
  
  dfReserve<-data.frame(month=as.Date(x$month,
                                      origin = "1899-12-30"), Reserved=x$Reserved, stringsAsFactors = T)
  dfReserve$month<-as.factor(format(dfReserve$month,"%d/%m/%Y" ))
  dfReserve.xts<<-as.xts(dfReserve[,-1], order.by = as.yearmon(as.character(dfReserve[,1]),
                                                               "%d/%m/%Y"))
  df.xts<<-dfReserve.xts
  colnames(dfReserve.xts)<-colnamesLabel
  df.ts <<- as.ts(dfReserve.xts, start=1900 + .indexyear(dfReserve.xts)[1] +
                    .indexmon(dfReserve.xts)[1] / frequency(dfReserve.xts), names=colnames(dfReserve.xts))
  forecast::ggseasonplot(df.ts, year.labels=TRUE, main=plotLabel)
}

#Show Scatter
getScatterplot<- function(data, dayValue, sumValue, scatterLabel){
  x<-data.frame(Date=as.Date(dayValue), Reserved=sumValue)
  x<- x %>% group_by(month=floor_date(Date, "20 days")) %>%
    summarize(Reserved=sum(Reserved))
  x %>%
    ggplot( aes(x=month, y=Reserved)) +
    geom_line( color="grey") +
    geom_point(shape=21, color="black", fill="#69b3a2", size=3) +
    theme_ipsum() +
    ggtitle(scatterLabel)
}

#Function to Plotting of Forecast Models with Confidence intervals.
Autoplot.Forecast <- function(x, y = NULL, testing=FALSE)
{
  # Load required parameters and packages
  stopifnot(is.forecast(x) | is.null(x) )
  if (!requireNamespace("forecast", quietly = TRUE)) 
    stop("forecast is needed for this function to work. Install it via install.packages(\"forecast\")", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
    stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")", call. = FALSE)
  
  if (length(x$method)>1)
    Title = "Ensemble of Prediction Models"
  else Title = x$method
  
  listing <- accuracy(x, y)     # In-sample & Out-sample Accuracy
  rownames(listing)[2] <- ifelse(testing, "Testing Set", "Validation set")
  cat("\nAccuracy of", Title, "\n")
  print(listing)
  
  p <- autoplot(x, main=Title, xlab="Year", ylab=colnames(df.xts)) + 
    geom_line(data = fortify(x$fitted, melt = FALSE), aes(x = x, y = y),
              na.rm = TRUE, col = "#0000AA", linetype = "dashed")
  print(p)
  
  if ( !is.null(y) ) {
    
    # Initialise ggplot object
    fcdata <- fortify(merge.xts(Forecast = as.xts(x$mean), Actual = y), melt = TRUE)
    Labels <- c("Forecast", "Measured")
    gptitle <- ifelse(testing, paste("Forecast of Testing Set by Final Method", Title), 
                      paste("Forecast of Validation Set by Method", Title))
    
    p <-  autoplot(x, include = 0)
    
    # Forecasted and Validationed points
    p <- p + geom_line(data = fcdata, aes_(x = ~Index, y = ~Value, color = ~Series, size = ~Series))
    
    p <- p + scale_colour_manual(name = "Data", labels = Labels,
                                 values=c("Forecast" = "#0000AA", "Actual" = "magenta")) +
      scale_size_manual(name = "Data", labels = Labels,
                        values = c("Forecast" = 0.5, "Actual" = 1)) +
      scale_x_yearmon(format = "%b-%Y", n = 5) +
      theme(legend.position = c(0, 1), legend.justification = c(0, 1)) +
      labs(title = gptitle, subtitle = Subtitle,
           x = "Months", y = colnames(df.xts)) +
      theme(plot.title = element_text(size = 12))
    print(p)
  }
}

# Analog of car::qqPlot for ggplot2
QQ.GGplot <- function(x, distribution = "norm", ..., line.estimate = NULL,
                      conf = 0.95, labels = names(x)){
  # from stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2/
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower,
                       labels[ord], "")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2], col = "red") +
    geom_ribbon(aes(ymin = lower, ymax=upper), fill="coral", alpha=0.15) +
    labs(title = "Quantiles of In-Sample & Normal Distribution",
         x = "Theoretical Quantiles", y = "In-Sample Quantiles")
  
  if(!is.null(labels)) p <- p + geom_text(aes(label = label),
                                          size = 3, vjust = 0, nudge_y = max(x)*0.03)
  print(p)
  # coef
}
#Function to Examination Residuals of time series forecasting.
Proof.Residuals <- function(x)
{
  # Load required parameters and packages
  stopifnot(is.forecast(x) | is.null(x))
  if (!requireNamespace("forecast", quietly = TRUE)) 
    stop("forecast is needed for this function to work. Install it via install.packages(\"forecast\")", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
    stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")", call. = FALSE)
  
  if (length(x$method)>1)
    Title = "Ensemble of Forecast Models"
  else Title = x$method
  
  # Plot some Charts of Residuals
  ggtsdisplay(residuals(x)[-1],main=paste("Residuals from",
                                          Title), ylab="", xlab="Months", plot.type="scatter")
  
  z <- as.numeric(residuals(x))
  names(z) <- as.yearmon(time(residuals(x)))
  QQ.GGplot(z[-1])
  
  # Test of Normality and AutoCorrelation of Residuals
  TestOfResiduals <- c(NormalityOfResiduals="",
                       AutoCorrelationOfResiduals="")
  print( shapiro.test(residuals(x)) )
  cat( "Normality of Residuals is", 
       shapiro.test(residuals(x))$p.value > 0.05, "\n" )
  
  TestOfResiduals$NormalityOfResiduals <-
    ifelse(shapiro.test(residuals(x))$p.value > 0.05, "Normality of Residuals is TRUE", "Normality of Residuals is FALSE")
  #"остатки распределены нормально","остатки не распределены нормально")
  
  print( Acf(residuals(x), lag.max=frequency(x$residuals)*2, plot=FALSE) )
  print( Pacf(residuals(x), lag.max=frequency(x$residuals)*2, plot=FALSE))
  
  # Has Autocorrelation of Residuals any Lags?
  print( Box.test(residuals(x), lag=24, type="Ljung") )
  cat( "AutoCorrelation is", Box.test(residuals(x),
                                      lag=frequency(x$residuals)*2, type="Ljung")$p.value < 0.05, "\n" )
  
  TestOfResiduals$AutoCorrelationOfResiduals <-
    ifelse(Box.test(residuals(x), lag=frequency(x$residuals)*2,
                    type="Ljung")$p.value < 0.05, "AutoCorrelation of Residuals is TRUE",
           "AutoCorrelation of Residuals is FALSE")
  return(TestOfResiduals)
}


#Analyzing  Reservation in Restaurant
getPlotByAllWeeks(dataReserve,
                  dataReserve$dayValue,
                  dataReserve$reservSum,
                  "Table Reserve Visitors, people",
                  "Reservation during 7 Days of Each Year ")

getPlotByAllYears(dataReserve,
                  dataReserve$dayValue,
                  dataReserve$reservSum,
                  "Table Reserve Visitors, people",
                  "Reservation during Months of All Years")

getScatterplot(dataReserve,
               dataReserve$dayValue,
               dataReserve$reservSum,
               "Evolution of reserving seats")


#Analyzing Attending Consumers in Restaurant
getPlotByAllWeeks(dataWholeAmount,
                  dataWholeAmount$dayValue,
                  dataWholeAmount$amountOfVisitors,
                  "Table Attending Consumers, people",
                  "Reservation during 7 Days of Each Year ")



getPlotByAllYears(dataWholeAmount,
                  dataWholeAmount$dayValue,
                  dataWholeAmount$amountOfVisitors,
                  "Table Attending Consumers, people",
                  "Reservation during Months of All Years"
                  )


getScatterplot(dataWholeAmount,
               dataWholeAmount$dayValue,
               dataWholeAmount$amountOfVisitors,
               "Evolution of attending consumers"
               )






#dividing data into 3 sets
train.ts <- window(df.ts, start=2016, end=2018 )
validation.ts <- window(df.ts, start=2017, end=2018 +4/frequency(df.ts))
test.ts <- window(df.ts, start=2018 , end=2018 +4/frequency(df.ts))
hf <- 15 # length(validation.ts)


Subtitle <- paste0("Building by ", round((tsp(train.ts)[2] -
                                            tsp(train.ts)[1]) * tsp(train.ts)[3] + 1, digits = 0), " months (",
                   zoo::yearmon(tsp(train.ts)[1]), " - ", zoo::yearmon(tsp(train.ts)[2]), ")")


#decomposing dataset into: data, trend, periodic, remainder
fit.stl = stl(train.ts, s.window="periodic", robust=TRUE)
autoplot(fit.stl)+ xlab("Year") +
  ggtitle("Classical periodical decomposition
    of reservation seats in restaurant")



#Forecasting train set
(fcast.stl <- forecast(fit.stl, h=hf, method="naive"))


Autoplot.Forecast(fcast.stl, validation.ts)




# Automatic selection of Box Cox transformation parameter for Forecasting
lam <- 0 

fit.ets <- ets(train.ts, lambda=lam, model="AAA")
fcast.ets <- forecast(fit.ets, h=hf, lambda=lam, level=c(80, 95)) 

res.ets <- Proof.Residuals(fcast.ets)



