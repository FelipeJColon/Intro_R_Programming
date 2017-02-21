
## ------------------------
## The Working Directory
## ------------------------
# Typically you set up a working directory before analysing any data

# Get working directory
getwd()

# Modify the working directory
setwd("C:/Users/zck07opu/Documents/")

# Explore the working directory
dir(); list.files()

# Create a new subdirectory
dir.create("RCourse")
list.files()

## ----------------------------------
## 2. Installing and Loading Packages
## ----------------------------------

# R Packages
# A package is a collection of functions for specific tasks

# Check for installed packages
installed.packages()

# Installing packages within the AD could be a real pain!
# so we will work a way around it

# Modify path to library
.libPaths(file.path(getwd(), "RCourse"))

# Install packages on the working directory
install.packages("dplyr", repos="http://cran.rstudio.com/")

# Load package 
require(dplyr)

# Load data
data("ToothGrowth")

# Unload package
detach("package:dplyr", unload=TRUE)

## ------------------------------------
## 3. Reading and Writing Tabular Data
## ------------------------------------

# Load a dummy dataset 

# Create string for url in two steps
url1  <- 'https://raw.githubusercontent.com/FelipeJColon/'
url2  <- 'Intro_R_Programming/master/mydata.csv'

# Concatenate strings
myUrl <- paste0(url1, url2)

# Read remote data
mydata <- read.csv(myUrl)


##--- Reading large sets 

# Load Dummy data 2 
url2  <- "Intro_R_Programming/master/traffic.csv"

# Concatenate strings
myUrl <- paste0(url1, url2)

# create sample
sample   <- read.csv(myUrl, nrows=100)

# Read classes
classes  <- sapply(sample, class)

# Read data
myData2  <- read.csv(myUrl, colClasses=classes)

# Writing data

dfr <- data.frame(foo=1:4, bar=letters[1:4]) # Create data.frame

write.csv(x=dfr, file="mydata1.csv")
write.csv(x=dfr, file="mydata2.csv", row.names=FALSE) # are they different?

## ---------------------
## 4. Dates and Times
## ---------------------     

# Dates are considered by R as the number of **days** since 1 Jan 1970
# Typically, dates are written using numbers as follows "YYYY-MM-DD" 
# (e.g. 1970-01-01)

# Dates may be coerced from character strings using `as.Date()` 
(mydate <- as.Date("2017-02-21"))

# get DoW
weekdays(mydate)

# Get month
months(mydate)

# Get quarter
quarters(mydate)

# Get year
format(mydate, "%Y")        # No generic function for Year


# Convert dates from other formats into R-friendly
as.Date("02/27/92", "%m/%d/%y")

as.Date("1jan1960", "%d%b%Y")

strptime("02/27/92", "%m/%d/%y") # Check ?strptime for details

strptime("December 31, 2015 23:40:59", "%B %d, %Y %H:%M:%S")

## ----------------------------------------
## Part III. Exploratory Data Analysis
## ----------------------------------------

# Let's work with the `mtcars` set (*Motor Trend Car Road Tests*) 
# from the `datasets` package. Load the data and look at the 
# dimensions and structure of the set. Explore the description 
# of the variables in the set

# Load packages
require(datasets)

# Load dataset
data(mtcars)

# Explore variables in the set
?mtcars


# You can explore the top and bottom six rows using `head()` and `tail()`
head(mtcars)
tail(mtcars)

head(mtcars, n=10) 

# Summary statistics
summary(mtcars)

# Frequency tables
table(mtcars$cyl)  # Number of cars by number of cylinders

with(mtcars, table(cyl, gear))  # Cars by cylinders and gears

# 3 way Frequency tables
ftable(mtcars[c("cyl", "am", "gear")])


# Frequency tables with missing entries
mtcars[1,2] <- NA
table(mtcars$cyl, useNA="always")

# Quantiles for specific probabilities 
quantile(mtcars$mpg) # Default
quantile(mtcars$mpg, probs=c(0.1, 0.5, 0.9))

## ---------------------
## Exploratory plots
## ---------------------

# Scatter plot
with(mtcars, plot(wt, mpg)) 

# Line plot
data(AirPassengers)
plot(AirPassengers)

# Customizing plots
with(mtcars, plot(wt, mpg, pch=19, col=factor(am))) # ?pch

with(mtcars, plot(wt, mpg, pch=19, col=factor(am), 
                  xlab="Weight", ylab="Miles per gallon"))


# Change the size of the points (equally or based on a condition)
with(mtcars, plot(wt, mpg, pch=19, col=factor(am), 
                  xlab="Weight", ylab="Miles per gallon",
                  cex=carb))


# Add a title 
with(mtcars, plot(wt, mpg, pch=19, col=factor(am), 
                  xlab="Weight", ylab="Miles per gallon",
                  cex=carb, main="My first plot"))


# And a legend!
with(mtcars, plot(wt, mpg, pch=19, col=factor(am), 
                  xlab="Weight", ylab="Miles per gallon",
                  cex=carb, main="My first plot"))
legend("topright", legend=unique(mtcars$am), pch=19,
       col=c("black", "red"))



# One-dimensional box and whisker 
data("ToothGrowth")
boxplot(ToothGrowth$len, col="blue")

# Two-dimensional box and whisker plots
with(ToothGrowth, boxplot(len ~ supp, col=c("orange","blue")))

# Histogram
x <- rnorm(1000)
hist(x)

# Histogram
hist(x, breaks=50, col="green")

# Bar plots
set.seed(234)
x <- rnorm(100)
barplot(x, main="Random numbers")

# Scatter plot matrices
pairs(mtcars[,3:6])

# Overlaying features
hist(x, breaks=15, col="steelblue4")
abline(v=median(x), col="tomato", lwd=6)

# Overlaying features - Horizontal dotted line
with(ToothGrowth, boxplot(len ~ supp, 
                          col=c("goldenrod","darkolivegreen")))
abline(h=28, col="red3", lwd=3, lty=3)

# Conditional colours
set.seed(234); x <- rnorm(100)
barplot(x, main="Random numbers", col=ifelse(x < 0, "blue", "red"))




# Multiple plots in one graphic device
par(mfrow=c(1, 2)) # One device, one row, two columns
with(mtcars, {
     plot(wt, mpg, main="Weight vs. miles per gallon")
     plot(hp, mpg, main="Horsepower vs. miles per gallon")
})

# Save plots to a file
pdf("myplot.pdf")    # Opens graphics device in working directory
with(mtcars, plot(wt, mpg, pch=19, col=cyl, cex=2))
dev.off()            # Closes graphics device

##---------------------
## Your turn!
##---------------------

# - Load the `mydata.csv` file
# - Compute **summary statistics** for both sulfate and nitrate 
# measurements (are there any missing values?)
# - Create a subset called `mySubset` with no missing observations 
# - Compute a frequency table for each ID level on your subset
# - Set the dates in the dataset as R-dates
# - Generate a 2 panel plot with the histograms of sulfate and nitrate
# - Create a subset where all negative sulfate values are equal to zero and 
# plot sulfate values with colours by ID


## To practice at your own pace

# - Install the interactive `swirl` package
# - Go through the `R programming` modules


