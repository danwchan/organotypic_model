#' ---
#' title: "Old import chunk"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---
#' 
#/*    An old piece of code for importing CFU data into R objects, made into a chunk*/
#/*    Copyright (C) 2016  Daniel Chan*/

#/*    This program is free software: you can redistribute it and/or modify*/
#/*    it under the terms of the GNU General Public License as published by*/
#/*    the Free Software Foundation, either version 3 of the License, or*/
#/*    (at your option) any later version.*/

#/*    This program is distributed in the hope that it will be useful,*/
#/*    but WITHOUT ANY WARRANTY; without even the implied warranty of*/
#/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*/
#/*    GNU General Public License for more details.*/

#/*    You should have received a copy of the GNU General Public License*/
#/*    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

# ---- import-process-data ----
eval(parse(text=readLines(tablename, n=1))) #extract the information in the old testnames comment
expdata <- read.table(tablename, header=TRUE, sep=',', skip=1) #read in the tab delimited
timepoints <- levels(as.factor(expdata$Time))
datacols <- names(expdata)
for (i in 1:length(testnames)) {
  count <- paste('test', as.character(i), sep='')
  dilution <- paste('test', as.character(i), '.D', sep='')
  expdata[[testnames[i]]] <- expdata[[count]] * expdata[[dilution]] * 50
}
melted_data <- melt(expdata[c('Time', testnames)], id.vars='Time', na.rm=TRUE, variable.name = "sample_id", value.name = "cfu") #reshape the data
