#/*    R function to take a numeric outpit in scientific notation and translate it into an expression for use as a geom_text label*/
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

numeric_to_label <- function (numeric, prefix = "",
                              precision = 3,
                              round_method = floor
                              ) {
  value <- as.integer(numeric)
  exponent <- round_method(log10(numeric))
  number <- signif(numeric * 10^(exponent * -1), digits = precision)
  label <- paste0(prefix, number, "%*%10^{", exponent, "}")
  return(label)
}