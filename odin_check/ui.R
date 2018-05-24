#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggmap)

library(readr)
library(RJSONIO)
library(curl)
library(base64enc)

pageWithSidebar(
  headerPanel('Summary of ODIN'),
  sidebarPanel('Not reporting ODIN',
    DT::dataTableOutput("table")
  ),
  mainPanel(
    plotOutput('plot1')
  )
)
