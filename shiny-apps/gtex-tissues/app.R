
#

library(shiny)
library(tidyverse)
library(data.table)


# Sample dataframe

gtex.meta.file <- '../../code/resources/GTEx/juncs/sampid-smts-smtsd-subjid.tsv'
gtex.meta <- fread(gtex.meta.file)

df <- gtex.meta[, .(Nsamp = uniqueN(SAMPID)), by = .(SMTS, SUBJID)][order(SMTS, -Nsamp)]



# Define UI
ui <- fluidPage(
  titlePanel("Tissue Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_tissue", "Select Tissue:", choices = unique(df$SMTS))
    ),
    mainPanel(
      tableOutput("table")
    )
  )
)

# Define server
server <- function(input, output) {
  output$table <- renderTable({
    filtered_df <- df %>%
      filter(SMTS == input$selected_tissue)
    return(filtered_df)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
