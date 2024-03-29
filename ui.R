library(shiny)
library(shinydashboard)

dashboardPage(
    dashboardHeader(title = "CNV-met Analysis Results"),
    
    dashboardSidebar(
        radioButtons(
            inputId="metOrCNV",label="Select the primary data source to build the classification model:",
            choices=c("CNV","DMR"), selected="met", width="90%"),
        conditionalPanel(
            condition = "input.metOrCNV=='DMR'",
            selectInput("metGene",
                        "Select a potential DMR-gene:",
                        choices=geneList
            ),
            selectInput("coadOnly",
                        "Select a classification model:",
                        choices=c("COAD vs. Normal","COAD tumor tissues vs. normal tissues")
            )
        ),
        conditionalPanel(
            condition = "input.metOrCNV=='CNV'",
            selectInput("segment",
                        "Select a CNV segment:",
                        choices=cnvList,
                        selected=NULL
            ),
            conditionalPanel(
                condition = "input.metOrCNV=='CNV'",
                sliderInput("UD","How many CNV down/upstream?",
                            label = div(style='width:200px;', 
                            div(style='float:left;', 'Downstream'), 
                            div(style='float:right;', 'Upstream')), 
                            min=-2,max=2,value=c(-2,2),step=1)
            )
        )
    ),
    
    dashboardBody(
        conditionalPanel(
            condition = "input.metOrCNV=='DMR'",
            h5("ROC Plot of Gradient Boosting Machine COAD classification models 
               on methylation data of the selected gene"),
            h6("COAD vs. Normal classifies all COAD from Normal samples."),
            h6("COAD tumor tissues vs. normal tissues only looks at COAD samples and 
               classifies tumor tissues from normal tissues."),
            plotOutput("rocPlot"),
            textOutput("auc")
        ),
        conditionalPanel(
            condition = "input.metOrCNV=='CNV'",
            h5("ROC Plot of Gradient Boosting Machine COAD tumor tissue classification models 
               on CNV-enriched data"),
            h6("CNV location based (chr_start_end): The CNVs were seleced through rank of importance 
            in COAD classification model on the 100 most common masked CNVs, 
               the methylation data at the corresponding location was extracted to build the model."),
            h6("CNV up/downstream: The methylation data at n CNV-length up/downstream regions of the CNVs 
               was used for model building."),
            h6("DMR-overlapped: The CNVs data was first overlapped with 18 DMR genes, then the top 
               50 most common CNVs were selected for modeling."),
            plotOutput("cnvROC"),
            tableOutput("cnvAUC")
        ),
       
    )
)
