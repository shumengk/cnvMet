library(shiny)
library(ggplot2)
library(dplyr)

# Define server logic required to output the ROC plots and corresponding AUC values
shinyServer(function(input, output) {
    
    ## Methylation branch
    # ROC Plots of (1) all COAD samples vs. normal samples or 
    # (2) COAD tumor tissue vs. COAD normal tissues
    output$rocPlot <- renderPlot({
        geneIn <- input$metGene
        if (input$coadOnly=="COAD vs. Normal"){
            rocPlot <- geneROCPlot 
        }
        else if (input$coadOnly=="COAD tumor tissues vs. normal tissues"){
            rocPlot <- coadOnlyROC
        }
        
        tmp <- filter(rocPlot, gene==geneIn)
        ggplot(tmp,aes(fpr,tpr))+geom_line() + 
            geom_segment(aes(x=0,y=0,xend = 1, yend = 1),
                         linetype = 2,col='grey')+theme_bw()+
            ggtitle(geneIn)
    })
    
    # Corresponding AUC values
    output$auc <- renderText({
        geneIn <- input$metGene
        if (input$coadOnly=="COAD vs. Normal"){
            aucValue <- geneAUC
        }
        else if (input$coadOnly=="COAD tumor tissues vs. normal tissues"){
            aucValue <- coadOnlyAUC 
        }
        gene_aucValue <- round(mean(unlist(filter(aucValue,gene==geneIn)[2:6])),5)
        paste0("The AUC value of the model is: ",gene_aucValue)
    })
    
    
    ## CNV Branch
    # ROC plots
    output$cnvROC <- renderPlot({
        cnvInput <- input$segment
        if (cnvInput=="DMR-overlapped"){
            cnvAUCValue <- "0.59682"
            cnvROC <- cnvROCPlot 
            
        }
        else {
            cnvROC <- filter(cnvIndROC,chrLoc==cnvInput)
        }
        
        ggplot(cnvROC,aes(fpr,tpr))+geom_line() + 
            geom_segment(aes(x=0,y=0,xend = 1, yend = 1),
                         linetype = 2,col='grey')+theme_bw()+
            ggtitle(cnvInput) 
    })
    
    output$cnvAUC <- renderText({
        cnvInput <- input$segment
        if (cnvInput=="DMR-overlapped"){
            cnvAUCValue <- "0.59682"
        }
        else {
            cnvAUCValue <- round(mean(unlist(filter(cnvAUC,chrLoc==cnvInput)[2:6])),5)
        }
        paste0("The AUC value of the CNV-based classification model is ",cnvAUCValue)
    })
    
    
})
