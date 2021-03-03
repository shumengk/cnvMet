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
            geom_segment(aes(x=0, y=0, xend = 1, yend = 1),
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
        gene_aucValue <- round(mean(as.numeric(filter(aucValue,gene==geneIn)[2:6])),5)
        paste0("The AUC value of the model is: ",gene_aucValue)
    })
    
    
    ## CNV Branch
    # ROC plots
    output$cnvROC <- renderPlot({
        cnvInput <- input$segment
        udInput <- input$UD
        if (cnvInput=="DMR-overlapped"){
            cnvAUCValue <- "0.59682"
            cnvROC <- cnvROCPlot 
            ggplot(cnvROC,aes(fpr,tpr))+geom_line() + 
                geom_segment(aes(x=0,y=0,xend = 1, yend = 1),
                             linetype = 2,col='grey')+theme_bw()+
                ggtitle(cnvInput) 
        }
        else {
            
            lowerBound <- min(udInput[1]:udInput[2])
            upperBound <- max(udInput[1]:udInput[2])
            tmp <- cnvIndROC[grep(cnvInput,cnvIndROC$chrLoc),]
            tmp1 <- filter(tmp,tmp$udLoc>=lowerBound & tmp$udLoc<=upperBound)
            rm(tmp)
            tmp1$udLoc <- factor(tmp1$udLoc,levels=c(-2,-1,0,1,2),
                            labels=c("2 CNVs downstream","1 CNV downstream",
                                    "CNV","1 CNV upstream","2 CNVs upstream"))
            ggplot(tmp1,aes(fpr,tpr,colour=udLoc))+geom_line() + 
                geom_segment(aes(x=0,y=0,xend = 1, yend = 1),
                                 linetype = 2,col='grey')+theme_bw()+
                ggtitle(cnvInput)+labs(colour="Up/downstream") 
        }
    })

    output$cnvAUC <- renderTable({
        cnvInput <- input$segment
        udInput <- input$UD
        if (cnvInput=="DMR-overlapped"){
            cnvAUCValue <- "0.59682"
            tmp1 <- data.frame(matrix(nrow=1,ncol=2))
            names(tmp1) <- c("Location","Mean AUC")
            tmp1[1,] <- c(cnvInput,cnvAUCValue)
        }
        else {
            tmp <- cnvAUC[grep(cnvInput,cnvAUC$chrLoc),]
            lowerBound <- min(udInput[1]:udInput[2])
            upperBound <- max(udInput[1]:udInput[2])
            tmp1 <- filter(tmp,tmp$udLoc>=lowerBound & tmp$udLoc<=upperBound)

            cnvAUCValue <- round(rowMeans(tmp1[,2:6]),5)
            tmp1$`Mean AUC` <- cnvAUCValue
            names(tmp1)[1] <- "Location"
        }
        #paste0("The AUC value of the CNV-based classification model is ",cnvAUCValue)
        return(tmp1[,c(1,ncol(tmp1))])
    })
})



