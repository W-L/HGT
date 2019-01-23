library(shiny)
library(ggplot2)
library(data.table)
library(tidyr)
library(ape)
library(fitdistrplus)

#setwd("/project/ngs_marsico/HGT/report/")

#data for oma full results
taxa<-fread("data/taxon_table",fill=T)
taxa_ngene<-taxa[order(-ngene),]
taxa_mean_dist<-taxa[order(mean),]

#data for oma - ordered by legionella genome
table_order<-fread("data/table_LEGPC_genome_order_add")
table_order$V5<-as.factor(table_order$V5)
table_order$V1<-as.factor(table_order$V1)

#data for oma standalone
taxa_info_50<-fread("data/50_taxa_table", header = F)
table_order50<-fread("data/table_LEGPC_genome_order_50")
table_order50$V5<-as.factor(table_order50$V5)
table_order50$V1<-as.factor(table_order50$V1)

taxa_50<-fread("data/taxon_table_50", fill=T)


#data for rbb plot
table_path<-"data/rbb_table_form"
#additional data for marking the oma results
orthos<-read.csv("data/oma50_1to1_ortho_mimi_legpc", header=F, sep=" ")

rbb<-fread(table_path, sep=' ')
rbb$V1<-as.factor(rbb$V1)
names(rbb)<-c("orth","score","eval","ident","pos","gap")

score_norm<-rbb$score/max(rbb$score)
pre_eval_norm<--(rbb$eval-max(rbb$eval))
eval_norm<-pre_eval_norm/max(pre_eval_norm)

rbb_new<-data.frame("orth"=rbb$orth,"score"=score_norm,"eval"=eval_norm)
rbb_long <- gather(rbb_new, metric, val, score:eval, factor_key=TRUE)
rbb_order<-rbb[order(-score),]

#create comparison table
comp<-data.frame()
idx_vec<-c()

for (i in orthos$V1){
  idx<-match(i,rbb_order$orth)
  prot<-rbb_order[idx,]
  idx_vec<-c(idx_vec,idx)
  comp<-rbind(comp,prot)
}
comp<-cbind(comp,idx_vec)
comp_order<-comp[order(comp$idx_vec),]

#add logic for labelling the rbb oma comparison
lab <- rbb_long$orth %in% comp_order$orth
rbb_long$val1<-rbb_long$val
rbb_long$val1<-as.character(rbb_long$val1)
rbb_long$val1[!lab] <- NA
rbb_long$val1[lab] <- "*"



#data for kolmogorov-smirnoff tests and protein info
dist_sort<-fread("data/distance_matrix_ks_testing",header =T)
dist_sort50<-fread("data/distance_matrix_ks_testing50",header =T)
proteins<-fread("data/oma_ncbi_LEGPC",header = F)

gen_links<- function(ID){
  link<-paste("http://www.genome.jp/dbget-bin/www_bget?lpc:",ID,sep = "")
  return(link)
}
kegg_links<-sapply(proteins$V2,gen_links)
proteins<-cbind(proteins,kegg_links)
colnames(proteins)<-c("OMA_ID","NCBI_ID","loc","pos","Desc.","KEGG-link")

#load trees from RAxML on the 50 species dataset
trees_50<-fread("data/RAxML_trees_50_species",header = F, sep='\t')

#data for model fit
c5<-fread("data/model_fit_closest_proteins",header = T)

#candidate
cand<-fread("data/cand_prot",header=T)
cand<-cand[,1:5]

# Define UI with navbarPage for top tabs
ui <- navbarPage(
  
  navbarPage("HGT - Internship project",
             tabPanel("Summary",includeMarkdown("mds/summary.md")),
             
             tabPanel("GC-content",includeMarkdown("mds/gc.md")),
             
             tabPanel("Mummer - dotplots",includeMarkdown("mds/mummer.md")),
             
             tabPanel("OMA full - analysis", includeMarkdown("mds/oma_full.md")),
             
             tabPanel("OMA full - results, total", includeMarkdown("mds/oma_full_res.md"),
                      
                      plotOutput("distPlot1", width = 2000, dblclick = "distPlot1_dblclick",click="distPlot1_click", brush = brushOpts(
                        id = "distPlot1_brush",resetOnNew = TRUE)),
                      
                      plotOutput("distPlot2", width = 2000, dblclick = "distPlot2_dblclick",click="distPlot2_click", brush = brushOpts(
                        id = "distPlot2_brush",resetOnNew = TRUE)),
                      
                      plotOutput("distPlot3", width = 2000, dblclick = "distPlot3_dblclick",click="distPlot3_click", brush = brushOpts(
                        id = "distPlot3_brush",resetOnNew = TRUE)),
                      
                      plotOutput("distPlot4", width = 2000, dblclick = "distPlot4_dblclick",click="distPlot4_click", brush = brushOpts(
                        id = "distPlot4_brush",resetOnNew = TRUE)),
                      
                      dataTableOutput('table')
             ),
             
             tabPanel("OMA full - results, per species", includeMarkdown("mds/oma_tax_select.md"),
                      
                      fluidRow(
                        
                        column(8,plotOutput("orderedPlot", dblclick = "dblclick",click="click", brush = brushOpts(
                          id = "brush",resetOnNew = TRUE))),
                        
                        column(4,plotOutput("taxPlot"))
                      ),
                      
                      selectInput("select", label = h3("Select taxon"),choices = levels(table_order$V1), selected = 1),
                      sliderInput("slider1", h3("Binsize"),min = 0, max = 1, value = 0.5)
             ),
             
             tabPanel("OMA standalone - analysis", includeMarkdown("mds/oma_standalone.md"),
                      dataTableOutput('table1')
             ),
             
             tabPanel("OMA standalone - results, total", includeMarkdown("mds/oma_standalone_res.md"),
                      dataTableOutput('table2')
             ),
             
             tabPanel("OMA standalone - results, per species", includeMarkdown("mds/oma_tax_select_50.md"),
                      
                      fluidRow(
                        
                        column(8,plotOutput("orderedPlot50", dblclick = "dblclick50",click="click50", brush = brushOpts(
                          id = "brush50",resetOnNew = TRUE))),
                        
                        column(4,plotOutput("taxPlot50"))
                      ),
                      
                      selectInput("select50", label = h3("Select taxon"),choices = levels(table_order50$V1), selected = 1),
                      sliderInput("slider50", h3("Binsize"),min = 0, max = 0.2, value = 0.1)
             ),
             
             tabPanel("Rbb - blasting", includeMarkdown("mds/rbb.md"),
                      
                      plotOutput("orderedPlot_rbb", width ="100%", height ="800px", dblclick = "dblclick_rbb",click="click_rbb", brush = brushOpts(
                        id = "brush_rbb",resetOnNew = TRUE)),
                      
                      dataTableOutput('table_rbb')
             ),
             
             tabPanel("Kolmogorov-Smirnov", includeMarkdown("mds/ks.md"),
                      
                      selectInput("ks_select",label=h3("Select Leg. Protein"),choices=colnames(dist_sort),selected = "LEGPC00001"),
                      dataTableOutput('protein_table'),
                      plotOutput("ks_plot",width = "100%",height="800px"),
                      dataTableOutput('ks_table')
             ),
             
             tabPanel("Distance trees",includeMarkdown("mds/trees.md"),
                      
                      selectInput("ks_select_tree",label=h3("Select Leg. Protein"),choices=colnames(dist_sort),selected = "LEGPC00001"),
                      dataTableOutput('protein_table_tree'),
                      plotOutput("tree50",width = "100%",height="800px"),
                      dataTableOutput('ks_table50')
                      ),
             
             tabPanel("Groups with 2 or 3 proteins", includeMarkdown("mds/small_groups.md")
                      
             ),
             
             tabPanel("Model fit - closest proteins", includeMarkdown("mds/model_fit.md"),
                      
                      selectInput("select_closest_species", label = h3("Select taxon"),choices = levels(table_order$V1), selected = 1),
                      dataTableOutput('table_closest_species'),
                      selectInput("select_closest_protein",label=h3("Select Leg. Protein"),choices=colnames(dist_sort),selected = "LEGPC00001"),
                      dataTableOutput('table_closest_protein')
             ),
             
             
             tabPanel("Interesting Genes",includeMarkdown("mds/candidate_genes.md"),
                      dataTableOutput('table_cand')
                      ),
             
             
             tabPanel("Species Info",
                      selectInput("select_info", label = h3("Select taxon"),choices = levels(table_order$V1), selected = 1),
                      dataTableOutput('table_info')
                      ),
             
             tabPanel("Papers, Software etc.", includeMarkdown("mds/papers_software.md"))
  )
  
  
)


# Define server logic
server<- function(input, output) {
  ranges <- reactiveValues(x = NULL, y = NULL)
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  ranges3 <- reactiveValues(x = NULL, y = NULL)
  ranges4 <- reactiveValues(x = NULL, y = NULL)
  
  output$distPlot1 <- renderPlot({
    ggplot(taxa_ngene, aes(x=taxid, y=ngene, fill=kingdom)) + 
      geom_bar(stat="identity") +
      scale_x_discrete(limits= taxa_ngene$taxid)+
      xlab("Taxon ID")+
      ylab("Number of Orthologs with Leg. pn.c.") +
      ggtitle("Number of Orthologs per taxon - ordered by number of orthologs") +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  
  output$distPlot2 <- renderPlot({
    ggplot(taxa_mean_dist, aes(x=taxid, y=ngene, fill=kingdom)) +
      geom_bar(stat="identity") +
      scale_x_discrete(limits= taxa_mean_dist$taxid)+
      xlab("Taxon ID")+
      ylab("Number of Orthologs with Leg. pn.c.") +
      ggtitle("Number of Orthologs per taxon - ordered by mean distance") +
      coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
  })
  
  output$distPlot3 <- renderPlot({
    ggplot(taxa_mean_dist, aes(x=taxid, y=mean, colour=kingdom)) +
      geom_point(size=.5)+
      scale_x_discrete(limits=taxa_mean_dist$taxid) +
      xlab("Taxon ID") +
      ylab("Mean distance to Leg. pn. c.") +
      ggtitle("Mean distance per taxon - ordered by mean distance") +
      coord_cartesian(xlim = ranges3$x, ylim = ranges3$y, expand = FALSE)
  })
  
  output$distPlot4 <- renderPlot({
    ggplot(taxa_ngene, aes(x=taxid, y=mean, colour=kingdom)) +
      geom_point(size=.5)+
      scale_x_discrete(limits=taxa_ngene$taxid) +
      xlab("Taxon ID") +
      ylab("Mean distance to Leg. pn. c.") +
      ggtitle("Mean distance per taxon - ordered by number of orthologs") +
      coord_cartesian(xlim = ranges4$x, ylim = ranges4$y, expand = FALSE)
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$distPlot1_dblclick, {
    brush <- input$distPlot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  observeEvent(input$distPlot2_dblclick, {
    brush <- input$distPlot2_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  observeEvent(input$distPlot3_dblclick, {
    brush <- input$distPlot3_brush
    if (!is.null(brush)) {
      ranges3$x <- c(brush$xmin, brush$xmax)
      ranges3$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges3$x <- NULL
      ranges3$y <- NULL
    }
  })
  
  observeEvent(input$distPlot4_dblclick, {
    brush <- input$distPlot4_brush
    if (!is.null(brush)) {
      ranges4$x <- c(brush$xmin, brush$xmax)
      ranges4$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges4$x <- NULL
      ranges4$y <- NULL
    }
  })
  
  output$table <- renderDataTable(taxa)
  
  #ordered plots
  
  dataset0 <- reactive({
    subset(table_order, V1==input$select)
  })
  
  ranges0 <- reactiveValues(x = NULL, y = NULL)
  
  output$orderedPlot <- renderPlot({
    ggplot(dataset0(), aes(x=V5, y=V4)) +
      geom_point(size=.5)+
      geom_hline(yintercept=mean(dataset0()$V4),lty="dashed",color="darkred") +
      xlab("Legionella ortholog (ordered)") +
      ylab("Distance between orthologs") +
      ggtitle(paste("Distances between all orthologs ordered by genome position on LEGPC")) +
      scale_colour_manual(values = "darkblue") +
      #theme(axis.text.x = element_text(angle = 90)) +
      coord_cartesian(xlim = ranges0$x, ylim = ranges0$y, expand = FALSE)
  })
  
  output$taxPlot <- renderPlot({
    ggplot(dataset0(), aes(x=V4)) + 
      geom_histogram(binwidth = input$slider1, color="black") +
      xlab("Distance inferred by protdist") +
      ylab("Freq.") +
      geom_vline(xintercept = mean(dataset0()$V4),lty="dashed",color="darkred") +
      ggtitle("Distribution of distances for selected species to LEGPC") +
      scale_x_continuous(position = "top") +
      coord_flip()
    
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$dblclick, {
    brush <- input$brush
    if (!is.null(brush)) {
      ranges0$x <- c(brush$xmin, brush$xmax)
      ranges0$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges0$x <- NULL
      ranges0$y <- NULL
    }
  })
  
  # oma standalone
  
  output$table1 <- renderDataTable(taxa_info_50)
  
  output$table2 <- renderDataTable(taxa_50)
  
  #standalone selector
  
  dataset50 <- reactive({
    subset(table_order50, V1==input$select50)
  })
  
  ranges50 <- reactiveValues(x = NULL, y = NULL)
  
  output$orderedPlot50 <- renderPlot({
    ggplot(dataset50(), aes(x=V5, y=V4)) +
      geom_point(size=.5)+
      geom_hline(yintercept=mean(dataset50()$V4),lty="dashed",color="darkred") +
      xlab("Legionella ortholog (ordered)") +
      ylab("Distance between orthologs") +
      ggtitle(paste("Distances between all orthologs ordered by genome position on LEGPC")) +
      scale_colour_manual(values = "darkblue") +
      #theme(axis.text.x = element_text(angle = 90)) +
      coord_cartesian(xlim = ranges50$x, ylim = ranges50$y, expand = FALSE)
  })
  
  output$taxPlot50 <- renderPlot({
    ggplot(dataset50(), aes(x=V4)) + 
      geom_histogram(binwidth = input$slider50, color="black") +
      xlab("Distance inferred by protdist") +
      ylab("Freq.") +
      geom_vline(xintercept = mean(dataset50()$V4),lty="dashed",color="darkred") +
      ggtitle("Distribution of distances for selected species to LEGPC") +
      scale_x_continuous(position = "top") +
      coord_flip()
    
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$dblclick50, {
    brush <- input$brush50
    if (!is.null(brush)) {
      ranges50$x <- c(brush$xmin, brush$xmax)
      ranges50$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges50$x <- NULL
      ranges50$y <- NULL
    }
  })
  
  
  #rbb plots
  
  ranges_rbb <- reactiveValues(x = NULL, y = NULL)
  
  output$orderedPlot_rbb <- renderPlot({
    ggplot(rbb_long, aes(x=orth,y=log10(val),color=metric)) +
      geom_point(size=1)+
      scale_x_discrete(limits= rbb_order$orth) +
      scale_color_manual("",values = c("darkred","darkblue","darkgreen","darkgoldenrod","brown")) +
      xlab("Reciprocal best blast pair") +
      ylab("(a.u.)") +
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=rel(1))) +
      geom_hline(yintercept=mean(log10(rbb_new$score)),lty="dashed",color="darkgreen") +
      geom_hline(yintercept=median(log10(rbb_new$score)),lty="dashed",color="darkgoldenrod") +
      geom_text(aes(label=val1),color="black",vjust=-1,size=7) +
      coord_cartesian(xlim = ranges_rbb$x, ylim = ranges_rbb$y, expand = FALSE)
  })
  
  dataset_rbb <- reactive({
    rbb_order[ranges_rbb$x[1]:ranges_rbb$x[2]]
  })
  
  output$table_rbb <- renderDataTable(dataset_rbb())
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$dblclick_rbb, {
    brush <- input$brush_rbb
    if (!is.null(brush)) {
      ranges_rbb$x <- c(brush$xmin, brush$xmax)
      ranges_rbb$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_rbb$x <- NULL
      ranges_rbb$y <- NULL
    }
  })
  
  #KS PLOTTING
  
  col_ks <- reactive({
    subset(dist_sort, select = input$ks_select)
  })
  
  dataset_ks <- reactive({
    na.omit(cbind(sp=dist_sort$species,means=dist_sort$species_mean,prot=col_ks()))
  })
  
  dataset_protein <- reactive({
    subset(proteins, OMA_ID==input$ks_select)
  })
  
  output$ks_table<- renderDataTable(dataset_ks())
  
  output$protein_table <- renderDataTable(dataset_protein())
  
  output$ks_plot <- renderPlot({
    ks<-ks.test(dataset_ks()$means,dataset_ks()$prot)
    n.ortho<-nrow(dataset_ks())
    ypos<-max(dataset_ks()$means,dataset_ks()$prot)
    
    plot(dataset_ks()$means,type = "l",col="darkred",ylab="Distance",xlab="Species - sorted by mean distance",
         ylim=c(min(dataset_ks()$means,dataset_ks()$prot),max(dataset_ks()$means,dataset_ks()$prot)),
         xlim=c(0,n.ortho))
    points(dataset_ks()$prot,pch=19,col="darkblue",cex=0.5)
    text(x=1,y=ypos-5,labels = paste(ks$method,"\n","D =",round(ks$statistic,2),"\n","p.val =",ks$p.value,
                                     "\n","n.ortho =",n.ortho),pos=4)
  })
  
  #Tree plotting
  
  newicktree<- reactive({
    subset(trees_50, V1==input$ks_select_tree)
  })
  
  output$tree50 <- renderPlot({
    phylo<-read.tree(text=newicktree()$V2)
    plot(phylo,adj=1)
    edgelabels(round(phylo$edge.length,2),adj = c(0.5, 0.5), bg="darkblue",col="white")
    nodelabels(phylo$node.label,bg="darkred",col="white")
    
  })
  
  col_ks50 <- reactive({
    subset(dist_sort50, select = input$ks_select_tree)
  })
  
  dataset_ks50 <- reactive({
    na.omit(cbind(sp=dist_sort50$species,means=dist_sort50$species_mean,prot=col_ks50()))
  })
  
  dataset_protein50 <- reactive({
    subset(proteins, OMA_ID==input$ks_select_tree)
  })
  
  output$ks_table50<- renderDataTable(dataset_ks50())
  
  output$protein_table_tree <- renderDataTable(dataset_protein50())
  
  
  #model fit - closest proteins
  
  dataset_closest_species <- reactive({
    c5[orthos %like% input$select_closest_species]
  })
  
  dataset_closest_protein <- reactive({
    subset(proteins, OMA_ID==input$select_closest_protein)
  })
  
  output$table_closest_species<- renderDataTable(dataset_closest_species())
  
  output$table_closest_protein<- renderDataTable(dataset_closest_protein())
  
  # taxa info
  
  dataset_info <- reactive({
    subset(taxa,taxid==input$select_info, select = c("taxid","kingdom","division","genus","species"))
  })
  
  output$table_info<- renderDataTable(dataset_info())
  
  #candidates
  
  output$table_cand<- renderDataTable(cand)
  
  
}

runApp(list(ui = ui, server = server),launch.browser = T)
