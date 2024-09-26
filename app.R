#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


# getwd()
library(shiny)
# install.packages('shinydashboard')
library(shinydashboard)
library(waiter)
# install.packages('shinycssloaders')
library(shinycssloaders)
# install.packages('markdown')
# install.packages("DT")
library(DT) # outputting tables
library(Seurat)
library(data.table) # functionality

# Plotting
library(ggrepel)
library(ggpubr)
library(ggExtra)

# Loading data

obj <- readRDS('./data/diet_Seurat_obj_full.Rds')
degs_ <- readRDS('./data/cluster_cond_degs.Rds')
all_markers <- readRDS('./data/cluster_find_all_markers.Rds')
centroids <- readRDS('./data/cluster_centroids.Rds')

marker_heatmap_plots <- readRDS('./data/marker_heatmap_plots.Rds')

lit_markers <- read.csv('./data/MarkerGenes091724.csv')

lit_markers <- lit_markers[lit_markers$Gene %in% rownames(obj),]

lit_markers <- lit_markers[!grepl('/', lit_markers$CellType),]

lit_markers$group <- paste0(lit_markers$Comparison,
                            ':',
                            lit_markers$Reference_Origin)


colors_ <- scales::hue_pal()(length(levels(obj$seurat_clusters)))
names(colors_) <- levels(obj$seurat_clusters)

# Wait screen
waiting_screen <- tagList(
  waiter::spin_loaders(12),
  h4("Loading...")
)

# Define UI for application that draws a histogram

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Project Summary", tabName = 'Home', icon = icon('house', lib='font-awesome')),
    menuItem("DimRed Plots", tabName = 'DimRed',icon = icon('map', lib='font-awesome')),
    menuItem("Quality Control", tabName= 'QC', icon = icon('ranking-star', lib='font-awesome')),
    menuItem("Cluster Abundance", tabName = "Abund", icon = icon('arrow-up-9-1', lib = 'font-awesome')),
    menuItem("Comp. Marker Genes", tabName = "Markers", icon = icon('location-dot', lib = 'font-awesome')),
    menuItem("Gene Expression", tabName = 'GeneExp', icon = icon('glass-water', lib='font-awesome')),
    menuItem("Two Gene Expression", tabName = 'TwoGene', icon = icon('dice-two', lib = 'font-awesome')),
    menuItem("Literature Markers", tabName = "Lit", icon = icon('book', lib = 'font-awesome')),
    menuItem("Condition DEGs", tabName = 'CondDEGs', icon = icon('magnifying-glass', lib='font-awesome'))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = 'Home',
            fluidPage(
              # tags$head(includeHTML(("www/google-analytics.html"))),
              # waiter::useWaiter(),
              # waiter::waiter_preloader(html=waiting_screen),
              includeMarkdown("./scripts/Intro.md"),
              br(),
              br()
              
            )),
    tabItem(tabName = "DimRed",
            h2('Dimensionality Reduction Plots'),
            fluidPage(
              fluidRow(column(6, withSpinner(plotOutput('Plot_PC2'), type = 6), style = 'padding:4px;'),
                       column(6, withSpinner(plotOutput('Plot_UMAP2'), type = 6), style = 'padding:4px;')),
              fluidRow(selectInput("Cluster2", label = "Cluster", choices = NULL)),
              fluidRow(column(6,withSpinner(plotOutput("PC_color_one"), type = 6), style = 'padding:4px;'),
                       column(6,withSpinner(plotOutput("UMAP_color_one"), type = 6), style = 'padding:4px;'))
            )),
    tabItem(tabName = 'QC',
            h2("Quality Control (QC) Figures"),
            fluidPage(fluidRow(selectInput("QC_metric", label = 'QC Metric', choices = NULL)),
                      fluidRow(column(width = 6, plotOutput('QC_featureplot'), style = 'padding:4px;'),
                               column(width = 6, plotOutput('QC_violinplot'), style = 'padding:4px;'),
                               column(width = 6, offset = 3, plotOutput('QC_dimplot'), style = 'padding:4px;'))
            )),
    tabItem(tabName = 'Abund',
            h2("Abundance of Samples"),
            fluidPage(
              fluidRow(column(width = 6, offset = 3,
                              plotOutput('Plot_UMAP_sample'))),
              fluidRow(column(width = 6, plotOutput('Plot_cluster_sample_table'), style = 'padding:4px;'),
                       column(width = 6, plotOutput('Plot_cluster_sample_bar'), style = 'padding:4px;'))
            )),
    tabItem(tabName = "Markers",
            h2("Computational Cluster Marker Genes"),
            fluidPage(
              fluidRow(selectInput("Cluster3", label = "Cluster of Interest", choices = NULL)),
              fluidRow(selectInput("Gene3", label = "Marker of Interest", choices = NULL)),
              fluidRow(column(width = 6,box(width = 12,"Marker Table", DT::dataTableOutput("Marker_table"), style = 'height:800px;overflow-y: scroll;overflow-x: scroll;')),
                       column(width = 6,
                              fluidRow(
                                column(width = 12, 
                                       withSpinner(plotOutput("Marker_feature"), type = 6),
                                       # style = "background-color:green", 
                                       # div(style = "height:200px;padding:4px;"),
                                       sylte = 'padding:4px;')),
                              fluidRow(
                                column(width = 12, 
                                       withSpinner(plotOutput("Marker_violin"),type = 6),
                                       # div(style = "height:200px;padding:4px;"),
                                       sylte = 'padding:4px;'))
                              ))
              )
            ),
    tabItem(tabName = "GeneExp",
            h2("Gene Expression"),
            fluidPage(
              fluidRow(selectInput("Gene1", label = "Gene of Interest", choices = NULL)),
              # fluidRow(column(6, withSpinner(plotOutput('Plot_PC'), type = 6), style = 'padding:4px;'),
              #          column(6, withSpinner(plotOutput('Plot_UMAP'), type = 6), style = 'padding:4px;')),
              fluidRow(column(6, withSpinner(plotOutput('AllFeatPlot.PC'), type=6), style='padding:4px;'),
                       column(6, withSpinner(plotOutput('AllFeatPlot.UMAP'), type=6), style='padding:4px;')),
              fluidRow(column(6, withSpinner(plotOutput('VlnPlot'), type=6), style='padding:4px;'),
                       column(6, withSpinner(plotOutput('VlnPlot.Split'), type=6), style='padding:4px;'))
              
            )),
    tabItem(tabName="TwoGene",
            h2('Co-expression of Genes'),
            fluidPage(
              fluidRow(column(4,selectInput("CoGene1", label = "Gene 1", choices = NULL)),
                       column(4,selectInput("CoGene2", label = "Gene 2", choices = NULL)),
                       column(1,actionButton("LoadCoGene", label = "GO"))),
              fluidRow(column(4, plotOutput("FeatPlot_CoGene1"), style = 'padding:4px;'),
                       column(4, plotOutput("FeatPlot_CoGene2"), style = 'padding:4px;'),
                       column(4, plotOutput("Scatter_CoGene"), style = 'padding:4px;')),
              fluidRow(plotOutput("FeatPlot_CoGene_CT"), style = 'padding:4px;'),
              fluidRow(selectInput("Cluster_CoGene", label = "Cluster", choices = NULL)),
              fluidRow(column(6, plotOutput("Scatter_CoGene_CT"), style = 'padding:4px;'))
            )),
    tabItem(tabName='Lit',
            h2("Collected Literature Marker Lists"),
            fluidPage(
              fluidRow(selectInput("MarkerList", label = "Marker List", choices = NULL)),
              fluidRow(selectInput("CellType", label = "Cell Type", choices = NULL)),
              fluidRow(column(6, plotOutput("Lit_heatmap"), style='padding:4px;'),
                       column(6, plotOutput("Lit_score_violin"), style = 'padding:4px;')),
              fluidRow(column(6,plotOutput("Lit_relabel_dim"), style = 'padding:4px;'),
                       column(6,plotOutput("Lit_dotplot"), style = 'padding:4px;'))
                       
            )),
    tabItem(tabName='CondDEGs',
            h2('DEGs by Condition Across Clusters'),
            fluidPage(
              fluidRow(selectInput("Cluster", label = "Cluster", choices = NULL)),
              fluidRow(column(6,withSpinner(plotOutput("DEG_VolcanoPlot",
                                                       click = "DEGvolcano_click"),
                                            type = 6), style = 'padding:4px;'),
                       column(6,withSpinner(plotOutput("DEG_ViolinPlot"), type = 6), style = 'padding:4px;')
                       )
            ),
            fluidRow(column(width = 12,
                            h4('Points Near Click'),
                            verbatimTextOutput('click_info'))))
  )
)

# ui dashboard
ui <- dashboardPage(
  dashboardHeader(title = 'scRNA Shiny Application Demo'),
  sidebar,
  body
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  
  ##############
  # For Dim Red
  
  ## Top Row
  Plot_PC <- DimPlot(obj, reduction = 'pca',raster = T, label = T, label.size = 5) +ggtitle("PCA")
  Plot_UMAP <- DimPlot(obj, reduction = 'umap',raster = T, label = T,label.size = 5) +ggtitle("UMAP")
  
  output$Plot_PC <- renderPlot({
    print(Plot_PC)
  })
  output$Plot_UMAP <- renderPlot({
    (Plot_UMAP)
  })
  
  ## Row Two
  
  clusters.to.select2 <- levels(obj$seurat_clusters)
  
  updateSelectizeInput(session, 'Cluster2', choices = clusters.to.select2, server = TRUE)
  Cluster2 <- reactive((input$Cluster2))
  
  
  
  output$PC_color_one <- renderPlot({
    colors_ <- scales::hue_pal()(length(levels(obj$seurat_clusters)))
    names(colors_) <- levels(obj$seurat_clusters)
    
    cols_ <- rep('lightgrey', length(levels(obj$seurat_clusters)))
    names(cols_) <- levels(obj$seurat_clusters)
    
    cols_[Cluster2()] <- colors_[Cluster2()]
    (DimPlot(obj, reduction = 'pca', cols = cols_, raster = T) + ggtitle("PCA"))
  })
  
  output$UMAP_color_one <- renderPlot({
    colors_ <- scales::hue_pal()(length(levels(obj$seurat_clusters)))
    names(colors_) <- levels(obj$seurat_clusters)
    
    cols_ <- rep('lightgrey', length(levels(obj$seurat_clusters)))
    names(cols_) <- levels(obj$seurat_clusters)
    
    cols_[Cluster2()] <- colors_[Cluster2()]
    (DimPlot(obj, reduction = 'umap', cols = cols_, raster = T) + ggtitle("UMAP"))
  })
  
  ##############
  # Quality Control
  
  colnames(obj@meta.data)
  
  qc_metrics <- c('log10_nCount_RNA', 'nFeature_RNA', 'percent.mt')
  updateSelectizeInput(session,'QC_metric', choices = qc_metrics,server = T)
  qc_metric <- reactive(input$QC_metric)
  
  output$QC_featureplot <- renderPlot({
    req(qc_metric() %in% colnames(obj@meta.data))
    FeaturePlot(obj, feature = qc_metric(), cols = c('white','blue'),raster = T)
  })
  
  output$QC_violinplot <- renderPlot({
    req(qc_metric() %in% colnames(obj@meta.data))
    VlnPlot(obj, feature = qc_metric(),pt.size = 0)
  })
  
  output$QC_dimplot <- renderPlot({
    DimPlot(obj, label = T, raster = T, label.size = 5)
  })
  
  
  ##############
  # Abundance
  
  output$Plot_UMAP_sample <- renderPlot({
    DimPlot(obj, group.by = 'orig.ident', raster = T)
  })
  
  
  cnt_table <- as.data.frame(table(obj$seurat_clusters))
  
  cnt_table$Var2 <- 'Total'
  
  cnt_table <- rbind(cnt_table[,c(1,3,2)],as.data.frame(table(obj$seurat_clusters, obj$orig.ident)))
  
  cnt_table$Var2 <- factor(cnt_table$Var2, levels = c('WT','MUT','Total'))
  
  cnt_table$Var1_total <- NA
  
  for(var_ in unique(cnt_table$Var1)){
    cnt_table[cnt_table$Var1 == var_,]$Var1_total <- cnt_table[cnt_table$Var2 == 'Total' & cnt_table$Var1 ==var_,]$Freq
  }
  
  cnt_table$fill <- ifelse(cnt_table$Var2 == 'Total', NA,
                           cnt_table$Freq / cnt_table$Var1_total)
  
  cnt_table$Var2_total <- NA
  
  for(var_ in unique(cnt_table$Var2)){
    cnt_table[cnt_table$Var2 == var_,]$Var2_total <- sum(cnt_table[cnt_table$Var2 == var_,]$Freq)
  }
  
  cnt_table$fill2 <- cnt_table$Freq / cnt_table$Var2_total
  output$Plot_cluster_sample_table <- renderPlot({
    ggplot(cnt_table, aes(x = Var2, y = Var1, label = Freq, fill = fill)) +
      geom_tile(color = 'black',na.rm = F) +
      labs(fill = 'Perc.\nCluster') +
      scale_fill_gradient(low = 'white', high = 'red') +
      geom_text() + 
      theme_classic() +
      xlab(NULL) +
      ylab('Cluster')
  })
  
  output$Plot_cluster_sample_bar <- renderPlot({
    ggplot(cnt_table[cnt_table$Var2 != 'Total',],
           aes(x = Var2, y = fill2, fill = Var1)) +
      geom_bar(stat = 'identity', colour = 'black',
               position = position_fill(reverse = T)) +
      theme_classic() +
      xlab(NULL) +
      ylab('Perc. Sample') +
      labs(fill = 'Clusters') +
      scale_color_manual(values = scales::hue_pal()(length(unique(cnt_table$Var1))))
  })
  
  ##############
  # Marker Genes
  
  clusters.to.select3 <- levels(obj$seurat_clusters)
  updateSelectizeInput(session, 'Cluster3', choices = clusters.to.select3,selected = "0", server = TRUE)
  Cluster3 <- reactive((input$Cluster3))
  
  observeEvent(input$Cluster3, {
    marker_table <- (all_markers[all_markers$cluster == Cluster3(),])
    
    genes.to.select3 <- marker_table$gene[1:100]
    updateSelectizeInput(session, 'Gene3', choices = genes.to.select3, server = TRUE)
    
  })
  
  Gene3 <- reactive((input$Gene3))
  
  output$Marker_table <- renderDataTable({
    req(Cluster3() %in% all_markers$cluster)
    marker_table <- all_markers[all_markers$cluster == Cluster3(),]
    rownames(marker_table) <- marker_table$gene
    datatable(marker_table[1:100,2:5], options = list(paging = F))
  })
  
  output$Marker_feature <- renderPlot({
    req(Gene3() %in% rownames(obj))
    
    FeaturePlot(obj, features = Gene3(),label = T, repel = T,
                cols = c('lightgrey','blue'),order = T)
  })
  
  output$Marker_violin <- renderPlot({
    req(Cluster3() %in% obj$seurat_clusters)
    req(Gene3() %in% rownames(obj))
    
    cluster_ <- Cluster3()
    p <- VlnPlot(obj, features = Gene3(), pt.size = 0)
    
    p$data$grp <- ifelse(p$data$ident == cluster_,'++', 'Other Clusters')
    
    p + facet_grid(cols = vars(grp),scales = 'free',space = 'free') +
      xlab(NULL)
  }) 
  
  ##############
  # For GeneExp
  
  output$Plot_PC2 <- renderPlot({
    print(Plot_PC)
  })
  output$Plot_UMAP2 <- renderPlot({
    print(Plot_UMAP)
  })
  
  genes.to.select <- rownames(obj)
  genes.to.select <- genes.to.select[order(genes.to.select)]
  updateSelectizeInput(session, 'Gene1', choices = genes.to.select, selected='Arg1', server = TRUE)
  Gene1 <- reactive((input$Gene1))
  
  output$AllFeatPlot.PC <- renderPlot({
    req(Gene1() %in% rownames(obj))
    print(FeaturePlot(obj, features = Gene1(), order = T, min.cutoff = 1,
                      label = T, repel = T,
                      reduction = 'pca'))
  })
  
  output$AllFeatPlot.UMAP <- renderPlot({
    req(Gene1() %in% rownames(obj))
    print(FeaturePlot(obj, features = Gene1(), order = T, min.cutoff = 1,
                      label = T, repel = T,
                      reduction = 'umap'))
  })
  
  output$VlnPlot <- renderPlot({
    req(Gene1() %in% rownames(obj))
    print(VlnPlot(obj, features = Gene1(), slot = 'data', pt.size = 0))
  })
  
  output$VlnPlot.Split <- renderPlot({
    req(Gene1() %in% rownames(obj))
    print(VlnPlot(obj, features = Gene1(), slot = 'data', pt.size = 0,split.by = 'orig.ident'))
  })
  
  ##############
  # Two Genes Co-expression
  
  updateSelectizeInput(session, 'CoGene1', choices = genes.to.select,selected = 'Chd3',server = TRUE)
  CoGene1 <- reactive((input$CoGene1))
  
  updateSelectizeInput(session, 'CoGene2', choices = genes.to.select,selected = 'Rps3', server = TRUE)
  CoGene2 <- reactive((input$CoGene2))
  
  clusters.to.select3 <- unique(Idents(obj))
  clusters.to.select3 <- clusters.to.select3[order(as.numeric(as.character(clusters.to.select3)))]
  updateSelectizeInput(session, 'Cluster_CoGene', choices = clusters.to.select3, selected=NULL, server = TRUE)
  Cluster_coGene <- reactive((input$Cluster_CoGene))
    

  output$FeatPlot_CoGene1 <- renderPlot({
    req(CoGene1() %in% rownames(obj))
    req(input$LoadCoGene)
    
    print(FeaturePlot(obj, features = CoGene1(), order = T, min.cutoff = 1,
                      label = T, repel = T,
                      reduction = 'umap'))
  })
  
  output$FeatPlot_CoGene2 <- renderPlot({
    req(CoGene2() %in% rownames(obj))
    req(input$LoadCoGene)
    print(FeaturePlot(obj, features = CoGene2(), order = T, min.cutoff = 1,
                      label = T, repel = T,
                      reduction = 'umap'))
  })
  
  md_for_coexpression <- eventReactive(inptu$LoadCoGene,{
    md <- obj@meta.data
    md$gene1 <- obj@assays$RNA$data[CoGene1(),]
    md$gene2 <- obj@assays$RNA$data[CoGene2(),]
  })
  
  output$Scatter_CoGene <- renderPlot({
    req(CoGene1() %in% rownames(obj))
    req(CoGene2() %in% rownames(obj))
    req(input$LoadCoGene)
    
    scatter_plot <- ggplot(md, aes(x = gene1, y = gene2, color = seurat_clusters)) +
      geom_point() +theme_bw() +
      theme(legend.position = 'bottom') +
      xlab(CoGene1()) +ylab(CoGene2())
    
    ggMarginal(scatter_plot,
               type = 'violin',
               margins = 'both', 
               groupColour = F,
               groupFill = F)
  })
  
  output$FeatPlot_CoGene_CT <- renderPlot({
    req(CoGene1() %in% rownames(obj))
    req(CoGene2() %in% rownames(obj))
    req(input$LoadCoGene)
    
    FeaturePlot(obj, features = c(CoGene1(), CoGene2()), blend = T, pt.size = 0.1)
  })
  
  output$Scatter_CoGene_CT <- renderPlot({
    req(CoGene1() %in% rownames(obj))
    req(CoGene2() %in% rownames(obj))
    req(input$LoadCoGene)
    req(Cluster_coGene() %in% md$seurat_clusters)
    
    scatter_plot <- ggplot(md[md$seurat_clusters == Cluster_coGene(),],
                           aes(x = gene1, y = gene2,
                               color = seurat_clusters)) + 
      geom_point(color = colors_[as.character(Cluster_coGene())]) +
      # geom_point() +
      theme_bw() +
      theme(legend.position = 'bottom') +
      xlab(CoGene1()) +ylab(CoGene2())
    
    # scatter_plot
    ggMarginal(scatter_plot,
               type = 'violin',
               margins = 'both')
  })
  
  ##############
  # Lit Markers
  
  marker_lists <- unique(lit_markers$group)
  updateSelectizeInput(session,"MarkerList", choices = marker_lists,server = T)
  MarkerList <- reactive(input$MarkerList)
  
  observeEvent(input$MarkerList, {
    temp <- lit_markers[lit_markers$group == MarkerList(),]
    
    cell_types <- unique(temp$CellType)
    cell_types <- cell_types[order(cell_types)]
    
    updateSelectizeInput(session, 'CellType', choices = cell_types, server = T)
  })
  
  CellType <- reactive(input$CellType)
  
  output$Lit_heatmap <- renderPlot({
    req(MarkerList() %in% names(marker_heatmap_plots))
    marker_heatmap_plots[[MarkerList()]]
  })
  
  output$Lit_score_violin <- renderPlot({
    req(MarkerList() %in% names(marker_heatmap_plots))
    req(CellType() %in% lit_markers$CellType)
    
    feat_ <- colnames(obj@meta.data)[grepl(MarkerList(), colnames(obj@meta.data)) &
                                             grepl(paste0(':',
                                                          CellType(),
                                                          '$'),
                                                   colnames(obj@meta.data))]
    
    obj$temp <- obj@meta.data[,feat_]
    
    VlnPlot(obj, features = 'temp', pt.size = 0) +
      ggtitle(paste0(MarkerList(),':', CellType()))
    
  })
  
  output$Lit_relabel_dim <- renderPlot({
    # req(MarkerList() %in% names(marker_heatmap_plots))
    temp <- colnames(obj@meta.data)[grepl(MarkerList(), colnames(obj@meta.data)) &
                                          grepl('LabeledClusters$',
                                                colnames(obj@meta.data))]
    
    obj$temp <- obj@meta.data[,temp]
    DimPlot(obj, group.by = 'temp', label = T) +
      ggtitle(paste0('Clusters relabeled by top match:', MarkerList()))
  })
  
  output$Lit_dotplot <- renderPlot({
    req(MarkerList() %in% names(marker_heatmap_plots))
    req(CellType() %in% lit_markers$CellType)
    
    temp_markers <- lit_markers[lit_markers$group == MarkerList() &
                                  lit_markers$CellType == CellType(),]$Gene
    
    DotPlot(obj, features = temp_markers, cols = c('white','blue')) +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste0(MarkerList(),':', CellType()))
  })
  
  
  ##############
  # For CondDEGs

  clusters.to.select <- names(degs_)
  
  updateSelectizeInput(session, 'Cluster', choices = clusters.to.select, selected='Cluster0', server = TRUE)
  Cluster <- reactive((input$Cluster))
  
  colors_volcano <- c('lightgrey','red','blue')
  names(colors_volcano) <- c('N.S.','Up\nin MUT', 'Up\nin WT')
  
  
  output$DEG_VolcanoPlot <- renderPlot({
    req(Cluster() %in% names(degs_))
    ggplot(degs_[[Cluster()]],aes(x = avg_log2FC, y = `-log10(p_val_adj)`,
                                  label = label, color = color)) +
      geom_point() + ggtitle(Cluster()) +
      scale_color_manual(values = colors_volcano[levels(as.factor(degs_[[Cluster()]]$color))]) + 
      geom_hline(yintercept = -log10(0.05), colour = 'black', linetype = 2) +
      ggrepel::geom_label_repel() +
      theme_bw()
      
  })
  
  selected_point <- reactive(nearPoints(degs_[[Cluster()]],input$DEGvolcano_click))
  
  output$click_info <- renderPrint({
    selected_point()
  })
  
  output$DEG_ViolinPlot <- renderPlot({
    req(selected_point() %in% rownames(obj))
    req(Cluster() %in% names(degs_))
    cluster_ <- data.table::tstrsplit(Cluster(),'r',keep = 2)[[1]]
    # VlnPlot(obj, features = selected_point()$gene, pt.size = 0,
    #         split.by = 'orig.ident') +
    #   scale_x_discrete(labels = function(x){
    #     ifelse(x == cluster_,paste0("**",x,"**"),x)
    #   })
    
    p <- VlnPlot(obj, features = selected_point()$gene, pt.size = 0,
                         split.by = 'orig.ident', cols = c('red','blue'))
    
    p$data$grp <- ifelse(p$data$ident == cluster_,'++', 'Other Clusters')
    
    p + facet_grid(cols = vars(grp),scales = 'free',space = 'free') +
      xlab(NULL)

  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
