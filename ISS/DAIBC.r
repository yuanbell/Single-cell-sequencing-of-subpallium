## DESCIRPTION ##
#   This Web App is used to appear the result of ISS Base Calling
#################

## AUTHOR ##
#   Yang Zhou (zhouyang@genomics.cn)
#   Hao Yu (yuhao@genomics.cn)
############

## VERSION ##
#   2019-04-04  v1.0    Initial release
#   2019-04-16  v1.1    To modify the large file (>5M) input
#   2019-12-09  v1.2    Add support to pre-defined color and shape
#############


library(shiny)
library(shinyWidgets)
library(shinythemes)
library(DT)
library(hash) 
library(tiff)
library(grid)
library(ggplot2)
library(Cairo)
library(rsconnect)


# Design the UI
ui <- fixedPage(
  theme=shinytheme('lumen'),
  
  titlePanel('DAIBC', 'Data Analysis after ISS Base Calling'),
  
  h5('Data Analysis after ISS Base Calling'),
  
  sidebarPanel(
    fixedRow(
      column(4,
             fileInput('bfile',
                       'Import Barcode Info',
                       accept=c('text/plain', 'text/csv'),
                       placeholder='Barcode Info File'
             ),
             
             fileInput('dfile',
                       'Import Base Calling Data',
                       accept=c('text/plain'),
                       placeholder='Base Calling Data File'
             ),
             
             fileInput('ifile',
                       'Import Background Image',
                       accept=c('image/tiff'),
                       placeholder='Background Image File'
             ),
             
             actionButton('parse', 'Parse Imported Data', width='100%', icon=icon('microchip')),
             conditionalPanel(condition='input.parse > 0', hr(), downloadButton('export', 'Export Result'))
      ),
      
      column(8,
             conditionalPanel(condition='input.parse > 0', 
                              
                              actionButton('all', 'Select ALL', width='30%', icon=icon('check-square')),
                              actionButton('non', 'Select Non', width='30%', icon=icon('square')),
                              downloadButton('downloadData', 'Export Data'),
                              br(),
                              br(),
                              
                              DT::dataTableOutput('binfo', height=400)
                              # verbatimTextOutput('debug')  # DEBUG
             )
      )
    ),
    
    width=14
  ),
  
  conditionalPanel(condition='input.parse > 0',
                   mainPanel(plotOutput('mplot', height=600), width=9),
                   
                   sidebarPanel(sliderInput('qua_thre',
                                            'Quality Threshold',
                                            min=34, max=74, step=2, value=34
                   ),
                   
                   hr(),
                   
                   sliderInput('sct_size',
                               'Scatter Size',
                               min=1, max=10, value=2
                   ),
                   
                   sliderInput('sct_lumi',
                               'Scatter Luminance',
                               min=0, max=1.0, step=0.1, value=0.4
                   ),
                   
                   hr(),
                   
                   materialSwitch('scope', 'Exp. Scope Option', width='100%', status='success'),
                   
                   conditionalPanel(condition='input.scope == true',
                                    sliderInput('win_length',
                                                'Window Length',
                                                min=20, max=100, step=2, value=60
                                    ),
                                    
                                    sliderInput('min_num',
                                                'Minimum Scatter Number',
                                                min=1, max=30, step=1, value=1
                                    )
                   ),
                   
                   width=3
                   )
  )
)




# Preprocess data
options(shiny.maxRequestSize=(128 * 1024 ^ 4))

server <- function(input, output, session) {
    observeEvent(input$parse,
                 {
                    barcodeInfoFile <- input$bfile$datapath
                    baseCallingDataFile <- input$dfile$datapath
                    bgImgFile <- input$ifile$datapath
        
                    if (is.null(barcodeInfoFile)) {barcodeInfoFile <- 'example/barcode_info.txt'}
                        
                    barcodeInfo <- read.table(barcodeInfoFile, sep='\t', comment.char='', quote='')
                    barcodeInfo <- barcodeInfo[order(barcodeInfo[, 1]),]
                    col_number = ncol(barcodeInfo)
                    if(col_number == 2){
                        barcodeInfo <- cbind(barcodeInfo, '')
                        barcodeInfo <- cbind(barcodeInfo, 17)
                    }
                    barcodeInfo <- cbind(barcodeInfo, 0)
                    rownames(barcodeInfo) <- barcodeInfo[, 1]
                    colnames(barcodeInfo) = c('Barcode', 'Gene Name', 'Color', 'Shape', 'Called Number')
                    
                    barcodeHash <- hash(keys=barcodeInfo$Barcode, values=barcodeInfo$Shape)
                    
                    if (is.null(baseCallingDataFile)) {baseCallingDataFile <- 'example/basecalling_data.txt'}
                    
                    baseCallingData <- read.table(baseCallingDataFile, sep='\t', comment.char='', quote='')
                    baseCallingData <- baseCallingData[grep('N', baseCallingData[, 2], invert=TRUE),]
                    baseCallingData <- baseCallingData[grep('!', baseCallingData[, 3], invert=TRUE),]
                    baseCallingData <- baseCallingData[order(baseCallingData[, 2]),]
                    colnames(baseCallingData) = c('ID', 'Barcode', 'Quality', 'Row', 'Column')
                    
                    qList <- c()
                    
                    for (q in baseCallingData$Quality) {qList <- c(qList, round(mean(strtoi(charToRaw(q), 16L))))}
                    
                    baseCallingData$Quality <- qList
                    
                    if (is.null(bgImgFile)) {bgImgFile <- 'example/bg_img.tif'}
                    
                    bgImg <- readTIFF(bgImgFile, as.is=TRUE)
                    max_row <- dim(bgImg)[1]
                    max_col <- dim(bgImg)[2]
                    bgImg <- bgImg / 255
                    bgImg <- rasterGrob(bgImg)

                    baseCallingData$Row <- max_row - baseCallingData$Row
                    
                    shapeVec <- rep(17, nrow(baseCallingData))
                    if(col_number == 4){
                        shapeVec = apply(baseCallingData, 1, function(x){
                            if(sum(barcodeInfo$Barcode == x[2]) == 0){
                                NA
                            }
                            else{
                                barcodeInfo$Shape[which(barcodeInfo$Barcode == x[2])]
                            }
                        })
                    }
                    shapeVec = as.vector(shapeVec)
                    baseCallingData = cbind(baseCallingData, shapeVec)
                    colnames(baseCallingData)[6] = 'Shape'
                    
                    combination = paste(baseCallingData$Barcode, baseCallingData$Shape, sep = '-')
                    geneName = as.character(apply(baseCallingData, 1, function(x){
                        if(sum(barcodeInfo$Barcode == x[2]) == 0){
                            'NA'
                        }
                        else{
                            barcodeInfo$`Gene Name`[which(barcodeInfo$Barcode == x[2])]
                        }
                    }))
                    baseCallingData = cbind(baseCallingData, combination, geneName)
                        
                    img <- ggplot() +
                        annotation_custom(bgImg, xmin=0, xmax=max_col, ymin=0, ymax=max_row) + 
                        scale_x_continuous(limits=c(0, max_col)) + scale_y_continuous(limits=c(0, max_row)) + 
                        theme(panel.background=element_blank(), 
                              panel.border=element_blank(),
                              panel.grid=element_blank(),
                              axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank(),
                              legend.position='right'
                        )
                    
                    barcodeNum <- hash(keys=baseCallingData[!duplicated(baseCallingData[, 2]), 2], 
                                       values=table(unlist(matrix(baseCallingData[, 2])))
                    )
                    
                    for (k in barcodeInfo[, 1]) {if(has.key(k, barcodeNum)) {barcodeInfo[k, 5] <- barcodeNum[[k]]}}
                    
                    if(col_number == 2){
                        color_list <- as.data.frame(col2rgb(colors()[seq(657, 1, -3)]))
                        color_list <- rgb(color_list[1,] / 255, color_list[2,] / 255, color_list[3,] / 255)
                        shape_list <- rep(17, nrow(barcodeInfo))
                    }
                    else{
                        color_list <- as.data.frame(col2rgb(barcodeInfo$Color))
                        color_list <- rgb(color_list[1,] / 255, color_list[2,] / 255, color_list[3,] / 255)
                        shape_list <- barcodeInfo$Shape
                    }
                    output$binfo <- renderDataTable(
                        {
                            dt <- DT::datatable(barcodeInfo,
                                                options=list(dom='ftipr', 
                                                             order=list(2, 'desc'),
                                                             columnDefs=list(list(className='dt-left', targets=2),
                                                                             list(width='10%', targets=3)
                                                             )
                                                ),
                                                
                                                rownames=FALSE,
                                                fillContainer=TRUE,
                                                style='bootstrap'
                            )
                            
                            formatStyle(dt,
                                        columns='Color',
                                        valueColumns='Barcode',
                                        backgroundColor=styleEqual(barcodeInfo[, 1], c(color_list[1:nrow(barcodeInfo)]))
                            )
                            
                        }
                    )
                    
                    
                    
                    proxy = DT::dataTableProxy('binfo')
            
                    observeEvent(input$all, {DT::selectRows(proxy, input$binfo_rows_all)})
                    observeEvent(input$non, {DT::selectRows(proxy, NULL)})
                    
                    output$downloadData <- downloadHandler(
                      filename = 'table.csv',
                      content = function(file) {
                        write.csv(barcodeInfo, file)
                      }
                    )
                    
                    
                    
                    output$mplot <- renderPlot(
                        {
                            if (! is.null(input$binfo_rows_selected)) {
                                used_data <- baseCallingData[baseCallingData[, 2] %in% barcodeInfo[sort(input$binfo_rows_selected), 1] & 
                                                             baseCallingData[, 3] >= input$qua_thre,]

                                r = sort(input$binfo_rows_selected)
                                r = r[which(barcodeInfo[r, 5] > 0)]
                                
                                f_barcode = factor(used_data$Barcode, levels = barcodeInfo[sort(r), 2])
                                used_data$Barcode = factor(used_data$Barcode, levels = barcodeInfo[sort(r), 2])
                                used_data$geneName = factor(used_data$geneName, levels = unique(used_data$geneName[order(f_barcode)]))
                                
                                img <- img + 
                                    geom_point(aes(x=Column, y=Row, color=geneName, shape=geneName),
                                               data=used_data, 
                                               size=input$sct_size, alpha=input$sct_lumi
                                    ) +
                                    scale_color_manual(values=color_list[sort(r)], labels = barcodeInfo[sort(r), 2]) + 
                                    scale_shape_manual(values=shape_list[sort(r)], labels = barcodeInfo[sort(r), 2])

                                if (input$scope == TRUE) {
                                    win_len <- input$win_length
                                    
                                    win_c_x <- seq(win_len / 2, max_col, win_len)
                                    win_c_y <- seq(win_len / 2, max_row, win_len)
                                    
                                    win_x_num <- length(win_c_x)
                                    win_y_num <- length(win_c_y)
                                    
                                    permut <- expand.grid(win_c_x, win_c_y)
                                    
                                    colnames(permut) <- c('x', 'y')
                                    
                                    tmp <- cbind(data.frame(num=rep(NA, win_x_num * win_y_num)), permut)
                                    idx <- paste(tmp$x, tmp$y, sep='-')
                                    tmp <- cbind(idx, tmp)
                                    
                                    xloc <- floor(used_data$Column / win_len) * win_len + win_len / 2
                                    yloc <- floor(used_data$Row / win_len) * win_len + win_len / 2
                                    loc <- paste(xloc, yloc, sep = '-')
                                    
                                    tab <- data.frame(table(loc))
                                    colnames(tab) <- c('idx', 'num')
                                    tab$num[tab$num < input$min_num] = NA
                                    tile <- rbind(tab, tmp[which(!(tmp$idx %in% tab$idx)), c(1, 2)])
                                    spt <- unlist(strsplit(as.character(tile$idx), '-'))
                                    
                                    xloc <- as.numeric(spt[seq(1, length(spt), by=2)])
                                    yloc <- as.numeric(spt[seq(2, length(spt), by=2)])
                                    
                                    tile_w <- rep(win_len, length(tile$idx))
                                    tile_h <- rep(win_len, length(tile$idx))
                                    tile <- cbind(tile, xloc, yloc, tile_w, tile_h)
                                    
                                    img <- img + 
                                        geom_raster(aes(x=xloc, y=yloc, width=tile_w, height=tile_h, fill=num), 
                                                  data=tile,
                                                  na.rm=TRUE,
                                                  alpha=0.3,
                                                  interpolate=TRUE
                                        ) +
                                        scale_fill_gradient(na.value=rgb(0, 0, 0, alpha=0))

                                }
                            }

                            img <- img + coord_equal()
                            
                            plot(img)

                            output$export <- downloadHandler(
                                filename=function() {return(paste('DAIBC.', format(Sys.time(), format='%Y%m%d%H%M%S'), '.pdf', sep=''))},
                                
                                content=function(file) {
                                    CairoPDF(file, width=15, height=15)
                                    plot(img)
                                    dev.off()
                                },
                                
                                contentType='application/pdf'
                            )
                        }
                    )
            
                    ## DEBUG ##
                    # output$debug <- renderPrint({
                    #     print(db)
                    # })
                    ###########
                 }
    )
}



# Run the application 
shinyApp(ui = ui, server = server)
