library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.refGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(shiny)
library(formattable)
setwd("E:/Users/cwang/DMGE/")
# setwd("/home/cwang/project/dmge")
ref <- readRDS("data/ref.rds")
tkRefseq <- readRDS("data/tkRefseqHg19.rds")
listTkIdeo <- readRDS("data/listTkIdeoHg19.rds")
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.refGene)
gr <- readRDS("result/roadmap/gr/gr_0.1_60.rds")
seqInfo <- seqinfo(gr)
tkGa <- GenomeAxisTrack(cex = 1.2)
displayPars(tkRefseq) <- list(cex.group = 1,
                              size = 0.3,
                              just.group = "left")
tkSeq <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg19)

eidColors <- setNames(ref$COLOR, ref$EID)
shortNameColors <- setNames(ref$COLOR, ref$SHORT_NAME)
markColors <- c("#F2EB17", "#EE1F23", "#0D803F", "#B8B9BC", "#8FCDDB")
mi <- -0.005748266
ma <- 1.031437
chromPx <- 150
refseqPx <- 43 * 0.3
csrePx <- 43 * 0.8

color_bar_my <- function(color = "lightgray",
                         fun = formattable::normalize, ...) {
  fun <- match.fun(fun)
  formatter("span",
            style = function(x)
              style(display = "inline-block",
                    direction = "ltr",
                    `border-radius` = "4px",
                    `padding-right` = "2px",
                    `background-color` = csscolor(color),
                    width = percent(fun(as.numeric(x), ...))),
            x ~ sprintf("%.2f", x))
}

server <- shinyServer(function(input, output, session) {

  rv <- reactiveValues(region = GRanges("chr22:22888118-22902174",
                                        seqinfo = seqInfo),
                       subRegion = GRanges("chr22:22888118-22902174",
                                           seqinfo = seqInfo))
  observeEvent(input$go, {
    gor <- input$Gene_OR_Region
    if (any(grepl("^chr", gor))) {
      region <- tryCatch(unstrand(trim(GRanges(gor, seqinfo = seqInfo))),
                         error = function(e) e)
      if (!inherits(region, "error")) {
        rv$region <- region
        rv$subRegion <- rv$region
      }
    } else {
      nameGene <- gor
      if (nameGene %in% keys(org.Hs.eg.db, "SYMBOL")) {
        nameRefseq <- mapIds(org.Hs.eg.db, nameGene, "REFSEQ", "SYMBOL")
        region <- unstrand(trim(tx[tx$tx_name == nameRefseq] + 2000))
        rv$region <- region
        rv$subRegion <- rv$region
      }
    }

  })

  observeEvent(input$left_95, {
    region <- rv$region
    wid <- width(region) * 0.95
    left <- start(region) - 1
    wid <- if (left >= wid) {
      wid
    } else {
      left
    }
    region <- shift(region, -wid)
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$left_47.5, {
    region <- rv$region
    wid <- width(region) * 0.475
    left <- start(region) - 1
    wid <- if (left >= wid) {
      wid
    } else {
      left
    }
    region <- shift(region, -wid)
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$left_10, {
    region <- rv$region
    wid <- width(region) * 0.1
    left <- start(region) - 1
    wid <- if (left >= wid) {
      wid
    } else {
      left
    }
    region <- shift(region, -wid)
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$right_10, {
    region <- rv$region
    wid <- width(region) * 0.1
    right <- seqlengths(region)[as.character(seqnames(region))] - end(region)
    wid <- if (right >= wid) {
      wid
    } else {
      right
    }
    region <- shift(region, wid)
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$right_47.5, {
    region <- rv$region
    wid <- width(region) * 0.475
    right <- seqlengths(region)[as.character(seqnames(region))] - end(region)
    wid <- if (right >= wid) {
      wid
    } else {
      right
    }
    region <- shift(region, wid)
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$right_95, {
    region <- rv$region
    wid <- width(region) * 0.95
    right <- seqlengths(region)[as.character(seqnames(region))] - end(region)
    wid <- if (right >= wid) {
      wid
    } else {
      right
    }
    region <- shift(region, wid)
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$in_1.5, {
    region <- rv$region
    regionNew <- region * 1.5
    if (width(regionNew) >= 50) {
      region <- regionNew
    }
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$in_3, {
    region <- rv$region
    regionNew <- region * 3
    if (width(regionNew) >= 50) {
      region <- regionNew
    }
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$in_10, {
    region <- rv$region
    regionNew <- region * 10
    if (width(regionNew) >= 50) {
      region <- regionNew
    }
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })

  observeEvent(input$out_1.5, {
    region <- rv$region
    region <- region * (1 / 1.5)
    lenChr <- seqlengths(region)[as.character(seqnames(region))]
    s <- start(region)
    e <- end(region)
    if (s < 1)
      region <- trim(shift(region, -s + 1))
    if (e > lenChr)
      region <- trim(shift(region, lenChr - e))
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$out_3, {
    region <- rv$region
    region <- region * (1 / 3)
    lenChr <- seqlengths(region)[as.character(seqnames(region))]
    s <- start(region)
    e <- end(region)
    if (s < 1)
      region <- trim(shift(region, -s + 1))
    if (e > lenChr)
      region <- trim(shift(region, lenChr - e))
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })
  observeEvent(input$out_10, {
    region <- rv$region
    region <- region * (1 / 10)
    lenChr <- seqlengths(region)[as.character(seqnames(region))]
    s <- start(region)
    e <- end(region)
    if (s < 1)
      region <- trim(shift(region, -s + 1))
    if (e > lenChr)
      region <- trim(shift(region, lenChr - e))
    updateTextInput(session, "Gene_OR_Region",
                    value = as.character(region))
    rv$region <- region
    rv$subRegion <- rv$region
  })

  hei <- eventReactive(rv$region, {
    region <- rv$region
    txInRegion <- subsetByOverlaps(tx, region)
    start(txInRegion) <- pmax(start(region), start(txInRegion))
    end(txInRegion) <- pmin(end(region), end(txInRegion))
    numRefseq <- if (length(txInRegion) == 0) {
      0
    } else {
      max(disjointBins(txInRegion))
    }
    grUsed <- subsetByOverlaps(gr, region)
    numCt <- length(unique(grUsed$cellType))
    height <- chromPx + refseqPx * numRefseq + csrePx * numCt
    list(height = height, numRefseq = numRefseq, numCt = numCt)
  })
  output$plotUI <- renderUI({
    plotOutput("plotCSRE", height = hei()$height,
               hover = hoverOpts("plot_hover", 50, "debounce"),
               dblclick = "plot_dblclick",
               brush = brushOpts(id = "plot_brush",
                                 delay = 100,
                                 delayType = "debounce",
                                 direction = "x",
                                 resetOnNew = TRUE))
  })

  observeEvent(rv$region, {
    output$plotCSRE <- renderPlot({
      region <- rv$region
      grUsed <- subsetByOverlaps(gr, region)
      chr <- as.character(seqnames(region))
      from <- start(region)
      to <- end(region)

      tkIdeo <- listTkIdeo[[chr]]
      displayPars(tkIdeo) <- list(cex = 1.2)

      #### stacked
      # tkCsre <- AnnotationTrack(grUsed,
      #                           feature = grUsed$cellType,
      #                           id = grUsed$cellType,
      #                           showFeatureId = TRUE,
      #                           col = NULL,
      #                           cex = 0.9,
      #                           fontcolor.item = "white",
      #                           group = grUsed$cellType,
      #                           # showId = TRUE,
      #                           # just.group = "below",
      #                           col.line = "white",
      #                           name = "CSRE")
      # displayPars(tkCsre) <- as.list(eidColors)

      #### each line for each cell type
      # if (length(grUsed) > 0) {
      #   ct <- unique(grUsed$cellType)
      #   w <- width(region)
      #   l <- min(start(grUsed), from)
      #   r <- max(start(grUsed), to)
      #   left <- GRanges(chr, rep(IRanges(l,
      #                                    l - 1),
      #                            length(ct)),
      #                   seqinfo = seqInfo,
      #                   cellType = Rle(ct))
      #   right <- GRanges(chr, rep(IRanges(r + 1,
      #                                     r),
      #                             length(ct)),
      #                    seqinfo = seqInfo,
      #                    cellType = Rle(ct))
      #   mcols(grUsed) <- mcols(grUsed)["cellType"]
      #   grPlot <- c(left, grUsed, right)
      #   tkCsre <- AnnotationTrack(grPlot,
      #                             feature = grPlot$cellType,
      #                             id = grPlot$cellType,
      #                             # showFeatureId = TRUE,
      #                             # cex = 0.9,
      #                             fontcolor.item = "white",
      #                             group = grPlot$cellType,
      #                             showId = TRUE,
      #                             cex.group = 1,
      #                             just.group = "above",
      #                             col.line = "white",
      #                             name = "CSRE",
      #                             cex.group = 1,
      #                             size = 0.8,
      #                             min.width = 0, ## do not show width 0 ranges
      #                             col = NULL, ## do not show border and width 0 ranges
      #                             stackHeight = 0.75,
      #                             rotation.group = 0)
      #   displayPars(tkCsre) <- as.list(eidColors)
      # } else {
      #   tkCsre <- NULL
      # }

      #### each line for each cell type
      if (length(grUsed) > 0) {
        ct <- unique(grUsed$cellType)
        w <- width(region)
        l <- min(start(grUsed), from)
        r <- max(start(grUsed), to)
        left <- GRanges(chr, rep(IRanges(l,
                                         l - 1),
                                 length(ct)),
                        seqinfo = seqInfo,
                        cellType = Rle(ct))
        right <- GRanges(chr, rep(IRanges(r + 1,
                                          r),
                                  length(ct)),
                         seqinfo = seqInfo,
                         cellType = Rle(ct))
        mcols(grUsed) <- mcols(grUsed)["cellType"]
        grPlot <- c(left, grUsed, right)
        grPlot$name <-
          Rle(as.character(ref$SHORT_NAME)[match(as.character(grPlot$cellType),
                                             ref$EID)])
        tkCsre <- AnnotationTrack(grPlot,
                                  feature = grPlot$name,
                                  # id = grPlot$name,
                                  # showFeatureId = TRUE,
                                  # cex = 0.9,
                                  fontcolor.item = "white",
                                  group = grPlot$cellType,
                                  groupAnnotation = "feature",
                                  showId = TRUE,
                                  cex.group = 1,
                                  just.group = "above",
                                  col.line = "white",
                                  name = "CSRE",
                                  cex.group = 1,
                                  size = 0.8,
                                  min.width = 0, ## do not show width 0 ranges
                                  col = NULL, ## do not show border and width 0 ranges
                                  stackHeight = 0.75,
                                  rotation.group = 0)
        displayPars(tkCsre) <- as.list(shortNameColors)
      } else {
        tkCsre <- NULL
      }

      tkRefseqToPlot <- tkRefseq
      if (length(subsetByOverlaps(tx, region)) == 0) tkRefseqToPlot <- NULL

      plotTracks(c(tkIdeo, tkGa, tkSeq, tkRefseqToPlot, tkCsre, tkGa),
                 chromosome = chr, from = from, to = to + 1,
                 col.title = "black", cex.title = 1)
    }, execOnResize = FALSE, height = hei()$height)
  })


  output$ft <- DT::renderDataTable({
    region <- rv$subRegion
    grUsed <- subsetByOverlaps(gr, region)
    tb <- data.frame(name = ref$SHORT_NAME[as.integer(match(grUsed$cellType, ref$EID))],
                     region = as.character(grUsed),
                     width = width(grUsed),
                     score = grUsed$score_pd_sh,
                     grUsed$profile_zscore,
                     stringsAsFactors = FALSE)
    color <- rep(markColors, each = nrow(tb))
    as.datatable(formattable(tb, list(
      area(col = c(H3K4me1:H3K9me3)) ~ color_bar_my(color, formattable::normalize),
      score = color_tile("white", "purple")
    )))
  })

  output$highlight <- renderPlot({
    # s <- input$ft_rows_selected
    # region <- rv$subRegion
    # par(mar = c(3, 2.3, 0, 0.6))
    # le <- start(region)
    # re <- end(region)
    # plot(c(le, re), c(0, 1), xaxs = "i", yaxt = "n", type = "n")
    # if (length(s)) {
    #   grHighlight <- subsetByOverlaps(gr, region)[s, ]
    #   lefts <- start(grHighlight)
    #   rights <- end(grHighlight)
    #   co <- eidColors[as.character(grHighlight$cellType)]
    #   rect(lefts, 0, rights, 1, col = co)
    # }

    s <- input$ft_rows_selected
    region <- rv$region
    subRegion <- rv$subRegion
    chr <- as.character(seqnames(region))
    from <- start(region)
    to <- end(region)
    grPlot <- subsetByOverlaps(gr, subRegion)[s, ]
    grPlot$name <-
      Rle(as.character(ref$SHORT_NAME)[match(as.character(grPlot$cellType),
                                             ref$EID)])
    tkHighlight <- AnnotationTrack(grPlot,
                                   feature = grPlot$name,
                                   # id = grPlot$name,
                                   # showFeatureId = TRUE,
                                   # cex = 0.9,
                                   fontcolor.item = "white",
                                   group = grPlot$cellType,
                                   groupAnnotation = "feature",
                                   showId = TRUE,
                                   cex.group = 1,
                                   just.group = "above",
                                   col.line = "white",
                                   name = "highlight",
                                   cex.group = 1,
                                   size = 0.8,
                                   min.width = 0, ## do not show width 0 ranges
                                   col = NULL, ## do not show border and width 0 ranges
                                   stackHeight = 0.75,
                                   rotation.group = 0)
    displayPars(tkHighlight) <- as.list(shortNameColors)
    plotTracks(tkHighlight,
               chromosome = chr, from = from, to = to + 1,
               col.title = "black", cex.title = 1)
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("CSRE_", as.character(unstrand(rv$subRegion)), ".txt")
    },
    content = function(file) {
      region <- rv$subRegion
      grUsed <- subsetByOverlaps(gr, region)
      tb <- data.frame(name = ref$SHORT_NAME[as.integer(match(grUsed$cellType, ref$EID))],
                       region = as.character(grUsed),
                       width = width(grUsed),
                       score = grUsed$score_pd_sh,
                       grUsed$profile_zscore,
                       stringsAsFactors = FALSE)
      write.table(tb, file,
                  sep = "\t",
                  row.names = FALSE, col.names = TRUE)
    },
    contentType = "text/plain"
  )

  output$region <- renderText({
    as.character(rv$region)
  })

  output$pos <- renderText({
    region <- rv$region
    x <- input$plot_hover$x
    pos <- floor(start(region) +
                   width(region) * (x - mi) / (ma - mi))
    if (is.null(x)) {
      NULL
    } else {
      paste0("Mouse position: ", seqnames(region), ":", pos)
    }
  })

  # output$pos <- renderPrint({
  #   region <- rv$region
  #   txInRegion <- subsetByOverlaps(tx, region)
  #   start(txInRegion) <- pmax(start(region), start(txInRegion))
  #   end(txInRegion) <- pmin(end(region), end(txInRegion))
  #   numRefseq <- max(disjointBins(txInRegion))
  #   grUsed <- subsetByOverlaps(gr, region)
  #   numCt <- length(unique(grUsed$cellType))
  #   height <- chromPx + refseqPx * numRefseq + csrePx * numCt
  #   print(numRefseq)
  #   print(numCt)
  #   print(height)
  #   print(input$plot_hover)
  #   print(outputOptions(output))
  # })

  output$subRegion <- renderText({
    as.character(rv$subRegion)
  })


  observeEvent(input$plot_brush, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      s <- start(rv$region)
      w <- width(rv$region)
      subs <- s + w * (brush$xmin - mi) / (ma - mi)
      sube <- s + w * (brush$xmax - mi) / (ma - mi)
      region <- trim(GRanges(seqnames(rv$region),
                             IRanges(subs, sube),
                             seqinfo = seqInfo))
      rv$subRegion <- region
    } else {
      rv$subRegion <- rv$region
    }
  })

  observeEvent(session$clientData$output_plotCSRE_width, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      s <- start(rv$region)
      w <- width(rv$region)
      subs <- s + w * (brush$xmin - mi) / (ma - mi)
      sube <- s + w * (brush$xmax - mi) / (ma - mi)
      region <- trim(GRanges(seqnames(rv$region),
                             IRanges(subs, sube),
                             seqinfo = seqInfo))
      rv$subRegion <- region
    } else {
      rv$subRegion <- rv$region
    }
  })

  observeEvent(input$plot_dblclick, {
    brush <- input$plot_brush
    if (!is.null(brush)) {
      s <- start(rv$region)
      w <- width(rv$region)
      subs <- s + w * (brush$xmin - mi) / (ma - mi)
      sube <- s + w * (brush$xmax - mi) / (ma - mi)
      region <- trim(GRanges(seqnames(rv$region),
                             IRanges(subs, sube),
                             seqinfo = seqInfo))
      rv$region <- region
      rv$subRegion <- rv$region
    }

  })

})
