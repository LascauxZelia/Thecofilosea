################################################################################
## Plotting a map of dot/icm systems in the Legionellales
################################################################################
## Libs
library(genoPlotR)

## Constants
artdb <- "/home/zelia/Desktop/Nanopore/T4BSS/artdb"

#### Read data
ids <- c("R64", "LegLon", "LegPhi", "Tara121", "LegOak", "AquLus", "Poke", "BerAqu", "Fisci", "RicIso", "CoxRSA93", "LegPhi", "PisSal", "FanHon", "AciFer", "R64")
labels <- c("R64 IncI plasmid (Salmonella)", "Legionella longbeachae", "Legionella pneumophila", "Tara121", "Legionella oakridgensis", "Aquicella lusitana", "Ca. Pokemonas kadabra", "Ca. Berkiella aquae", "Ca. Fiscibacter pecunius", "Rickettsiella isopodorum", "Coxiella burnetii RSA93", "Legionella pneumophila", "Piscirickettsia salmonis", "Fangia hongkongensis", "Acidithiobacillus ferrivorans", "R64 IncI plasmid (Salmonella)")

## DNA segs
dna_segs <- list()

annots <- list()
for (i in 1:length(ids)){
  path <- paste(artdb, "/", ids[i], ".gbk", sep = "")
  dna_segs[[i]] <- read_dna_seg_from_file(path, gene_type = "arrows",
                                          fill = grey(0.5), col = "black")
  annot <- NULL
  if (ids[i] %in% c("R64", "LegPhi", "CoxRSA93")){
    annot <- annotation(middle(dna_segs[[i]]), text=dna_segs[[i]]$gene,
                        rot = 45)
    annot <- annot[annot$text != "-",]
  }
  annots[[i]] <- annot
  ##annots[[i]] <- auto_annotate(dna_segs[[i]], dna_segs[[i]]$gene, rot = 30)
}

## Comparisons
comps <- list()
for (i in 2:length(ids)){
  path <- paste(artdb, "/", ids[i-1], "_vs_", ids[i], ".tblastx.parsed",
                sep = "")
  comps[[i-1]] <- read_comparison_from_blast(path, length = 20)
}

## Offsets
xtraoffset <- 500
xlimsdf <- read.table("/home/zelia/Desktop/Nanopore/T4BSS/offsets_Fisci.tab", h=T, sep = "\t", stringsAsFactors=F, fill=TRUE)
xlims <- list()
for (i in 1:length(ids)){
  xs <- numeric(0)
  xlimsdf_sub <- xlimsdf[xlimsdf$Genome == ids[i],]
  for (j in 1:nrow(xlimsdf_sub)){
    x1 <- xlimsdf_sub$Start[j]-xtraoffset*xlimsdf_sub$Strand[j]
    if (x1 < 0)
      x1 <- 0
    x2 <- xlimsdf_sub$End[j]+xtraoffset*xlimsdf_sub$Strand[j]
    xs <- c(xs, x1, x2)
  }
  xlims[[i]] <- xs
}

## gene names
icmannots <- read.table("/home/zelia/Desktop/Nanopore/T4BSS/geneNames.tab", h=T, sep = "\t", stringsAsFactors=F)
icmannots$color <- rainbow(nrow(icmannots))
for (i in 1:length(ids)){
  for (j in grep(ids[i], names(icmannots))){
    for (k in 1:nrow(icmannots)){
      locustag <- icmannots[k,j]
      if (locustag != "-"){
        dna_segs[[i]]$fill[dna_segs[[i]]$synonym == locustag] <-
          icmannots$color[k]
      }
    }
  }
}

cairo_pdf("T4BSS_map_230214.pdf", w=14, h=9)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1, 2, widths = unit(c(1,0.15),
                                                             rep("null", 2))),
                      name="overall_vp"))
pushViewport(viewport(layout.pos.col=2, name="legend"))
grid.legend(labels = icmannots$genome, pch = 22,
            gp = gpar(col = "black", fill = icmannots$color),
            vgap = unit(0.3, "lines"))
upViewport()
pushViewport(viewport(layout.pos.col=1, name="genoplotr"))
plot_gene_map(dna_segs, comparisons=comps, xlims=xlims, annotations=annots,
              annotation_height = 1.5, annotation_cex = 0.6,
              dna_seg_labels = labels, plot_new = F)
upViewport(0)
dev.off()
