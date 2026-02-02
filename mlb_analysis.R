#! /usr/bin/Rscript

library(magrittr)
library(stringr)
library(dplyr)
library(tibble)

gname <- 'ERVFRD-1'
gid <- 'ENSG00000244476'
tename <- 'HERVFRD_6p24.2'

calcCPM <- function(df) {
    m <- as.matrix(df)
    cpmdf <- data.frame(t(t(m) / colSums(m)) * 1e6)
    names(cpmdf) <- names(df)
    cpmdf
}


# gtex sample meta
if(!file.exists('meta.gtex.rds')) {
    saveRDS(
        (
            meta.gtex <- read.delim(
                "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",
                sep='\t',
                header=TRUE
            )
        ),
        file="meta.gtex.rds"
    )
} else {
    meta.gtex <- readRDS('meta.gtex.rds')
}

# gtex counts
counts.gtex <- arrow::read_parquet(
        'GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_reads.parquet'
    )
fmeta.gtex <- counts.gtex %>%
    select(gene_id = Name, gene_name=Description) %>%
    mutate(gene_id_prefix=str_remove(gene_id, '\\.\\d+$'))

stopifnot(!any(duplicated(fmeta.gtex$gene_id_prefix)))

mat.gtex <- counts.gtex %>%
    mutate(rowname=str_remove(Name, '\\.\\d+$')) %>%
    select(-c(Name,Description)) %>%
    column_to_rownames() %>%
    as.matrix()

stopifnot(
    all(rownames(mat.gtex) == fmeta.gtex$gene_id_prefix)
)

rm(counts.gtex)
gc()

sampid.gtex <- colnames(mat.gtex)
stopifnot(
    all(sampid.gtex %in% meta.gtex$SAMPID),
    gname %in% fmeta.gtex$gene_name
)

libsize.gtex <- colSums(mat.gtex)
cpm.gtex <- t((t(mat.gtex) / libsize.gtex) * 1e6)
rm(mat.gtex)
gc()

### load this study data
local({
    load('1.count_data.RData')

    meta.cm <<- metadata %>%
        rownames_to_column() %>%
        mutate(response_sample = rowname %in% rownames(response_samples)) %>%
        mutate(
            response = case_when(
                !response_sample ~ NA,
                treatment_best_response %in% c('Clinical Progressive Disease','Stable Disease') ~ 'non-responder',
                treatment_best_response %in% c('Complete Response', 'Partial Response') ~ 'responder',
                .default = 'err'
            )
        )

    fmeta.cm <<- gene_annot %>%
        rownames_to_column() %>%
        mutate(gene_id_prefix=str_remove(gene_id, '\\.\\d+$')) %>%
        select(gene_id, gene_name, gene_type, display, gene_id_prefix)

    counts.cm <<- raw_gene %>%
        rownames_to_column() %>%
        mutate(new_rowname= str_remove(rowname, '\\.\\d+$')) %>%
        column_to_rownames('new_rowname')
})

stopifnot(
    gname %in% fmeta.cm$gene_name,
    all(rownames(counts.cm) %in% fmeta.cm$gene_id_prefix)
)
mat.cm <- counts.cm %>% select(-rowname) %>% as.matrix()
libsize.cm <- colSums(mat.cm)
cpm.cm <- t((t(mat.cm) / libsize.cm) * 1e6)

rm(counts.cm, mat.cm)
gc()




### categories
cat0.levels <- c("TCGA-SKCM", "GTEX")
cat1.levels <- c(
    "Primary Tumor",
    "Metastatic",
    "Blood",
    "Brain",
    "Adipose Tissue",
    "Muscle",
    "Blood Vessel",
    "Heart",
    "Thyroid",
    "Kidney",
    "Uterus",
    "Vagina",
    "Breast",
    "Skin",
    "Salivary Gland",
    "Adrenal Gland",
    "Lung",
    "Spleen",
    "Pancreas",
    "Esophagus",
    "Stomach",
    "Colon",
    "Small Intestine",
    "Prostate",
    "Testis",
    "Nerve",
    "Liver",
    "Pituitary",
    "Ovary",
    "Bladder",
    "Cervix Uteri",
    "Fallopian Tube"
)


### cm samples
samples.cm <- tibble(name = colnames(cpm.cm)) %>%
    left_join(meta.cm, by=join_by(name == rowname)) %>%
    select(
        name,
        cat0 = Project.ID,
        cat1 = Sample.Type,
        cat2 = REGIONAL_VS_PRIMARY,
        cat3 = response
    ) %>%
    mutate(cat1 = ifelse(cat1=="Additional Metastatic", "Metastatic", cat1)) %>%
    mutate(cat2 = case_when(
        cat1=="Metastatic" & is.na(cat2) ~ "Metastatic",
        cat1=="Primary Tumor" & is.na(cat2) ~ "Primary",
        cat1=="Primary Tumor" & cat2=="Primary_Disease" ~ "Primary",
        .default = cat2
    )) %>%
    mutate(
        cat3 = factor(cat3, levels=c('responder','non-responder'))
    )

# > ggsci::pal_material("amber")(10)
# [1] "#FFF8E0FF" "#FFEBB2FF" "#FFDF81FF" "#FFD44EFF" "#FFCA27FF" "#FFC006FF" "#FFB200FF" "#FF9F00FF" "#FF8E00FF" "#FF6E00FF"
# > ggsci::pal_material("deep-orange")(10)
# [1] "#FAE9E6FF" "#FFCCBBFF" "#FFAB91FF" "#FF8A65FF" "#FF7043FF" "#FF5721FF" "#F3511EFF" "#E54A19FF" "#D84314FF" "#BF350CFF"

pal.tcga <- c(
    `Primary Tumor` = "#FFD44EFF",
    `Metastatic` = "#FF5721FF",
    `Primary` = "#FFD44EFF",
    `Distant_Metastases` = "#FF5721FF",
    `Regional_Lymph_Node` = "#BF350CFF",
    `Regional_Skin_or_Soft_Tissue` = "#E54A19FF"
)

cat.cm <- samples.cm %>%
    count(cat1, cat2) %>%
    mutate(
        cat1.col = pal.tcga[cat1],
        cat2.col = pal.tcga[cat2]
    ) %>%
    mutate(
        cat0="TCGA-SKCM",
        cat1=factor(cat1, levels=cat1.levels)
    ) %>%
    arrange(cat1, cat2) %>%
    select(cat0,cat1,cat1.col,cat2,cat2.col)

# gtex samples
samples.gtex <- tibble(SAMPID = colnames(cpm.gtex)) %>%
    left_join(meta.gtex, by='SAMPID') %>%
    mutate(cat0 = "GTEX") %>%
    select(
        name = SAMPID,
        cat0,
        cat1 = SMTS,
        cat2 = SMTSD
    ) %>%
    mutate(
        cat2 = case_match(
            cat2,
            "Cells - Cultured fibroblasts" ~ "Cells - Transformed fibroblasts",
            .default = cat2
        )
    ) %>%
    mutate(cat3 = factor(NA, levels=c('responder','non-responder')))


pal.gtexvis <- sapply(
    jsonlite::fromJSON("https://github.com/broadinstitute/gtex-viz/raw/f9d18299b2bcb69e4a919a4afa359c99a33fbc3b/boxplot/dev/colors.json"),
    function(jobj) str_c('#', jobj$tissue_color_hex)
)

cat.gtex <- samples.gtex %>%
    count(cat1,cat2) %>%
    mutate(
        cat2alt = str_replace(cat2, " - [\\w ]+$", "")
    ) %>%
    mutate(
        cat1.col.orig = pal.gtexvis[cat1],
        cat2.col = ifelse(
            is.na(pal.gtexvis[cat2]),
            pal.gtexvis[cat2alt],
            pal.gtexvis[cat2]
        )
    ) %>%
    select(-cat2alt)

stopifnot(
    cat.gtex %>% filter(is.na(cat2.col)) %>% nrow() == 0
)

tmppal.cat1 <- cat.gtex %>%
    group_by(cat1,cat2) %>%
    summarize(n=sum(n)) %>%
    slice_max(n, n=1) %>%
    left_join(
        cat.gtex %>% select(cat2, cat1.col=cat2.col),
        by='cat2'
    ) %>%
    pull(cat1.col, name=cat1)

cat.gtex %<>%
    mutate(
        cat1.col = tmppal.cat1[cat1]
    )

stopifnot(
    cat.gtex %>% filter(!is.na(cat1.col.orig) & (cat1.col.orig != cat1.col)) %>% nrow() == 0
)

cat.gtex %<>%
    mutate(cat0="GTEX") %>%
    mutate(cat1=factor(cat1, levels=cat1.levels)) %>%
    arrange(cat1, cat2) %>%
    select(cat0,cat1,cat1.col,cat2,cat2.col)

rm(tmppal.cat1)

### combine
cat.all <- rbind(cat.cm, cat.gtex)
samples.all <- rbind(samples.cm, samples.gtex) %>%
    mutate(
        cat0=factor(cat0,levels=cat0.levels),
        cat1=factor(cat1,levels=cat1.levels),
        cat2=factor(cat2,levels=cat.all$cat2)
    )


# d1 <- setdiff(rownames(cpm.cm), rownames(cpm.gtex))
# d2 <- setdiff(rownames(cpm.gtex), rownames(cpm.cm))
i1 <- intersect(rownames(cpm.cm), rownames(cpm.gtex))

# Keep only genes in both
rm(mat.gtex, mat.cm)
gc()
cpm.both <- cbind(cpm.cm[i1, ], cpm.gtex[i1, ])

stopifnot(
    all(rownames(cpm.both) == i1),
    all(colnames(cpm.both) == c(colnames(cpm.cm),colnames(cpm.gtex))),
    all(samples.all$name == colnames(cpm.both))
)


# sanity check: all values for `gid` are aligned correctly
stopifnot(
    all(cpm.both[gid,] == c(cpm.cm[gid,], cpm.gtex[gid,]))
)

rm(cpm.cm, cpm.gtex)
gc()

### stat
sumstats <- function(gdf){
    gdf %>%
    summarise(
        n = n(),
        mean = mean(value),
        var = var(value),
        p1 = quantile(value, 0.01),
        p5 = quantile(value, 0.05),
        q1 = quantile(value, 0.25),
        median = median(value),
        q3 = quantile(value, 0.75),
        p95 = quantile(value, 0.95),
        p99 = quantile(value, 0.99)
    )
}

summary.cat0 <- tibble::enframe(cpm.both[gid, ]) %>%
    left_join(samples.all, by="name") %>%
    group_by(cat0) %>%
    sumstats()

summary.cat0

summary.cat1 <- tibble::enframe(cpm.both[gid, ]) %>%
    left_join(samples.all, by="name") %>%
    group_by(cat1) %>%
    sumstats() %>%
    arrange(-mean)

summary.cat1

summary.cat2 <- tibble::enframe(cpm.both[gid, ]) %>%
    left_join(samples.all, by="name") %>%
    group_by(cat2, cat1) %>%
    sumstats() %>%
    arrange(-mean)

summary.cat2

summary.cat3 <- tibble::enframe(cpm.both[gid, ]) %>%
    left_join(samples.all, by="name") %>%
    group_by(cat3) %>%
    sumstats() %>%
    arrange(cat3)

summary.cat3




#### Plots
library(ggplot2)

stopifnot(
    fmeta.cm[fmeta.cm$gene_name==gname, 'gene_id_prefix'] == fmeta.gtex[fmeta.gtex$gene_name==gname, 'gene_id_prefix']
)

pal.cat <- list(
    `cat1` = cat.all %>% pull(cat1.col, name=cat1),
    `cat2` = cat.all %>% pull(cat2.col, name=cat2)
)

the.gid <- fmeta.cm[fmeta.cm$gene_name==gname, 'gene_id_prefix']
the.gname <- fmeta.cm[fmeta.cm$gene_name==gname, 'gene_name']

plt.list <- list()

## Plot 1: all tissues, cat1
source("lib.VlnPltPlus.R")
padj.th <- 1e-3

plt.list$all.cat1 <- VlnPlotPlus(
    cpm.both,
    the.gid,
    col.meta = samples.all,
    group_var = "cat1",
    pal.colors = pal.cat[["cat1"]],
    show.boxplot = "mean_cl_boot",
    show.summary = "mean",
    ylab="CPM", title=the.gname
)

pdf("ERVFRD-1.all_tissues.pdf", width=11, height=4)
    plt.list$all.cat1 +
        coord_cartesian(ylim = c(0,3)) +
        VlnPlotPlus.theme() +
        theme(
            axis.text.x = element_text(
                size = rel(0.8),
                angle = 50,
                vjust = 1,
                hjust=1
            )
        )
dev.off()

pdf("ERVFRD-1.all_tissues.notrim.pdf", width=11, height=4)
plt.list$all.cat1 +
    VlnPlotPlus.theme() +
    theme(
        axis.text.x = element_text(
            size = rel(0.8),
            angle = 50,
            vjust = 1,
            hjust=1
        )
    )
dev.off()

## Plot 2: select tissues, cat1
selected.tissues <- c("Primary Tumor", "Metastatic", "Skin", "Liver")
selected.tissues <- c(selected.tissues,
    summary.cat1 %>%
        filter(!(cat1 %in% selected.tissues)) %>%
        slice_max(mean, n=3) %>%
        mutate(cat1char=as.character(cat1)) %>%
        pull(cat1char)
)
plt.list$sel.cat1 <- VlnPlotPlus(
    cpm.both,
    the.gid,
    col.meta = samples.all,
    group_var = "cat1",
    filter_val = selected.tissues,
    pal.colors = pal.cat[["cat1"]],
    show.summary = "mean",
    show.boxplot = "mean_cl_boot",
    ylab="CPM", title=the.gname
) %>%
    VlnPlotPlus.test(
        padj.th = padj.th,
        test.comparisons=list(c('Metastatic'), c("Primary"))
    )

pdf("ERVFRD-1.select_tissues.pdf", width=11, height=5)
    plt.list$sel.cat1 +
        coord_cartesian(ylim = c(0,3)) +
        VlnPlotPlus.theme() +
        theme(
            axis.text.x = element_text(
                size = rel(0.8),
                angle = 50,
                vjust = 1,
                hjust=1
            )
        )
dev.off()
pdf("ERVFRD-1.select_tissues.notrim.pdf", width=11, height=5)
plt.list$sel.cat1 +
    coord_cartesian(ylim = c(0,3)) +
    VlnPlotPlus.theme() +
    theme(
        axis.text.x = element_text(
            size = rel(0.8),
            angle = 50,
            vjust = 1,
            hjust=1
        )
    )
dev.off()




## Plot 3: select tissues, cat2
selected.tissues <- c("Primary Tumor", "Metastatic", "Skin")
plt.list$sel.cat2 <- VlnPlotPlus(
    cpm.both,
    the.gid,
    col.meta = samples.all,
    group_var = "cat2",
    filter_var = "cat1",
    filter_val = selected.tissues,
    pal.colors = pal.cat[["cat2"]],
    show.summary = "mean",
    show.boxplot = "mean_cl_boot",
    ylab="CPM", title=the.gname
) %>%
    VlnPlotPlus.test(
        padj.th = padj.th,
        test.comparisons = (
            summary.cat2 %>%
                filter(cat1 %in% c("Primary Tumor", "Metastatic")) %>%
                mutate(cat2chr = as.character(cat2)) %>%
                pull(cat2chr) %>%
                as.list()
        )
    )

pdf("ERVFRD-1.select_detail.pdf", width=11, height=5)
    plt.list$sel.cat2 +
        coord_cartesian(ylim = c(0,6)) +
        VlnPlotPlus.theme() +
        theme(
            axis.text.x = element_text(
                size = rel(0.8),
                angle = 50,
                vjust = 1,
                hjust=1
            )
        )

    # legend.only(plt.list$sel.cat2)
dev.off()

pdf("ERVFRD-1.select_detail.notrim.pdf", width=11, height=5)
plt.list$sel.cat2 +
    VlnPlotPlus.theme() +
    theme(
        axis.text.x = element_text(
            size = rel(0.8),
            angle = 50,
            vjust = 1,
            hjust=1
        )
    )

# legend.only(plt.list$sel.cat2)
dev.off()


## Plot 3.1: select tissues, select cat2
selected.tissues <- c(
    "Regional_Lymph_Node",
    "Distant_Metastases",
    "Regional_Skin_or_Soft_Tissue",
    "Skin - Sun Exposed (Lower leg)",
    "Skin - Not Sun Exposed (Suprapubic)"
    # "Cells - Transformed fibroblasts",
)
plt.list$specific <- VlnPlotPlus(
    cpm.both,
    the.gid,
    col.meta = samples.all,
    group_var = "cat2",
    filter_val = selected.tissues,
    pal.colors = pal.cat[["cat2"]],
    show.summary = "mean",
    show.boxplot = "mean_cl_boot",
    ylab="CPM", title=the.gname
) %>%
    VlnPlotPlus.test(
        padj.th = padj.th
    )

pdf("ERVFRD-1.specific.pdf", width=11, height=5)
    plt.list$specific +
        coord_cartesian(ylim = c(0,6)) +
        VlnPlotPlus.theme() +
        theme(
            axis.text.x = element_text(
                size = rel(0.8),
                angle = 50,
                vjust = 1,
                hjust=1
            )
        )
    # legend.only(plt.list$specific)
dev.off()

pdf("ERVFRD-1.specific.notrim.pdf", width=11, height=5)
plt.list$specific +
    VlnPlotPlus.theme() +
    theme(
        axis.text.x = element_text(
            size = rel(0.8),
            angle = 50,
            vjust = 1,
            hjust=1
        )
    )
# legend.only(plt.list$specific)
dev.off()


## Plot 4
plt.list$cat3 <- VlnPlotPlus(
    cpm.both,
    the.gid,
    col.meta = samples.all,
    group_var = "cat3",
    filter_var = "cat3",
    filter_val = c("responder", "non-responder"),
    show.summary = "mean",
    show.boxplot = "mean_cl_boot",
    ylab="CPM", title=the.gname
) %>%
    VlnPlotPlus.test(
        padj.th = padj.th,
        show.ns = TRUE
    )

pdf("ERVFRD-1.response.pdf", width=7, height=5)
    plt.list$cat3 +
        coord_cartesian(ylim = c(0,10)) +
        VlnPlotPlus.theme() +
        theme(
            axis.text.x = element_text(
                size = rel(0.8),
                angle = 50,
                vjust = 1,
                hjust=1
            )
        )
    # legend.only(plt.list$cat3)
dev.off()

pdf("ERVFRD-1.response.notrim.pdf", width=7, height=5)
plt.list$cat3 +
    VlnPlotPlus.theme() +
    theme(
        axis.text.x = element_text(
            size = rel(0.8),
            angle = 50,
            vjust = 1,
            hjust=1
        )
    )
# legend.only(plt.list$cat3)
dev.off()



