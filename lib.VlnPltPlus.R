#!/usr/bin/env Rscript

require(dplyr)
require(magrittr)
require(ggplot2)
require(ggpubr)
require(see)
require(ggpp)

VlnPlotPlus <- function(
    mat,         # data matrix
    rowname,     # row name to plot
    col.meta,    # metadata annotating column
    group_var,   #
    filter_var = NULL,
    filter_val = NULL,
    pal.colors = NULL,
    show.summary = FALSE,
    show.summary.color  = c('#DC0000B2'),
    show.boxplot = FALSE,
    show.boxplot.color = c('#DC0000B2'),
    ylab = "value",
    xlab = group_var,
    title = rowname
) {
    .df <- tibble::enframe(mat[rowname,]) %>%
        dplyr::left_join(samples.all, by="name")

    if (!is.null(filter_val)) {
        if(is.null(filter_var)) filter_var <- group_var
        .df %<>%
            filter(.data[[filter_var]] %in% filter_val) %>%
            mutate({{filter_var}} := factor(.data[[filter_var]], levels=filter_val))
    }

    plt <- (
        ggplot(
            .df,
            aes(.data[[group_var]], value, fill=.data[[group_var]])
        ) +
        see::geom_violinhalf(
            position = position_nudge(x=0.01, y=0),
            scale='width'
        ) +
        geom_point(
            alpha=0.05,
            size=0.4,
            position = ggpp::position_jitternudge(
                width=0.15,
                x=-0.15,
                nudge.from = 'jittered'
            )
        )
    )

    if(!is.null(pal.colors)) {
        plt <- plt + scale_fill_manual(values = pal.colors)
    }

    if(is.logical(show.summary) && show.summary)
        show.summary <- c("median")
    if(is.character(show.summary)) {
        show.summary.color <- rep(show.summary.color, length(show.summary))
        for(i in 1:length(show.summary)) {
            plt <- plt + stat_summary(
                fun = show.summary[i],
                geom = "crossbar",
                width = 0.5,
                color = show.summary.color[i], fill='#00000000',
                position=position_nudge(x=0, y=0),
                show.legend = FALSE
            )
        }
    }

    if(is.logical(show.boxplot) && show.boxplot)
        show.boxplot <- c("median_hilow")
    if(is.character(show.boxplot)) {
        show.boxplot.color <- rep(show.boxplot.color, length(show.boxplot))
        for(i in 1:length(show.boxplot)) {
            plt <- plt + stat_summary(
                fun.data = show.boxplot[i],
                geom="errorbar",
                width=0.35,
                color=show.boxplot.color[i],
                position=position_nudge(x=0, y=0),
                show.legend = FALSE
            )
        }
    }

    plt + labs(x=xlab, y=ylab, title=title)
}

VlnPlotPlus.theme <- function(...) {
    list(
        ggpubr::theme_pubclean(...),
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=rel(0.8)),
            legend.title = element_blank()
        ),
        theme(
            legend.position = 'none',
            validate = TRUE
        ),
        theme(
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    )
}

VlnPlotPlus.test <- function(ggplt, padj.th=1e-3, test.comparisons=NULL, show.ns=FALSE, ...) {
    fm <- as.formula(str_c(as_label(ggplt$mapping$y), ' ~ ', as_label(ggplt$mapping$x)))
    res.kruskal <- ggpubr::compare_means(
        fm,
        data=ggplt$data,
        method="kruskal.test"
    )
    if(res.kruskal$p.adj[1] < padj.th){
        ggplt <- ggplt + ggpubr::stat_compare_means(
            method='kruskal.test',
            label.y.npc='top',
            label.x.npc = 'left',
            size = rel(3)
        )
        res.wilcox <- ggpubr::compare_means(
            fm,
            data=ggplt$data,
            method = "wilcox.test",
            p.adjust.method = "holm"
        )
        if(is.null(test.comparisons)) {
            if(show.ns) {
                # get all comps
                comps <- res.wilcox %>%
                    rowwise() %>%
                    mutate(x=list(c(group1,group2))) %>%
                    pull(x)
            } else {
                # get all < padj.th
                comps <- res.wilcox %>%
                    filter(p.adj < padj.th) %>%
                    rowwise() %>%
                    mutate(x=list(c(group1,group2))) %>%
                    pull(x)
            }
        } else {
            sel.wilcox <- lapply(test.comparisons, function(v) {
                if(length(v)==2) {
                    res.wilcox %>%
                        filter((group1==v[1] & group2==v[2]) | (group1==v[2] & group2==v[1]))
                } else if (length(v)==1) {
                    res.wilcox %>%
                        filter((group1==v[1]) | (group2==v[1]))
                }
            }) %>%
                bind_rows %>%
                distinct()
            if(show.ns) {
                # get all comps
                comps <- sel.wilcox %>%
                    rowwise() %>%
                    mutate(x=list(c(group1,group2))) %>%
                    pull(x)
            } else {
                # get all < padj.th
                comps <- sel.wilcox %>%
                    filter(p.adj < padj.th) %>%
                    rowwise() %>%
                    mutate(x=list(c(group1,group2))) %>%
                    pull(x)
            }
        }
        if(length(comps)>0) {
            ggplt <- ggplt + ggpubr::stat_compare_means(
                comparisons = comps,
                label.y.npc='middle',
                label.x.npc='middle',
                # tip.length = 0.03, # height that bar goes down
                bracket.size = 0.2, # width of bracket lines
                step.increase = rel(0.2),
                size=rel(2)
            )
        }
    } else {
        if(show.ns) {
            ggplt <- ggplt + ggpubr::stat_compare_means(
                method='kruskal.test',
                label.y.npc='top',
                label.x.npc = 'left',
                size = rel(3)
            )
        }
    }
    ggplt
}

legend.only <- function(p) {
    cowplot::ggdraw(suppressWarnings({
        cowplot::get_plot_component(
            (p + guides(colour = guide_legend(override.aes = list(size=3))) +
                 theme(
                     legend.key.size=unit(0.1, 'points'),
                     legend.text=element_text(size=7)
                 )
            ),
            pattern = "guide-box",
            return_all=FALSE
        )
    }))
}

