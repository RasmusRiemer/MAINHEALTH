##################
## load colors ##
#################
#Coulour pallettes

library(RColorBrewer)

col_fil <- brewer.pal(8, "Dark2")#pal_jco("default")(10)

col_scale <- scale_fill_distiller(palette = "Dark2")

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
###########################################
## Selection from original Rmd.R-script ##
###########################################


# Save/load function for 'heavy' objects
sl <- function(name, ..., overwrite = FALSE, dir_path = here::here("results", "RData", subdir_name)) {
  # Possibility to add name as name or literal character string
  name <- as.character(substitute(name))
  assign(name, 
    if(file.exists(glue::glue("{dir_path}/{name}.Rds")) && !overwrite) {
     readRDS(glue::glue("{dir_path}/{name}.Rds"))
    }
    else { 
     dir.create(dir_path, showWarnings = F, recursive = T)
     saveRDS(..., file=glue::glue("{dir_path}/{name}.Rds"))
     readRDS(glue::glue("{dir_path}/{name}.Rds"))
    }, envir=.GlobalEnv)
}

# Save/load function for 'heavy' MicroEco objects
sl_meco <- function(name, ..., overwrite = FALSE, dir_path = here::here("results", "RData", subdir_name)) {
  # Possibility to add name as name or literal character string
  name <- as.character(substitute(name))
  (...)
  assign(name, 
         if(file.exists(glue::glue("{dir_path}/{name}.Rdata")) && !overwrite) {
           load(glue::glue("{dir_path}/{name}.Rdata"))
         }
         else { 
           dir.create(dir_path, showWarnings = F, recursive = T)
           save(object, file=glue::glue("{dir_path}/{name}.Rdata"))
           load(glue::glue("{dir_path}/{name}.Rdata"))
         }, envir=.GlobalEnv)
}


# Format printing numbers with comma separating thousands
f <- function(x, ...) {
  format(x = x, big.mark=",", ...)
}

# Round
r <- function(x, ...) {
  round(x = x, digits = 3, ...)
}

# Keep trailing zero's
# Significant figures: # http://stackoverflow.com/questions/3245862/
s <- function(vec, digits = 3, format = "f", flag = ""){
  return(gsub("\\.$", "", 
              formatC(vec, 
                      digits = digits, 
                      # use "fg" for significant digits
                      format = format, 
                      flag = flag)))
}

pval_lab <- function(p.value) {
  if_else(p.value < 0.001,
          paste("list(~italic(p)==", gsub("e", " %*% 10^", scales::scientific_format(digits = 2)(p.value)), ")", sep = ""), 
          paste("list(~italic(p)==", format(round(p.value, 2), nsmall = 2), ")", sep = ""))
}

###########
## mb ##
###########

# Conversion functions

to_RA <- function(ps) {
  phyloseq::transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
}

pres_abund_filter <- function(ps, pres = 2, abund = 0.001, verbose = TRUE) { # Subramanian filter
  is_raw <- max(phyloseq::otu_table(ps)) > 1 # detect is ps has raw reads
  
  if(is_raw) { # if raw reads; convert to RA for filtering
    ps_raw <- ps
    ps <- ps_raw %>% to_RA()
  }
  
  ps_filt <- phyloseq::filter_taxa(ps, function(x) sum(x > abund) >= pres, TRUE)
  
  if(verbose) {
    message(glue::glue("A total of {phyloseq::ntaxa(ps_filt)} ASVs were found to be present at or above a level of confident detection ({(abund * 100)}% relative abundance) in at least {pres} samples (n = {phyloseq::ntaxa(ps) - phyloseq::ntaxa(ps_filt)} ASVs excluded)."))
  }
  
  if(!is_raw) {
    return(ps_filt)
  } else {
    return(phyloseq::prune_taxa(phyloseq::taxa_names(ps_filt), ps_raw))
  }
}

prev_abund_filter <- function(ps, prev = 0.05, abund = 0.0001) {
  n <- floor(nsamples(ps) * prev) # prevalence at least 5% across time points
  print(n)
  
  ps_filt <- ps %>%
    pres_abund_filter(pres = n, abund = abund) 
  return(ps_filt)
}

mean_abund_filter <- function(ps, mean_abund = 0.001) {
  filter_taxa(ps, function(x) mean(x) > mean_abund, TRUE) } 

asinsqrt <- function(otu_table) {
  asin(sqrt(otu_table))
}

top_n_mean <- function (ps, n = ntaxa(ps)) {
  names(sort(rowMeans(as(otu_table(ps), "matrix")), decreasing = TRUE)[seq_len(n)])
}

get_topn <- function(ps, n = 15, residuals = TRUE) {
  otu_tab <- phyloseq::otu_table(ps)
  otu_tab_n <- otu_tab[order(rowMeans(otu_tab), decreasing = TRUE)[1:n], ]
  otu_tab_res <- otu_tab[-order(rowMeans(otu_tab), decreasing = TRUE)[1:n], ]
  
  if(residuals) {
    otu_tab_n <- otu_tab_n %>%
      t %>%
      data.frame(residuals = colSums(otu_tab_res), check.names = FALSE) %>%
      t
  }
  return(phyloseq::phyloseq(phyloseq::otu_table(otu_tab_n, taxa_are_rows = TRUE),
                            phyloseq::sample_data(ps)))
}

prep_bar <- function(ps, n, residuals = TRUE) {
  excl_cols <- c("sample_id", colnames(phyloseq::sample_data(ps)))
  
  df_topn <- ps %>%
    get_topn(n = n, residuals = residuals) %>%
    ps_to_df(sample_name = "sample_id") %>%
    tidyr::pivot_longer(-dplyr::all_of(excl_cols), names_to = "OTU", values_to = "value") %>%
    dplyr::mutate(OTU = format_OTU(.data$OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev()) %>%
    dplyr::arrange(.data$sample_id) %>%
    dplyr::mutate(sample_id = forcats::fct_inorder(.data$sample_id))
  
  return(df_topn)
}

otu_tab_to_df <- function(ps, sample_name = "sample_id") {
  otu_tab <- phyloseq::otu_table(ps)
  df <- otu_tab %>%
    t %>%
    data.frame(check.names = F) %>%
    tibble::rownames_to_column(sample_name)
  return(df)
}

meta_to_df <- function(ps, sample_name = "sample_id") {
  meta <- phyloseq::sample_data(ps)
  df <- meta %>%
    data.frame(check.names = F) %>%
    tibble::rownames_to_column(sample_name)
  return(df)
}

ps_to_df <- function(ps, sample_name = "sample_id") {
  df_meta <- meta_to_df(ps, sample_name)
  df_otu <- otu_tab_to_df(ps, sample_name)
  
  dplyr::left_join(df_meta, df_otu, by = sample_name) %>% tibble::as_tibble()
}

log10_px <- function(ps, pseudocount = 1) {
  ps %>% transform_sample_counts(., function(x) log10(x + pseudocount)) }

format_OTU <- function(OTU_names, short = F, parse = F) {
  
  OTU_names_tb <- tibble::tibble(OTU_names)
  
  OTU_names_num <- OTU_names_tb %>%
    dplyr::mutate(num = purrr::map_chr(OTU_names, ~stringr::str_extract(., "(?<=_)[0-9]+$")))
  
  OTU_names_tb_format <- OTU_names_num %>%
    dplyr::filter(!is.na(.data$num) & !duplicated(OTU_names)) %>%
    dplyr::mutate(
      name1 = unlist(purrr::map(stringr::str_split(OTU_names, "__|_"), ~ head(.x, -1) %>% glue::glue_collapse(., " "))), #everything except last
      name2 = unlist(purrr::map(stringr::str_split(OTU_names, "__|_"), ~ head(.x, 1))), #last
      long_name = as.character(glue::glue("*{name1}* ({num})")),
      short_name = as.character(glue::glue("*{name2}* ({num})")),
      parse_name = as.character(glue::glue("italic('{name2}')~({num})")))
  
  OTU_names_nonum <- OTU_names_num %>%
    dplyr::filter(is.na(.data$num) & !duplicated(OTU_names))
  
  OTU_names_final <- dplyr::left_join(
    OTU_names_tb,
    dplyr::bind_rows(OTU_names_tb_format, OTU_names_nonum) %>%
      dplyr::mutate(dplyr::across(dplyr::ends_with("_name"), ~dplyr::if_else(is.na(.data$num), stringr::str_to_title(OTU_names), .x)))
    , by = "OTU_names")
  
  if(!short & !parse) {
    return(OTU_names_final %>% dplyr::pull(.data$long_name))
  } else if(parse) {
    return(OTU_names_final %>% dplyr::pull(.data$parse_name))
  } else {
    return(OTU_names_final %>% dplyr::pull(.data$short_name)) }
}

# Plotting functions

prep_bar <- function(ps, n, residuals = TRUE, rename_OTU = TRUE) {
  excl_cols <- c("sample_id", colnames(phyloseq::sample_data(ps)))
  
  OTU_namer <- function(OTU, rename_OTU = T) {
    if(rename_OTU) {
      OTU_adj <- format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev() 
    } else {
      OTU_adj <- OTU %>% forcats::fct_inorder() %>% forcats::fct_rev()
    }
    return(OTU_adj)
  }
  
  df_topn <- ps %>%
    get_topn(n = n, residuals = residuals) %>%
    ps_to_df(sample_name = "sample_id") %>%
    tidyr::pivot_longer(-dplyr::all_of(excl_cols), names_to = "OTU", values_to = "value") %>%
    dplyr::mutate(OTU = OTU_namer(OTU, rename_OTU)) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::mutate(sample_id = forcats::fct_inorder(sample_id))
  
  return(df_topn)
}

make_color_scheme <- function(name, n) {
  max_n <-  RColorBrewer::brewer.pal.info %>%
    tibble::rownames_to_column("name_pal") %>%
    dplyr::filter(.data$name_pal == name) %>%
    dplyr::pull(.data$maxcolors)
  
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(
    dplyr::case_when(n < 3 ~ 3, n < max_n ~ n, TRUE ~ as.numeric(max_n)), name))(n)
}


create_bar <- function(ps = NULL, df_topn = NULL, id = "sample_id", y = "value", n = 15, ncol_legend = 3, name_legend = "OTU", RA = TRUE, colour = "white", residuals = T, rename_OTU = T) {
  
  # accepts either ps (running prep_bar) or a dataframe already prepared with prep_bar
  if(any(class(ps)=="phyloseq")) { df_topn <- prep_bar(ps = ps, n = n, residuals = residuals, rename_OTU = rename_OTU) }
  
  if(residuals) { 
    cols <- c("grey90", rev(make_color_scheme("Paired", n))) 
  } else {
    cols <- rev(make_color_scheme("Paired", n))
  }
  
  bar <- df_topn %>%
    ggplot2::ggplot(ggplot2::aes_string(x = id, y = y, fill = "OTU")) +
    ggplot2::geom_bar(stat = "identity") + #colour = colour
    ggplot2::scale_fill_manual(name = name_legend, values = cols) +
    ggplot2::labs(x = "" , y = "Number of reads") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom", legend.text = ggtext::element_markdown(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_legend))
  
  if(RA) {
    bar <- bar +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
      ggplot2::ylab("Relative abundance")
  }
  return(bar)
  #https://bookdown.org/rdpeng/RProgDA/non-standard-evaluation.html
  
  #TODO: add y = value to plot mean values 
}

hc_order <- function(ps) {
  ps.m <- ps %>%
    pres_abund_filter(pres = 1, verbose = F) %>% # prefilter for faster BC-calculation
    otu_tab_to_df() %>%
    column_to_rownames("sample_id") %>%
    as.matrix()
  
  bc <- vegan::vegdist(ps.m, method = "bray")
  bc.m <-as.matrix(bc)
  hc <- hclust(bc, method="average") 
  
  ps_reorder_samples <- function(ps, order_by) { # order_by is numeric
    sample_data(ps) <- sample_data(ps)[order_by, ]
    otu_table(ps) <- otu_table(ps)[, order_by]
    return(ps)
  }
  
  ps_ord <- ps_reorder_samples(ps, hc$order)
  sample_data(ps_ord)$hc_sample_id <- fct_inorder(rownames(sample_data(ps_ord)))
  return(ps_ord)
}

create_dendro_gg <- function(hc, hang_height=0.05) {
  
  library(dendextend)
  library(ggdendro)
  
  dendro_data_bl <- hc %>% as.dendrogram %>% hang.dendrogram(hang_height=hang_height) %>% dendro_data
  dendro_data_bl$segments$yend[dendro_data_bl$segments$yend<0] <- 0
  
  dendro_data_gr <- hc %>% as.dendrogram %>% dendro_data
  hc_order <- dendro_data_gr$labels$label
  
  plot <- ggplot(segment(dendro_data_gr)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend),colour="grey75") + 
    geom_segment(data=segment(dendro_data_bl),aes(x=x, y=y, xend=xend, yend=yend)) +
    scale_x_continuous(expand = rep(1/length(hc_order)/2, 2)) + 
    scale_y_continuous(expand=c(0,0.02)) +
    #theme(plot.margin=unit(c(0,0,0,0),"lines")) +
    theme_void()
  return(list(plot=plot, hc_order=hc_order)) 
}

meta_bar <- function(data, var, name = NULL, legend_name = NULL, color_scale = NULL, ...) {
  data <- data %>% mutate(index=1:nrow(.))
  
  p <- ggplot(data, aes_string(x = "index", y = 1, fill = var)) + 
    geom_tile() + 
    scale_y_continuous(expand = c(0,0), breaks = 1, labels = name) + 
    scale_x_continuous(expand = c(0,0)) + 
    theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.key.size = unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0, 0, 0, 0), "lines")) +
    labs(fill = legend_name)
  
  if(!is.null(color_scale)) { 
    p <- p + color_scale
  } else {
    n_col <- length(levels(data[[var]])) %>% as.numeric()
    p <- p + scale_fill_manual(values = c(make_color_scheme("BrBG", n_col)[-1], "grey70"))
  }
  return(p)
}

create_age_group <- function(var, name = "virus")  {
  case_when(var <= 45 ~ "≤m1",
            var > 45 & var <= 110 ~ "m1-m3",
            var > 110 & var <= 200 ~ "m3-m6",
            ((var > 200) | is.na(var)) ~ str_c(">m6/no ", name)) %>%
    factor(levels = c("≤m1", "m1-m3", ">m3", "m3-m6", str_c(">m6/no ", name)))
}

tax_glom_rename <- function(ps, taxrank = "genus") {
  ps_glom <- ps %>%
    tax_glom(taxrank = taxrank)
  
  taxa_names(ps_glom) <- tax_table(ps_glom)[,taxrank]
  
  return(ps_glom)
}

effect.size.variables <- c("Batch",
                           "DNA.conc",
                           "Days",
                           "Maternal.age",
                           "Siblings",
                           "ppBMI",
                           "Gestational.age",
                           "C.section",
                           "GDM",
                           #"Induction",
                           "Infant.sex",
                           "Birth.weight",
                           #"Blood.type",
                           #"Placental.weight",
                           #"Maternal.AB.pre.post.birth",
                           #"HMO.group",
                           #"delta.haz",
                           "delta.whz",
                           #"delta.waz",
                           "AB.mother",
                           #"AB.mother.since_birth",
                           "AB.infant",
                           #"AB.infant_since_birth",
                           "Breast.infection",
                           #"Breast.infection.since.birth",
                           #"Formula",                        
                           "Breastfeeding.issues",
                           #"BMI.delta.preg",
                           "BMI.delta.3m",
                           "delta.waz",
                           "delta.haz",
                           "delta.whz",
                           "Secretor",
                           "Lewis"
                           )

get_tax_rank_best <- function(ps = "PSB",n_show = 20){

ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

ps.rel@sam_data <- sample_data(ps.rel)[,1]

#Create melted dataframe
df <- psmelt(ps.rel)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df_best <- df

df_best$best_species <- apply(df, 1, function(x) tail(na.omit(x), 1))

df_best$best_genus <- df %>% select(-"Species") %>% apply(1, function(x) tail(na.omit(x), 1))

df_best$best_family <-  df %>% select(-c("Species","Genus")) %>% apply(1, function(x) tail(na.omit(x), 1))

for (taxa in c("best_species","best_genus","best_family")){
  print(noquote(taxa))
}

for (taxa in c("best_species","best_genus","best_family")){
  #Arrange samples by mean abundance
  top <- df_best %>%
    group_by(!!ensym(taxa)) %>%
    dplyr::summarize(Mean = mean(Abundance)) %>%
    arrange(-Mean)

top_n <- top[taxa][1:n_show,]

#taxa_other <- paste(taxa,"other",sep="_")
df_best <- df_best %>%
  mutate("{taxa}_other" := fct_other(!!ensym(taxa), c(as.matrix(top_n))))
}

return(df_best)

}

# Return taxonomy dataframe with highest level taxonomy per entry
tax_best <- function(ps = "phyloseq object"){
  tax <- as.data.frame(tax_table(ps))
  tax[tax==""] <- NA
  
  tax <- tax %>%
    rowwise() %>%
    mutate(best = tail(na.omit(c_across(cols = everything())), 1)) %>%
    ungroup() %>%
    dplyr::mutate(best = ifelse(is.na(Genus), paste(best, "unclassified"), best)) %>%
    group_by(best) %>%
    mutate(duplicate_count = row_number()) %>%
    ungroup() %>%
    mutate(best = if_else(duplicate_count > 1, 
                                 paste(best, duplicate_count, sep = "_"), 
                          best)) %>%
  
  return(tax)
}
