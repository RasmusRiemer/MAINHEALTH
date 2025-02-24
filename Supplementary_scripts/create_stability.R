create_stability <- function(dist, meta, pair = "infant - mother") {
  
  dist %>%
    as(., "matrix") %>%
    data.frame(check.names = F) %>%
    rownames_to_column("id_from") %>%
    pivot_longer(-id_from, names_to = "id_to", values_to = "dist") %>%
    left_join(., meta  %>% 
                select(sample_id, subject, time, niche = niche_toNP, mother_inf) %>%
                rename_with(~str_c(., "_from")), by = c("id_from" = "sample_id_from")) %>%
    left_join(., meta %>% 
                select(sample_id, subject, time, niche, mother_inf) %>%
                rename_with(~str_c(., "_to")), by = c("id_to" = "sample_id_to")) %>%
    mutate(type = if_else(subject_from == subject_to, "related", "unrelated")) %>%
    unite("time_to_from", time_to, time_from, sep = " - ", remove = F) %>%
    unite("niche_to_from", niche_to, niche_from, sep = " - ", remove = F) %>%
    unite("mi_to_from", mother_inf_to, mother_inf_from, sep = " - ", remove = F) %>%
    filter(mi_to_from == pair) %>%
    mutate(subject_from = fct_inorder(subject_from))
}
