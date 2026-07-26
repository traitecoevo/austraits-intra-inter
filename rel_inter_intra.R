source("R/fit_model2d_bayesian.R")

traits <- unique(data_traits_climate$trait_name)
traits <- traits[-c(6,8,10)]
trait_data <- list()
for(variable in c("log_annual_precipitation","average_temperature")){
  for (trait in traits) {
    
    trait_data[[trait]] <- readRDS(paste0("outputs/model_fits/",trait[1],"_",variable,".RDS"))[[2]]
    
    trait_data[[trait]] <- posterior::as_draws_matrix(trait_data[[trait]]) %>%
      as_tibble() %>%
      select(1:3) %>%
      mutate(trait_name = trait) %>%
      mutate(environment = variable)
  }
  
  if(trait_name == "leaf_mass_per_area"){trait_label = expression(paste(LMA," (g ",m^{-2},")"))}
  if(trait_name == "huber_value"){trait_label = expression(paste(SA:LA))}
  if(trait_name == "leaf_area"){trait_label = expression(paste(Leaf~area," (",mm^{-2},")"))}
  if(trait_name == "leaf_delta13C"){trait_label = expression(paste(delta^{13},C))}
  if(trait_name == "leaf_N_per_area"){trait_label = expression(paste(N[area]," (g ",m^{-2},")"))}
  if(trait_name == "wood_density"){trait_label = expression(paste(Wood~density," (g ",cm^{-3},")"))}
  if(trait_name == "leaf_P_per_area"){trait_label = expression(paste(P[area]," (g ",m^{-2},")"))}
  
  trait_data %>%
    bind_rows() %>%
    select(2:5) %>%
    rename(`Mn` = `beta[2]`,`C` = `beta[3]`) %>%
    mutate(`Mn - C` = Mn - C) %>%
    pivot_longer(cols = c(Mn, C, `Mn - C`)) %>%
    mutate(trait_name = factor(trait_name, levels = c("leaf_mass_per_area","huber_value","leaf_area",
                                                      "leaf_delta13C","leaf_N_per_area","wood_density","leaf_P_per_area"))) %>%
    ggplot(aes(x = value, y = trait_name, group = trait_name)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggridges::geom_density_ridges(aes(colour = name, fill = name),alpha=0.5,scale = 0.5,size = 1, parse = TRUE) +
    facet_wrap(~name) +
    theme_classic() +
    scale_y_discrete(labels = c(
      bquote(LMA),
      bquote(SA:LA),
      bquote(LA),
      bquote(delta),
      bquote(N),
      bquote(Wood~density),
      bquote(P)
    )) +
    scale_colour_manual(values = c("aquamarine","orange","purple"))+
    scale_fill_manual(values = c("aquamarine","orange","purple"))+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          # axis.ticks.y = element_blank(),
          strip.border = element_rect(color = NA),
          legend.position = "none",
          text = element_text(size = 12),
          axis.text = element_text(size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 12)) 
    
}

leaf_area_lme <- readRDS("outputs/model_fits/leaf_area_log_annual_precipitationPFT.RDS")[[1]]




leaf_area <- readRDS("outputs/model_fits/leaf_area_log_annual_precipitationPFT.RDS")[[2]]
l_n_area <- readRDS("outputs/model_fits/leaf_N_per_area_log_annual_precipitationPFT.RDS")[[2]]
leaf_mass_per_area <- readRDS("outputs/model_fits/_log_annual_precipitationPFT.RDS")[[2]]

leaf_mass_per_area <- posterior::as_draws_matrix(leaf_mass_per_area) %>%
  as_tibble() %>%
  select(1:15) %>%
  mutate(trait_name = "leaf_mass_per_area")

leaf_area_df <- posterior::as_draws_matrix(leaf_area) %>%
  as_tibble() %>%
  select(1:15) %>%
  mutate(trait_name = "leaf_area")

l_n_area_df <- posterior::as_draws_matrix(l_n_area) %>%
  as_tibble() %>%
  select(1:15) %>%
  mutate(trait_name = "l_n_area")

leaf_area_df %>%
  bind_rows(l_n_area_df) %>%
  bind_rows(leaf_mass_per_area)-> trait_df

names(trait_df) <- c(names(coef(leaf_area_lme)$cond$taxon_name), "trait_name")

trait_df %>%
  select(trait_name, contains("mn")) %>%
  pivot_longer(cols = c(-trait_name, -log_annual_precipitation_mn)) %>%
  mutate(real_value = log_annual_precipitation_mn + value) %>%
  select(-value) %>%
  pivot_wider(names_from = "name", values_from = "real_value") %>%
  pivot_longer(-trait_name) %>%
  mutate(name = factor(x = name, levels = c("log_annual_precipitation_mn","log_annual_precipitation_mn:growth_form_simplergraminoid","log_annual_precipitation_mn:growth_form_simplerherb","log_annual_precipitation_mn:growth_form_simplershrub","log_annual_precipitation_mn:growth_form_simplertree"), labels = c("Climber","Graminoid","Herb","Shrub","Tree"))) %>%
  mutate(type = "mn") -> mn

trait_df %>%
  select(trait_name, contains("_c")) %>%
  pivot_longer(cols = c(-trait_name, -log_annual_precipitation_c)) %>%
  mutate(real_value = log_annual_precipitation_c + value) %>%
  select(-value) %>%
  pivot_wider(names_from = "name", values_from = "real_value") %>%
  pivot_longer(-trait_name) %>%
  mutate(name = factor(x = name, levels = c("log_annual_precipitation_c","log_annual_precipitation_c:growth_form_simplergraminoid","log_annual_precipitation_c:growth_form_simplerherb","log_annual_precipitation_c:growth_form_simplershrub","log_annual_precipitation_c:growth_form_simplertree"), labels = c("Climber","Graminoid","Herb","Shrub","Tree"))) %>%
  mutate(type = "c") -> c_data

c_data %>%
  rename("c" = "value") %>%
  select(-type)-> c_data

mn %>%
  rename("mn" = "value") %>% 
  select(-type) -> mn


c_data %>%
  mutate(mn = mn$mn) %>%
  rename() %>%
  mutate("mn - c" = c - mn)  |>
  mutate(trait_name = factor(trait_name, levels = c("l_n_area", "leaf_area", "leaf_mass_per_area"), labels = c("Leaf N per area","LA","LMA"))) %>%
  pivot_longer(cols = c(c,mn,"mn - c"), names_to = "type") %>%
  group_by(trait_name,name,type) %>% 
  summarise(median = median(value), lower_quant = quantile(value, 0.05), upper_quant = quantile(value, 0.95)) %>%
  ungroup()-> quants

quants %>%
  rename(type = name, name =type) -> quants

c_data %>%
  mutate(mn = mn$mn) %>%
  rename() %>%
  mutate(diff = c - mn) %>%
  rename(type = name) %>%
  pivot_longer(cols = c(c,mn,diff)) %>%
  mutate(trait_name = factor(trait_name, levels = c("l_n_area", "leaf_area", "leaf_mass_per_area"), labels = c("Leaf N per area","LA","LMA"))) %>%
  mutate(name = factor(name, levels = c("diff", "c", "mn"), labels = c("mn - c","c","mn"))) %>%
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(data = quants, aes(y = name, yend = name, x = lower_quant, xend = upper_quant), size = 0.5) +
  geom_point(data = quants, aes(y = name, x = median), size = 1.5) +
  ggridges::geom_density_ridges(aes(y = name, x = as.numeric(value), group = name, colour = name, fill = name),alpha=0.5,scale = 0.5,size = 1) +
  facet_grid(type~trait_name) +
  theme_classic() +
  scale_colour_manual(values = c("aquamarine","orange","purple"))+
  scale_fill_manual(values = c("aquamarine","orange","purple"))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        # axis.ticks.y = element_blank(),
        strip.border = element_rect(color = NA),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) 





