# libs ----------------------------------------------------------------------- #

library(gamlss)
library(gamlss.ggplots)
library(dplyr)
library(readr)
library(ggspatial)
library(sf)
library(ggplot2)
library(summarytools)
library(tidyr)
library(caret)
library(patchwork)
library(colorspace)

theme_set(theme_bw())

# loading data --------------------------------------------------------------- #

data <- read_csv('data/data.csv')

## variables used for modeling

pxrf_vars <- data %>%
  dplyr::select(starts_with('pXRF_')) %>%
  colnames()
dry_color_vars <- c('LAB_L_dry', 'LAB_A_dry', 'LAB_B_dry')
pca_vars <- c('PC1', 'PC2')
cat_vars <- c('Depth', 'PM')
aux_vars <- c(
  'Profile', 'Lat', 'Long', 'Depth', 'Symbol',
  'SiBCS', 'Soil taxonomy', 'Soil order', 'PM'
)
target_vars <- c('Lab_SOM')

# pca ------------------------------------------------------------------------ #

pxrf_pca <- prcomp(data %>% select(all_of(c(pxrf_vars))), scale = T)
pca_imp <- summary(pxrf_pca)$importance[2,]

scale_factor <- 25
x_loadings <- pxrf_pca$rotation[,1] * scale_factor
y_loadings <- pxrf_pca$rotation[,2] * scale_factor
x_lab <- paste0('PC1 (', round(pca_imp[1]*100, 1), '% variance explained)')
y_lab <- paste0('PC2 (', round(pca_imp[2]*100, 1), '% variance explained)')
remove_pXRF <- function(x) stringr::str_extract(x, '(?<=pXRF_).*')

ggplot() +
  xlab(x_lab) + ylab(y_lab) +
  geom_point(
    data = pxrf_pca$x,
    aes(
      x = PC1, y = PC2,
      color = data$`Soil order`, shape = data$PM
    ),
    size = 3,
    alpha = 0.6
  ) +
  scale_x_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
  scale_y_continuous(limits = function(x) c(-max(abs(x)), max(abs(x)))) +
  geom_segment(
    aes(x = 0, y = 0, xend = x_loadings, yend = y_loadings),
    arrow = arrow(length = unit(2, 'mm'), type = 'closed')
  ) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_text(
    aes(x = x_loadings, y = y_loadings, label = remove_pXRF(names(x_loadings))),
    hjust = ifelse(x_loadings > 0, 0, 1),
    vjust = 0.5,
    angle = atan(y_loadings / x_loadings) * (360 / (2*pi)),
    size = 3
  ) +
  guides(
    color = guide_legend(title = 'Soil order'),
    shape = guide_legend(title = 'PM')
  )

ggsave('figures/pca.png', units = 'mm', dpi = 300, width = 150, height = 130)
knitr::plot_crop('figures/pca.png')

x_loading_ranks <- data.frame(x = pxrf_pca$rotation[,1]) %>%
  mutate(sign = ifelse(x < 0, 'Negative', 'Positive')) %>%
  mutate(names = rownames(.), x = abs(x))
y_loading_ranks <- data.frame(x = pxrf_pca$rotation[,2]) %>%
  mutate(sign = ifelse(x < 0, 'Negative', 'Positive')) %>%
  mutate(names = rownames(.), x = abs(x))

pc1_loadings <- ggplot(data = x_loading_ranks) +
  xlab('PC1') + ylab('Elements') +
  geom_bar(
    aes(x = x, y = reorder(remove_pXRF(names), x), fill = sign),
    stat = 'identity'
  ) +
  guides(fill = guide_legend(title = ''))
pc2_loadings <- ggplot(data = y_loading_ranks) +
  xlab('PC2') + ylab('') +
  geom_bar(
    aes(x = x, y = reorder(remove_pXRF(names), x), fill = sign),
    stat = 'identity'
  ) +
  guides(fill = guide_legend(title = ''))
pc1_loadings + pc2_loadings + plot_layout(guides = 'collect')

ggsave('figures/loadings.png', units = 'mm', dpi = 300, width = 180, height = 120)
knitr::plot_crop('figures/loadings.png')

## creating the modeling dataset
## using results from pca

model_data <- data %>%
  mutate(PC1 = pxrf_pca$x[,1], PC2 = pxrf_pca$x[,2]) %>%
  select(all_of(c(target_vars, cat_vars, dry_color_vars, pca_vars))) %>%
  dplyr::mutate(color_id = factor(seq(1, nrow(data), 1), ordered = T))

# descriptive stats ---------------------------------------------------------- #

stats <- c(
  'mean', 'sd', 'min', 'q1',
  'med','q3', 'max', 'iqr', 'cv',
  'skewness', 'kurtosis'
)
desc_stats <- descr(model_data, order = 'p', stats = stats) %>%
  as.data.frame() %>%
  mutate(stat = rownames(.))

soil_table <- data %>%
  group_by(Symbol, SiBCS, `Soil taxonomy`, `Soil order`, PM) %>%
  summarize(n = n())

# color analysis ------------------------------------------------------------- #

hex_colors <- hex(
  LAB(model_data$LAB_L_dry, model_data$LAB_A_dry, model_data$LAB_B_dry),
)

LA_plot <- ggplot(data = model_data) +
  xlab('a* (Redness)') + ylab('L* (Brightness)') +
  geom_point(
    aes(x = LAB_A_dry, y = LAB_L_dry, color = color_id, shape = Depth, size = Lab_SOM),
  ) +
  scale_color_manual(values = hex_colors, guide = 'none') +
  guides(
    shape = guide_legend(title = 'Depth (cm)'),
    size = guide_legend(title = 'SOC (%)')
  )

LB_plot <- ggplot(data = model_data) +
  xlab('b* (Yellowness)') + ylab('') +
  geom_point(
    aes(x = LAB_B_dry, y = LAB_L_dry, color = color_id, shape = Depth, size = Lab_SOM),
  ) +
  scale_color_manual(values = hex_colors, guide = 'none') +
  guides(
    shape = guide_legend(title = 'Depth (cm)'),
    size = guide_legend(title = 'SOC (%)')
  )
LA_plot + LB_plot + plot_layout(guides = 'collect')

ggsave('figures/colors.png', units = 'mm', dpi = 300, width = 200, height = 120)
knitr::plot_crop('figures/colors.png')

# study area map ------------------------------------------------------------- #

area_bound <- sf::read_sf('data/study_area_bound.shp')

plot_data <- model_data %>%
  bind_cols(data[c('Lat', 'Long', 'Soil order')]) %>%
  group_by(Depth) %>%
  mutate(depth_count = n()) %>% 
  mutate(Depth = paste0(Depth, ' cm (n=', depth_count,')'))

ggplot(plot_data) +
  xlab('') + ylab('') +
  geom_sf(data = area_bound) +
  coord_sf(datum = sf::st_crs(32723)) +
  geom_point(aes(x = Long, y = Lat, color = `Soil order`, shape = PM)) +
  facet_wrap(~Depth, ) +
  annotation_scale() +
  annotation_north_arrow(
    height = unit(8, 'mm'),
    width = unit(8, 'mm'),
    pad_x = unit(6, 'mm'),
    pad_y = unit(6, 'mm'),
    style = north_arrow_fancy_orienteering()
  ) +
  scale_size(range = c(0.5, 4)) +
  scale_x_continuous(breaks = seq(501500, 504000, by = 1000)) +
  guides(size = guide_legend(title = 'SOC (%)'))

ggsave('figures/map.png', units = 'mm', dpi = 300, width = 230, height = 180)
knitr::plot_crop('figures/map.png')

# dry color ------------------------------------------------------------------ #

## stepwise algorithm -------------------------------------------------------- #

## using the stepwise algorithm to select variables relevant to
## each distributional parameter (takes some time to run)


dists <- fitDist(log(model_data$Lab_SOM), type = 'realplus', k = 2)
dists$fits

params_form <- formula(
  ~ Depth + pb(LAB_L_dry) + pb(LAB_A_dry) + pb(LAB_B_dry) + pb(PC1) + pb(PC2)
)

set.seed(200)
m_dry <- gamlss(
  Lab_SOM ~ 1,
  data = model_data,
  family = BCPEo,
  method = RS(100)
)

nc <- detectCores()
set.seed(200)
m_dry_step <- stepGAICAll.A(
  m_dry,
  scope = list(
    lower = ~ 1,
    upper = params_form
  ),
  k = 2,
  parallel = 'snow',
  ncpus = nc
)

summary(m_dry_step)
resid_plots(m_dry_step, value = 3)
resid_wp(m_dry_step)
moment_bucket(m_dry_step)

par(mfrow = c(3,3))
term.plot(m_dry_step, pages = 1)
term.plot(m_dry_step, what = 'sigma', pages = 1)
term.plot(m_dry_step, what = 'nu', pages = 1)
term.plot(m_dry_step, what = 'tau', pages = 1)

## final model --------------------------------------------------------------- #

## creating the definitive model based on what was learned
## from the stepwise analysis (it should run pretty fast)

set.seed(200)
m_dry <- gamlss(
  Lab_SOM ~ Depth + LAB_L_dry + LAB_A_dry + pb(LAB_B_dry) + PC1,
  sigma.fo = ~ Depth +  LAB_A_dry + LAB_B_dry + PC1,
  nu.fo = ~ PC2,
  tau.fo = ~ pb(PC2),
  data = model_data,
  family = BCPEo,
  method = RS(500)
)

summary(m_dry)
resid_plots(m_dry, value = 3)
resid_wp(m_dry, ylim = 1)
moment_bucket(m_dry)

term.plot(m_dry, pages = 1)
term.plot(m_dry, what = 'sigma', pages = 1)
term.plot(m_dry, what = 'nu', pages = 1)
term.plot(m_dry, what = 'tau', pages = 1)

# model plots ---------------------------------------------------------------- #

## plotting smoothing functions

b_smo <- lpred(m_dry, what = 'mu', type = 'terms', terms = 4, se.fit = T)
b_smo$se.low <- b_smo$fit - b_smo$se.fit
b_smo$se.up <- b_smo$fit + b_smo$se.fit

pc2_smo <- lpred(m_dry, what = 'tau', type = 'terms', terms = 1, se.fit = T)
pc2_smo$se.low <- pc2_smo$fit - pc2_smo$se.fit
pc2_smo$se.up <- pc2_smo$fit + pc2_smo$se.fit

b_label <- 'pb(b*), fitted for the median (\u03BC)'
pc2_label <- 'pb(PC2), fitted for kurtosis (\u03C4)'
smo_funs <- data.frame(
  x = c(model_data$LAB_B_dry, model_data$PC2),
  term = c(rep(b_label, length(model_data$LAB_B_dry)), rep(pc2_label, length(model_data$PC2))),
  smo = c(b_smo$fit, pc2_smo$fit),
  smo_low = c(b_smo$se.low, pc2_smo$se.low),
  smo_up = c(b_smo$se.up, pc2_smo$se.up)
)

ggplot(smo_funs) +
  xlab('') + ylab('Partial function') +
  geom_rug(aes(y = 0, x = x), sides = 'b', alpha = 0.3) +
  geom_line(aes(x = x, y = smo)) +
  geom_ribbon(aes(x = x, ymin = smo_low, ymax = smo_up), alpha = 0.5) +
  facet_wrap(~term, scales = 'free')

ggsave('smoothing_funs.png', units = 'mm', dpi = 300, width = 150, height = 90)
knitr::plot_crop('smoothing_funs.png')

## plotting wormp lots and bucket plots

worm_p <- resid_wp(
  m_dry,
  ylim = 1,
  title = 'b) Worm-plot'
)
bucket_p <- moment_bucket(m_dry, text_to_show = '') +
  ggtitle('c) Moment bucket') +
  xlab('Transformed moment skewness') +
  ylab('Transformed moment excess kurtosis')
resid_p <- resid_mu(m_dry, value = 3, title = 'd) Residuals vs. fitted') +
  xlab('Fitted values (SOC, %)') +
  ylab('Quantile residuals')
qq_p <- resid_qqplot(m_dry, value = 3, title = 'a) QQ-plot')
(qq_p + worm_p) / (bucket_p + resid_p)

ggsave('figures/residuals.png', units = 'mm', dpi = 300, width = 180, height = 180)
knitr::plot_crop('figures/residuals.png')
