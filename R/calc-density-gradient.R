#' calculate density gradient
#'
#' Calculates a modifier for the distance based adult movement matrix based
#' on biomass
#'
#' @param pop
#' @param y
#' @param num_patches
#' @param density_modifier
#'
#' @return a movement matrix modifier
#' @export
#'
calc_density_gradient <- function(pop, y, num_patches, density_modifier, b0_buffer = 1.25) {

adults <- pop %>%
  filter(year == y) %>%
  mutate(biomass = numbers * weight_at_age) %>%
  group_by(patch) %>%
  summarise(b = sum(biomass))

adults_zero <- pop %>%
  filter(year == min(pop$year)) %>%
  mutate(biomass = numbers * weight_at_age) %>%
  group_by(patch) %>%
  summarise(b0 = b0_buffer * sum(biomass))

adults <- adults %>%
  left_join(adults_zero, by = "patch") %>%
  mutate(density = pmin(1,(b / b0) * density_modifier)) %>%
  select(patch, density)

density_gradient <-
  expand.grid(
    source = 1:num_patches,
    sink = 1:num_patches
  ) %>%
  left_join(adults %>% rename(source_density = density), by = c("source" = "patch")) %>%
  left_join(adults %>% rename(sink_density = density), by = c("sink" = "patch")) %>%
  mutate(gradient = source_density - sink_density  + 1) %>%
  select(source, sink, gradient) %>%
  spread(sink, gradient) %>%
  select(-source) %>%
  as.matrix()


return(density_gradient)

}
