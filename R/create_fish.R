#' create_fish creates a fish list object with all the life history goodies
#'
#' @param common_name
#' @param scientific_name
#' @param linf
#' @param vbk
#' @param t0
#' @param max_age
#' @param weight_a
#' @param weight_b
#' @param length_50_mature
#' @param length_95_mature
#' @param age_50_mature
#' @param age_95_mature
#' @param age_mature
#' @param length_mature
#' @param m
#' @param steepness
#' @param density_dependence_form
#' @param adult_movement
#' @param larval_movement
#' @param query_fishbase
#' @param lmat_to_linf_ratio
#' @param r0
#'
#' @return a fish list object
#' @export
#'
#' @examples white_seabass = create_fish(scientific_name = "Atractoscion nobilis")
#'
create_fish <- function(common_name = 'white seabass',
                        scientific_name = "Atractoscion nobilis",
                        linf = NA,
                        vbk = NA,
                        t0 = NA,
                        length_units = 'cm',
                        max_age = 20,
                        weight_a = NA,
                        weight_b = NA,
                        weight_units = 'kg',
                        length_50_mature = NA,
                        length_95_mature = NA,
                        age_50_mature = NA,
                        age_95_mature = NA,
                        age_mature = NA,
                        length_mature = NA,
                        m = 0.2,
                        steepness = 0.8,
                        r0 = 1000,
                        density_dependence_form = 1,
                        adult_movement = 2,
                        larval_movement = 2,
                        query_fishbase = T,
                        lmat_to_linf_ratio = 0.6) {
  fish <- list()


  # check fishbase -------------

  if (is.na(scientific_name) == F & query_fishbase == T) {
    # scientific_name <- tolower(scientific_name)

    taxonomy <-
      str_split(scientific_name, pattern = ' ', simplify = T) %>%
      as_data_frame() %>%
      setNames(c('genus', 'species')) %>%
      mutate(sci_name = paste(genus, species))

    fishbase_matches <- fishbase %>%
      setNames(nm = tolower(colnames(.))) %>%
      mutate(genus = (genus),
             species = (species)) %>%
      mutate(
        genus_match = str_detect(genus, taxonomy$genus),
        species_match = str_detect(species, taxonomy$species),
        complete_match = genus_match == T & species_match == T,
        full_name = paste(genus, species)
      ) %>%
      filter(genus_match == T)

    if (any(fishbase_matches$complete_match) == T) {
      fishbase_matches <- fishbase_matches %>%
        filter(complete_match == T)

    }

    fb_life_history <- popgrowth(fishbase_matches$full_name,
                                 fields = c('Loo', 'K', 'to', 'M'))

    # a= poplw(fishbase_matches$full_name)
    fb_length_weight <-
      poplw(fishbase_matches$full_name, fields = c('a', 'b'))

    fb_length_weight$a <-  fb_length_weight$a / 1000 # convert to kg
  }


  # process lengths ---------------------------------------------------------

  if (is.na(linf)) {
    linf <- fb_life_history$Loo

  }

  if (is.na(vbk)) {
    vbk <- fb_life_history$K

  }
  if (is.na(t0)) {
    t0 <- fb_life_history$to

  }

  if (any(c(is.na(linf), is.na(vbk), is.na(t0)))) {
    stop("Not enough Von Bert Data")

  }

  fish$length_at_age <- linf * (1 - exp(-vbk * ((1:max_age - t0))))

  # process weight

  if (is.na(weight_a)) {
    weight_a <- fb_length_weight$a

  }
  if (is.na(weight_b)) {
    weight_b <- fb_length_weight$b

  }

  fish$weight_at_age <- weight_a * fish$length_at_age ^ weight_b

  # process maturity
  if ((is.na(age_50_mature) |
       is.na(age_95_mature)) & is.na(age_mature) == F) {
    age_50_mature <- age_mature

    age_95_mature <- 1.01 * age_50_mature

  } else if (is.na(age_mature)) {
    if (is.na(length_mature)) {
      length_mature <-  linf * lmat_to_linf_ratio
    }

    age_mature <- (log(1 - length_mature / linf) / -vbk) + t0

    age_50_mature <- age_mature

    age_95_mature <- 1.01 * age_50_mature
  }

  fish$maturity_at_age <-
    ((1 / (1 + exp(-log(
      19
    ) * (((1:max_age) - age_50_mature) / (age_95_mature - age_50_mature)
    )))))

  fish$ssb_at_age <- fish$maturity_at_age * fish$weight_at_age

  fish$linf <- linf
  fish$vbk  <-  vbk
  fish$t0 <-  t0
  fish$max_age <-  max_age
  fish$weight_a <-  weight_a
  fish$weight_b <-  weight_b
  fish$length_50_mature <-  length_50_mature
  fish$length_95_mature <-  length_95_mature
  fish$age_50_mature <-  age_50_mature
  fish$age_95_mature <-  age_95_mature
  fish$age_mature <-  age_mature
  fish$length_mature <-  length_mature
  fish$m <-  m
  fish$steepness <- steepness
  fish$r0 <- r0
  fish$density_dependence_form = density_dependence_form
  fish$adult_movement <-  adult_movement
  fish$larval_movement <-  larval_movement
  fish$lmat_to_linf_ratio <-  lmat_to_linf_ratio
  fish$length_units <-  length_units
  fish$weight_units <-  weight_units

  return(fish)
}