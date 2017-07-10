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
#' @param cv_len
#' @param length_units
#' @param min_age
#' @param time_step
#' @param weight_units
#' @param delta_mature
#' @param lhi_type
#' @param price
#' @param sigma_r
#' @param rec_ac
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
                        t0 = 0,
                        cv_len = 0.1,
                        length_units = 'cm',
                        min_age = 0,
                        max_age = 20,
                        time_step = 1,
                        weight_a = NA,
                        weight_b = NA,
                        weight_units = 'kg',
                        length_50_mature = NA,
                        length_95_mature = NA,
                        delta_mature = .1,
                        age_50_mature = NA,
                        age_95_mature = NA,
                        age_mature = NA,
                        length_mature = NA,
                        m = 0.2,
                        steepness = 0.8,
                        r0 = 10000,
                        density_dependence_form = 1,
                        adult_movement = 2,
                        larval_movement = 2,
                        query_fishbase = F,
                        lhi_type = 1,
                        lmat_to_linf_ratio = NA,
                        price = 1,
                        sigma_r = 0,
                        rec_ac = 0) {

  lhi_groups <- lhi %>%
    group_by(type) %>%
    summarise(mean_m_v_k = mean(m_v_k, na.rm = T),
              mean_lmat_v_linf = mean(lmat_v_linf, na.rm = T),
              mean_m_by_tm = mean(m_times_tmat, na.rm = T),
              mean_wa = mean(lw_a, na.rm = T),
              mean_wb = mean(lw_b, na.rm = T))

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
    if (dim(fb_life_history)[1] == 0){

      fb_life_history <-  data_frame(Loo = NA, K = NA, to = NA, M = NA)

    }

    fb_life_history <- fb_life_history %>%
      # group_by(sciname) %>%
      summarise(Loo = median(Loo),
                K = median(K),
                to = median(to),
                M = median(M))


    # a= poplw(fishbase_matches$full_name)
    fb_length_weight <-
      poplw(fishbase_matches$full_name, fields = c('a', 'b'))

    if (dim(fb_length_weight)[1] == 0){

      fb_length_weight <-  data_frame(a = NA, b = NA)

    }

    fb_length_weight <- fb_length_weight %>%
      summarise(a = median(a),
                b = median(b),
                a = a / 1000)

    # fb_length_weight$a <-  fb_length_weight$a / 1000 # convert to kg


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

    if (is.na(weight_a)) {
      weight_a <- fb_length_weight$a

    }
    if (is.na(weight_b)) {
      weight_b <- fb_length_weight$b

    }

} #close fishbase query

  # if (any(c(is.na(linf), is.na(vbk), is.na(t0)))) {
  #   warning("Not enough Von Bert Data")
  # }
  max_age <- ((-log(0.01)/m)) %>% floor()

  if (is.na(vbk)){

    vbk <- m / (lhi_groups$mean_m_v_k[lhi_groups$type == lhi_type])

  }
  if (is.na(weight_a)){

    weight_a <-lhi_groups$mean_wa[lhi_groups$type == lhi_type]

    weight_b <-lhi_groups$mean_wb[lhi_groups$type == lhi_type]

  }

  fish$length_at_age <- linf * (1 - exp(-vbk * (seq(min_age,max_age, by = time_step) - t0)))

  # process weight

  fish$weight_at_age <- weight_a * fish$length_at_age ^ weight_b

 if (is.na(lmat_to_linf_ratio)) {

   lmat_to_linf_ratio <- lhi_groups$mean_lmat_v_linf[lhi_groups$type == lhi_type]

 }

  # process maturity
  if ((is.na(age_50_mature) |
       is.na(age_95_mature)) & is.na(age_mature) == F) {
    age_50_mature <- age_mature

    age_95_mature <-  age_50_mature + delta_mature

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
    ) * ((seq(min_age,max_age, by = time_step) - age_50_mature) / (age_95_mature - age_50_mature)
    )))))

  if (is.na(length_50_mature)){

    length_50_mature <- length_mature

    length_95_mature <- length_50_mature + delta_mature

  }

  fish$scientific_name <- scientific_name
  fish$common_name <- common_name
  fish$ssb_at_age <- fish$maturity_at_age * fish$weight_at_age
  fish$linf <- linf
  fish$vbk  <-  vbk
  fish$t0 <-  t0
  fish$cv_len <- cv_len
  fish$max_age <-  max_age
  fish$min_age <- min_age
  fish$weight_a <-  weight_a
  fish$weight_b <-  weight_b
  fish$length_50_mature <-  length_50_mature
  fish$length_95_mature <-  length_95_mature
  fish$age_50_mature <-  age_50_mature
  fish$age_95_mature <-  age_95_mature
  fish$delta_mature <- delta_mature
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
  fish$price <- price
  fish$sigma_r <- sigma_r
  fish$rec_ac <- rec_ac
  fish$time_step <- time_step

  return(fish)
}