#' length_to_age
#'
#' @param length_samples
#' @param cv
#' @param k
#' @param linf
#' @param t0
#' @param max_age
#'
#' @return a data frame of estimated numbers at age
#' @export
#'
length_to_age <-
  function(length_samples, cv, k, linf, t0, max_age, min_age = 1) {
    lengths <- length_samples %>% {
      map2(.$length_bin, .$numbers, ~ rep(.x, .y))
    } %>%
      unlist()

    mean_length_at_age <- linf * (1 - exp(-k * ((min_age:max_age) - t0)))

    length_at_age_vars <- data_frame(
      age = min_age:max_age,
      mean_length_at_age = mean_length_at_age,
      sigma_at_age = cv * mean_length_at_age
    ) #calculate standard deviation of length at age for each age bin

    # now calculate the probability of being in each length bin at each age

    p_length_at_age <-
      expand.grid(age = min_age:max_age, length_bin = 0:(1.5 * linf)) %>%
      as_data_frame() %>%
      left_join(length_at_age_vars, by = 'age') %>%
      arrange(age, length_bin)

    p_length_at_age <- p_length_at_age %>%
      group_by(age) %>%
      mutate(next_length_bin = lead(length_bin, 1)) %>%
      mutate(p_bin = ifelse(
        is.na(next_length_bin) == F,
        pnorm(next_length_bin, mean_length_at_age, sigma_at_age),
        1
      ) -
        pnorm(length_bin, mean_length_at_age, sigma_at_age))

    #rescale probabilities by the probability of being in an age bin at a given length``
    p_length_at_age <- p_length_at_age %>%
      group_by(length_bin) %>%
      mutate(p_age_at_length = p_bin / sum(p_bin, na.rm = T))

    # p_length_at_age %>%
    #   ggplot(aes(age, p_age_at_length)) +
    #   geom_col() +
    #   facet_wrap(~length_bin)
    #

    # Bin lengths into age bins, and then join the probability of being in each age as a function
    # of length size
    lengths_to_ages <- data_frame(lengths = lengths) %>%
      mutate(length_bin = floor(lengths)) %>%
      group_by(length_bin) %>%
      summarise(samples = length(lengths)) %>%
      left_join(p_length_at_age %>% select(length_bin, age, p_age_at_length),
                by = 'length_bin')


    # Assign lengths to ages in proportion to the probability of age at length

    age_comp <- lengths_to_ages %>%
      nest(-length_bin) %>%
      mutate(numbers = map(data,  ~ data_frame(
        age = .x$age,
        numbers = rmultinom(1, unique(.x$samples), .x$p_age_at_length) %>% as.numeric()
      ))) %>%
      select(numbers) %>%
      unnest() %>%
      group_by(age) %>%
      summarise(numbers = sum(numbers, na.rm = T)) %>%
      mutate(age = age)

    return(age_comp)

  }