#' sample_lenghts
#'
#' @param n_at_age
#' @param cv
#' @param k
#' @param linf
#' @param t0
#' @param sample_type
#' @param percent_sampled
#'
#' @return a data frame of sampled length frequencies
#' @export
#'
sample_lengths <-
  function(n_at_age,
           cv,
           k,
           linf,
           t0,
           time_step = 1,
           sample_type = 'catch',
           percent_sampled = .1,
           sample_col = NA,
           linf_buffer = 10) {
    if (sample_type == 'catch') {
      sample_col <- quo(numbers_caught)

    } else if (sample_type == 'population') {
      sample_col <- quo(numbers)
    }

    length_comp_samples <- n_at_age %>%
      summarise(samps = percent_sampled * (sum(!!sample_col))) %>%  {
        .$samps
      }

    min_age <- min(n_at_age$age)

    max_age <- max(n_at_age$age)

    mean_length_at_age <-
      linf * (1 - exp(-k * (seq(
        min_age, max_age, by = time_step
      ) - t0)))
    p_n_at_age <- n_at_age %>%
      group_by(age) %>%
      summarise(numbers = sum(!!sample_col)) %>%
      ungroup() %>%
      mutate(p_sampled_at_age = numbers / sum(numbers))

    length_at_age_vars <- data_frame(
      age = seq(min_age, max_age, by = time_step),
      mean_length_at_age = mean_length_at_age,
      sigma_at_age = cv * mean_length_at_age
    ) #calculate standard deviation of length at age for each age bin

    # now calculate the probability of being in each length bin at each age
    p_length_at_age <-
      expand.grid(
        age = seq(min_age, max_age, by = time_step),
        length_bin = 0:(linf_buffer * linf)
      ) %>%
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

    p_length_at_age <- p_length_at_age %>%
      left_join(p_n_at_age, by = 'age')

    p_sampling_length_bin <- p_length_at_age %>%
      group_by(length_bin) %>%
      summarise(prob_sampled = sum(p_bin * p_sampled_at_age))

    if (length_comp_samples > 0) {
      length_comps <-
        rmultinom(1, size = length_comp_samples, prob = p_sampling_length_bin$prob_sampled) %>% as.numeric()

      # p_sampling_length_bin %>%
      #   ggplot(aes(length_bin, prob_sampled)) +
      #   geom_point() +
      #   geom_vline(data = data_frame(lata = fish$length_at_age %>% floor()), aes(xintercept =lata))
      length_comps <-
        data_frame(length_bin = unique(p_length_at_age$length_bin),
                   numbers = length_comps)

    } else {
      length_comps <-
        data_frame(length_bin = unique(p_length_at_age$length_bin),
                   numbers = 0)

    }
    return(length_comps)

  }