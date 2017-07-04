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
  function(n_at_age, cv, k, linf, t0,sample_type = 'catch',
           percent_sampled = .1) {

    if (sample_type == 'catch'){

      sample_col <- quo(numbers_caught)

    } else if (sample_type == 'population'){

      sample_col <- quo(numbers)
    }

    length_comp_samples <- n_at_age %>%
      summarise(samps = percent_sampled * (sum(!!sample_col))) %>%  {
        .$samps
      }

    max_age <- max(n_at_age$age) - 1

    mean_length_at_age <- linf * (1 - exp(-k * ((0:max_age) - t0)))

    n_at_age$age <- n_at_age$age - 1

    p_n_at_age <- n_at_age %>%
      group_by(age) %>%
      summarise(numbers = sum(!!sample_col)) %>%
      ungroup() %>%
      mutate(p_sampled_at_age = numbers / sum(numbers))

    length_at_age_vars <- data_frame(
      age = 0:max_age,
      mean_length_at_age = mean_length_at_age,
      sigma_at_age = cv * mean_length_at_age
    ) #calculate standard deviation of length at age for each age bin

    # now calculate the probability of being in each length bin at each age

    p_length_at_age <-
      expand.grid(age = 0:max_age, length_bin = 0:(1.5*linf)) %>%
      as_data_frame() %>%
      left_join(length_at_age_vars, by = 'age') %>%
      arrange(age, length_bin)

    p_length_at_age <- p_length_at_age %>%
      group_by(age) %>%
      mutate(next_length_bin = lead(length_bin, 1)) %>%
      mutate(
        p_bin = ifelse(is.na(next_length_bin) == F, pnorm(next_length_bin, mean_length_at_age, sigma_at_age), 1) -
          pnorm(length_bin, mean_length_at_age, sigma_at_age)
      )

    p_length_at_age %>%
      ggplot(aes(length_bin, p_bin, fill = factor(age))) +
      geom_density(stat ='identity', alpha = 0.75)

    p_length_at_age <- p_length_at_age %>%
      left_join(p_n_at_age, by = 'age')

    p_sampling_length_bin <- p_length_at_age %>%
      group_by(length_bin) %>%
      summarise(prob_sampled = sum(p_bin * p_sampled_at_age))

    p_sampling_length_bin %>%
      ggplot(aes(length_bin, prob_sampled)) +
      geom_col()

    length_comps <- rmultinom(1, size = length_comp_samples, prob = p_sampling_length_bin$prob_sampled) %>% as.numeric()

    length_comps <- data_frame(length_bin = unique(p_length_at_age$length_bin),
                               numbers = length_comps)

    return(length_comps)

  }