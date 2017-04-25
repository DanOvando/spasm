# library(tidyverse)

# lhi <- read_csv("~/Box Sync/Databases/prince-lhi/prince-lhi.csv")

num_foo <- function(x) {

  if (any(is.na(as.numeric(x)) == F)){

    x <- as.numeric(x)

  } else {
    x <- x
  }

  return(x)
}

lhi <- lhi %>%
  map_df(num_foo) %>%
  mutate(taxa = map_chr(taxa, ~str_split(.x, ' ', simplify = T)[1])) %>%
  separate(species,c('genus','species'), sep = ' ') %>%
  mutate(taxa = paste(genus, species)) %>%
  mutate(m_v_k = m / k,
         lmat_v_linf = l50 / linf,
         tmat = (log(1 - l50/linf) / -k) + 0,
         m_times_tmat = m * tmat
  )


devtools::


