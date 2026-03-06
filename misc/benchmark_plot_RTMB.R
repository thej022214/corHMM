library(tidyverse)
theme_set(theme_bw())

## exploring differences ...
dd <- readRDS("benchmark2.rds")
if (!is.data.frame(dd)) dd <- do.call(rbind, dd)
print(nrow(dd))
tail(dd)

with(dd, summary(RTMB_loglik-orig_loglik))
dd2 <- dd |>
  mutate(ldiff = RTMB_loglik-orig_loglik) |>
  arrange(ldiff) |>
  mutate(n = seq(n()))

## filter(dd2, ldiff<0) |> View()

ggplot(dd2, aes(n, ldiff)) + geom_point()

## will have to hack corHMM further: store convergence codes from optimizer?
## explore: which cases are bad?
## worst for smallest/overparameterized cases (3 traits/20 taxa, 2 traits/41 taxa)

## for benchmark-1, 'model' was ignored
ddt <- (dd
  |> select(c(seed, ntrait, ntaxa, model, ends_with("time")))
  |> pivot_longer(ends_with("time"))
  |> mutate(across(name, ~ stringr::str_remove(., "\\.time$")))
  |> separate(name, into = c("method", "type"), sep = "_")
)

dd_bad <- (dd |>
           mutate(ldiff = RTMB_loglik-orig_loglik,
                  bad = case_when(abs(ldiff)<(1e-2) ~ "OK",
                                  ldiff < 0 ~ "RTMB",
                                  ldiff > 0 ~ "orig")) |>
           select(c(seed, ntrait, ntaxa, model, bad, ldiff))
)

ddtb <- full_join(ddt, dd_bad,
                  by = c("seed", "ntrait", "ntaxa", "model"))

ggplot(ddt, aes(ntaxa, value, colour = method)) +
  geom_point(aes(shape = factor(ntrait)), size = 8) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~type) +
  geom_smooth(aes(linetype = factor(ntrait)))
## ntrait = 3 is *faster* on average?? hmmm ...


ggplot(ddt, aes(ntaxa, value, colour = interaction(model, ntrait))) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~type) +
  geom_smooth(aes(group=interaction(method, model, ntrait),
                  fill=interaction(model, ntrait)))


ggplot(filter(ddtb, type == "opt"), aes(ntaxa, value)) +
  geom_point(aes(colour = bad, size = 0.2 + abs(ldiff))) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(adjustcolor("grey", alpha.f = 0.2),
                                 palette()[c(4, 2)])) +
  scale_size(range = c(3, 10))


