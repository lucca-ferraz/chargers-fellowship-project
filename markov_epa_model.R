library(nflfastR)
library(tidyverse)
# install.packages("markovchain")
library(markovchain)

# loads all play-by-play data since 1999
pbp_data <- load_pbp(TRUE)

pbp_data <- pbp_data |> 
  filter(!is.na(fixed_drive_result), !is.na(down), !is.na(yardline_100), !is.na(ydstogo),
         !is.na(posteam_timeouts_remaining), !is.na(half_seconds_remaining)) |> 
  filter(posteam_timeouts_remaining != -1) |> 
  # bin continuous variables for states
  mutate(
         yardline_100_group = case_when(
           # Backed up (own 10)
           yardline_100 >= 95 ~ "96-99",
           yardline_100 >= 90 ~ "90-95",
           # Red zone (25 and in) — finer bins
           yardline_100 <= 5  ~ "1-5",
           yardline_100 <= 10 ~ "6-10",
           yardline_100 <= 15 ~ "11-15",
           yardline_100 <= 20 ~ "16-20",
           # Middle of field — 10-yard bins
           TRUE ~ paste0(
             floor((yardline_100 - 1) / 10) * 10 + 1, "-",
             floor((yardline_100 - 1) / 10) * 10 + 10)),
         ydstogo_group = case_when(
           down == 1 & ydstogo > 10 ~ "10+", # first and long
           down == 1 & ydstogo < 10 & yardline_100 > 10 ~ "<10", # first and short
           ydstogo <= 3 ~ "1-3", # short
           ydstogo <= 7 ~ "4-7", # medium
           ydstogo <= 10 ~ "8-10", # standard
           ydstogo <= 15 ~ "11-15", # long
           ydstogo <= 20 ~ "16-20", # longer
           ydstogo > 20 ~ ">20" # very long
         ), 
         time_group = case_when(
           half_seconds_remaining > 900 ~ "900+", #first/third quarters
           half_seconds_remaining > 300 ~ "300-900", # 15:00-5:00
           half_seconds_remaining > 120 ~ "120-300", # under 5 minutes
           half_seconds_remaining > 40  ~ "40-120", # 2 minute drill
           half_seconds_remaining > 10  ~ "10-40", # under 40 seconds
           TRUE                         ~ "<=10" # under 10 seconds
         ),
         timeout = ifelse(posteam_timeouts_remaining >= 1, 1, 0), # binary indicator
         state = paste0(down, "/", ydstogo_group, "/", yardline_100_group, "/", 
                        time_group, "/", timeout))

condensed_data <- pbp_data |> 
  select(play_id, game_id, drive, state, fixed_drive_result)

absorptions <- condensed_data |> 
  group_by(game_id, drive) |> 
  summarize(state = last(fixed_drive_result))

condensed_data <- condensed_data |> mutate(is_absorption = 0)
absorptions <- absorptions |> mutate(is_absorption = 1)

condensed_data <- condensed_data |> select(-fixed_drive_result)

full_data <- rbind(condensed_data |> select(-play_id), absorptions) |> 
  group_by(game_id, drive) |> 
  arrange(game_id, drive, is_absorption) |> 
  select(-is_absorption) |> 
  mutate(next_state = lead(state)) |> 
  ungroup() 

unique_states <- full_data |> 
  count(state)

# compute transition matrix
transition_matrix <- matrix(0, nrow = length(unique_states$state), 
                            ncol = length(unique_states$state), 
                            dimnames = list(unique_states$state, unique_states$state))

for (i in 1:nrow(full_data)) {
  current_state <- full_data$state[i]
  next_state <- full_data$next_state[i]
  
  if (!is.na(next_state)) {
    ci <- match(current_state, unique_states$state)
    ni <- match(next_state, unique_states$state)
    transition_matrix[ci, ni] <- transition_matrix[ci, ni] + 1
  }
}

for (i in 1:nrow(transition_matrix)) {
  transition_matrix[i, ] <- transition_matrix[i, ] / sum(transition_matrix[i, ])
}

transition_matrix[is.nan(transition_matrix)] <- 0

absorption_states <- unique_states$state[(nrow(unique_states) - 8):nrow(unique_states)]

for (state in absorption_states) {
  transition_matrix[state, state] <- 1
}

# create markov chain model
mc <- new("markovchain", states = unique_states$state, transitionMatrix = transition_matrix)

# calculate absorption probabilities
absorption_probs <- absorptionProbabilities(mc) 
absorption_probs_df <- as.data.frame(absorption_probs)
absorption_probs_df$state <- rownames(absorption_probs)

# calculate expected points
absorption_probs_df <- absorption_probs_df |> 
  mutate(markov_ep = 7 * Touchdown + 3 * `Field goal` - 2 * Safety - 7 * `Opp touchdown`)

absorption_data <- absorption_probs_df |> select(state, markov_ep)

# merge with play by play data
markov_pbp <- pbp_data |> 
  left_join(absorption_data, by = "state") |> 
  left_join(absorptions |> select(-is_absorption) |> rename(end_state = state), by = join_by(game_id, drive)) |>
  group_by(game_id, drive) |> 
  mutate(epa = ep - lag(ep), markov_epa = markov_ep - lag(markov_ep), abs_diff = abs(epa - markov_epa)) |> 
  select(play_id, game_id, posteam, defteam, desc, state, ep, markov_ep, epa, markov_epa, abs_diff, end_state)

mean(markov_pbp$abs_diff, na.rm = TRUE)
mean(abs(markov_pbp$ep - markov_pbp$markov_ep), na.rm = TRUE)

summary(lm(ep ~ markov_ep, data = markov_pbp))
summary(lm(epa ~ markov_epa, data = markov_pbp))

# filter to second and shorts
report_data <- pbp_data |> 
  left_join(absorption_data, by = "state") |> 
  group_by(game_id, drive) |> 
  mutate(epa = ep - lag(ep), markov_epa = markov_ep - lag(markov_ep), abs_diff = abs(epa - markov_epa)) |> 
  select(play_id, game_id, down, ydstogo, yardline_100, half_seconds_remaining, 
         shotgun, no_huddle, rush_attempt, pass_attempt, touchdown, 
         play_type, run_location, run_gap, pass_location, pass_length,
         markov_ep, epa, markov_epa, abs_diff) |> 
  ungroup() |> 
  filter(down == 2, ydstogo <= 3, play_type %in% c("run", "pass"))

# examine run-pass EPA splits
report_data |> 
  group_by(play_type) |> 
  summarise(count = n(), total_epa = sum(markov_epa, na.rm = TRUE),
            avg_epa = total_epa / count)

report_data |>
  group_by(play_type) |>
  summarise(
    mean_epa = mean(markov_epa, na.rm = TRUE),
    sd_epa = sd(markov_epa, na.rm = TRUE),
    p10 = quantile(markov_epa, 0.10, na.rm = TRUE),
    p25 = quantile(markov_epa, 0.25, na.rm = TRUE),
    median = quantile(markov_epa, 0.50, na.rm = TRUE),
    p75 = quantile(markov_epa, 0.75, na.rm = TRUE),
    p90 = quantile(markov_epa, 0.90, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

report_data |>
  filter(yardline_100 <= 5) |> 
  group_by(play_type) |>
  summarise(
    mean_epa = mean(markov_epa, na.rm = TRUE),
    sd_epa = sd(markov_epa, na.rm = TRUE),
    p10 = quantile(markov_epa, 0.10, na.rm = TRUE),
    p25 = quantile(markov_epa, 0.25, na.rm = TRUE),
    median = quantile(markov_epa, 0.50, na.rm = TRUE),
    p75 = quantile(markov_epa, 0.75, na.rm = TRUE),
    p90 = quantile(markov_epa, 0.90, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

report_data |>
  filter(yardline_100 > 5 & yardline_100 <= 25) |> 
  group_by(play_type) |>
  summarise(
    mean_epa = mean(markov_epa, na.rm = TRUE),
    sd_epa = sd(markov_epa, na.rm = TRUE),
    p10 = quantile(markov_epa, 0.10, na.rm = TRUE),
    p25 = quantile(markov_epa, 0.25, na.rm = TRUE),
    median = quantile(markov_epa, 0.50, na.rm = TRUE),
    p75 = quantile(markov_epa, 0.75, na.rm = TRUE),
    p90 = quantile(markov_epa, 0.90, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

report_data |>
  filter(yardline_100 > 90) |> 
  group_by(play_type) |>
  summarise(
    mean_epa = mean(markov_epa, na.rm = TRUE),
    sd_epa = sd(markov_epa, na.rm = TRUE),
    p10 = quantile(markov_epa, 0.10, na.rm = TRUE),
    p25 = quantile(markov_epa, 0.25, na.rm = TRUE),
    median = quantile(markov_epa, 0.50, na.rm = TRUE),
    p75 = quantile(markov_epa, 0.75, na.rm = TRUE),
    p90 = quantile(markov_epa, 0.90, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

report_data |>
  filter(yardline_100 > 40 & yardline_100 < 60) |> 
  group_by(play_type) |>
  summarise(
    mean_epa = mean(markov_epa, na.rm = TRUE),
    sd_epa = sd(markov_epa, na.rm = TRUE),
    p10 = quantile(markov_epa, 0.10, na.rm = TRUE),
    p25 = quantile(markov_epa, 0.25, na.rm = TRUE),
    median = quantile(markov_epa, 0.50, na.rm = TRUE),
    p75 = quantile(markov_epa, 0.75, na.rm = TRUE),
    p90 = quantile(markov_epa, 0.90, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
